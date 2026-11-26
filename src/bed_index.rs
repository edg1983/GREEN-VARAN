use anyhow::{Context, Result};
use rust_htslib::tbx::{self, Read as TbxRead};
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::path::Path;

use crate::cli::ConnectionMode;
use crate::config::DbSchema;
use crate::intervals::{region_overlapped_by_sv, SvInterval};

/// Compact contiguous string pool: stores unique strings in one flat Vec<u8>.
#[derive(Debug, Default)]
pub struct StringPool {
    bytes: Vec<u8>,
    /// (byte offset, byte len)
    spans: Vec<(u32, u16)>,
}

impl StringPool {
    pub fn with_capacity(bytes_cap: usize, count_cap: usize) -> Self {
        Self {
            bytes: Vec::with_capacity(bytes_cap),
            spans: Vec::with_capacity(count_cap),
        }
    }

    #[inline]
    pub fn push(&mut self, s: &[u8]) -> u32 {
        let id = self.spans.len() as u32;
        let off = self.bytes.len() as u32;
        let len = s.len().min(u16::MAX as usize) as u16;
        self.bytes.extend_from_slice(&s[..len as usize]);
        self.spans.push((off, len));
        id
    }

    #[inline]
    pub fn get(&self, id: u32) -> &str {
        let (off, len) = self.spans[id as usize];
        let slice = &self.bytes[off as usize..(off as usize + len as usize)];
        // Safety: all inserted slices are validated UTF-8 (or lossy-converted) on load.
        unsafe { std::str::from_utf8_unchecked(slice) }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.spans.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.spans.is_empty()
    }
}

struct Interner {
    pool: StringPool,
    map: FxHashMap<Vec<u8>, u32>,
}

impl Interner {
    fn new(bytes_cap: usize, count_cap: usize) -> Self {
        Self {
            pool: StringPool::with_capacity(bytes_cap, count_cap),
            map: FxHashMap::default(),
        }
    }

    #[inline]
    fn intern(&mut self, s: &[u8]) -> u32 {
        if let Some(&id) = self.map.get(s) {
            return id;
        }
        let id = self.pool.push(s);
        self.map.insert(s.to_vec(), id);
        id
    }

    fn into_pool(self) -> StringPool {
        self.pool
    }
}

/// Compact Structure-of-Arrays index for a single chromosome of GREEN-DB (~12-15 MB max for chr1).
#[derive(Debug)]
pub struct ChromIndex {
    pub starts: Vec<i64>,
    pub stops: Vec<i64>,
    /// Monotonically non-decreasing prefix max of `stops` for O(log N) range/chunk seeding.
    pub prefix_max_stops: Vec<i64>,
    /// `f64::NAN` denotes missing ("NA").
    pub constraints: Vec<f64>,
    pub stdtypes: Vec<u8>,
    pub stdtype_pool: StringPool,

    /// Contiguous region ID bytes + (offset, len) per region.
    id_bytes: Vec<u8>,
    id_spans: Vec<(u32, u16)>,

    /// Interned DB sources flattened across regions.
    pub dbsource_pool: StringPool,
    dbsource_ids: Vec<u16>,
    dbsource_spans: Vec<(u32, u16)>,

    /// Interned genes flattened across regions (only columns required by `ConnectionMode`).
    pub gene_pool: StringPool,
    gene_ids: Vec<u32>,
    gene_spans: Vec<(u32, u16)>,

    /// Precomputed per-region VOI flag (O(1) check in hot loop).
    pub region_has_voi: Vec<bool>,
}

impl ChromIndex {
    pub fn empty() -> Self {
        Self {
            starts: Vec::new(),
            stops: Vec::new(),
            prefix_max_stops: Vec::new(),
            constraints: Vec::new(),
            stdtypes: Vec::new(),
            stdtype_pool: StringPool::default(),
            id_bytes: Vec::new(),
            id_spans: Vec::new(),
            dbsource_pool: StringPool::default(),
            dbsource_ids: Vec::new(),
            dbsource_spans: Vec::new(),
            gene_pool: StringPool::default(),
            gene_ids: Vec::new(),
            gene_spans: Vec::new(),
            region_has_voi: Vec::new(),
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.starts.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    #[inline]
    pub fn region_id(&self, idx: usize) -> &str {
        let (off, len) = self.id_spans[idx];
        let slice = &self.id_bytes[off as usize..(off as usize + len as usize)];
        unsafe { std::str::from_utf8_unchecked(slice) }
    }

    #[inline]
    pub fn region_stdtype_id(&self, idx: usize) -> u8 {
        self.stdtypes[idx]
    }

    #[inline]
    pub fn stdtype_str(&self, id: u8) -> &str {
        self.stdtype_pool.get(id as u32)
    }

    #[inline]
    pub fn region_dbsource_ids(&self, idx: usize) -> &[u16] {
        let (off, len) = self.dbsource_spans[idx];
        &self.dbsource_ids[off as usize..(off as usize + len as usize)]
    }

    #[inline]
    pub fn dbsource_str(&self, id: u16) -> &str {
        self.dbsource_pool.get(id as u32)
    }

    #[inline]
    pub fn region_gene_ids(&self, idx: usize) -> &[u32] {
        let (off, len) = self.gene_spans[idx];
        &self.gene_ids[off as usize..(off as usize + len as usize)]
    }

    #[inline]
    pub fn gene_str(&self, id: u32) -> &str {
        self.gene_pool.get(id)
    }

    /// Load a single chromosome from an indexed BGZF GREEN-DB BED file.
    pub fn load_chromosome(
        db_path: &Path,
        chrom: &str,
        schema: DbSchema,
        connection: ConnectionMode,
        voi_genes: &FxHashSet<String>,
        threads: usize,
    ) -> Result<Self> {
        let mut reader = tbx::Reader::from_path(db_path)
            .with_context(|| format!("Failed to open indexed GREEN-DB BED: {}", db_path.display()))?;
        if threads > 1 {
            let _ = reader.set_threads(threads.min(4));
        }

        let tid = match reader.tid(chrom) {
            Ok(t) => t,
            Err(_) => return Ok(Self::empty()),
        };

        // htslib's hts_pos_t is signed i64 (HTS_POS_MAX = i64::MAX >> 1); u64::MAX wraps to -1.
        reader
            .fetch(tid, 0, (i64::MAX >> 1) as u64)
            .with_context(|| format!("Failed to fetch chromosome {chrom} from GREEN-DB"))?;

        let mut builder = ChromIndexBuilder::new(schema, connection, voi_genes);
        let mut line_buf: Vec<u8> = Vec::with_capacity(512);
        while reader.read(&mut line_buf)? {
            builder.parse_line(&line_buf);
        }
        Ok(builder.finish())
    }

    /// Build a `ChromIndex` directly from BED lines in memory (useful for unit tests).
    pub fn from_lines<'a, I: IntoIterator<Item = &'a [u8]>>(
        lines: I,
        schema: DbSchema,
        connection: ConnectionMode,
        voi_genes: &FxHashSet<String>,
    ) -> Self {
        let mut builder = ChromIndexBuilder::new(schema, connection, voi_genes);
        for line in lines {
            builder.parse_line(line);
        }
        builder.finish()
    }

    /// Pruned range scan for SV intervals: finds all GREEN-DB regions passing
    /// `region_overlapped_by_sv` in BED-row order.
    pub fn query_sv(
        &self,
        sv: SvInterval,
        min_overlap: f64,
        min_bp: i64,
        hits: &mut Vec<u32>,
    ) {
        hits.clear();
        if self.is_empty() {
            return;
        }
        // First index where any region at or before it could reach past `sv.start`.
        let first = self.prefix_max_stops.partition_point(|&m| m <= sv.start);
        let n = self.starts.len();
        for i in first..n {
            let r_start = self.starts[i];
            if r_start >= sv.stop {
                break;
            }
            let r_stop = self.stops[i];
            if r_stop > sv.start
                && region_overlapped_by_sv(r_start, r_stop, sv, min_overlap, min_bp)
            {
                hits.push(i as u32);
            }
        }
    }
}

struct ChromIndexBuilder<'a> {
    schema: DbSchema,
    max_col: usize,
    connection: ConnectionMode,
    voi_genes: &'a FxHashSet<String>,
    has_voi: bool,

    starts: Vec<i64>,
    stops: Vec<i64>,
    constraints: Vec<f64>,
    stdtypes: Vec<u8>,
    stdtype_interner: Interner,

    id_bytes: Vec<u8>,
    id_spans: Vec<(u32, u16)>,

    dbsource_interner: Interner,
    dbsource_ids: Vec<u16>,
    dbsource_spans: Vec<(u32, u16)>,

    gene_interner: Interner,
    gene_is_voi: Vec<bool>,
    gene_ids: Vec<u32>,
    gene_spans: Vec<(u32, u16)>,

    region_has_voi: Vec<bool>,
    is_sorted: bool,
    prev_start: i64,
}

impl<'a> ChromIndexBuilder<'a> {
    fn new(
        schema: DbSchema,
        connection: ConnectionMode,
        voi_genes: &'a FxHashSet<String>,
    ) -> Self {
        Self {
            schema,
            max_col: schema.max_col(),
            connection,
            voi_genes,
            has_voi: !voi_genes.is_empty(),
            starts: Vec::with_capacity(65_536),
            stops: Vec::with_capacity(65_536),
            constraints: Vec::with_capacity(65_536),
            stdtypes: Vec::with_capacity(65_536),
            stdtype_interner: Interner::new(256, 16),
            id_bytes: Vec::with_capacity(1_048_576),
            id_spans: Vec::with_capacity(65_536),
            dbsource_interner: Interner::new(1024, 64),
            dbsource_ids: Vec::with_capacity(131_072),
            dbsource_spans: Vec::with_capacity(65_536),
            gene_interner: Interner::new(512 * 1024, 32_768),
            gene_is_voi: Vec::with_capacity(32_768),
            gene_ids: Vec::with_capacity(262_144),
            gene_spans: Vec::with_capacity(65_536),
            region_has_voi: Vec::with_capacity(65_536),
            is_sorted: true,
            prev_start: i64::MIN,
        }
    }

    #[inline]
    fn intern_gene(&mut self, bytes: &[u8]) -> (u32, bool) {
        let prev_len = self.gene_interner.pool.len();
        let id = self.gene_interner.intern(bytes);
        if self.gene_interner.pool.len() > prev_len {
            let is_v = if self.has_voi {
                std::str::from_utf8(bytes)
                    .map(|s| self.voi_genes.contains(s))
                    .unwrap_or(false)
            } else {
                false
            };
            self.gene_is_voi.push(is_v);
            (id, is_v)
        } else {
            (id, self.gene_is_voi[id as usize])
        }
    }

    fn parse_line(&mut self, raw: &[u8]) {
        // Strip trailing CR/LF
        let mut end = raw.len();
        while end > 0 && (raw[end - 1] == b'\n' || raw[end - 1] == b'\r') {
            end -= 1;
        }
        if end == 0 || raw[0] == b'#' {
            return;
        }
        let line = &raw[..end];

        // Split tabs into stack array for the common <= 32-column case
        let mut cols_arr: [&[u8]; 32] = [b""; 32];
        let mut cols_vec: Vec<&[u8]> = Vec::new();
        let mut ncol = 0usize;
        let mut field_start = 0usize;
        for (i, &b) in line.iter().enumerate() {
            if b == b'\t' {
                let slice = &line[field_start..i];
                if ncol < 32 {
                    cols_arr[ncol] = slice;
                } else {
                    if cols_vec.is_empty() {
                        cols_vec.extend_from_slice(&cols_arr);
                    }
                    cols_vec.push(slice);
                }
                ncol += 1;
                field_start = i + 1;
            }
        }
        let last = &line[field_start..];
        if ncol < 32 {
            cols_arr[ncol] = last;
        } else {
            if cols_vec.is_empty() {
                cols_vec.extend_from_slice(&cols_arr);
            }
            cols_vec.push(last);
        }
        ncol += 1;

        if ncol <= self.max_col || ncol < 3 {
            return;
        }
        let cols: &[&[u8]] = if cols_vec.is_empty() {
            &cols_arr[..ncol]
        } else {
            &cols_vec[..]
        };

        let start = match parse_i64_ascii(cols[1]) {
            Some(v) => v,
            None => return,
        };
        let stop = match parse_i64_ascii(cols[2]) {
            Some(v) => v,
            None => return,
        };
        if start > stop {
            return;
        }

        if start < self.prev_start {
            self.is_sorted = false;
        }
        self.prev_start = start;

        self.starts.push(start);
        self.stops.push(stop);

        // Constraint
        let c_col = cols[self.schema.constraint];
        let constraint = if c_col == b"NA" || c_col.is_empty() {
            f64::NAN
        } else {
            std::str::from_utf8(c_col)
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(f64::NAN)
        };
        self.constraints.push(constraint);

        // Region ID
        let id_col = ensure_utf8_or_empty(cols[self.schema.id]);
        let id_off = self.id_bytes.len() as u32;
        let id_len = id_col.len().min(u16::MAX as usize) as u16;
        self.id_bytes.extend_from_slice(&id_col[..id_len as usize]);
        self.id_spans.push((id_off, id_len));

        // stdtype
        let st_col = ensure_utf8_or_empty(cols[self.schema.stdtype]);
        let st_id = self.stdtype_interner.intern(st_col) as u8;
        self.stdtypes.push(st_id);

        // dbsource (comma-separated)
        let ds_off = self.dbsource_ids.len() as u32;
        let ds_col = ensure_utf8_or_empty(cols[self.schema.dbsource]);
        for part in ds_col.split(|&b| b == b',') {
            let d_id = self.dbsource_interner.intern(part) as u16;
            self.dbsource_ids.push(d_id);
        }
        let ds_len = (self.dbsource_ids.len() as u32 - ds_off).min(u16::MAX as u32) as u16;
        self.dbsource_spans.push((ds_off, ds_len));

        // genes according to ConnectionMode (matching exact Nim order in gdbinfo.nim)
        let g_off = self.gene_ids.len() as u32;
        let mut r_voi = false;
        match self.connection {
            ConnectionMode::All => {
                let g_col = ensure_utf8_or_empty(cols[self.schema.genes]);
                for part in g_col.split(|&b| b == b',') {
                    if !part.is_empty() {
                        let (gid, is_v) = self.intern_gene(part);
                        self.gene_ids.push(gid);
                        r_voi |= is_v;
                    }
                }
                let cl_col = ensure_utf8_or_empty(cols[self.schema.closest]);
                let (gid1, is_v1) = self.intern_gene(cl_col);
                self.gene_ids.push(gid1);
                r_voi |= is_v1;

                let cp_col = ensure_utf8_or_empty(cols[self.schema.closest_prot]);
                let (gid2, is_v2) = self.intern_gene(cp_col);
                self.gene_ids.push(gid2);
                r_voi |= is_v2;
            }
            ConnectionMode::Closest => {
                let cl_col = ensure_utf8_or_empty(cols[self.schema.closest]);
                let (gid1, is_v1) = self.intern_gene(cl_col);
                self.gene_ids.push(gid1);
                r_voi |= is_v1;

                let cp_col = ensure_utf8_or_empty(cols[self.schema.closest_prot]);
                let (gid2, is_v2) = self.intern_gene(cp_col);
                self.gene_ids.push(gid2);
                r_voi |= is_v2;
            }
            ConnectionMode::Annotated => {
                let g_col = ensure_utf8_or_empty(cols[self.schema.genes]);
                for part in g_col.split(|&b| b == b',') {
                    if !part.is_empty() {
                        let (gid, is_v) = self.intern_gene(part);
                        self.gene_ids.push(gid);
                        r_voi |= is_v;
                    }
                }
            }
        }
        let g_len = (self.gene_ids.len() as u32 - g_off).min(u16::MAX as u32) as u16;
        self.gene_spans.push((g_off, g_len));
        self.region_has_voi.push(r_voi);
    }

    fn finish(mut self) -> ChromIndex {
        let n = self.starts.len();
        if !self.is_sorted && n > 1 {
            // Stable sort by start ascending so sweep-line and prefix_max_stops are always valid.
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by_key(|&i| self.starts[i]);

            let reorder = |v: &Vec<i64>| -> Vec<i64> { order.iter().map(|&i| v[i]).collect() };
            self.starts = reorder(&self.starts);
            self.stops = reorder(&self.stops);
            self.constraints = order.iter().map(|&i| self.constraints[i]).collect();
            self.stdtypes = order.iter().map(|&i| self.stdtypes[i]).collect();
            self.id_spans = order.iter().map(|&i| self.id_spans[i]).collect();
            self.dbsource_spans = order.iter().map(|&i| self.dbsource_spans[i]).collect();
            self.gene_spans = order.iter().map(|&i| self.gene_spans[i]).collect();
            self.region_has_voi = order.iter().map(|&i| self.region_has_voi[i]).collect();
        }

        let mut prefix_max_stops = Vec::with_capacity(n);
        let mut running_max = i64::MIN;
        for &stop in &self.stops {
            if stop > running_max {
                running_max = stop;
            }
            prefix_max_stops.push(running_max);
        }

        ChromIndex {
            starts: self.starts,
            stops: self.stops,
            prefix_max_stops,
            constraints: self.constraints,
            stdtypes: self.stdtypes,
            stdtype_pool: self.stdtype_interner.into_pool(),
            id_bytes: self.id_bytes,
            id_spans: self.id_spans,
            dbsource_pool: self.dbsource_interner.into_pool(),
            dbsource_ids: self.dbsource_ids,
            dbsource_spans: self.dbsource_spans,
            gene_pool: self.gene_interner.into_pool(),
            gene_ids: self.gene_ids,
            gene_spans: self.gene_spans,
            region_has_voi: self.region_has_voi,
        }
    }
}

/// O(1)-amortized sweep-line cursor for sorted point queries in `smallvars` mode.
/// Exact boundary semantics matching Nim's `greendb.query(chrom, v.POS, v.POS+1)`
/// (where `v.POS` is 1-based and htslib `tbx_itr_queryi(tid, beg, end)` matches BED
/// `[start, stop)` iff `start < end && stop > beg`, i.e., `start <= pos_1based && stop > pos_1based`).
#[derive(Debug, Default)]
pub struct SmallVarCursor {
    cursor: usize,
    active: SmallVec<[u32; 32]>,
    last_pos: i64,
    initialized: bool,
}

impl SmallVarCursor {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self) {
        self.cursor = 0;
        self.active.clear();
        self.last_pos = i64::MIN;
        self.initialized = false;
    }

    fn seek(&mut self, index: &ChromIndex, pos_1based: i64) {
        self.active.clear();
        if index.is_empty() {
            self.cursor = 0;
            self.last_pos = pos_1based;
            self.initialized = true;
            return;
        }
        let first = index.prefix_max_stops.partition_point(|&m| m <= pos_1based);
        let n = index.starts.len();
        let mut i = first;
        while i < n && index.starts[i] <= pos_1based {
            if index.stops[i] > pos_1based {
                self.active.push(i as u32);
            }
            i += 1;
        }
        self.cursor = i;
        self.last_pos = pos_1based;
        self.initialized = true;
    }

    /// Query all GREEN-DB region indices overlapping 1-based `pos_1based`.
    /// Returned slice is in exact BED-row order and requires zero allocations.
    #[inline]
    pub fn query_point(&mut self, index: &ChromIndex, pos_1based: i64) -> &[u32] {
        if !self.initialized || pos_1based < self.last_pos {
            self.seek(index, pos_1based);
            return &self.active;
        }
        if pos_1based == self.last_pos {
            return &self.active;
        }

        // Remove regions that ended at or before `pos_1based` (preserve order).
        if !self.active.is_empty() {
            let stops = &index.stops;
            self.active.retain(|idx| stops[*idx as usize] > pos_1based);
        }

        // Advance cursor to include newly starting regions.
        let starts = &index.starts;
        let stops = &index.stops;
        let n = starts.len();
        let mut c = self.cursor;
        while c < n && starts[c] <= pos_1based {
            if stops[c] > pos_1based {
                self.active.push(c as u32);
            }
            c += 1;
        }
        self.cursor = c;
        self.last_pos = pos_1based;
        &self.active
    }
}

#[inline]
fn ensure_utf8_or_empty(bytes: &[u8]) -> &[u8] {
    if std::str::from_utf8(bytes).is_ok() {
        bytes
    } else {
        b""
    }
}

#[inline]
fn parse_i64_ascii(bytes: &[u8]) -> Option<i64> {
    if bytes.is_empty() {
        return None;
    }
    let mut neg = false;
    let mut i = 0;
    if bytes[0] == b'-' {
        neg = true;
        i = 1;
        if i == bytes.len() {
            return None;
        }
    }
    let mut acc: i64 = 0;
    while i < bytes.len() {
        let b = bytes[i];
        if !b.is_ascii_digit() {
            return None;
        }
        acc = acc.checked_mul(10)?.checked_add((b - b'0') as i64)?;
        i += 1;
    }
    Some(if neg { -acc } else { acc })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_schema() -> DbSchema {
        DbSchema {
            id: 3,
            stdtype: 4,
            dbsource: 5,
            constraint: 8,
            closest: 10,
            closest_prot: 12,
            genes: 14,
        }
    }

    #[test]
    fn test_smallvar_cursor_boundaries_and_seek() {
        let lines: Vec<&[u8]> = vec![
            b"#chromosome\tstart\tstop\tregionID\tstd_type\tDB_source\tN_s\tN_m\tconstrain_pct\tphylo\tclosest\td1\tclosestProt\td2\tgenes",
            b"chr1\t10000\t10600\t0_enh\tenhancer\tDECRES\t1\t1\tNA\tNA\tDDX11L1\t10\tOR4F5\t50\tGENE_A,GENE_B",
            b"chr1\t10244\t10357\t502_pro\tpromoter\tBENGI,ENCODE\t1\t0\t0.85\tNA\tDDX11L1\t10\tOR4F5\t50\t",
            b"chr1\t10600\t11000\t999_sil\tsilencer\tBENGI\t1\t0\t0.5\tNA\tG3\t0\tG4\t0\tG5",
        ];
        let mut voi = FxHashSet::default();
        voi.insert("GENE_B".to_string());
        let idx = ChromIndex::from_lines(lines, test_schema(), ConnectionMode::All, &voi);
        assert_eq!(idx.len(), 3);
        assert!(idx.region_has_voi[0]);
        assert!(!idx.region_has_voi[1]);

        let mut cursor = SmallVarCursor::new();
        // Before first region
        assert!(cursor.query_point(&idx, 9999).is_empty());
        // Exact start boundary (start <= pos && stop > pos)
        assert_eq!(cursor.query_point(&idx, 10000), &[0]);
        // Overlapping both 0 and 1
        assert_eq!(cursor.query_point(&idx, 10244), &[0, 1]);
        assert_eq!(cursor.query_point(&idx, 10356), &[0, 1]);
        // Exact stop boundary of region 1 (10357 excluded from region 1)
        assert_eq!(cursor.query_point(&idx, 10357), &[0]);
        // Exact boundary between region 0 [10000,10600) and region 2 [10600,11000)
        assert_eq!(cursor.query_point(&idx, 10600), &[2]);
        // Re-seek backward (e.g. out-of-order query)
        assert_eq!(cursor.query_point(&idx, 10300), &[0, 1]);
    }

    #[test]
    fn test_connection_modes() {
        let lines: Vec<&[u8]> = vec![
            b"chr1\t100\t200\tr1\tenhancer\tSRC\t1\t1\t0.9\t0\tCLOSEST_G\t0\tPROT_G\t0\tANN1,ANN2",
        ];
        let voi = FxHashSet::default();

        let idx_all = ChromIndex::from_lines(lines.clone(), test_schema(), ConnectionMode::All, &voi);
        let g_all: Vec<&str> = idx_all.region_gene_ids(0).iter().map(|&id| idx_all.gene_str(id)).collect();
        assert_eq!(g_all, vec!["ANN1", "ANN2", "CLOSEST_G", "PROT_G"]);

        let idx_cl = ChromIndex::from_lines(lines.clone(), test_schema(), ConnectionMode::Closest, &voi);
        let g_cl: Vec<&str> = idx_cl.region_gene_ids(0).iter().map(|&id| idx_cl.gene_str(id)).collect();
        assert_eq!(g_cl, vec!["CLOSEST_G", "PROT_G"]);

        let idx_an = ChromIndex::from_lines(lines, test_schema(), ConnectionMode::Annotated, &voi);
        let g_an: Vec<&str> = idx_an.region_gene_ids(0).iter().map(|&id| idx_an.gene_str(id)).collect();
        assert_eq!(g_an, vec!["ANN1", "ANN2"]);
    }
}
