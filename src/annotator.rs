use anyhow::{Context, Result};
use rust_htslib::bcf::{header::HeaderView, record::Numeric, Header, Record};
use rustc_hash::FxHashSet;
use smallvec::SmallVec;
use std::collections::HashSet;

use crate::bed_index::ChromIndex;
use crate::cli::{Impact, PrioritizationStrategy};
use crate::config::PrioritizeConfig;
use crate::csq::{format_ann_entries, CsqSchema, DESC_ANN};

pub const DESC_ID: &str = "GREENDB region IDs";
pub const DESC_STDTYPE: &str = "Regulatory region type as reported in GREENDB";
pub const DESC_DBSOURCE: &str = "Sources of this GREENDB annotation";
pub const DESC_GENES: &str = "Genes potentially controlled by regulatory regions according to GREENDB";
pub const DESC_CONSTRAINT: &str = "Maximum constraint value for GREENDB regions overlapping this variant";
pub const DESC_LEVEL: &str = "GREENDB prioritization level. 0:overlap, 1:0+rare, 2:1+TFBS-DNase-UCNE, 3:2+high score,4:3+high constraint";
pub const DESC_MORE_SUPPORT: &str = "GREENDB annotation";
pub const DESC_VOI: &str = "Variant affecting a region connected to a gene of interest in GREENDB";

/// Add GREEN-VARAN INFO definitions to the output VCF/BCF header.
pub fn update_header(
    header: &mut Header,
    csq_schema: Option<&CsqSchema>,
    skip_level: bool,
    has_voi_genes: bool,
) {
    if let Some(csq) = csq_schema {
        if csq.create_ann {
            let rec = format!(
                "##INFO=<ID={},Number=1,Type=String,Description=\"{}\">",
                csq.csq_field_name, DESC_ANN
            );
            header.push_record(rec.as_bytes());
        }
    }

    header.push_record(
        format!("##INFO=<ID=greendb_id,Number=1,Type=String,Description=\"{DESC_ID}\">").as_bytes(),
    );
    header.push_record(
        format!("##INFO=<ID=greendb_stdtype,Number=1,Type=String,Description=\"{DESC_STDTYPE}\">")
            .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID=greendb_dbsource,Number=1,Type=String,Description=\"{DESC_DBSOURCE}\">"
        )
        .as_bytes(),
    );
    header.push_record(
        format!("##INFO=<ID=greendb_genes,Number=1,Type=String,Description=\"{DESC_GENES}\">")
            .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID=greendb_constraint,Number=1,Type=Float,Description=\"{DESC_CONSTRAINT}\">"
        )
        .as_bytes(),
    );
    if !skip_level {
        header.push_record(
            format!("##INFO=<ID=greendb_level,Number=1,Type=Integer,Description=\"{DESC_LEVEL}\">")
                .as_bytes(),
        );
    }
    header.push_record(
        format!(
            "##INFO=<ID=greendb_more_support,Number=1,Type=Integer,Description=\"{DESC_MORE_SUPPORT}\">"
        )
        .as_bytes(),
    );

    if has_voi_genes {
        header.push_record(
            format!("##INFO=<ID=greendb_VOI,Number=0,Type=Flag,Description=\"{DESC_VOI}\">")
                .as_bytes(),
        );
    }
}

/// Pre-encoded byte tags for prioritization fields so no Strings are allocated in the hot loop.
#[derive(Debug, Clone)]
pub struct CompiledPrioritize {
    pub af: Vec<Vec<u8>>,
    pub maxaf: f32,
    pub regions: Vec<Vec<u8>>,
    pub scores: Vec<(Vec<u8>, f32)>,
    pub constraint: f64,
    pub more_regions: Vec<Vec<u8>>,
    pub more_values: Vec<(Vec<u8>, f32)>,
    pub pileup: bool,
}

impl CompiledPrioritize {
    pub fn from_config(cfg: &PrioritizeConfig, strategy: PrioritizationStrategy) -> Self {
        Self {
            af: cfg.af.iter().map(|s| s.as_bytes().to_vec()).collect(),
            maxaf: cfg.maxaf,
            regions: cfg.regions.iter().map(|s| s.as_bytes().to_vec()).collect(),
            scores: cfg
                .scores
                .iter()
                .map(|(k, &v)| (k.as_bytes().to_vec(), v))
                .collect(),
            constraint: cfg.constraint as f64,
            more_regions: cfg
                .more_regions
                .iter()
                .map(|s| s.as_bytes().to_vec())
                .collect(),
            more_values: cfg
                .more_values
                .iter()
                .map(|(k, &v)| (k.as_bytes().to_vec(), v))
                .collect(),
            pileup: strategy == PrioritizationStrategy::Pileup,
        }
    }

    /// Check that all required fields exist in header; returns list of missing field names.
    pub fn check_header(&self, header: &HeaderView) -> Vec<String> {
        let mut missing = Vec::new();
        for tag in &self.af {
            if header.info_type(tag).is_err() {
                missing.push(String::from_utf8_lossy(tag).into_owned());
            }
        }
        for tag in &self.regions {
            if header.info_type(tag).is_err() {
                missing.push(String::from_utf8_lossy(tag).into_owned());
            }
        }
        for (tag, _) in &self.scores {
            if header.info_type(tag).is_err() {
                missing.push(String::from_utf8_lossy(tag).into_owned());
            }
        }
        missing
    }
}

/// Reusable scratch space to avoid per-variant heap allocations.
#[derive(Default)]
pub struct VariantScratch {
    seen_stdtypes: SmallVec<[u8; 8]>,
    seen_dbsources: SmallVec<[u16; 16]>,
    seen_genes: SmallVec<[u32; 32]>,
    seen_genes_hash: FxHashSet<u32>,
    seen_id_strs: FxHashSet<String>,

    pub id_buf: String,
    pub stdtype_buf: String,
    pub dbsource_buf: String,
    pub genes_buf: String,
    pub ann_pairs: Vec<(u32, u8)>,
    pub new_ann_buf: String,
    pub combined_ann_buf: String,
    entry_buf: String,
    seen_ann_spans: Vec<(usize, usize)>,
    seen_ann_set: HashSet<String>,
}

/// Aggregate GREEN-DB hits into `scratch`, evaluate VOI/filter/prioritization,
/// and apply INFO + consequence updates directly onto `rec`.
///
/// Returns `(annotated, is_voi, should_write)`:
/// - `annotated`: true iff `hits` is non-empty
/// - `is_voi`: true iff overlapping regions hit at least one gene of interest
/// - `should_write`: false when `filter_mode` drops this variant
#[allow(clippy::too_many_arguments)]
pub fn annotate_record_in_place(
    rec: &mut Record,
    hits: &[u32],
    chrom_index: &ChromIndex,
    filter_mode: bool,
    has_voi_genes: bool,
    prioritize: Option<&CompiledPrioritize>,
    skip_level: bool,
    csq_schema: Option<&CsqSchema>,
    ann_tag_bytes: &[u8],
    impact: Impact,
    scratch: &mut VariantScratch,
) -> Result<(bool, bool, bool)> {
    if hits.is_empty() {
        return Ok((false, false, !filter_mode));
    }

    // Fast O(1) VOI check across overlapping regions (precomputed per region!)
    let mut is_voi = false;
    if has_voi_genes {
        for &ridx in hits {
            if chrom_index.region_has_voi[ridx as usize] {
                is_voi = true;
                break;
            }
        }
        if filter_mode && !is_voi {
            return Ok((true, false, false));
        }
    }

    // Aggregate GREEN-DB fields into scratch buffers
    let max_constraint = aggregate_hits(hits, chrom_index, csq_schema.is_some(), scratch);

    // Compute prioritization level (smallvars only)
    let (level, more_support) = if skip_level {
        (0, 0)
    } else if let Some(prio) = prioritize {
        compute_level(rec, prio, max_constraint)
    } else {
        (0, 0)
    };

    // Push greendb_* INFO tags (matching exact Nim field order)
    if !scratch.id_buf.is_empty() {
        rec.push_info_string(b"greendb_id", &[scratch.id_buf.as_bytes()])
            .context("Failed to set INFO/greendb_id")?;
    }
    if !scratch.stdtype_buf.is_empty() {
        rec.push_info_string(b"greendb_stdtype", &[scratch.stdtype_buf.as_bytes()])
            .context("Failed to set INFO/greendb_stdtype")?;
    }
    if !scratch.dbsource_buf.is_empty() {
        rec.push_info_string(b"greendb_dbsource", &[scratch.dbsource_buf.as_bytes()])
            .context("Failed to set INFO/greendb_dbsource")?;
    }
    if !scratch.genes_buf.is_empty() {
        rec.push_info_string(b"greendb_genes", &[scratch.genes_buf.as_bytes()])
            .context("Failed to set INFO/greendb_genes")?;
    }
    if let Some(mc) = max_constraint {
        rec.push_info_float(b"greendb_constraint", &[mc as f32])
            .context("Failed to set INFO/greendb_constraint")?;
    }
    if !skip_level {
        rec.push_info_integer(b"greendb_level", &[level])
            .context("Failed to set INFO/greendb_level")?;
    }
    rec.push_info_integer(b"greendb_more_support", &[more_support])
        .context("Failed to set INFO/greendb_more_support")?;

    if is_voi {
        rec.push_info_flag(b"greendb_VOI")
            .context("Failed to set INFO/greendb_VOI")?;
    }

    // Update ANN / CSQ / BCSQ consequence field
    if let Some(schema) = csq_schema {
        update_csq_in_place(rec, chrom_index, schema, ann_tag_bytes, impact, scratch)?;
    }

    Ok((true, is_voi, true))
}

/// Pure-read aggregation + prioritization computation (used by parallel batch workers
/// if applying updates on the writer thread).
pub fn aggregate_hits(
    hits: &[u32],
    chrom_index: &ChromIndex,
    collect_ann_pairs: bool,
    scratch: &mut VariantScratch,
) -> Option<f64> {
    scratch.id_buf.clear();
    scratch.stdtype_buf.clear();
    scratch.dbsource_buf.clear();
    scratch.genes_buf.clear();
    scratch.seen_stdtypes.clear();
    scratch.seen_dbsources.clear();
    scratch.seen_genes.clear();
    scratch.ann_pairs.clear();

    let mut max_constraint: Option<f64> = None;
    let use_id_hash = hits.len() > 16;
    if use_id_hash {
        scratch.seen_id_strs.clear();
    }
    let use_gene_hash = hits.len() > 8;
    if use_gene_hash {
        scratch.seen_genes_hash.clear();
    }

    for (h_pos, &ridx_u32) in hits.iter().enumerate() {
        let ridx = ridx_u32 as usize;

        // Region ID (deduplicated preserving BED order)
        let rid = chrom_index.region_id(ridx);
        let id_is_new = if use_id_hash {
            scratch.seen_id_strs.insert(rid.to_string())
        } else {
            // Check earlier hits directly (typically 1-5 regions -> super fast!)
            !hits[..h_pos]
                .iter()
                .any(|&prev| chrom_index.region_id(prev as usize) == rid)
        };
        if id_is_new {
            if !scratch.id_buf.is_empty() {
                scratch.id_buf.push(',');
            }
            scratch.id_buf.push_str(rid);
        }

        // stdtype (at most ~5 distinct values -> tiny linear scan)
        let st_id = chrom_index.region_stdtype_id(ridx);
        if !scratch.seen_stdtypes.contains(&st_id) {
            scratch.seen_stdtypes.push(st_id);
            if !scratch.stdtype_buf.is_empty() {
                scratch.stdtype_buf.push(',');
            }
            scratch.stdtype_buf.push_str(chrom_index.stdtype_str(st_id));
        }

        // constraint
        let c = chrom_index.constraints[ridx];
        if !c.is_nan() {
            max_constraint = Some(match max_constraint {
                Some(prev) => prev.max(c),
                None => c,
            });
        }

        // dbsources (at most ~15 distinct values -> tiny linear scan)
        for &ds_id in chrom_index.region_dbsource_ids(ridx) {
            if !scratch.seen_dbsources.contains(&ds_id) {
                scratch.seen_dbsources.push(ds_id);
                if !scratch.dbsource_buf.is_empty() {
                    scratch.dbsource_buf.push(',');
                }
                scratch
                    .dbsource_buf
                    .push_str(chrom_index.dbsource_str(ds_id));
            }
        }

        // genes + ann pairs
        for &g_id in chrom_index.region_gene_ids(ridx) {
            let is_new = if use_gene_hash {
                scratch.seen_genes_hash.insert(g_id)
            } else if !scratch.seen_genes.contains(&g_id) {
                scratch.seen_genes.push(g_id);
                true
            } else {
                false
            };
            if is_new {
                if !scratch.genes_buf.is_empty() {
                    scratch.genes_buf.push(',');
                }
                scratch.genes_buf.push_str(chrom_index.gene_str(g_id));
            }
            if collect_ann_pairs {
                scratch.ann_pairs.push((g_id, st_id));
            }
        }
    }

    max_constraint
}

pub fn compute_level(
    rec: &Record,
    prio: &CompiledPrioritize,
    max_constraint: Option<f64>,
) -> (i32, i32) {
    let mut eval_seq = [0i32; 4];
    let mut more_support = 0i32;

    // 0: AF check
    let mut has_af = false;
    let mut max_af = f32::NEG_INFINITY;
    for tag in &prio.af {
        if let Ok(Some(vals)) = rec.info(tag).float() {
            if let Some(&v) = vals.first() {
                if !v.is_missing() && !v.is_nan() {
                    has_af = true;
                    if v > max_af {
                        max_af = v;
                    }
                }
            }
        }
    }
    if !has_af || max_af < prio.maxaf {
        eval_seq[0] = 1;
    }

    // 1: Functional region overlap flags
    for tag in &prio.regions {
        if rec.info(tag).flag().unwrap_or(false) {
            eval_seq[1] = 1;
            break;
        }
    }

    // 2: Prediction scores
    for (tag, thresh) in &prio.scores {
        if let Ok(Some(vals)) = rec.info(tag).float() {
            if let Some(&v) = vals.first() {
                if !v.is_missing() && !v.is_nan() && v >= *thresh {
                    eval_seq[2] = 1;
                    break;
                }
            }
        }
    }

    // 3: Region constraint
    if let Some(mc) = max_constraint {
        if mc >= prio.constraint {
            eval_seq[3] = 1;
        }
    }

    // More regions
    for tag in &prio.more_regions {
        if rec.info(tag).flag().unwrap_or(false) {
            more_support += 1;
        }
    }

    // More values
    for (tag, thresh) in &prio.more_values {
        if let Ok(Some(vals)) = rec.info(tag).float() {
            if let Some(&v) = vals.first() {
                if !v.is_missing() && !v.is_nan() && v >= *thresh {
                    more_support += 1;
                }
            }
        }
    }

    let level = if prio.pileup {
        eval_seq[0] + eval_seq[1] + eval_seq[2] + eval_seq[3]
    } else {
        let mut lvl = 0;
        for v in eval_seq {
            if v > 0 {
                lvl += v;
            } else {
                break;
            }
        }
        lvl
    };

    (level, more_support)
}

fn update_csq_in_place(
    rec: &mut Record,
    chrom_index: &ChromIndex,
    schema: &CsqSchema,
    ann_tag_bytes: &[u8],
    impact: Impact,
    scratch: &mut VariantScratch,
) -> Result<()> {
    if scratch.ann_pairs.is_empty() {
        return Ok(());
    }

    // Format new ANN entries in an inner block so the immutable borrow of `rec` ends before `push_info_string`.
    let generated = {
        let all_alleles = rec.alleles();
        let alt_alleles: SmallVec<[&[u8]; 4]> = if all_alleles.len() > 1 {
            all_alleles[1..].iter().copied().collect()
        } else {
            SmallVec::new()
        };

        let pairs: SmallVec<[(&str, &str); 32]> = scratch
            .ann_pairs
            .iter()
            .map(|&(gid, stid)| (chrom_index.gene_str(gid), chrom_index.stdtype_str(stid)))
            .collect();

        scratch.new_ann_buf.clear();
        format_ann_entries(
            &alt_alleles,
            &pairs,
            impact.as_str(),
            schema,
            &mut scratch.entry_buf,
            &mut scratch.seen_ann_spans,
            &mut scratch.seen_ann_set,
            &mut scratch.new_ann_buf,
        )
    };

    if !generated || scratch.new_ann_buf.is_empty() {
        return Ok(());
    }

    scratch.combined_ann_buf.clear();
    if let Ok(Some(existing)) = rec.info(ann_tag_bytes).string() {
        for (i, part) in existing.iter().enumerate() {
            if i > 0 {
                scratch.combined_ann_buf.push(',');
            }
            if let Ok(s) = std::str::from_utf8(part) {
                scratch.combined_ann_buf.push_str(s);
            }
        }
        if !scratch.combined_ann_buf.is_empty() {
            scratch.combined_ann_buf.push(',');
        }
    }
    scratch.combined_ann_buf.push_str(&scratch.new_ann_buf);

    if !scratch.combined_ann_buf.is_empty() {
        rec.push_info_string(ann_tag_bytes, &[scratch.combined_ann_buf.as_bytes()])
            .with_context(|| format!("Failed to update INFO/{}", schema.csq_field_name))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::ConnectionMode;
    use crate::config::DbSchema;

    #[test]
    fn test_aggregate_hits_and_max_constraint() {
        let schema = DbSchema {
            id: 3,
            stdtype: 4,
            dbsource: 5,
            constraint: 8,
            closest: 10,
            closest_prot: 12,
            genes: 14,
        };
        let lines: Vec<&[u8]> = vec![
            b"chr1\t100\t200\tr1\tenhancer\tSRC1,SRC2\t1\t1\t0.4\t0\tG1\t0\tG2\t0\tG3",
            b"chr1\t150\t250\tr2\tpromoter\tSRC2,SRC3\t1\t1\t0.85\t0\tG2\t0\tG4\t0\t",
            b"chr1\t160\t260\tr3\tenhancer\tSRC1\t1\t1\tNA\t0\tG1\t0\tG2\t0\tG3",
        ];
        let voi = FxHashSet::default();
        let idx = ChromIndex::from_lines(lines, schema, ConnectionMode::All, &voi);
        let mut scratch = VariantScratch::default();
        let mc = aggregate_hits(&[0, 1, 2], &idx, true, &mut scratch);
        assert_eq!(mc, Some(0.85));
        assert_eq!(scratch.id_buf, "r1,r2,r3");
        assert_eq!(scratch.stdtype_buf, "enhancer,promoter");
        assert_eq!(scratch.dbsource_buf, "SRC1,SRC2,SRC3");
        assert_eq!(scratch.genes_buf, "G3,G1,G2,G4");
    }
}
