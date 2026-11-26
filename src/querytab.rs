use anyhow::{bail, Context, Result};
use clap::Args;
use rust_htslib::bcf::{
    header::HeaderView, record::Numeric, IndexedReader, Read as BcfRead, Reader, Record,
};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use crate::cli::resolve_threads;
use crate::items::get_items;
use crate::logging::{format_elapsed, DualLogger};

#[derive(Debug, Args)]
#[command(about = "Generate a table for greendb_query from VCF")]
pub struct QuerytabArgs {
    /// path to input VCF/BCF
    #[arg(short = 'i', long = "vcf")]
    pub vcf: Option<PathBuf>,

    /// output tsv file
    #[arg(short = 'o', long = "out")]
    pub out: Option<PathBuf>,

    /// min accepted greenvaran level
    #[arg(short = 'l', long = "minlevel")]
    pub minlevel: Option<i32>,

    /// Extract only variants of interest marked as greendb_VOI
    #[arg(short = 'x', long = "voi")]
    pub voi: bool,

    /// bed file or comma-separated list of regions as chr:start-end
    #[arg(short = 'r', long = "regions")]
    pub regions: Option<String>,

    /// Number of threads for BGZF decompression (0 = auto)
    #[arg(long = "threads", default_value_t = 1)]
    pub threads: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

/// Parse a region string either as `chr:start-end` or BED tab-separated `chr\tstart\tend...`.
pub fn parse_region(raw: &str) -> Result<ParsedRegion> {
    let s = raw.trim();
    if let Some((chrom, coords)) = s.split_once(':') {
        if let Some((st, en)) = coords.split_once('-') {
            if let (Ok(start), Ok(end)) = (st.parse::<u64>(), en.parse::<u64>()) {
                return Ok(ParsedRegion {
                    chrom: chrom.to_string(),
                    start,
                    end,
                });
            }
        }
    }

    let fields: Vec<&str> = s.split('\t').collect();
    if fields.len() >= 3 {
        if let (Ok(start), Ok(end)) = (fields[1].parse::<u64>(), fields[2].parse::<u64>()) {
            return Ok(ParsedRegion {
                chrom: fields[0].to_string(),
                start,
                end,
            });
        }
    }

    bail!("Failed to parse region {raw}")
}

fn check_headers(header: &HeaderView, logger: &DualLogger) -> Result<bool> {
    let use_end = header.info_type(b"END").is_ok();
    logger.info(&format!(
        "End position for variant from {}",
        if use_end { "END" } else { "POS" }
    ));

    for tag in ["greendb_id", "greendb_level"] {
        if header.info_type(tag.as_bytes()).is_err() {
            logger.fatal(&format!(
                "{tag} annotation missing in header. Is your VCF annotated with GREEN-VARAN?"
            ));
            std::process::exit(1);
        }
    }
    Ok(use_end)
}

fn extract_variant<W: Write>(
    rec: &Record,
    header: &HeaderView,
    use_end: bool,
    minlevel: i32,
    voi: bool,
    out: &mut W,
) -> Result<usize> {
    let regions_vec = match rec.info(b"greendb_id").string()? {
        Some(v) if !v.is_empty() => v,
        _ => return Ok(0),
    };

    if minlevel > 0 {
        if let Ok(Some(lvl)) = rec.info(b"greendb_level").integer() {
            if let Some(&v) = lvl.first() {
                if !v.is_missing() && v < minlevel {
                    return Ok(0);
                }
            }
        }
    }

    if voi && !rec.info(b"greendb_VOI").flag().unwrap_or(false) {
        return Ok(0);
    }

    let pos_1based = rec.pos() + 1;
    let vend: i64 = if use_end {
        match rec.info(b"END").integer() {
            Ok(Some(vals)) if !vals.is_empty() && !vals[0].is_missing() => vals[0] as i64,
            _ => pos_1based,
        }
    } else {
        pos_1based
    };

    let chrom = match rec.rid() {
        Some(rid) => std::str::from_utf8(header.rid2name(rid)?).unwrap_or("."),
        None => ".",
    };
    let alleles = rec.alleles();
    let ref_allele = alleles
        .first()
        .and_then(|a| std::str::from_utf8(a).ok())
        .unwrap_or(".");
    let alt0 = alleles
        .get(1)
        .and_then(|a| std::str::from_utf8(a).ok())
        .unwrap_or(".");

    write!(out, "{chrom}_{pos_1based}_{vend}_{ref_allele}_{alt0}\t")?;
    for (i, part) in regions_vec.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        out.write_all(part)?;
    }
    out.write_all(b"\n")?;
    Ok(1)
}

pub fn run_querytab(args: QuerytabArgs) -> Result<()> {
    let logger = DualLogger::console_only();

    let (vcf_path, out_path) = match (&args.vcf, &args.out) {
        (Some(v), Some(o)) => (v, o),
        _ => {
            logger.fatal("One of --vcf or --out is missing");
            std::process::exit(1);
        }
    };

    let mut regions: Vec<ParsedRegion> = Vec::new();
    if let Some(ref reg_arg) = args.regions {
        if !reg_arg.is_empty() {
            let raw_items = get_items(Some(reg_arg), "regions", false);
            for item in raw_items {
                match parse_region(&item) {
                    Ok(r) => regions.push(r),
                    Err(e) => {
                        logger.fatal(&e.to_string());
                        std::process::exit(1);
                    }
                }
            }
            logger.info(&format!("{} regions parsed from input", regions.len()));
        }
    }

    let mut starttime = Instant::now();
    let minlevel = args.minlevel.unwrap_or(0);
    let threads = resolve_threads(args.threads);

    logger.info("=== Start processing VCF ===");

    let out_file = File::create(out_path)
        .with_context(|| format!("Failed to create output file: {}", out_path.display()))?;
    let mut out = BufWriter::with_capacity(256 * 1024, out_file);

    let mut n: usize = 0;
    let mut nwritten: usize = 0;
    let mut t0 = Instant::now();

    if !regions.is_empty() {
        let mut reader = IndexedReader::from_path(vcf_path)
            .with_context(|| format!("Failed to open indexed VCF/BCF: {}", vcf_path.display()))?;
        if threads > 1 {
            let _ = reader.set_threads(threads.min(4));
        }
        let use_end = check_headers(reader.header(), &logger)?;
        let header_view = reader.header().clone();
        let mut rec = reader.empty_record();

        for r in &regions {
            let rid = match header_view.name2rid(r.chrom.as_bytes()) {
                Ok(id) => id,
                Err(_) => continue,
            };
            // htslib fetch takes 0-based start, exclusive end; Nim passes `chr:start-end` (1-based inclusive)
            let start_0 = r.start.saturating_sub(1);
            if reader.fetch(rid, start_0, Some(r.end)).is_err() {
                continue;
            }
            while let Some(res) = reader.read(&mut rec) {
                res.context("Failed to read VCF/BCF record")?;
                n += 1;
                if n % 10_000 == 0 {
                    logger.info(&format!(
                        "{n} vars processed, last batch in {}",
                        format_elapsed(&mut t0)
                    ));
                }
                nwritten += extract_variant(&rec, &header_view, use_end, minlevel, args.voi, &mut out)?;
            }
        }
    } else {
        let mut reader = Reader::from_path(vcf_path)
            .with_context(|| format!("Failed to open VCF/BCF: {}", vcf_path.display()))?;
        if threads > 1 {
            let _ = reader.set_threads(threads.min(4));
        }
        let use_end = check_headers(reader.header(), &logger)?;
        let header_view = reader.header().clone();
        let mut rec = reader.empty_record();

        while let Some(res) = reader.read(&mut rec) {
            res.context("Failed to read VCF/BCF record")?;
            n += 1;
            if n % 10_000 == 0 {
                logger.info(&format!(
                    "{n} vars processed, last batch in {}",
                    format_elapsed(&mut t0)
                ));
            }
            nwritten += extract_variant(&rec, &header_view, use_end, minlevel, args.voi, &mut out)?;
        }
    }

    out.flush()?;
    logger.info(&format!("{nwritten} variants written to output"));
    logger.info(&format!(
        "All done - Completed in {}",
        format_elapsed(&mut starttime)
    ));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region_formats() {
        assert_eq!(
            parse_region("chr1:100-200").unwrap(),
            ParsedRegion {
                chrom: "chr1".into(),
                start: 100,
                end: 200
            }
        );
        assert_eq!(
            parse_region("22\t1000\t2000\textra").unwrap(),
            ParsedRegion {
                chrom: "22".into(),
                start: 1000,
                end: 2000
            }
        );
        assert!(parse_region("invalid_region").is_err());
    }
}
