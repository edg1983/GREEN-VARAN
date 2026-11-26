use anyhow::{Context, Result};
use clap::Args;
use rust_htslib::bcf::{Header, IndexedReader, Read as BcfRead};
use rustc_hash::FxHashSet;
use std::path::PathBuf;
use std::time::Instant;

use crate::annotator::{annotate_record_in_place, update_header, VariantScratch};
use crate::bed_index::ChromIndex;
use crate::cli::{resolve_threads, ConnectionMode, Impact};
use crate::config::DbSchema;
use crate::csq::{parse_csq_schema, CsqSchema};
use crate::intervals::make_sv_interval;
use crate::items::get_items;
use crate::logging::{format_elapsed, DualLogger};
use crate::vcf_io::{open_vcf_writer, split_io_threads};

#[derive(Debug, Args)]
#[command(about = "Annotate structural variants VCF")]
pub struct SvArgs {
    /// path to input VCF/BCF
    #[arg(short = 'i', long = "invcf")]
    pub invcf: Option<PathBuf>,

    /// output VCF file
    #[arg(short = 'o', long = "outvcf")]
    pub outvcf: Option<PathBuf>,

    /// GREEN-DB bed.gz file
    #[arg(short = 'd', long = "db")]
    pub db: Option<PathBuf>,

    /// json file containig greendb column mapping
    #[arg(short = 's', long = "dbschema")]
    pub dbschema: Option<PathBuf>,

    /// Which impact to assign when updating snpEff field
    #[arg(short = 'm', long = "impact", value_enum, default_value = "MODIFIER")]
    pub impact: Impact,

    /// Consequence field containing gene consequences
    #[arg(short = 'q', long = "csq_field")]
    pub csq_field: Option<String>,

    /// do not update consequence field
    #[arg(short = 'u', long = "noupdate")]
    pub noupdate: bool,

    /// Filter instead of annotate. Only vars with greendb overlap and eventually gene of interest will be written
    #[arg(short = 'f', long = "filter")]
    pub filter: bool,

    /// Chromosome names in input file do not have chr prefix
    #[arg(short = 'n', long = "nochr")]
    pub nochr: bool,

    /// Annotate only for specific chromosome
    #[arg(long = "chrom")]
    pub chrom: Option<String>,

    /// Genes of interest, variants connected to those will be flagged with greendb_VOI
    #[arg(short = 'g', long = "genes")]
    pub genes: Option<String>,

    /// Region-gene connections accepted for annotation
    #[arg(long = "connection", value_enum, default_value = "all")]
    pub connection: ConnectionMode,

    /// Value to add on each side of BND/INS, this override the CIPOS when set
    #[arg(short = 'p', long = "padding")]
    pub padding: Option<i64>,

    /// INFO field listing the confidence interval around POS, it is expected to have 2 comma-separated values
    #[arg(long = "cipos", default_value = "CIPOS")]
    pub cipos: String,

    /// INFO field listing the confidence interval around END, it is expected to have 2 comma-separated values
    #[arg(long = "ciend", default_value = "CIEND")]
    pub ciend: String,

    /// Max length of SV to be annotated
    #[arg(long = "maxlen", default_value_t = 2_000_000)]
    pub maxlen: i64,

    /// Min fraction of GREENDB region to be overlapped by a SV
    #[arg(short = 't', long = "minoverlap", default_value_t = 0.000001)]
    pub minoverlap: f64,

    /// Min number of bases of GREENDB region to be overlapped by a SV
    #[arg(short = 'b', long = "minbp", default_value_t = 1)]
    pub minbp: i64,

    /// Log file. Default is greenvaran_[now].log
    #[arg(long = "log")]
    pub log: Option<String>,

    /// Number of threads for BGZF decompression/compression (0 = auto)
    #[arg(long = "threads", default_value_t = 1)]
    pub threads: usize,
}

pub fn run_sv(args: SvArgs) -> Result<()> {
    let out_str = args
        .outvcf
        .as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_default();
    let logger = DualLogger::new(args.log.as_deref(), &out_str)?;

    let (invcf, db_path, dbschema_path) = match (&args.invcf, &args.db, &args.dbschema) {
        (Some(i), Some(d), Some(s)) => (i, d, s),
        _ => {
            logger.fatal("One of --vcf, --db or --dbschema is missing");
            std::process::exit(1);
        }
    };

    let threads = resolve_threads(args.threads);
    let (read_threads, write_threads) = split_io_threads(threads);
    let mut starttime = Instant::now();

    let dbschema = DbSchema::from_file(dbschema_path)?;

    let chromosomes = get_items(args.chrom.as_deref(), "chromosome", args.nochr);
    let genes_vec = get_items(args.genes.as_deref(), "gene", args.nochr);
    let genes: FxHashSet<String> = genes_vec.into_iter().collect();
    let has_voi_genes = !genes.is_empty();

    logger.info(&format!("N selected chromosomes: {}", chromosomes.len()));
    logger.info(&format!("N selected genes: {}", genes.len()));
    logger.debug(&format!("Selected chromosomes: {chromosomes:?}"));
    logger.debug(&format!("Selected genes: {genes:?}"));

    logger.info(&format!(
        "Min overlap - fraction: {}; bp: {}",
        args.minoverlap, args.minbp
    ));
    logger.info(&format!("Max length for SVs: {}", args.maxlen));

    if let Some(pad) = args.padding {
        logger.info(&format!("Fixed padding for INS/BND: {pad}"));
    } else {
        logger.info(&format!(
            "Padding for INS/BND based on INFO field: {}",
            args.cipos
        ));
    }

    let update_ann = !args.noupdate;
    logger.info(&format!("Update existing gene annotations: {update_ann}"));
    if let Some(ref o) = args.outvcf {
        logger.info(&format!("Output to: {}", o.display()));
    }
    logger.info(&format!("Filter mode active: {}", args.filter));
    if threads > 1 {
        logger.info(&format!("Threads: {threads}"));
    }

    logger.info("=== Start processing VCF ===");

    let mut reader = IndexedReader::from_path(invcf)
        .with_context(|| format!("Failed to open indexed VCF/BCF: {}", invcf.display()))?;
    if read_threads > 1 {
        let _ = reader.set_threads(read_threads);
    }

    if reader.header().info_type(b"END").is_err() {
        logger.warn(
            "No END field defined in VCF header. SV intervals will be built from SVLEN and padding only",
        );
    }
    if args.padding.is_none() {
        logger.warn(
            "No padding value set for BND/INS. Intervals will be built based on CIPOS only.",
        );
    }
    if reader.header().info_type(args.cipos.as_bytes()).is_err() {
        logger.warn("The configured CIPOS field is not in header. Zero will be used as value");
    }
    if reader.header().info_type(args.ciend.as_bytes()).is_err() {
        logger.warn("The configured CIEND field is not in header. Zero will be used as value");
    }

    let csq_schema: Option<CsqSchema> = if update_ann {
        match parse_csq_schema(
            reader.header(),
            &invcf.to_string_lossy(),
            args.csq_field.as_deref(),
            &logger,
        ) {
            Ok(s) => Some(s),
            Err(e) => {
                logger.fatal(&e.to_string());
                std::process::exit(1);
            }
        }
    } else {
        None
    };

    let ann_tag_bytes: Vec<u8> = csq_schema
        .as_ref()
        .map(|s| s.csq_field_name.as_bytes().to_vec())
        .unwrap_or_default();

    let mut out_header = Header::from_template(reader.header());
    update_header(&mut out_header, csq_schema.as_ref(), true, has_voi_genes);
    let mut writer = open_vcf_writer(args.outvcf.as_deref(), &out_header, write_threads)?;

    let cipos_bytes = args.cipos.as_bytes();
    let ciend_bytes = args.ciend.as_bytes();

    let mut total_n: usize = 0;
    let mut total_written: usize = 0;
    let mut total_annotated: usize = 0;
    let mut total_voi: usize = 0;
    let mut warning = false;
    let mut t0 = Instant::now();

    let mut rec = reader.empty_record();
    let mut scratch = VariantScratch::default();
    let mut hits_buf: Vec<u32> = Vec::with_capacity(64);

    for chrom in &chromosomes {
        logger.debug(&format!("Processing chromosome {chrom}"));
        let rid = match reader.header().name2rid(chrom.as_bytes()) {
            Ok(r) => r,
            Err(_) => {
                logger.debug(&format!(
                    "chromosome {chrom}: 0 processed; 0 annotated, 0 vars of interest, 0 written"
                ));
                continue;
            }
        };

        if reader.fetch(rid, 0, None).is_err() {
            continue;
        }

        match reader.read(&mut rec) {
            None => {
                logger.debug(&format!(
                    "chromosome {chrom}: 0 processed; 0 annotated, 0 vars of interest, 0 written"
                ));
                continue;
            }
            Some(Err(e)) => return Err(e).context("Failed to read VCF/BCF record"),
            Some(Ok(())) => {}
        }

        let db_chrom = if args.nochr {
            format!("chr{chrom}")
        } else {
            chrom.clone()
        };
        let chrom_index = ChromIndex::load_chromosome(
            db_path,
            &db_chrom,
            dbschema,
            args.connection,
            &genes,
            threads,
        )?;

        let mut nchrom: usize = 0;
        let mut nannotchrom: usize = 0;
        let mut nvoichrom: usize = 0;
        let mut nwrittenchrom: usize = 0;

        loop {
            nchrom += 1;
            total_n += 1;
            if total_n % 5_000 == 0 {
                logger.info(&format!(
                    "{total_n} vars analyzed, last batch in {}",
                    format_elapsed(&mut t0)
                ));
            }

            let vinterval = match make_sv_interval(&rec, cipos_bytes, ciend_bytes, args.padding) {
                Some(iv) => iv,
                None => {
                    warning = true;
                    logger.debug(&format!(
                        "Unable to parse an interval for variant at {chrom}:{}",
                        rec.pos() + 1
                    ));
                    writer.translate(&mut rec);
                    writer.write(&rec)?;
                    nwrittenchrom += 1;
                    total_written += 1;
                    match reader.read(&mut rec) {
                        None => break,
                        Some(Err(e)) => return Err(e).context("Failed to read VCF/BCF record"),
                        Some(Ok(())) => continue,
                    }
                }
            };

            if vinterval.len() > args.maxlen {
                warning = true;
                logger.debug(&format!(
                    "Skipping variant at {chrom}:{} with interval length > {}",
                    rec.pos() + 1,
                    args.maxlen
                ));
                writer.translate(&mut rec);
                writer.write(&rec)?;
                nwrittenchrom += 1;
                total_written += 1;
                match reader.read(&mut rec) {
                    None => break,
                    Some(Err(e)) => return Err(e).context("Failed to read VCF/BCF record"),
                    Some(Ok(())) => continue,
                }
            }

            chrom_index.query_sv(vinterval, args.minoverlap, args.minbp, &mut hits_buf);

            if hits_buf.is_empty() {
                if !args.filter {
                    writer.translate(&mut rec);
                    writer.write(&rec)?;
                    nwrittenchrom += 1;
                    total_written += 1;
                }
            } else {
                writer.translate(&mut rec);
                let (annotated, is_voi, keep) = annotate_record_in_place(
                    &mut rec,
                    &hits_buf,
                    &chrom_index,
                    args.filter,
                    has_voi_genes,
                    None,
                    true,
                    csq_schema.as_ref(),
                    &ann_tag_bytes,
                    args.impact,
                    &mut scratch,
                )?;
                if annotated {
                    nannotchrom += 1;
                    total_annotated += 1;
                }
                if is_voi {
                    nvoichrom += 1;
                    total_voi += 1;
                }
                if keep {
                    writer.write(&rec)?;
                    nwrittenchrom += 1;
                    total_written += 1;
                }
            }

            match reader.read(&mut rec) {
                None => break,
                Some(Err(e)) => return Err(e).context("Failed to read VCF/BCF record"),
                Some(Ok(())) => {}
            }
        }

        logger.debug(&format!(
            "chromosome {chrom}: {nchrom} processed; {nannotchrom} annotated, {nvoichrom} vars of interest, {nwrittenchrom} written"
        ));
    }

    drop(writer);
    drop(reader);

    logger.info(&format!(
        "{total_annotated} variants annotated with greendb information"
    ));
    logger.info(&format!(
        "{total_voi} vars of interest based on the input gene list if any"
    ));
    logger.info(&format!("{total_written} variants written to output"));
    if warning {
        logger.error(
            "Some variants were not annotated due to failed interval parsing of lenght limit. See log",
        );
    }
    logger.info(&format!(
        "All done - Completed in {}",
        format_elapsed(&mut starttime)
    ));
    Ok(())
}
