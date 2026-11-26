use anyhow::{Context, Result};
use clap::Args;
use rust_htslib::bcf::{Header, IndexedReader, Read as BcfRead};
use rustc_hash::FxHashSet;
use std::path::PathBuf;
use std::time::Instant;

use crate::annotator::{
    annotate_record_in_place, update_header, CompiledPrioritize, VariantScratch,
};
use crate::bed_index::{ChromIndex, SmallVarCursor};
use crate::cli::{resolve_threads, ConnectionMode, Impact, PrioritizationStrategy};
use crate::config::{DbSchema, PrioritizeConfig};
use crate::csq::{parse_csq_schema, CsqSchema};
use crate::items::get_items;
use crate::logging::{format_elapsed, DualLogger};
use crate::vcf_io::{open_vcf_writer, split_io_threads};

#[derive(Debug, Args)]
#[command(about = "Annotate small variants VCF")]
pub struct SmallvarsArgs {
    /// path to input VCF/BCF
    #[arg(short = 'i', long = "invcf")]
    pub invcf: Option<PathBuf>,

    /// output VCF file
    #[arg(short = 'o', long = "outvcf")]
    pub outvcf: Option<PathBuf>,

    /// json config file for prioritization
    #[arg(short = 'c', long = "config")]
    pub config: Option<PathBuf>,

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

    /// Perform prioritization even if one of INFO fields required by prioritization config is missing
    #[arg(short = 'p', long = "permissive")]
    pub permissive: bool,

    /// Chromosome names in input file do not have chr prefix
    #[arg(short = 'n', long = "nochr")]
    pub nochr: bool,

    /// Annotate only for specific chromosome
    #[arg(long = "chrom")]
    pub chrom: Option<String>,

    /// Genes of interest, variants connected to those will be flagged with greendb_VOI
    #[arg(short = 'g', long = "genes")]
    pub genes: Option<String>,

    /// Set the behaviour to compute prioritization levels
    #[arg(long = "prioritization_strategy", value_enum, default_value = "levels")]
    pub prioritization_strategy: PrioritizationStrategy,

    /// Region-gene connections accepted for annotation
    #[arg(long = "connection", value_enum, default_value = "all")]
    pub connection: ConnectionMode,

    /// Log file. Default is greenvaran_[now].log
    #[arg(long = "log")]
    pub log: Option<String>,

    /// Number of threads for BGZF decompression/compression (0 = auto)
    #[arg(long = "threads", default_value_t = 1)]
    pub threads: usize,
}

pub fn run_smallvars(args: SmallvarsArgs) -> Result<()> {
    let out_str = args
        .outvcf
        .as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_default();
    let logger = DualLogger::new(args.log.as_deref(), &out_str)?;

    let (invcf, outvcf, config_path, db_path, dbschema_path) = match (
        &args.invcf,
        &args.outvcf,
        &args.config,
        &args.db,
        &args.dbschema,
    ) {
        (Some(i), Some(o), Some(c), Some(d), Some(s)) => (i, o, c, d, s),
        _ => {
            logger.fatal("One of --vcf, --config, --db or --dbschema is missing");
            std::process::exit(1);
        }
    };

    let threads = resolve_threads(args.threads);
    let (read_threads, write_threads) = split_io_threads(threads);
    let mut starttime = Instant::now();

    logger.info(&format!(
        "Reading config from file: {}",
        config_path.display()
    ));
    let config = PrioritizeConfig::from_file(config_path)?;
    logger.debug(&format!("config: {config:?}"));

    let dbschema = DbSchema::from_file(dbschema_path)?;

    let chromosomes = get_items(args.chrom.as_deref(), "chromosome", args.nochr);
    let genes_vec = get_items(args.genes.as_deref(), "gene", args.nochr);
    let genes: FxHashSet<String> = genes_vec.into_iter().collect();
    let has_voi_genes = !genes.is_empty();

    logger.info(&format!("N selected chromosomes: {}", chromosomes.len()));
    logger.info(&format!("N selected genes: {}", genes.len()));
    logger.debug(&format!("Selected chromosomes: {chromosomes:?}"));
    logger.debug(&format!("Selected genes: {genes:?}"));

    let update_ann = !args.noupdate;
    logger.info(&format!("Update existing gene annotations: {update_ann}"));
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

    let compiled_prio = CompiledPrioritize::from_config(&config, args.prioritization_strategy);
    let mut prioritize = true;
    for missing in compiled_prio.check_header(reader.header()) {
        logger.warn(&format!(
            "{missing} field defined in config not present in the VCF header"
        ));
        if !args.permissive {
            prioritize = false;
        }
    }
    if !prioritize {
        logger.warn("Prioritize is not active and all variants will get level zero");
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
    update_header(&mut out_header, csq_schema.as_ref(), false, has_voi_genes);
    let mut writer = open_vcf_writer(Some(outvcf), &out_header, write_threads)?;

    let prio_opt = if prioritize {
        Some(&compiled_prio)
    } else {
        None
    };

    let mut total_n: usize = 0;
    let mut total_written: usize = 0;
    let mut total_annotated: usize = 0;
    let mut total_voi: usize = 0;
    let mut t0 = Instant::now();

    let mut rec = reader.empty_record();
    let mut scratch = VariantScratch::default();
    let mut cursor = SmallVarCursor::new();

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

        // Peek first record: if 0 variants on this chrom, skip loading GREEN-DB BED!
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

        cursor.reset();
        loop {
            nchrom += 1;
            total_n += 1;
            if total_n % 10_000 == 0 {
                logger.info(&format!(
                    "{total_n} vars analyzed, last batch in {}",
                    format_elapsed(&mut t0)
                ));
            }

            let pos_1based = rec.pos() + 1;
            let hits = cursor.query_point(&chrom_index, pos_1based);

            if hits.is_empty() {
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
                    hits,
                    &chrom_index,
                    args.filter,
                    has_voi_genes,
                    prio_opt,
                    false,
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
    logger.info(&format!(
        "All done - Completed in {}",
        format_elapsed(&mut starttime)
    ));
    Ok(())
}
