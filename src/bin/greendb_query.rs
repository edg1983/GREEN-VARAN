use anyhow::{Context, Result};
use clap::Parser;
use greenvaran::cli::GenomeBuild;
use greenvaran::items::get_items;
use greenvaran::logging::{format_elapsed, DualLogger};
use rusqlite::{params_from_iter, types::ValueRef, Connection, OpenFlags};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

const GREENDB_VERSIONS: &[&str] = &["v2.0", "v2.5"];

const TABLES: &[(&str, &[&str])] = &[
    (
        "regions",
        &[
            "regionID",
            "chrom",
            "start",
            "stop",
            "type",
            "std_type",
            "DB_source",
            "detection_method",
            "PhyloP100_median",
            "constraint_pct",
            "controlled_genes",
            "closestGene_symbol",
            "closestGene_dist",
            "closestProtGene_symbol",
            "closestProtGene_dist",
            "cell_or_tissues",
            "phenotype",
        ],
    ),
    (
        "gene_details",
        &[
            "regionID",
            "chrom",
            "start",
            "stop",
            "std_type",
            "controlled_gene",
            "same_TAD",
            "detection_method",
            "tissue_of_interaction",
        ],
    ),
    (
        "pheno_details",
        &[
            "regionID",
            "chrom",
            "start",
            "stop",
            "std_type",
            "phenotype",
            "detection_method",
            "DB_source",
        ],
    ),
    (
        "DNase",
        &[
            "regionID",
            "DNase_chrom",
            "DNase_start",
            "DNase_stop",
            "DNase_ID",
            "DNase_cell_or_tissue",
        ],
    ),
    (
        "dbSuper",
        &[
            "regionID",
            "dbSuper_chrom",
            "dbSuper_start",
            "dbSuper_stop",
            "dbSuper_ID",
            "dbSuper_cell_or_tissue",
        ],
    ),
    (
        "TFBS",
        &[
            "regionID",
            "TFBS_chrom",
            "TFBS_start",
            "TFBS_stop",
            "TF_name",
            "TFBS_cell_or_tissue",
        ],
    ),
    (
        "UCNE",
        &[
            "regionID",
            "UCNE_chrom",
            "UCNE_start",
            "UCNE_stop",
            "UCNE_ID",
        ],
    ),
];

#[derive(Debug, Parser)]
#[command(name = "greendb_query", about = "Query GREEN-DB SQLite database")]
struct Args {
    /// GREEN-DB SQlite file .db
    #[arg(short = 'd', long = "db")]
    db: Option<PathBuf>,

    /// output prefix. A file prefix for --regions or a directory for --tab
    #[arg(short = 'o', long = "out")]
    out: Option<String>,

    /// input tab-separated file [col1:varid, col2:regiond_IDs]
    #[arg(short = 't', long = "tab")]
    tab: Option<PathBuf>,

    /// list of region IDs. either comma-separated list of txt file
    #[arg(short = 'r', long = "regids")]
    regids: Option<String>,

    /// list of variants (chr_pos_ref_alt). either comma-separated list of txt file
    #[arg(short = 'v', long = "variants")]
    variants: Option<String>,

    /// genome build for query
    #[arg(short = 'g', long = "genome", value_enum, default_value = "GRCh38")]
    genome: GenomeBuild,

    /// Log file. Default is greenvaran_[now].log
    #[arg(long = "log")]
    log: Option<String>,
}

#[derive(Debug, Clone)]
struct QueryInterval {
    chrom: String,
    start: i64,
    stop: i64,
}

fn parse_variant_id(v: &str) -> Option<(String, i64, i64)> {
    let fields: Vec<&str> = v.split('_').collect();
    if fields.len() == 4 {
        let start = fields[1].parse::<i64>().ok()?;
        Some((fields[0].to_string(), start, start))
    } else if fields.len() == 5 {
        let start = fields[1].parse::<i64>().ok()?;
        let stop = fields[2].parse::<i64>().ok()?;
        Some((fields[0].to_string(), start, stop))
    } else {
        None
    }
}

fn make_interval_from_row(row: &[String], idx: [usize; 3]) -> Option<QueryInterval> {
    let chrom = row.get(idx[0])?.replace("chr", "");
    let start = row.get(idx[1])?.parse::<i64>().ok()?;
    let stop = row.get(idx[2])?.parse::<i64>().ok()?;
    if start > stop {
        return None;
    }
    Some(QueryInterval { chrom, start, stop })
}

/// Match Nim's `vinterval.overlap(rinterval)` overlap_n > 0 check cleanly without NaN division.
#[inline]
fn intervals_overlap_n(x: &QueryInterval, y: &QueryInterval) -> i64 {
    if x.chrom != y.chrom {
        return 0;
    }
    if x.start <= y.start {
        x.stop - y.start
    } else {
        y.stop - x.start
    }
}

fn build_query_sql(table: &str, build: &str, n_ids: usize) -> String {
    let placeholders = vec!["?"; n_ids].join(",");
    match table {
        "DNase" | "TFBS" => format!(
            "SELECT r.regionID AS regionID, \
            t.chromosome AS {table}_chrom, \
            t.start AS {table}_start, \
            t.stop AS {table}_stop, \
            t.name AS {table}_name, \
            tt.cell_or_tissue AS {table}_tissue \
            FROM {build}_Regions AS r \
            INNER JOIN {build}_regionID_to_{table} AS link ON r.regionID = link.regionID \
            INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID \
            INNER JOIN {build}_{table}_tissues AS tt ON t.regionID = tt.regionID \
            WHERE r.regionID IN ({placeholders})"
        ),
        "dbSuper" => format!(
            "SELECT r.regionID AS regionID, \
            t.chromosome AS {table}_chrom, \
            t.start AS {table}_start, \
            t.stop AS {table}_stop, \
            t.regionID AS {table}_name, \
            tt.cell_or_tissue AS {table}_tissue \
            FROM {build}_Regions AS r \
            INNER JOIN regionID_to_{table} AS link ON r.regionID = link.regionID \
            INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID \
            INNER JOIN {table}_tissues AS tt ON t.regionID = tt.regionID \
            WHERE r.regionID IN ({placeholders})"
        ),
        "UCNE" => format!(
            "SELECT r.regionID AS regionID, \
            t.chromosome AS {table}_chrom, \
            t.start AS {table}_start, \
            t.stop AS {table}_stop, \
            t.regionID AS {table}_name \
            FROM {build}_Regions AS r \
            INNER JOIN regionID_to_{table} AS link ON r.regionID = link.regionID \
            INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID \
            WHERE r.regionID IN ({placeholders})"
        ),
        "regions" => format!(
            "SELECT r.regionID, r.chromosome, r.start, r.stop, \
            r.type, r.std_type, r.DB_source, GROUP_CONCAT(DISTINCT m.method), \
            r.PhyloP100_median, r.constrain_pct, \
            GROUP_CONCAT(DISTINCT g.gene_symbol), \
            r.closestGene_symbol, r.closestGene_dist, \
            r.closestProt_symbol, r.closestProt_dist, \
            GROUP_CONCAT(DISTINCT t.cell_or_tissue), \
            GROUP_CONCAT(DISTINCT p.phenotype) \
            FROM {build}_Regions AS r \
            LEFT JOIN tissues AS t ON r.regionID = t.regionID \
            LEFT JOIN genes AS g ON r.regionID = g.regionID \
            LEFT JOIN methods AS m ON r.regionID = m.regionID \
            LEFT JOIN phenos AS p ON r.regionID = p.regionID \
            WHERE r.regionID IN ({placeholders}) \
            GROUP BY r.regionID"
        ),
        "gene_details" => format!(
            "SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, \
            g.gene_symbol, g.sameTAD, m.method, t.cell_or_tissue \
            FROM {build}_Regions AS r \
            LEFT JOIN genes AS g ON r.regionID = g.regionID \
            LEFT JOIN tissues AS t ON g.interactionID = t.regionID \
            LEFT JOIN methods AS m ON g.interactionID = m.regionID \
            WHERE r.regionID IN ({placeholders})"
        ),
        "pheno_details" => format!(
            "SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, \
            p.phenotype, p.method, p.DB_source \
            FROM {build}_Regions AS r \
            LEFT JOIN phenos AS p ON r.regionID = p.regionID \
            WHERE r.regionID IN ({placeholders})"
        ),
        _ => unreachable!(),
    }
}

fn value_to_string(v: ValueRef<'_>) -> String {
    match v {
        ValueRef::Null => String::new(),
        ValueRef::Integer(i) => i.to_string(),
        ValueRef::Real(f) => {
            if f.is_finite() && f.fract() == 0.0 && f.abs() < 1e15 {
                format!("{f:.1}")
            } else {
                format!("{f}")
            }
        }
        ValueRef::Text(t) => String::from_utf8_lossy(t).into_owned(),
        ValueRef::Blob(b) => String::from_utf8_lossy(b).into_owned(),
    }
}

fn make_tables(
    conn: &Connection,
    regions: &[String],
    build: &str,
    out_prefix: &str,
    variant: Option<&str>,
    logger: &DualLogger,
) -> Result<bool> {
    let mut run_ok = true;
    let parsed_var: Option<Option<QueryInterval>> = variant.map(|v| {
        parse_variant_id(v).and_then(|(c, st, en)| {
            if st > en {
                None
            } else {
                Some(QueryInterval {
                    chrom: c.replace("chr", ""),
                    start: st,
                    stop: en,
                })
            }
        })
    });

    for &(table, header) in TABLES {
        let out_path = format!("{out_prefix}.{table}.tsv");
        let file = File::create(&out_path)
            .with_context(|| format!("Failed to create output table: {out_path}"))?;
        let mut w = BufWriter::with_capacity(64 * 1024, file);

        if variant.is_some() {
            writeln!(w, "{}\tVariantID", header.join("\t"))?;
        } else {
            writeln!(w, "{}", header.join("\t"))?;
        }

        if regions.is_empty() {
            w.flush()?;
            continue;
        }

        // Chunk IN (...) parameters safely under SQLite variable limit (32766)
        for chunk in regions.chunks(10_000) {
            let sql = build_query_sql(table, build, chunk.len());
            let mut stmt = conn.prepare(&sql)?;
            let col_count = stmt.column_count();
            let mut rows = stmt.query(params_from_iter(chunk.iter()))?;

            while let Some(row) = rows.next()? {
                let mut cols: Vec<String> = Vec::with_capacity(col_count);
                for c in 0..col_count {
                    cols.push(value_to_string(row.get_ref(c)?));
                }

                if let Some(var_str) = variant {
                    let varout = match parsed_var.as_ref().unwrap() {
                        Some(v_iv) => match make_interval_from_row(&cols, [1, 2, 3]) {
                            Some(r_iv) => {
                                if intervals_overlap_n(v_iv, &r_iv) > 0 {
                                    var_str
                                } else {
                                    "NOVARIANT"
                                }
                            }
                            None => {
                                logger.debug(&format!(
                                    "Error occurred parsing interval for region {} and variant {var_str}",
                                    cols.first().map(|s| s.as_str()).unwrap_or("")
                                ));
                                run_ok = false;
                                "NA"
                            }
                        },
                        None => {
                            logger.debug(&format!("Failed to parse variant {var_str}"));
                            run_ok = false;
                            "NA"
                        }
                    };
                    writeln!(w, "{}\t{varout}", cols.join("\t"))?;
                } else {
                    writeln!(w, "{}", cols.join("\t"))?;
                }
            }
        }
        w.flush()?;
    }

    Ok(run_ok)
}

fn run() -> Result<()> {
    let args = Args::parse();
    let out_str = args.out.as_deref().unwrap_or("");
    let logger = DualLogger::new(args.log.as_deref(), out_str)?;

    let mut start_time = Instant::now();
    logger.info("=== GREEN-DB QUERY ===");
    logger.info(&format!(
        "Compatible with GREEN-DB ver: {GREENDB_VERSIONS:?}"
    ));
    logger.info(&format!("Genome build: {}", args.genome.as_str()));
    if let Some(ref db) = args.db {
        logger.info(&format!("Open database from: {}", db.display()));
    }

    let (db_path, out_prefix) = match (&args.db, &args.out) {
        (Some(d), Some(o)) if !o.is_empty() => (d, o.as_str()),
        _ => {
            logger.fatal("One of --db, --genome or --out is missing");
            std::process::exit(1);
        }
    };

    let mode_count = args.regids.is_some() as u8
        + args.tab.is_some() as u8
        + args.variants.is_some() as u8;
    if mode_count == 0 {
        logger.fatal("One of --regids, --tab or --variants must be provided");
        std::process::exit(1);
    }
    if mode_count > 1 {
        logger.fatal("Only one of --regids, --tab or --variants can be used");
        std::process::exit(1);
    }

    let conn = Connection::open_with_flags(
        db_path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .with_context(|| format!("Failed to open SQLite database: {}", db_path.display()))?;
    let _ = conn.execute_batch("PRAGMA cache_size = -65536; PRAGMA temp_store = MEMORY;");

    let build = args.genome.as_str();
    let mut warning = false;

    if let Some(ref regids_arg) = args.regids {
        let regions = get_items(Some(regids_arg), "regions", false);
        logger.info(&format!("Regions in input: {}", regions.len()));
        make_tables(&conn, &regions, build, out_prefix, None, &logger)?;
    } else if let Some(ref tab_path) = args.tab {
        logger.info(&format!("Reading from: {}", tab_path.display()));
        let file = File::open(tab_path)
            .with_context(|| format!("Failed to open tab file: {}", tab_path.display()))?;
        let reader = BufReader::new(file);
        let mut n = 0usize;
        let mut nregions = 0usize;
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            n += 1;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                continue;
            }
            let var_id = fields[0];
            let dir = format!("{out_prefix}/{var_id}");
            fs::create_dir_all(&dir)?;
            let sub_prefix = format!("{dir}/{var_id}");
            let regions: Vec<String> = fields[1]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .collect();
            nregions += regions.len();
            if !make_tables(&conn, &regions, build, &sub_prefix, Some(var_id), &logger)? {
                warning = true;
            }
        }
        logger.info(&format!(
            "{n} variants and {nregions} regions processed from input"
        ));
    } else if let Some(ref var_arg) = args.variants {
        let variants = get_items(Some(var_arg), "variants", false);
        logger.info(&format!("Variants in input: {}", variants.len()));
        let mut nregions = 0usize;
        let query_sql = format!(
            "SELECT regionID FROM {build}_Regions WHERE chromosome = ? AND start <= ? AND stop >= ?"
        );
        let mut stmt = conn.prepare(&query_sql)?;

        for v in &variants {
            let dir = format!("{out_prefix}/{v}");
            fs::create_dir_all(&dir)?;
            let sub_prefix = format!("{dir}/{v}");

            let mut regions = Vec::new();
            match parse_variant_id(v) {
                Some((chrom, start, stop)) => {
                    let mut rows = stmt.query(rusqlite::params![chrom, start, stop])?;
                    while let Some(row) = rows.next()? {
                        regions.push(row.get::<_, String>(0)?);
                    }
                }
                None => {
                    logger.debug(&format!("Failed to parse variant {v}"));
                    warning = true;
                }
            }
            nregions += regions.len();
            logger.debug(&format!("{} regions extracted for variant {v}", regions.len()));
            if !make_tables(&conn, &regions, build, &sub_prefix, Some(v), &logger)? {
                warning = true;
            }
        }
        logger.info(&format!(
            "Total of {nregions} regions extracted for input variants"
        ));
    }

    if warning {
        logger.error("Some regions of variants failed to parse, see log for more information");
    }
    logger.info(&format!(
        "All done - Completed in {}",
        format_elapsed(&mut start_time)
    ));
    Ok(())
}

fn main() {
    if let Err(err) = run() {
        eprintln!("FATAL: {err:#}");
        std::process::exit(1);
    }
}
