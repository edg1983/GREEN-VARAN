use anyhow::{Context, Result};
use rust_htslib::bcf::{Format, Header, Writer};
use std::path::Path;

/// Split requested `--threads` budget between htslib BGZF reader and BGZF writer.
/// BGZF compression is ~3x heavier than decompression, so most threads go to the writer.
pub fn split_io_threads(threads: usize) -> (usize, usize) {
    match threads {
        0 | 1 => (1, 1),
        2 => (1, 1),
        3 => (1, 2),
        n => {
            let read_t = (n / 4).max(1);
            let write_t = (n - read_t).max(1);
            (read_t, write_t)
        }
    }
}

/// Open a BCF/VCF/VCF.GZ writer matching hts-nim's extension-based format detection:
/// - `.bcf` -> compressed BCF
/// - `.vcf.gz` / `.gz` / `.bgz` -> BGZF compressed VCF
/// - `-` or empty -> uncompressed VCF to stdout
/// - otherwise -> uncompressed VCF
pub fn open_vcf_writer(
    out_path: Option<&Path>,
    header: &Header,
    write_threads: usize,
) -> Result<Writer> {
    let mut writer = match out_path {
        None => Writer::from_stdout(header, true, Format::Vcf)
            .context("Failed to open stdout VCF writer")?,
        Some(p) if p.as_os_str().is_empty() || p.as_os_str() == "-" => {
            Writer::from_stdout(header, true, Format::Vcf)
                .context("Failed to open stdout VCF writer")?
        }
        Some(p) => {
            let s = p.to_string_lossy();
            let (uncompressed, format) = if s.ends_with(".bcf") {
                (false, Format::Bcf)
            } else if s.ends_with(".gz") || s.ends_with(".bgz") {
                (false, Format::Vcf)
            } else {
                (true, Format::Vcf)
            };
            Writer::from_path(p, header, uncompressed, format)
                .with_context(|| format!("Failed to open output VCF/BCF: {}", p.display()))?
        }
    };

    if write_threads > 1 {
        let _ = writer.set_threads(write_threads);
    }

    Ok(writer)
}
