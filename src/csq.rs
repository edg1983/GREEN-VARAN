use anyhow::{bail, Result};
use rust_htslib::bcf::{header::HeaderView, HeaderRecord};
use std::collections::HashSet;

pub use crate::cli::DESC_ANN;
use crate::logging::DualLogger;

#[derive(Debug, Clone)]
pub struct CsqSchema {
    pub allele: Option<usize>,
    pub gene_symbol: usize,
    pub consequence: usize,
    pub impact: Option<usize>,
    pub csq_len: usize,
    pub csq_field_name: String,
    pub create_ann: bool,
}

fn find_info_description(header: &HeaderView, tag: &str) -> Option<String> {
    if tag.is_empty() {
        return None;
    }
    for rec in header.header_records() {
        if let HeaderRecord::Info { values, .. } = rec {
            if values.get("ID").map(|s| s.as_str()) == Some(tag) {
                if let Some(desc) = values.get("Description") {
                    // Strip surrounding quotes that htslib header_records may keep
                    let trimmed = desc.trim_matches('"');
                    return Some(trimmed.to_string());
                }
            }
        }
    }
    None
}

/// Read VCF header and determine CSQ/ANN/BCSQ column indices matching Nim's `parse_csq_schema`.
pub fn parse_csq_schema(
    header: &HeaderView,
    vcf_name: &str,
    requested_field: Option<&str>,
    logger: &DualLogger,
) -> Result<CsqSchema> {
    let req = requested_field.unwrap_or("");
    let possible_fields = [req, "ANN", "CSQ", "BCSQ"];

    let mut desc = String::new();
    let mut field_name = String::new();
    let mut create_ann = false;

    for &try_field in &possible_fields {
        if let Some(d) = find_info_description(header, try_field) {
            desc = d;
            field_name = try_field.to_string();
            break;
        } else if try_field == req && !req.is_empty() {
            logger.warn(&format!(
                "Didn't find requested {req} field in header in {vcf_name} trying other fields"
            ));
        }
    }

    if desc.is_empty() {
        logger.warn(&format!(
            "None of the possible consequence fields {:?} found. The ANN field will be added.",
            possible_fields
        ));
        field_name = "ANN".to_string();
        create_ann = true;
        desc = DESC_ANN.to_string();
    }

    logger.info(&format!("Writing/updating consequences to {field_name} field"));

    parse_csq_description(&desc, field_name, create_ann)
}

pub fn parse_csq_description(
    desc: &str,
    csq_field_name: String,
    create_ann: bool,
) -> Result<CsqSchema> {
    let spl = if desc.contains("Format: '") {
        "Format: '"
    } else if desc.contains("Format: ") {
        "Format: "
    } else {
        ": '"
    };

    let after_spl = match desc.split_once(spl) {
        Some((_, rhs)) => rhs,
        None => {
            bail!(
                "Error parsing annotation structure frm {csq_field_name}. Check the format is properly defined in the header"
            );
        }
    };

    let raw_fmt = after_spl
        .split('\'')
        .next()
        .unwrap_or("")
        .trim()
        .trim_matches(|c| c == '"' || c == '\'');

    let cleaned: String = raw_fmt
        .chars()
        .filter(|&c| c != '[' && c != ']' && c != '\'' && c != '*')
        .collect();

    let adesc: Vec<String> = cleaned
        .split('|')
        .map(|s| s.trim().to_ascii_uppercase())
        .collect();

    let csq_len = adesc.len();
    let allele = adesc.iter().position(|s| s == "ALLELE");
    let impact = adesc
        .iter()
        .position(|s| s == "IMPACT" || s == "ANNOTATION_IMPACT");
    let gene_symbol = adesc
        .iter()
        .position(|s| s == "SYMBOL" || s == "GENE" || s == "GENE_NAME");
    let consequence = adesc
        .iter()
        .position(|s| s == "CONSEQUENCE" || s == "ANNOTATION");

    if csq_field_name == "ANN" || csq_field_name == "CSQ" {
        if allele.is_none() {
            bail!("unable to find expected Allele field in {csq_field_name} description");
        }
        if impact.is_none() {
            bail!("unable to find expected impact field in {csq_field_name} description");
        }
    }
    let consequence = match consequence {
        Some(idx) => idx,
        None => bail!("unable to find expected consequence field in {csq_field_name} description"),
    };
    let gene_symbol = match gene_symbol {
        Some(idx) => idx,
        None => bail!("unable to find expected gene symbol field in {csq_field_name} description"),
    };

    Ok(CsqSchema {
        allele,
        gene_symbol,
        consequence,
        impact,
        csq_len,
        csq_field_name,
        create_ann,
    })
}

/// Format deduplicated consequence entries into `out` preserving first-seen order
/// (exact equivalent of Nim's `makeAnnField` with zero steady-state heap allocations).
#[allow(clippy::too_many_arguments)]
pub fn format_ann_entries(
    alt_alleles: &[&[u8]],
    gene_consequences: &[(&str, &str)],
    impact_str: &str,
    schema: &CsqSchema,
    entry_buf: &mut String,
    seen_spans: &mut Vec<(usize, usize)>,
    seen_set: &mut HashSet<String>,
    out: &mut String,
) -> bool {
    seen_spans.clear();
    let use_hash = gene_consequences.len() * alt_alleles.len().max(1) > 16;
    if use_hash {
        seen_set.clear();
    }

    let emit_one = |allele_str: Option<&str>,
                    gene: &str,
                    csq: &str,
                    entry_buf: &mut String,
                    seen_spans: &mut Vec<(usize, usize)>,
                    seen_set: &mut HashSet<String>,
                    out: &mut String| {
        entry_buf.clear();
        for col in 0..schema.csq_len {
            if col > 0 {
                entry_buf.push('|');
            }
            if schema.allele == Some(col) {
                if let Some(a) = allele_str {
                    entry_buf.push_str(a);
                }
            } else if col == schema.consequence {
                entry_buf.push_str(csq);
            } else if schema.impact == Some(col) {
                entry_buf.push_str(impact_str);
            } else if col == schema.gene_symbol {
                entry_buf.push_str(gene);
            }
        }

        let is_new = if use_hash {
            if seen_set.contains(entry_buf.as_str()) {
                false
            } else {
                seen_set.insert(entry_buf.clone());
                true
            }
        } else {
            let needle = entry_buf.as_str();
            !seen_spans
                .iter()
                .any(|&(s, e)| &out[s..e] == needle)
        };

        if is_new {
            if !out.is_empty() {
                out.push(',');
            }
            let start = out.len();
            out.push_str(entry_buf);
            if !use_hash {
                seen_spans.push((start, out.len()));
            }
        }
    };

    for &(gene, csq) in gene_consequences {
        if schema.allele.is_some() {
            for &allele_bytes in alt_alleles {
                let allele_str = std::str::from_utf8(allele_bytes).unwrap_or("");
                emit_one(
                    Some(allele_str),
                    gene,
                    csq,
                    entry_buf,
                    seen_spans,
                    seen_set,
                    out,
                );
            }
        } else {
            emit_one(None, gene, csq, entry_buf, seen_spans, seen_set, out);
        }
    }

    !out.is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{DESC_BCSQ, DESC_CSQ};

    #[test]
    fn test_parse_ann_csq_bcsq() {
        let ann = parse_csq_description(DESC_ANN, "ANN".into(), true).unwrap();
        assert_eq!(ann.allele, Some(0));
        assert_eq!(ann.consequence, 1);
        assert_eq!(ann.impact, Some(2));
        assert_eq!(ann.gene_symbol, 3);
        assert_eq!(ann.csq_len, 16);

        let bcsq = parse_csq_description(DESC_BCSQ, "BCSQ".into(), false).unwrap();
        assert_eq!(bcsq.allele, None);
        assert_eq!(bcsq.consequence, 0);
        assert_eq!(bcsq.impact, None);
        assert_eq!(bcsq.gene_symbol, 1);
        assert_eq!(bcsq.csq_len, 7);

        let csq = parse_csq_description(DESC_CSQ, "CSQ".into(), false).unwrap();
        assert_eq!(csq.allele, Some(0));
        assert_eq!(csq.consequence, 1);
        assert_eq!(csq.impact, Some(2));
        assert_eq!(csq.gene_symbol, 3);
        assert_eq!(csq.csq_len, 15);
    }

    #[test]
    fn test_format_ann_entries_dedup_and_multiallelic() {
        let ann = parse_csq_description(DESC_ANN, "ANN".into(), true).unwrap();
        let mut entry = String::new();
        let mut seen_spans = Vec::new();
        let mut seen_s = HashSet::new();
        let mut out = String::new();

        let alts: [&[u8]; 2] = [b"A", b"T"];
        let pairs = [("Gene1", "enhancer"), ("Gene1", "enhancer"), ("Gene2", "promoter")];
        format_ann_entries(
            &alts,
            &pairs,
            "MODIFIER",
            &ann,
            &mut entry,
            &mut seen_spans,
            &mut seen_s,
            &mut out,
        );
        assert_eq!(
            out,
            "A|enhancer|MODIFIER|Gene1||||||||||||,T|enhancer|MODIFIER|Gene1||||||||||||,A|promoter|MODIFIER|Gene2||||||||||||,T|promoter|MODIFIER|Gene2||||||||||||"
        );
    }
}
