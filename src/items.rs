use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::cli::STD_CHROMS;

/// Parse a CLI list argument matching Nim's `getItems`:
/// - If `raw` is non-empty and points to an existing readable file, read one item per line.
/// - Otherwise, if `raw` is non-empty, split by `,`.
/// - If resulting list is empty and `class == "chromosome"`, return `STD_CHROMS`
///   (stripping `"chr"` prefix when `nochr == true`).
pub fn get_items(raw: Option<&str>, class: &str, nochr: bool) -> Vec<String> {
    let mut result = Vec::new();
    if let Some(s) = raw {
        if !s.is_empty() {
            if let Ok(file) = File::open(s) {
                let reader = BufReader::new(file);
                for line in reader.lines().map_while(Result::ok) {
                    let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
                    result.push(trimmed.to_string());
                }
            } else {
                result.extend(s.split(',').map(|p| p.to_string()));
            }
        }
    }

    if result.is_empty() && class == "chromosome" {
        if nochr {
            result.extend(STD_CHROMS.iter().map(|c| c.replace("chr", "")));
        } else {
            result.extend(STD_CHROMS.iter().map(|&c| c.to_string()));
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_items_default_chroms() {
        let with_chr = get_items(None, "chromosome", false);
        assert_eq!(with_chr.len(), 25);
        assert_eq!(with_chr[0], "chr1");
        assert_eq!(with_chr[24], "chrM");

        let no_chr = get_items(None, "chromosome", true);
        assert_eq!(no_chr[0], "1");
        assert_eq!(no_chr[24], "M");
    }

    #[test]
    fn test_get_items_comma_list() {
        let items = get_items(Some("chr1,chr22"), "chromosome", false);
        assert_eq!(items, vec!["chr1", "chr22"]);
    }
}
