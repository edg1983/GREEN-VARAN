use rust_htslib::bcf::{record::Numeric, Record};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SvInterval {
    pub start: i64,
    pub stop: i64,
    pub is_ins_bnd: bool,
}

impl SvInterval {
    #[inline]
    pub fn len(&self) -> i64 {
        self.stop - self.start
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.stop <= self.start
    }
}

/// Build an SV interval from a BCF/VCF record.
/// Includes the requested fix: reads `INFO/SVTYPE` so `BND` and `INS` properly set
/// `is_ins_bnd = true` and apply `--padding` when provided (overriding `CIPOS`).
pub fn make_sv_interval(
    rec: &Record,
    cipos_tag: &[u8],
    ciend_tag: &[u8],
    padding: Option<i64>,
) -> Option<SvInterval> {
    let pos_1based = rec.pos() + 1;

    // 1. Determine SV end from INFO/END or INFO/SVLEN
    let svend: i64 = if let Ok(Some(end_vals)) = rec.info(b"END").integer() {
        let v = *end_vals.first()?;
        if v.is_missing() {
            return None;
        }
        v as i64
    } else if let Ok(Some(len_vals)) = rec.info(b"SVLEN").integer() {
        let v = *len_vals.first()?;
        if v.is_missing() {
            return None;
        }
        pos_1based + (v as i64).abs()
    } else {
        return None;
    };

    // 2. Read SVTYPE to detect BND / INS (fixing Nim's unread svtype bug)
    let is_ins_bnd = if let Ok(Some(vals)) = rec.info(b"SVTYPE").string() {
        vals.first()
            .map(|&s| s == b"BND" || s == b"INS")
            .unwrap_or(false)
    } else {
        false
    };

    // 3. Apply CIPOS and CIEND
    let cipos0: i64 = match rec.info(cipos_tag).integer() {
        Ok(Some(vals)) if !vals.is_empty() && !vals[0].is_missing() => (vals[0] as i64).abs(),
        _ => 0,
    };
    let ciend1: i64 = match rec.info(ciend_tag).integer() {
        Ok(Some(vals)) if vals.len() >= 2 && !vals[1].is_missing() => vals[1] as i64,
        _ => 0,
    };

    let mut start = pos_1based - cipos0;
    let mut stop = svend + ciend1;

    // 4. For BND/INS, override with fixed --padding when set
    if is_ins_bnd {
        if let Some(pad) = padding {
            start = pos_1based - pad;
            stop = svend + pad;
        }
    }

    if start > stop {
        return None;
    }

    Some(SvInterval {
        start,
        stop,
        is_ins_bnd,
    })
}

/// Evaluate whether GREEN-DB region `[r_start, r_stop)` is sufficiently overlapped by `sv`.
/// Includes the requested fix: uses exact clamped intersection
/// `min(r_stop, sv.stop) - max(r_start, sv.start)`.
#[inline]
pub fn region_overlapped_by_sv(
    r_start: i64,
    r_stop: i64,
    sv: SvInterval,
    min_overlap: f64,
    min_bp: i64,
) -> bool {
    let overlap_n = r_stop.min(sv.stop) - r_start.max(sv.start);
    if overlap_n <= 0 {
        return false;
    }
    if sv.is_ins_bnd {
        return true;
    }
    let r_len = r_stop - r_start;
    if r_len <= 0 {
        return false;
    }
    let overlap_pct = (overlap_n as f64) / (r_len as f64);
    overlap_n >= min_bp && overlap_pct >= min_overlap
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clamped_overlap_contained_sv() {
        // Region [1000, 2000) (len 1000), SV [1200, 1250) (len 50, contained inside region)
        // True overlap = 50 bp = 5% of region. (Old Nim bug computed 2000-1200 = 800 bp = 80%!)
        let sv = SvInterval {
            start: 1200,
            stop: 1250,
            is_ins_bnd: false,
        };
        assert!(region_overlapped_by_sv(1000, 2000, sv, 0.05, 50));
        assert!(!region_overlapped_by_sv(1000, 2000, sv, 0.10, 1));
        assert!(!region_overlapped_by_sv(1000, 2000, sv, 0.01, 100));
    }

    #[test]
    fn test_clamped_overlap_contained_region() {
        // Region [1200, 1400) (len 200), SV [1000, 5000) (completely covers region)
        let sv = SvInterval {
            start: 1000,
            stop: 5000,
            is_ins_bnd: false,
        };
        assert!(region_overlapped_by_sv(1200, 1400, sv, 1.0, 200));
        assert!(!region_overlapped_by_sv(1200, 1400, sv, 1.0, 201));
    }

    #[test]
    fn test_ins_bnd_ignores_thresholds() {
        let sv = SvInterval {
            start: 1500,
            stop: 1501,
            is_ins_bnd: true,
        };
        assert!(region_overlapped_by_sv(1000, 2000, sv, 0.9, 500));
        assert!(!region_overlapped_by_sv(1000, 1500, sv, 0.0, 1));
    }
}
