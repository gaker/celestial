use super::record::{EopFlags, EopQuality, EopRecord, EopSource};
use crate::CoordResult;

pub fn load_bundled_c04() -> CoordResult<Vec<EopRecord>> {
    let entries = celestial_eop_data::c04_data();
    convert_entries(entries, EopSource::IersC04)
}

pub fn load_bundled_combined() -> CoordResult<Vec<EopRecord>> {
    let c04 = celestial_eop_data::c04_data();
    let finals = celestial_eop_data::finals_data();
    let c04_max = c04.last().map(|e| e.mjd).unwrap_or(0.0);

    let mut records = convert_entries(c04, EopSource::IersC04)?;
    let finals_ext = finals.iter().filter(|e| e.mjd > c04_max);
    for entry in finals_ext {
        records.push(convert_entry(entry, EopSource::IersFinals)?);
    }
    Ok(records)
}

pub fn bundled_time_span() -> (f64, f64) {
    celestial_eop_data::data_time_span()
}

pub fn bundled_timestamp() -> &'static str {
    celestial_eop_data::data_timestamp()
}

fn convert_entries(
    entries: &[celestial_eop_data::EopEntry],
    source: EopSource,
) -> CoordResult<Vec<EopRecord>> {
    let mut records = Vec::with_capacity(entries.len());
    for entry in entries {
        records.push(convert_entry(entry, source)?);
    }
    Ok(records)
}

fn convert_entry(
    entry: &celestial_eop_data::EopEntry,
    source: EopSource,
) -> CoordResult<EopRecord> {
    let mut record = EopRecord::new(entry.mjd, entry.x_p, entry.y_p, entry.ut1_utc, entry.lod)?;

    let has_cip = entry.dx != 0.0 || entry.dy != 0.0;
    if has_cip {
        record = record.with_cip_offsets(entry.dx, entry.dy)?;
    }

    let flags = EopFlags {
        source,
        quality: EopQuality::HighPrecision,
        has_polar_motion: true,
        has_ut1_utc: true,
        has_cip_offsets: has_cip,
        has_pole_rates: false,
    };
    record = record.with_flags(flags);
    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_bundled_c04() {
        let records = load_bundled_c04().unwrap();
        assert!(records.len() > 20000);

        for window in records.windows(2) {
            assert!(window[0].mjd < window[1].mjd);
        }

        let first = &records[0];
        assert_eq!(first.flags.source, EopSource::IersC04);
        assert_eq!(first.flags.quality, EopQuality::HighPrecision);
    }

    #[test]
    fn test_load_bundled_combined() {
        let combined = load_bundled_combined().unwrap();
        let c04_only = load_bundled_c04().unwrap();
        assert!(combined.len() > c04_only.len());

        for window in combined.windows(2) {
            assert!(window[0].mjd < window[1].mjd);
        }
    }

    #[test]
    fn test_bundled_time_span() {
        let (start, end) = bundled_time_span();
        assert!(start <= 37665.0); // 1962-01-01
        assert!(end >= 60000.0); // ~2023
    }

    #[test]
    fn test_bundled_timestamp() {
        let ts = bundled_timestamp();
        assert_eq!(ts.len(), 10); // YYYY-MM-DD
        assert!(ts.starts_with("20"));
    }

    #[test]
    fn test_precision_preservation() {
        let records = load_bundled_c04().unwrap();
        for record in records.iter().take(100) {
            let params = record.to_parameters();
            assert!(params.x_p.abs() < 1.0);
            assert!(params.y_p.abs() < 1.0);
            assert!(params.ut1_utc.abs() < 1.0);
            assert!(params.lod.abs() < 0.01);
        }
    }
}
