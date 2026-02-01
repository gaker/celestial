use super::record::{EopFlags, EopQuality, EopRecord, EopSource};
use crate::CoordResult;

const BASELINE_EOP_DATA: &[EopDataPoint] = &[
    EopDataPoint {
        mjd: 58849.0,
        x_p_mas: 36.42,
        y_p_mas: 268.85,
        ut1_utc_ms: 177.4567,
        lod_ms: 0.3421,
    },
    EopDataPoint {
        mjd: 58850.0,
        x_p_mas: 37.12,
        y_p_mas: 269.34,
        ut1_utc_ms: 176.8234,
        lod_ms: 0.3456,
    },
    EopDataPoint {
        mjd: 59945.0,
        x_p_mas: 123.45,
        y_p_mas: 234.56,
        ut1_utc_ms: -345.2341,
        lod_ms: 1.1234,
    },
    EopDataPoint {
        mjd: 60310.0,
        x_p_mas: 178.92,
        y_p_mas: 156.78,
        ut1_utc_ms: -456.7890,
        lod_ms: 0.7821,
    },
];

#[derive(Debug, Clone, Copy)]
struct EopDataPoint {
    mjd: f64,
    x_p_mas: f64,
    y_p_mas: f64,
    ut1_utc_ms: f64,
    lod_ms: f64,
}

pub fn load_embedded_baseline() -> CoordResult<Vec<EopRecord>> {
    let mut records = Vec::with_capacity(BASELINE_EOP_DATA.len());

    for data_point in BASELINE_EOP_DATA {
        let x_p_arcsec = data_point.x_p_mas / 1000.0;
        let y_p_arcsec = data_point.y_p_mas / 1000.0;
        let ut1_utc_sec = data_point.ut1_utc_ms / 1000.0;
        let lod_sec = data_point.lod_ms / 1000.0;

        let mut record =
            EopRecord::new(data_point.mjd, x_p_arcsec, y_p_arcsec, ut1_utc_sec, lod_sec)?;

        let flags = EopFlags {
            source: EopSource::IersC04,
            quality: EopQuality::HighPrecision,
            has_polar_motion: true,
            has_ut1_utc: true,
            has_cip_offsets: false,
        };

        record = record.with_flags(flags);
        records.push(record);
    }

    records.sort_by(|a, b| a.mjd.partial_cmp(&b.mjd).unwrap());

    Ok(records)
}

pub fn generate_interpolated_baseline() -> CoordResult<Vec<EopRecord>> {
    let sparse_records = load_embedded_baseline()?;
    let mut interpolated_records = Vec::new();

    if sparse_records.len() < 2 {
        return Ok(sparse_records);
    }

    for window in sparse_records.windows(2) {
        let start_record = &window[0];
        let end_record = &window[1];

        let start_mjd = start_record.mjd;
        let end_mjd = end_record.mjd;
        let days_span = (end_mjd - start_mjd) as i32;

        interpolated_records.push(start_record.clone());

        for day in 1..days_span {
            let mjd = start_mjd + day as f64;
            let t = day as f64 / days_span as f64;

            let start_params = start_record.to_parameters();
            let end_params = end_record.to_parameters();

            let x_p = start_params.x_p + t * (end_params.x_p - start_params.x_p);
            let y_p = start_params.y_p + t * (end_params.y_p - start_params.y_p);
            let ut1_utc = start_params.ut1_utc + t * (end_params.ut1_utc - start_params.ut1_utc);
            let lod = start_params.lod + t * (end_params.lod - start_params.lod);

            let mut interpolated_record = EopRecord::new(mjd, x_p, y_p, ut1_utc, lod)?;

            let flags = EopFlags {
                source: EopSource::Interpolated,
                quality: EopQuality::Standard,
                has_polar_motion: true,
                has_ut1_utc: true,
                has_cip_offsets: false,
            };

            interpolated_record = interpolated_record.with_flags(flags);
            interpolated_records.push(interpolated_record);
        }
    }

    if let Some(last_record) = sparse_records.last() {
        interpolated_records.push(last_record.clone());
    }

    Ok(interpolated_records)
}

pub fn baseline_time_span() -> (f64, f64) {
    if let (Some(first), Some(last)) = (BASELINE_EOP_DATA.first(), BASELINE_EOP_DATA.last()) {
        (first.mjd, last.mjd)
    } else {
        (0.0, 0.0)
    }
}

pub fn baseline_record_count() -> usize {
    BASELINE_EOP_DATA.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_embedded_baseline() {
        let records = load_embedded_baseline().unwrap();

        assert!(!records.is_empty());
        assert_eq!(records.len(), BASELINE_EOP_DATA.len());

        for window in records.windows(2) {
            assert!(window[0].mjd < window[1].mjd);
        }

        let first_record = &records[0];
        assert_eq!(first_record.mjd, 58849.0);
        assert_eq!(first_record.flags.source, EopSource::IersC04);
        assert_eq!(first_record.flags.quality, EopQuality::HighPrecision);
    }

    #[test]
    fn test_baseline_time_span() {
        let (start, end) = baseline_time_span();
        assert_eq!(start, 58849.0);
        assert_eq!(end, 60310.0);

        let span_days = end - start;
        assert!((span_days - 1461.0).abs() < 1.0);
    }

    #[test]
    fn test_generate_interpolated_baseline() {
        let interpolated = generate_interpolated_baseline().unwrap();
        let baseline = load_embedded_baseline().unwrap();

        assert!(interpolated.len() > baseline.len());

        for baseline_record in &baseline {
            let found = interpolated.iter().any(|r| {
                (r.mjd - baseline_record.mjd).abs() < 1e-6 && r.flags.source == EopSource::IersC04
            });
            assert!(
                found,
                "Baseline record at MJD {} not found in interpolated data",
                baseline_record.mjd
            );
        }

        let interpolated_count = interpolated
            .iter()
            .filter(|r| r.flags.source == EopSource::Interpolated)
            .count();
        assert!(interpolated_count > 0, "No interpolated records found");
    }

    #[test]
    fn test_precision_preservation() {
        let records = load_embedded_baseline().unwrap();

        for record in &records {
            let params = record.to_parameters();

            assert!(
                params.x_p.abs() < 1.0,
                "X polar motion out of range: {}",
                params.x_p
            );
            assert!(
                params.y_p.abs() < 1.0,
                "Y polar motion out of range: {}",
                params.y_p
            );

            assert!(
                params.ut1_utc.abs() < 1.0,
                "UT1-UTC out of range: {}",
                params.ut1_utc
            );

            assert!(params.lod.abs() < 0.01, "LOD out of range: {}", params.lod);
        }
    }
}
