use super::record::{EopFlags, EopQuality, EopRecord, EopSource};
use crate::CoordResult;

pub fn parse_finals(content: &str) -> CoordResult<Vec<EopRecord>> {
    let mut records = Vec::new();

    for line in content.lines() {
        if let Some(record) = parse_finals_line(line) {
            records.push(record);
        }
    }

    if records.is_empty() {
        return Err(crate::CoordError::parsing_error(
            "No valid records found in finals2000A data",
        ));
    }

    records.sort_by(|a, b| a.mjd.partial_cmp(&b.mjd).unwrap());
    Ok(records)
}

pub fn parse_finals_line(line: &str) -> Option<EopRecord> {
    if line.len() < 79 {
        return None;
    }

    let mjd = parse_field(line, 7, 15)?;
    let xp = parse_field(line, 18, 27)?;
    let yp = parse_field(line, 37, 46)?;
    let ut1_utc = parse_field(line, 58, 68)?;
    let lod = parse_field(line, 79, 86).unwrap_or(0.0) * 0.001;
    let dx = parse_field(line, 97, 106).unwrap_or(0.0);
    let dy = parse_field(line, 116, 125).unwrap_or(0.0);

    let mut record = EopRecord::new(mjd, xp, yp, ut1_utc, lod).ok()?;

    let has_cip = dx != 0.0 || dy != 0.0;
    if has_cip {
        record = record.with_cip_offsets(dx, dy).ok()?;
    }

    let flags = EopFlags {
        source: EopSource::IersFinals,
        quality: EopQuality::HighPrecision,
        has_polar_motion: true,
        has_ut1_utc: true,
        has_cip_offsets: has_cip,
        has_pole_rates: false,
    };
    record = record.with_flags(flags);

    Some(record)
}

fn parse_field(line: &str, start: usize, end: usize) -> Option<f64> {
    let s = line.get(start..end)?.trim();
    if s.is_empty() {
        return None;
    }
    s.parse::<f64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_finals_line() -> String {
        let mut line = vec![b' '; 188];

        let mjd = b"60000.00";
        line[7..15].copy_from_slice(mjd);

        let xp = b"  0.10000";
        line[18..27].copy_from_slice(xp);

        let yp = b"  0.25000";
        line[37..46].copy_from_slice(yp);

        let ut1 = b" -0.050000";
        line[58..68].copy_from_slice(ut1);

        // LOD in milliseconds
        let lod = b"  1.500";
        line[79..86].copy_from_slice(lod);

        let dx = b"   0.2000";
        line[97..106].copy_from_slice(dx);

        let dy = b"  -0.1000";
        line[116..125].copy_from_slice(dy);

        String::from_utf8(line).unwrap()
    }

    #[test]
    fn test_parse_single_line() {
        let line = sample_finals_line();
        let record = parse_finals_line(&line).unwrap();
        let params = record.to_parameters();

        assert_eq!(params.mjd, 60000.0);
        assert!((params.x_p - 0.1).abs() < 1e-6);
        assert!((params.y_p - 0.25).abs() < 1e-6);
        assert!((params.ut1_utc - (-0.05)).abs() < 1e-6);
        assert!((params.lod - 0.0015).abs() < 1e-7);
        assert_eq!(params.dx, Some(0.2));
        assert_eq!(params.dy, Some(-0.1));
        assert_eq!(params.flags.source, EopSource::IersFinals);
        assert!(params.flags.has_cip_offsets);
    }

    #[test]
    fn test_parse_line_too_short() {
        assert!(parse_finals_line("short line").is_none());
    }

    #[test]
    fn test_parse_line_missing_required() {
        let line = " ".repeat(188);
        assert!(parse_finals_line(&line).is_none());
    }

    fn sample_finals_line_at(mjd: &[u8]) -> String {
        let mut line = vec![b' '; 188];
        line[7..7 + mjd.len()].copy_from_slice(mjd);
        let xp = b"  0.10000";
        line[18..27].copy_from_slice(xp);
        let yp = b"  0.25000";
        line[37..46].copy_from_slice(yp);
        let ut1 = b" -0.050000";
        line[58..68].copy_from_slice(ut1);
        let lod = b"  1.500";
        line[79..86].copy_from_slice(lod);
        String::from_utf8(line).unwrap()
    }

    #[test]
    fn test_parse_finals_multi_line() {
        let line1 = sample_finals_line_at(b"60000.00");
        let line2 = sample_finals_line_at(b"60001.00");

        let content = format!("{}\n{}\n", line1, line2);
        let records = parse_finals(&content).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].mjd, 60000.0);
        assert_eq!(records[1].mjd, 60001.0);
    }

    #[test]
    fn test_parse_finals_skips_bad_lines() {
        let good = sample_finals_line();
        let content = format!("bad line\n{}\nalso bad\n", good);
        let records = parse_finals(&content).unwrap();
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn test_parse_finals_empty_errors() {
        let result = parse_finals("bad\nlines\nonly\n");
        assert!(result.is_err());
    }

    #[test]
    fn test_lod_zero_when_missing() {
        let mut line = vec![b' '; 188];

        let mjd = b"60000.00";
        line[7..15].copy_from_slice(mjd);
        let xp = b"  0.10000";
        line[18..27].copy_from_slice(xp);
        let yp = b"  0.25000";
        line[37..46].copy_from_slice(yp);
        let ut1 = b" -0.050000";
        line[58..68].copy_from_slice(ut1);
        // LOD columns left blank

        let line = String::from_utf8(line).unwrap();
        let record = parse_finals_line(&line).unwrap();
        let params = record.to_parameters();
        assert_eq!(params.lod, 0.0);
    }

    #[test]
    fn test_no_cip_when_zero() {
        let mut line = vec![b' '; 188];

        let mjd = b"60000.00";
        line[7..15].copy_from_slice(mjd);
        let xp = b"  0.10000";
        line[18..27].copy_from_slice(xp);
        let yp = b"  0.25000";
        line[37..46].copy_from_slice(yp);
        let ut1 = b" -0.050000";
        line[58..68].copy_from_slice(ut1);
        // dX/dY columns left blank

        let line = String::from_utf8(line).unwrap();
        let record = parse_finals_line(&line).unwrap();
        assert!(!record.flags.has_cip_offsets);
    }
}
