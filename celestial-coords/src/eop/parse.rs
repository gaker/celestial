use super::record::{EopFlags, EopQuality, EopRecord, EopSource};
use crate::errors::{CoordError, CoordResult};

pub(super) fn parse_finals(content: &str) -> CoordResult<Vec<EopRecord>> {
    let mut records = Vec::new();

    for line in content.lines() {
        if let Some(record) = parse_finals_line(line)? {
            records.push(record);
        }
    }

    if records.is_empty() {
        return Err(CoordError::parsing_error(
            "No valid records found in finals2000A data",
        ));
    }

    records.sort_by(|a, b| a.mjd.total_cmp(&b.mjd));
    Ok(records)
}

// Finals files list future days with the polar motion and UT1-UTC columns
// blank, so those lines are Ok(None) rather than errors.
fn parse_finals_line(line: &str) -> CoordResult<Option<EopRecord>> {
    let (Some(mjd), Some(xp), Some(yp), Some(ut1_utc)) = (
        parse_field(line, 7, 15)?,
        parse_field(line, 18, 27)?,
        parse_field(line, 37, 46)?,
        parse_field(line, 58, 68)?,
    ) else {
        return Ok(None);
    };
    let mut record = EopRecord::new(mjd, xp, yp, ut1_utc)?;
    if let Some(lod_ms) = parse_field(line, 79, 86)? {
        record = record.with_lod(lod_ms * 0.001)?;
    }

    let cip = (parse_field(line, 97, 106)?, parse_field(line, 116, 125)?);
    let has_cip = matches!(cip, (Some(_), Some(_)));
    if let (Some(dx), Some(dy)) = cip {
        record = record.with_cip_offsets(dx, dy)?;
    }
    Ok(Some(
        record.with_flags(finals_flags(has_cip, finals_quality(line))),
    ))
}

// Columns 17, 58 and 96 mark the polar motion, UT1-UTC and nutation values as
// IERS (I) or predicted (P).
fn finals_quality(line: &str) -> EopQuality {
    let bytes = line.as_bytes();
    if [16, 57, 95].iter().any(|&i| bytes.get(i) == Some(&b'P')) {
        EopQuality::Predicted
    } else {
        EopQuality::HighPrecision
    }
}

fn finals_flags(has_cip_offsets: bool, quality: EopQuality) -> EopFlags {
    EopFlags {
        source: EopSource::IersFinals,
        quality,
        has_polar_motion: true,
        has_ut1_utc: true,
        has_cip_offsets,
        has_pole_rates: false,
    }
}

fn parse_field(line: &str, start: usize, end: usize) -> CoordResult<Option<f64>> {
    let Some(text) = line.get(start..end).map(str::trim) else {
        return Ok(None);
    };
    if text.is_empty() {
        return Ok(None);
    }
    match text.parse::<f64>() {
        Ok(value) if value.is_finite() => Ok(Some(value)),
        _ => Err(CoordError::parsing_error(format!(
            "finals2000A columns {}-{} hold {:?}, not a finite number",
            start + 1,
            end,
            text
        ))),
    }
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
        let record = parse_finals_line(&line).unwrap().unwrap();
        let params = record.to_parameters();

        assert_eq!(params.mjd, 60000.0);
        assert_eq!(params.x_p, 0.1);
        assert_eq!(params.y_p, 0.25);
        assert_eq!(params.ut1_utc, -0.05);
        assert_eq!(params.lod, Some(0.0015));
        assert_eq!(params.dx, Some(0.2));
        assert_eq!(params.dy, Some(-0.1));
        assert_eq!(params.flags.source, EopSource::IersFinals);
        assert!(params.flags.has_cip_offsets);
    }

    fn quality_with_flags(polar_motion: u8, ut1: u8) -> EopQuality {
        let mut line = sample_finals_line().into_bytes();
        line[16] = polar_motion;
        line[57] = ut1;
        let line = String::from_utf8(line).unwrap();
        let record = parse_finals_line(&line).unwrap().unwrap();
        record.to_parameters().flags.quality
    }

    #[test]
    fn test_prediction_flags_set_quality() {
        assert_eq!(quality_with_flags(b'I', b'I'), EopQuality::HighPrecision);
        assert_eq!(quality_with_flags(b'P', b'I'), EopQuality::Predicted);
        assert_eq!(quality_with_flags(b'I', b'P'), EopQuality::Predicted);
    }

    #[test]
    fn test_parse_line_too_short() {
        assert_eq!(parse_finals_line("short line").unwrap(), None);
    }

    #[test]
    fn test_parse_line_missing_required() {
        let line = " ".repeat(188);
        assert_eq!(parse_finals_line(&line).unwrap(), None);
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
    fn test_parse_finals_skips_lines_without_data() {
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

    fn line_with_lod(lod: &[u8; 7]) -> String {
        let mut line = vec![b' '; 188];
        line[7..15].copy_from_slice(b"60000.00");
        line[18..27].copy_from_slice(b"  0.10000");
        line[37..46].copy_from_slice(b"  0.25000");
        line[58..68].copy_from_slice(b" -0.050000");
        line[79..86].copy_from_slice(lod);
        String::from_utf8(line).unwrap()
    }

    #[test]
    fn test_blank_lod_is_absent_not_zero() {
        let blank = parse_finals_line(&line_with_lod(b"       "))
            .unwrap()
            .unwrap();
        let zero = parse_finals_line(&line_with_lod(b"  0.000"))
            .unwrap()
            .unwrap();
        assert_eq!(blank.to_parameters().lod, None);
        assert_eq!(zero.to_parameters().lod, Some(0.0));
    }

    #[test]
    fn test_no_cip_when_blank() {
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
        let record = parse_finals_line(&line).unwrap().unwrap();
        assert!(!record.flags.has_cip_offsets);
        assert_eq!(record.dx_encoded, None);
    }

    #[test]
    fn zero_cip_offsets_are_data() {
        let mut line = sample_finals_line().into_bytes();
        line[97..106].copy_from_slice(b"   0.0000");
        line[116..125].copy_from_slice(b"   0.0000");
        let line = String::from_utf8(line).unwrap();
        let params = parse_finals_line(&line).unwrap().unwrap().to_parameters();
        assert_eq!((params.dx, params.dy), (Some(0.0), Some(0.0)));
        assert!(params.flags.has_cip_offsets);
    }

    #[test]
    fn unreadable_column_is_an_error() {
        let mut line = sample_finals_line().into_bytes();
        line[18..27].copy_from_slice(b"  0.1x000");
        let err = parse_finals_line(&String::from_utf8(line).unwrap()).unwrap_err();
        assert!(err.to_string().contains("columns 19-27"));
    }

    #[test]
    fn non_finite_column_is_an_error() {
        let mut line = sample_finals_line().into_bytes();
        line[7..15].copy_from_slice(b"     NaN");
        assert!(parse_finals_line(&String::from_utf8(line).unwrap()).is_err());
    }

    #[test]
    fn out_of_range_value_is_an_error() {
        let mut line = sample_finals_line().into_bytes();
        line[58..68].copy_from_slice(b"  1.500000");
        let content = String::from_utf8(line).unwrap();
        assert!(parse_finals(&content).is_err());
    }
}
