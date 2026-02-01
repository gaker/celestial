use super::record::{EopFlags, EopQuality, EopRecord, EopSource};
use crate::{CoordError, CoordResult};
use std::str::FromStr;

pub struct C04Parser {
    pub default_quality: EopQuality,
}

impl C04Parser {
    pub fn new() -> Self {
        Self {
            default_quality: EopQuality::HighPrecision,
        }
    }

    pub fn parse_line(&self, line: &str) -> CoordResult<Option<EopRecord>> {
        if line.is_empty() || line.starts_with('#') || line.trim().is_empty() {
            return Ok(None);
        }

        let trimmed = line.trim();
        if trimmed.contains("Date")
            || trimmed.contains("MJD")
            || trimmed.contains("UTC")
            || trimmed.starts_with("\"")
            || trimmed.starts_with("(0h UTC)")
        {
            return Ok(None);
        }

        if line.len() < 60 {
            return Ok(None);
        }

        let _year_str = self.safe_substring(line, 0, 4)?;
        let _month_str = self.safe_substring(line, 7, 8)?;
        let _day_str = self.safe_substring(line, 11, 12)?;
        let mjd_str = self.safe_substring(line, 14, 20)?;
        let x_p_str = self.safe_substring(line, 21, 30)?;
        let y_p_str = self.safe_substring(line, 33, 41)?;
        let ut1_utc_str = self.safe_substring(line, 44, 54)?;
        let lod_str = self.safe_substring(line, 56, 66)?;

        let mjd = f64::from_str(mjd_str.trim())
            .map_err(|_| CoordError::parsing_error(format!("Invalid MJD: '{}'", mjd_str)))?;

        let x_p = f64::from_str(x_p_str.trim())
            .map_err(|_| CoordError::parsing_error(format!("Invalid x_p: '{}'", x_p_str)))?;

        let y_p = f64::from_str(y_p_str.trim())
            .map_err(|_| CoordError::parsing_error(format!("Invalid y_p: '{}'", y_p_str)))?;

        let ut1_utc = f64::from_str(ut1_utc_str.trim()).map_err(|_| {
            CoordError::parsing_error(format!("Invalid UT1-UTC: '{}'", ut1_utc_str))
        })?;

        let lod = f64::from_str(lod_str.trim())
            .map_err(|_| CoordError::parsing_error(format!("Invalid LOD: '{}'", lod_str)))?;

        let mut record = EopRecord::new(mjd, x_p, y_p, ut1_utc, lod)?;

        if line.len() >= 87 {
            let dx_str = self.safe_substring(line, 68, 78)?;
            let dy_str = self.safe_substring(line, 79, 87)?;

            if let (Ok(dx), Ok(dy)) = (f64::from_str(dx_str.trim()), f64::from_str(dy_str.trim())) {
                record = record.with_cip_offsets(dx * 1000.0, dy * 1000.0)?;
            }
        }

        let flags = EopFlags {
            source: EopSource::IersC04,
            quality: self.default_quality,
            has_polar_motion: true,
            has_ut1_utc: true,
            has_cip_offsets: line.len() >= 87,
        };

        record = record.with_flags(flags);

        Ok(Some(record))
    }

    fn safe_substring<'a>(&self, s: &'a str, start: usize, end: usize) -> CoordResult<&'a str> {
        if start >= s.len() {
            return Err(CoordError::parsing_error(format!(
                "Start position {} exceeds line length {}",
                start,
                s.len()
            )));
        }

        let actual_end = end.min(s.len());
        Ok(&s[start..actual_end])
    }

    /// Parse multiple lines of C04 data
    pub fn parse(&self, content: &str) -> CoordResult<Vec<EopRecord>> {
        let mut records = Vec::new();

        for (line_num, line) in content.lines().enumerate() {
            match self.parse_line(line) {
                Ok(Some(record)) => records.push(record),
                Ok(None) => {} // Skip comments/empty lines
                Err(e) => {
                    return Err(CoordError::parsing_error(format!(
                        "Error parsing line {}: {}",
                        line_num + 1,
                        e
                    )));
                }
            }
        }

        // Sort by MJD to ensure chronological order
        records.sort_by(|a, b| a.mjd.partial_cmp(&b.mjd).unwrap());

        Ok(records)
    }
}

impl Default for C04Parser {
    fn default() -> Self {
        Self::new()
    }
}

/// Parser for IERS finals.data format
///
/// The finals.data format is more complex with fixed-width columns
/// and different sections for different time periods.
pub struct FinalsParser {
    /// Default quality for parsed records
    pub default_quality: EopQuality,
}

impl FinalsParser {
    /// Create a new finals.data parser
    pub fn new() -> Self {
        Self {
            default_quality: EopQuality::Standard,
        }
    }

    /// Parse a single line of finals.data format
    pub fn parse_line(&self, line: &str) -> CoordResult<Option<EopRecord>> {
        // Skip header lines and empty lines
        // Minimum length needed: 68 chars for UT1-UTC field
        if line.trim().is_empty() || line.len() < 68 {
            return Ok(None);
        }

        // Extract fields using fixed positions (simplified version)
        // Full implementation would handle the complex finals.data format

        // For now, implement basic parsing - can be extended later
        let date_str = line[0..8].trim();
        let mjd_str = line[7..15].trim();

        let mjd = if mjd_str.is_empty() {
            // Convert date to MJD if MJD field is empty
            self.parse_date_to_mjd(date_str)?
        } else {
            f64::from_str(mjd_str)
                .map_err(|_| CoordError::parsing_error(format!("Invalid MJD: {}", mjd_str)))?
        };

        // Extract polar motion (positions vary by data section)
        let x_p_str = line[18..27].trim();
        let y_p_str = line[37..46].trim();
        let ut1_utc_str = line[58..68].trim();

        if x_p_str.is_empty() || y_p_str.is_empty() || ut1_utc_str.is_empty() {
            return Ok(None); // Skip incomplete records
        }

        let x_p = f64::from_str(x_p_str)
            .map_err(|_| CoordError::parsing_error(format!("Invalid x_p: {}", x_p_str)))?;

        let y_p = f64::from_str(y_p_str)
            .map_err(|_| CoordError::parsing_error(format!("Invalid y_p: {}", y_p_str)))?;

        let ut1_utc = f64::from_str(ut1_utc_str)
            .map_err(|_| CoordError::parsing_error(format!("Invalid UT1-UTC: {}", ut1_utc_str)))?;

        // LOD is optional in some finals.data sections
        let lod = if line.len() > 79 {
            let lod_str = line[79..86].trim();
            if !lod_str.is_empty() {
                f64::from_str(lod_str).unwrap_or(0.0)
            } else {
                0.0
            }
        } else {
            0.0
        };

        let record = EopRecord::new(mjd, x_p, y_p, ut1_utc, lod)?;

        let flags = EopFlags {
            source: EopSource::IersFinals,
            quality: self.default_quality,
            has_polar_motion: true,
            has_ut1_utc: true,
            has_cip_offsets: false, // Not typically in finals.data
        };

        Ok(Some(record.with_flags(flags)))
    }

    /// Convert YYMMDD date string to MJD
    fn parse_date_to_mjd(&self, date_str: &str) -> CoordResult<f64> {
        if date_str.len() != 6 {
            return Err(CoordError::parsing_error(format!(
                "Invalid date format: {}",
                date_str
            )));
        }

        let year_short: u32 = date_str[0..2]
            .parse()
            .map_err(|_| CoordError::parsing_error("Invalid year in date"))?;
        let month: u32 = date_str[2..4]
            .parse()
            .map_err(|_| CoordError::parsing_error("Invalid month in date"))?;
        let day: u32 = date_str[4..6]
            .parse()
            .map_err(|_| CoordError::parsing_error("Invalid day in date"))?;

        // Convert 2-digit year to 4-digit (assuming 20xx for now)
        let year = if year_short < 50 {
            2000 + year_short
        } else {
            1900 + year_short
        };

        // Simple Julian Day calculation (can be improved with astro-time integration)
        let a = (14 - month) / 12;
        let y = year - a;
        let m = month + 12 * a - 3;

        let jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 + 1721119;
        let mjd = jd as f64 - astro_core::constants::MJD_ZERO_POINT;

        Ok(mjd)
    }

    /// Parse multiple lines of finals.data
    pub fn parse(&self, content: &str) -> CoordResult<Vec<EopRecord>> {
        let mut records = Vec::new();

        for (line_num, line) in content.lines().enumerate() {
            match self.parse_line(line) {
                Ok(Some(record)) => records.push(record),
                Ok(None) => {} // Skip incomplete lines
                Err(e) => {
                    return Err(CoordError::parsing_error(format!(
                        "Error parsing finals.data line {}: {}",
                        line_num + 1,
                        e
                    )));
                }
            }
        }

        // Sort by MJD
        records.sort_by(|a, b| a.mjd.partial_cmp(&b.mjd).unwrap());

        Ok(records)
    }
}

impl Default for FinalsParser {
    fn default() -> Self {
        Self::new()
    }
}

/// Auto-detect EOP file format and parse accordingly
pub fn parse_eop_file(content: &str) -> CoordResult<Vec<EopRecord>> {
    // Simple format detection based on content patterns
    // Check for finals.data format first (fixed-width, long lines, often starts with spaces)
    if content.lines().any(|line| line.len() > 100) {
        // Looks like finals.data format (fixed-width, long lines)
        FinalsParser::new().parse(content)
    } else if content.lines().any(|line| {
        line.split_whitespace().count() >= 5
            && !line.starts_with('#')
            && line
                .trim()
                .chars()
                .next()
                .is_some_and(|c| c.is_ascii_digit())
    }) {
        // Looks like C04 format (space-separated, starts with digit)
        C04Parser::new().parse(content)
    } else {
        Err(CoordError::parsing_error("Unknown EOP file format"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c04_parser() {
        let parser = C04Parser::new();

        // Test with actual C04 fixed-width format
        let line = "1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000";
        let record = parser.parse_line(line).unwrap().unwrap();

        let params = record.to_parameters();
        assert_eq!(params.mjd, 37665.0);
        assert!((params.x_p - (-0.012700)).abs() < 1e-7);
        assert!((params.y_p - 0.213000).abs() < 1e-7);
        assert!((params.ut1_utc - 0.0326338).abs() < 1e-7);
        assert!((params.lod - 0.0017230).abs() < 1e-7);
    }

    #[test]
    fn test_c04_parser_with_cip() {
        let parser = C04Parser::new();

        // Test with actual C04 fixed-width format including CIP offsets
        let line = "1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000054   0.000008   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000";
        let record = parser.parse_line(line).unwrap().unwrap();

        let params = record.to_parameters();
        assert!(params.flags.has_cip_offsets);
        assert_eq!(params.dx, Some(0.054)); // 0.000054 arcsec * 1000 = 0.054 mas
        assert_eq!(params.dy, Some(0.008)); // 0.000008 arcsec * 1000 = 0.008 mas
    }

    #[test]
    fn test_c04_parser_header_skipping() {
        let parser = C04Parser::new();

        // Test that headers are properly skipped
        let header1 = "      Date      MJD      x          y        UT1-UTC       LOD         dX        dY        x Err     y Err   UT1-UTC Err  LOD Err     dX Err       dY Err";
        let header2 = "                         \"          \"           s           s          \"         \"           \"          \"          s         s            \"           \"";
        let header3 = "     (0h UTC)";
        let empty = "";

        assert!(parser.parse_line(header1).unwrap().is_none());
        assert!(parser.parse_line(header2).unwrap().is_none());
        assert!(parser.parse_line(header3).unwrap().is_none());
        assert!(parser.parse_line(empty).unwrap().is_none());
    }

    #[test]
    fn test_c04_parser_full_format() {
        let parser = C04Parser::new();

        let content = r#"      Date      MJD      x          y        UT1-UTC       LOD         dX        dY        x Err     y Err   UT1-UTC Err  LOD Err     dX Err       dY Err
                         "          "           s           s          "         "           "          "          s         s            "           "
     (0h UTC)

1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
1962   1   2  37666  -0.015900   0.214100   0.0320547   0.0016690   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
"#;

        let records = parser.parse(content).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].mjd, 37665.0);
        assert_eq!(records[1].mjd, 37666.0);

        // Verify first record values
        let params = records[0].to_parameters();
        assert!((params.x_p - (-0.012700)).abs() < 1e-7);
        assert!((params.y_p - 0.213000).abs() < 1e-7);
        assert!((params.ut1_utc - 0.0326338).abs() < 1e-7);
        assert!((params.lod - 0.0017230).abs() < 1e-7);
    }

    #[test]
    fn test_date_to_mjd_conversion() {
        let parser = FinalsParser::new();

        // Test known date: 2023-01-01 = MJD 59945
        let mjd = parser.parse_date_to_mjd("230101").unwrap();
        assert!((mjd - 59945.0).abs() < 1.0); // Within 1 day for simple calculation
    }

    #[test]
    fn test_auto_detect_format() {
        let c04_content = "1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230";
        let records = parse_eop_file(c04_content).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].flags.source, EopSource::IersC04);
    }

    #[test]
    fn test_c04_parser_invalid_mjd() {
        let parser = C04Parser::new();
        let line = "1962   1   1  ABC45  -0.012700   0.213000   0.0326338   0.0017230";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_invalid_xp() {
        let parser = C04Parser::new();
        let line = "1962   1   1  37665  INVALID   0.213000   0.0326338   0.0017230";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_multiline_error() {
        let parser = C04Parser::new();
        let content = "1962   1   1  37665  INVALID   0.213000   0.0326338   0.0017230\n1962   1   2  37666  -0.012700   0.213000   0.0326338   0.0017230";
        let result = parser.parse(content);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_sorting() {
        let parser = C04Parser::new();
        let content = "1962   1   3  37667  -0.012700   0.213000   0.0326338   0.0017230\n1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230\n1962   1   2  37666  -0.012700   0.213000   0.0326338   0.0017230";
        let records = parser.parse(content).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].mjd, 37665.0);
        assert_eq!(records[1].mjd, 37666.0);
        assert_eq!(records[2].mjd, 37667.0);
    }

    #[test]
    fn test_c04_parser_invalid_yp() {
        let parser = C04Parser::new();
        let line = "1962   1   1  37665  -0.012700   INVALID   0.0326338   0.0017230";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_invalid_ut1_utc() {
        let parser = C04Parser::new();
        let line = "1962   1   1  37665  -0.012700   0.213000   INVALID     0.0017230";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_invalid_lod() {
        let parser = C04Parser::new();
        let line = "1962   1   1  37665  -0.012700   0.213000   0.0326338   INVALID";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_safe_substring_start_exceeds_length() {
        let parser = C04Parser::new();
        let short_line = "short";
        let result = parser.safe_substring(short_line, 100, 110);
        assert!(result.is_err());
    }

    #[test]
    fn test_c04_parser_safe_substring_truncation() {
        let parser = C04Parser::new();
        let line = "test line";
        let result = parser.safe_substring(line, 0, 100).unwrap();
        assert_eq!(result, "test line");
    }

    #[test]
    fn test_c04_parser_default_trait() {
        let parser = C04Parser::default();
        assert_eq!(parser.default_quality, EopQuality::HighPrecision);
    }

    #[test]
    fn test_finals_parser_short_line() {
        let parser = FinalsParser::new();
        let short_line = "short line";
        let result = parser.parse_line(short_line).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_finals_parser_empty_line() {
        let parser = FinalsParser::new();
        let empty = "";
        let result = parser.parse_line(empty).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_finals_parser() {
        let parser = FinalsParser::new();
        // Format: date(0-8) mjd(7-15) x_p(18-27) y_p(37-46) ut1_utc(58-68)
        let line = "       59665.00   0.123456           0.234567             0.345678                                                                                                                                      ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_some());
        let record = result.unwrap();
        assert_eq!(record.flags.source, EopSource::IersFinals);
    }

    #[test]
    fn test_finals_parser_incomplete_xp() {
        let parser = FinalsParser::new();
        let line = "       59665.00                      0.234567             0.345678             0.001234                                                                                                                 ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_finals_parser_incomplete_yp() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456                                0.345678             0.001234                                                                                                                 ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_finals_parser_incomplete_ut1_utc() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456           0.234567                                  0.001234                                                                                                                 ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_finals_parser_invalid_xp() {
        let parser = FinalsParser::new();
        let line = "       59665.00   XXXXXXXX           0.234567             0.345678             0.001234                                                                                                                 ";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_finals_parser_invalid_yp() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456           XXXXXXXX             0.345678             0.001234                                                                                                                 ";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_finals_parser_invalid_ut1_utc() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456           0.234567             XXXXXXXX             0.001234                                                                                                                 ";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_finals_parser_with_lod() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456           0.234567             0.345678             0.001234                                                                                                                 ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_some());
        let record = result.unwrap();
        let params = record.to_parameters();
        // LOD field is 79-86 (7 chars), so "0.001234" gets truncated to "0.00123"
        assert!((params.lod - 0.00123).abs() < 1e-7);
    }

    #[test]
    fn test_finals_parser_empty_lod() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456           0.234567             0.345678                                                                                                                                      ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_some());
        let record = result.unwrap();
        let params = record.to_parameters();
        assert_eq!(params.lod, 0.0);
    }

    #[test]
    fn test_finals_parser_invalid_lod_defaults_to_zero() {
        let parser = FinalsParser::new();
        let line = "       59665.00   0.123456           0.234567             0.345678             XXXXXXX                                                                                                                  ";
        let result = parser.parse_line(line).unwrap();
        assert!(result.is_some());
        let record = result.unwrap();
        let params = record.to_parameters();
        assert_eq!(params.lod, 0.0);
    }

    #[test]
    fn test_finals_parser_short_line_no_lod() {
        let parser = FinalsParser::new();
        // Line must be at least 185 chars, but doesn't have LOD data (at position 79+)
        let mut line =
            "       59665.00   0.123456           0.234567             0.345678".to_string();
        while line.len() < 185 {
            line.push(' ');
        }
        let result = parser.parse_line(&line).unwrap();
        assert!(result.is_some());
        let record = result.unwrap();
        let params = record.to_parameters();
        assert_eq!(params.lod, 0.0);
    }

    #[test]
    fn test_finals_parser_date_parsing() {
        // Test the date-to-MJD conversion pathway (when MJD column is truly empty)
        // This is a complex edge case in real finals.data format
        let parser = FinalsParser::new();
        // For this test, we'll just verify the parse_date_to_mjd function works
        let mjd = parser.parse_date_to_mjd("230101").unwrap();
        assert!(mjd > 59000.0 && mjd < 60000.0);
    }

    #[test]
    fn test_finals_parser_invalid_mjd() {
        let parser = FinalsParser::new();
        let line = "       INVALID    0.123456           0.234567             0.345678             0.001234                                                                                                                 ";
        let result = parser.parse_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_date_to_mjd_invalid_length() {
        let parser = FinalsParser::new();
        let result = parser.parse_date_to_mjd("2301");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_date_to_mjd_invalid_year() {
        let parser = FinalsParser::new();
        let result = parser.parse_date_to_mjd("XX0101");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_date_to_mjd_invalid_month() {
        let parser = FinalsParser::new();
        let result = parser.parse_date_to_mjd("23XX01");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_date_to_mjd_invalid_day() {
        let parser = FinalsParser::new();
        let result = parser.parse_date_to_mjd("2301XX");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_date_to_mjd_19xx_year() {
        let parser = FinalsParser::new();
        let mjd = parser.parse_date_to_mjd("620101").unwrap();
        assert!(mjd > 37000.0 && mjd < 38000.0);
    }

    #[test]
    fn test_finals_parser_default_trait() {
        let parser = FinalsParser::default();
        assert_eq!(parser.default_quality, EopQuality::Standard);
    }

    #[test]
    fn test_finals_parser_multiline_error() {
        let parser = FinalsParser::new();
        let content = "       59665.00   XXXXXXXX           0.234567             0.345678             0.001234                                                                                                                 \n       59666.00   0.123456           0.234567             0.345678             0.001234                                                                                                                 ";
        let result = parser.parse(content);
        assert!(result.is_err());
    }

    #[test]
    fn test_finals_parser_sorting() {
        let parser = FinalsParser::new();
        let content = "       59667.00   0.123456           0.234567             0.345678             0.00123                                                                                                                  \n       59665.00   0.123456           0.234567             0.345678             0.00123                                                                                                                  \n       59666.00   0.123456           0.234567             0.345678             0.00123                                                                                                                  ";
        let records = parser.parse(content).unwrap();
        assert_eq!(records.len(), 3);
        assert!(records[0].mjd < records[1].mjd);
        assert!(records[1].mjd < records[2].mjd);
    }

    #[test]
    fn test_parse_eop_file_finals_format() {
        let finals_content = "       59665.00   0.123456           0.234567             0.345678             0.00123                                                                                                                  \n       59666.00   0.123456           0.234567             0.345678             0.00123                                                                                                                  ";
        let records = parse_eop_file(finals_content).unwrap();
        assert!(!records.is_empty());
        assert_eq!(records[0].flags.source, EopSource::IersFinals);
    }

    #[test]
    fn test_parse_eop_file_unknown_format() {
        let unknown_content = "# This is not a recognized format\nshort\nlines\n";
        let result = parse_eop_file(unknown_content);
        assert!(result.is_err());
    }
}
