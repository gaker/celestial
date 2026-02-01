use crate::{JulianDate, TimeError, TimeResult};

#[derive(Debug, Clone)]
pub struct ParsedDateTime {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: f64,
}

impl ParsedDateTime {
    pub fn to_julian_date(&self) -> JulianDate {
        JulianDate::from_calendar(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
        )
    }
}

pub fn parse_iso8601(s: &str) -> TimeResult<ParsedDateTime> {
    let s = s.trim();

    const MAX_ISO8601_LENGTH: usize = 32;
    if s.len() > MAX_ISO8601_LENGTH {
        return Err(TimeError::ParseError("Input too long".to_string()));
    }

    let s = s.strip_suffix('Z').unwrap_or(s);

    let separator_pos = s.find('T').or_else(|| s.find(' ')).ok_or_else(|| {
        TimeError::ParseError(format!(
            "Invalid datetime format: '{}'. Expected YYYY-MM-DDTHH:MM:SS",
            s
        ))
    })?;

    let (date_part, time_part_with_sep) = s.split_at(separator_pos);
    let time_part = &time_part_with_sep[1..];

    let date_components: Vec<&str> = date_part.split('-').collect();
    if date_components.len() != 3 {
        return Err(TimeError::ParseError(format!(
            "Invalid date format: '{}'. Expected YYYY-MM-DD",
            date_part
        )));
    }

    let year = if date_components[0].len() == 4 {
        let bytes = date_components[0].as_bytes();
        if bytes.iter().all(|&b| b.is_ascii_digit()) {
            (bytes[0] - b'0') as i32 * 1000
                + (bytes[1] - b'0') as i32 * 100
                + (bytes[2] - b'0') as i32 * 10
                + (bytes[3] - b'0') as i32
        } else {
            return Err(TimeError::ParseError(format!(
                "Invalid year: '{}'",
                date_components[0]
            )));
        }
    } else {
        return Err(TimeError::ParseError(format!(
            "Invalid year format: '{}'",
            date_components[0]
        )));
    };

    let month = match date_components[1].len() {
        1 => {
            let b = date_components[1].as_bytes()[0];
            if b.is_ascii_digit() {
                b - b'0'
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid month: '{}'",
                    date_components[1]
                )));
            }
        }
        2 => {
            let bytes = date_components[1].as_bytes();
            if bytes.iter().all(|&b| b.is_ascii_digit()) {
                (bytes[0] - b'0') * 10 + (bytes[1] - b'0')
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid month: '{}'",
                    date_components[1]
                )));
            }
        }
        _ => {
            return Err(TimeError::ParseError(format!(
                "Invalid month format: '{}'",
                date_components[1]
            )))
        }
    };

    let day = match date_components[2].len() {
        1 => {
            let b = date_components[2].as_bytes()[0];
            if b.is_ascii_digit() {
                b - b'0'
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid day: '{}'",
                    date_components[2]
                )));
            }
        }
        2 => {
            let bytes = date_components[2].as_bytes();
            if bytes.iter().all(|&b| b.is_ascii_digit()) {
                (bytes[0] - b'0') * 10 + (bytes[1] - b'0')
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid day: '{}'",
                    date_components[2]
                )));
            }
        }
        _ => {
            return Err(TimeError::ParseError(format!(
                "Invalid day format: '{}'",
                date_components[2]
            )))
        }
    };

    if !(1..=12).contains(&month) {
        return Err(TimeError::ParseError(format!(
            "Month out of range: {}",
            month
        )));
    }
    if !(1..=31).contains(&day) {
        return Err(TimeError::ParseError(format!("Day out of range: {}", day)));
    }

    let time_components: Vec<&str> = time_part.split(':').collect();
    if time_components.len() != 3 {
        return Err(TimeError::ParseError(format!(
            "Invalid time format: '{}'. Expected HH:MM:SS",
            time_part
        )));
    }

    let hour = match time_components[0].len() {
        1 => {
            let b = time_components[0].as_bytes()[0];
            if b.is_ascii_digit() {
                b - b'0'
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid hour: '{}'",
                    time_components[0]
                )));
            }
        }
        2 => {
            let bytes = time_components[0].as_bytes();
            if bytes.iter().all(|&b| b.is_ascii_digit()) {
                (bytes[0] - b'0') * 10 + (bytes[1] - b'0')
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid hour: '{}'",
                    time_components[0]
                )));
            }
        }
        _ => {
            return Err(TimeError::ParseError(format!(
                "Invalid hour format: '{}'",
                time_components[0]
            )))
        }
    };

    let minute = match time_components[1].len() {
        1 => {
            let b = time_components[1].as_bytes()[0];
            if b.is_ascii_digit() {
                b - b'0'
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid minute: '{}'",
                    time_components[1]
                )));
            }
        }
        2 => {
            let bytes = time_components[1].as_bytes();
            if bytes.iter().all(|&b| b.is_ascii_digit()) {
                (bytes[0] - b'0') * 10 + (bytes[1] - b'0')
            } else {
                return Err(TimeError::ParseError(format!(
                    "Invalid minute: '{}'",
                    time_components[1]
                )));
            }
        }
        _ => {
            return Err(TimeError::ParseError(format!(
                "Invalid minute format: '{}'",
                time_components[1]
            )))
        }
    };

    let second = time_components[2]
        .parse::<f64>()
        .map_err(|_| TimeError::ParseError(format!("Invalid second: '{}'", time_components[2])))?;

    if hour > 23 {
        return Err(TimeError::ParseError(format!(
            "Hour out of range: {}",
            hour
        )));
    }
    if minute > 59 {
        return Err(TimeError::ParseError(format!(
            "Minute out of range: {}",
            minute
        )));
    }
    if second >= 60.0 {
        return Err(TimeError::ParseError(format!(
            "Second out of range: {}",
            second
        )));
    }

    Ok(ParsedDateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iso8601() {
        let dt = parse_iso8601("2000-01-01T12:00:00").unwrap();
        assert_eq!(dt.year, 2000);
        assert_eq!(dt.month, 1);
        assert_eq!(dt.day, 1);
        assert_eq!(dt.hour, 12);
        assert_eq!(dt.minute, 0);
        assert_eq!(dt.second, 0.0);
    }

    #[test]
    fn test_iso8601_with_fractional_seconds() {
        let dt = parse_iso8601("2000-01-01T12:00:00.123").unwrap();
        assert_eq!(dt.second, 0.123);
    }

    #[test]
    fn test_iso8601_with_z_suffix() {
        let dt = parse_iso8601("2000-01-01T12:00:00Z").unwrap();
        assert_eq!(dt.year, 2000);
        assert_eq!(dt.hour, 12);
    }

    #[test]
    fn test_iso8601_space_separator() {
        let dt = parse_iso8601("2000-01-01 12:00:00").unwrap();
        assert_eq!(dt.year, 2000);
        assert_eq!(dt.hour, 12);
    }

    #[test]
    fn test_invalid_format() {
        assert!(parse_iso8601("not-a-date").is_err());
        assert!(parse_iso8601("2000-01-01").is_err());
        assert!(parse_iso8601("12:00:00").is_err());
    }

    #[test]
    fn test_invalid_ranges() {
        assert!(parse_iso8601("2000-13-01T12:00:00").is_err());
        assert!(parse_iso8601("2000-01-32T12:00:00").is_err());
        assert!(parse_iso8601("2000-01-01T25:00:00").is_err());
        assert!(parse_iso8601("2000-01-01T12:60:00").is_err());
        assert!(parse_iso8601("2000-01-01T12:00:60").is_err());
    }

    #[test]
    fn test_to_julian_date() {
        let dt = parse_iso8601("2000-01-01T12:00:00").unwrap();
        let jd = dt.to_julian_date();
        assert_eq!(jd.to_f64(), celestial_core::constants::J2000_JD);
    }

    #[test]
    fn test_input_too_long() {
        let long_input = "2000-01-01T12:00:00.".repeat(10);
        assert!(parse_iso8601(&long_input).is_err());
        if let Err(TimeError::ParseError(msg)) = parse_iso8601(&long_input) {
            assert_eq!(msg, "Input too long");
        } else {
            panic!("Expected ParseError with 'Input too long'");
        }
    }

    #[test]
    fn test_invalid_date_component_counts() {
        assert!(parse_iso8601("2000T12:00:00").is_err());
        assert!(parse_iso8601("2000-01T12:00:00").is_err());
        assert!(parse_iso8601("2000-01-01-01T12:00:00").is_err());
    }

    #[test]
    fn test_invalid_year_formats() {
        assert!(parse_iso8601("20a0-01-01T12:00:00").is_err());
        assert!(parse_iso8601("200-01-01T12:00:00").is_err());
        assert!(parse_iso8601("20000-01-01T12:00:00").is_err());
    }

    #[test]
    fn test_invalid_month_formats() {
        assert!(parse_iso8601("2000-a-01T12:00:00").is_err());
        assert!(parse_iso8601("2000-ab-01T12:00:00").is_err());
        assert!(parse_iso8601("2000-123-01T12:00:00").is_err());
    }

    #[test]
    fn test_invalid_day_formats() {
        assert!(parse_iso8601("2000-01-aT12:00:00").is_err());
        assert!(parse_iso8601("2000-01-abT12:00:00").is_err());
        assert!(parse_iso8601("2000-01-123T12:00:00").is_err());
    }

    #[test]
    fn test_invalid_time_component_counts() {
        assert!(parse_iso8601("2000-01-01T12").is_err());
        assert!(parse_iso8601("2000-01-01T12:00").is_err());
        assert!(parse_iso8601("2000-01-01T12:00:00:00").is_err());
    }

    #[test]
    fn test_invalid_hour_formats() {
        assert!(parse_iso8601("2000-01-01Ta:00:00").is_err());
        assert!(parse_iso8601("2000-01-01Tab:00:00").is_err());
        assert!(parse_iso8601("2000-01-01T123:00:00").is_err());
    }

    #[test]
    fn test_invalid_minute_formats() {
        assert!(parse_iso8601("2000-01-01T12:a:00").is_err());
        assert!(parse_iso8601("2000-01-01T12:ab:00").is_err());
        assert!(parse_iso8601("2000-01-01T12:123:00").is_err());
    }

    #[test]
    fn test_invalid_second_format() {
        assert!(parse_iso8601("2000-01-01T12:00:ab").is_err());
        assert!(parse_iso8601("2000-01-01T12:00:").is_err());
    }

    #[test]
    fn test_single_digit_components() {
        let dt = parse_iso8601("2000-1-1T1:1:1").unwrap();
        assert_eq!(dt.year, 2000);
        assert_eq!(dt.month, 1);
        assert_eq!(dt.day, 1);
        assert_eq!(dt.hour, 1);
        assert_eq!(dt.minute, 1);
        assert_eq!(dt.second, 1.0);
    }

    #[test]
    fn test_edge_case_ranges() {
        assert!(parse_iso8601("2000-00-01T12:00:00").is_err());
        assert!(parse_iso8601("2000-01-00T12:00:00").is_err());
        assert!(parse_iso8601("2000-12-31T23:59:59.999").is_ok());
    }

    #[test]
    fn test_whitespace_handling() {
        let dt = parse_iso8601("  2000-01-01T12:00:00  ").unwrap();
        assert_eq!(dt.year, 2000);
        assert_eq!(dt.hour, 12);
    }

    #[test]
    fn test_z_suffix_with_fractional_seconds() {
        let dt = parse_iso8601("2000-01-01T12:00:00.123Z").unwrap();
        assert_eq!(dt.second, 0.123);
    }
}
