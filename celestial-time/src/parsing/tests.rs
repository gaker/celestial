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
fn test_z_suffix_is_rejected() {
    for input in ["2000-01-01T12:00:00Z", " 2000-01-01T12:00:00.123Z "] {
        assert_eq!(
            parse_iso8601(input).unwrap_err(),
            TimeError::ParseError(format!("Only UTC takes a Z suffix: '{}'", input.trim()))
        );
    }
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
    assert_eq!(
        parse_iso8601("2000-01-01T12:00:61").unwrap_err(),
        TimeError::ParseError("Second out of range: 61".to_string())
    );
}

#[test]
fn test_leap_second_is_left_to_the_scale() {
    let dt = parse_iso8601("2016-12-31T23:59:60.5").unwrap();
    assert_eq!(dt.second, 60.5);
    assert_eq!(
        dt.to_julian_date(),
        Err(TimeError::InvalidDate(
            "time 23:59:60.5 is out of range".to_string()
        ))
    );
}

#[test]
fn test_to_julian_date() {
    let dt = parse_iso8601("2000-01-01T12:00:00").unwrap();
    let jd = dt.to_julian_date().unwrap();
    assert_eq!(jd.to_f64(), celestial_core::constants::J2000_JD);
}

#[test]
fn test_long_fractions_are_accepted() {
    for (input, expected) in [
        ("2000-01-01T12:00:00.123456789012", 0.123456789012),
        ("2000-01-01T12:00:59.99999999999999", 59.99999999999999),
        ("2000-01-01T12:00:00.000000000000000000001", 1e-21),
    ] {
        assert_eq!(parse_iso8601(input).unwrap().second, expected, "{}", input);
    }
    assert!(parse_iso8601(&"2000-01-01T12:00:00.".repeat(10)).is_err());
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
fn test_seconds_accept_only_digits_and_one_decimal_point() {
    for seconds in [
        "NaN", "nan", "inf", "-inf", "infinity", "-1", "+5", "-0", "1e1", "5.", ".5", "123",
        "5.5.5", "0x1", "5,5", "1_0", "5", "5.25", "059",
    ] {
        let input = format!("2000-01-01T12:00:{}", seconds);
        match parse_iso8601(&input) {
            Err(TimeError::ParseError(msg)) => {
                assert_eq!(msg, format!("Invalid second: '{}'", seconds))
            }
            other => panic!("{}: expected ParseError, got {:?}", input, other),
        }
    }
}

#[test]
fn test_seconds_values() {
    for (seconds, expected) in [
        ("05", 5.0),
        ("05.25", 5.25),
        ("59.999999999", 59.999999999),
        ("00.000000000001", 1e-12),
    ] {
        let input = format!("2000-01-01T12:00:{}", seconds);
        assert_eq!(parse_iso8601(&input).unwrap().second, expected, "{}", input);
    }
}

#[test]
fn test_fields_must_be_zero_padded() {
    for (input, message) in [
        ("2000-1-01T00:00:00", "Invalid month: '1'"),
        ("2000-01-1T00:00:00", "Invalid day: '1'"),
        ("2000-01-01T0:00:00", "Invalid hour: '0'"),
        ("2000-01-01T00:0:00", "Invalid minute: '0'"),
        ("2000-01-01T00:00:0", "Invalid second: '0'"),
        ("2000-1-1T0:0:0", "Invalid month: '1'"),
    ] {
        assert_eq!(
            parse_iso8601(input).unwrap_err(),
            TimeError::ParseError(message.to_string()),
            "{}",
            input
        );
    }
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
