use super::*;

#[test]
fn test_utc_string_parsing() {
    assert_eq!(
        UTC::from_str("2000-01-01T12:00:00")
            .unwrap()
            .to_julian_date()
            .parts(),
        (2451544.5, 0.5)
    );
    assert_eq!(
        UTC::from_str("2000-01-01T12:00:00.123"),
        utc_from_calendar(2000, 1, 1, 12, 0, 0.123)
    );
    assert_eq!(
        UTC::from_str("invalid-date"),
        Err(TimeError::ParseError(
            "Invalid datetime format: 'invalid-date'. Expected YYYY-MM-DDTHH:MM:SS".into()
        ))
    );
}

#[test]
fn test_new_matches_erfa_dtf2d() {
    // eraDtf2d on the Unix time's date and time of day.
    let cases = [
        (
            (1_576_800_000, 123_456_789),
            (2458837.5, 1.4288980208333333e-6),
        ),
        ((0, 0), (2440587.5, 0.0)),
        ((1_700_000_000, 1), (2460262.5, 0.9259259259259376)),
        // 2016-12-31 is 86,401 s long.
        ((1_483_185_600, 0), (2457753.5, 0.4999942130299418)),
        (
            (1_483_228_799, 500_000_000),
            (2457753.5, 0.9999826390898253),
        ),
        ((1_483_228_800, 0), (2457754.5, 0.0)),
        // The pre-1972 steps: +0.107758 s on 1971-12-31, -0.1 s on 1968-01-31.
        ((63_071_999, 250_000_000), (2441316.5, 0.9999900722577523)),
        ((-60_523_200, 0), (2439886.5, 0.5000005787043735)),
        ((-2, 500_000_000), (2440586.5, 0.9999826388888889)),
    ];
    for ((seconds, nanos), (jd1, jd2)) in cases {
        assert_eq!(
            UTC::new(seconds, nanos).unwrap().to_julian_date().parts(),
            (jd1, jd2),
            "{seconds} s {nanos} ns"
        );
    }
}

#[test]
fn test_new_rejects_out_of_range_input() {
    assert_eq!(
        UTC::new(0, 1_000_000_000),
        Err(TimeError::ConversionError(
            "nanos must be less than 1000000000, got 1000000000".into()
        ))
    );
    assert_eq!(
        UTC::new(i64::MAX, 0),
        Err(TimeError::ConversionError(
            "Julian Date 106751993607887.5 out of valid range [-68569.5, 1000000000]".into()
        ))
    );
}

#[test]
fn test_tai_utc_offset_edge_cases() {
    assert_eq!(get_tai_utc_offset(1950, 6, 15, 0.5), Ok(0.0));
    assert_eq!(get_tai_utc_offset(1960, 1, 1, 0.0), Ok(0.9434819999999999));
}

#[test]
fn test_utc_from_calendar_matches_erfa_dtf2d() {
    let cases = [
        ((2016, 12, 31, 23, 59, 60.0), 0.9999884260598836),
        ((2016, 12, 31, 23, 59, 60.5), 0.9999942130299417),
        ((2016, 12, 31, 23, 59, 60.999), 0.9999999884260599),
        ((2016, 12, 31, 12, 0, 0.0), 0.4999942130299418),
        ((1961, 7, 31, 23, 59, 59.94), 0.9999998842591924),
        ((2000, 1, 1, 0, 0, -0.0), 0.0),
    ];
    for ((y, mo, d, h, mi, s), expected_jd2) in cases {
        let jd = utc_from_calendar(y, mo, d, h, mi, s)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd.jd2(), expected_jd2, "{y}-{mo:02}-{d:02} {h}:{mi}:{s}");
    }
}

#[test]
fn test_from_str_takes_an_optional_z_suffix() {
    for input in ["2000-01-01T12:00:00.123Z", "  2000-01-01T12:00:00.123Z  "] {
        assert_eq!(
            UTC::from_str(input).unwrap(),
            UTC::from_str("2000-01-01T12:00:00.123").unwrap()
        );
    }
    assert!(UTC::from_str("2000-01-01T12:00:00ZZ").is_err());
}

#[test]
fn test_from_system_time_before_and_after_1970() {
    use std::time::Duration;
    let at = |offset_ms: i64| {
        let step = Duration::from_millis(offset_ms.unsigned_abs());
        let time = if offset_ms < 0 {
            UNIX_EPOCH.checked_sub(step).unwrap()
        } else {
            UNIX_EPOCH.checked_add(step).unwrap()
        };
        UTC::from_system_time(time).unwrap()
    };
    assert_eq!(at(1500), UTC::new(1, 500_000_000).unwrap());
    assert_eq!(at(0), UTC::new(0, 0).unwrap());
    assert_eq!(at(-3000), UTC::new(-3, 0).unwrap());
    assert_eq!(at(-1500), UTC::new(-2, 500_000_000).unwrap());
    assert_eq!(
        at(-1500).to_julian_date().parts(),
        (UNIX_EPOCH_JD - 1.0, 0.9999826388888889)
    );
}

#[test]
fn test_from_str_matches_erfa_dtf2d_on_leap_second_days() {
    let cases = [
        ("2016-12-31T12:00:00", (2457753.5, 0.4999942130299418)),
        ("2016-12-31T23:59:59.5", (2457753.5, 0.9999826390898253)),
        ("2016-12-31T23:59:60", (2457753.5, 0.9999884260598836)),
        ("2016-12-31T23:59:60.5Z", (2457753.5, 0.9999942130299417)),
        ("1972-06-30T23:59:60.999", (2441498.5, 0.9999999884260599)),
        ("1961-07-31T23:59:59.94", (2437511.5, 0.9999998842591924)),
        ("1963-10-31T23:59:60.05", (2438333.5, 0.9999994212969661)),
        ("2016-12-30T23:59:59.5", (2457752.5, 0.999994212962963)),
    ];
    for (text, (jd1, jd2)) in cases {
        let jd = UTC::from_str(text).unwrap().to_julian_date();
        assert_eq!((jd.jd1(), jd.jd2()), (jd1, jd2), "{text}");
    }
}

#[test]
fn test_from_str_rejects_seconds_the_day_does_not_have() {
    let out_of_range = |t: &str| TimeError::InvalidDate(format!("time {t} is out of range"));
    let cases = [
        ("2016-12-30T23:59:60", out_of_range("23:59:60")),
        ("2016-12-31T12:00:60", out_of_range("12:00:60")),
        ("1961-07-31T23:59:59.96", out_of_range("23:59:59.96")),
        (
            "2016-12-31T23:59:61",
            TimeError::ParseError("Second out of range: 61".to_string()),
        ),
    ];
    for (text, expected) in cases {
        assert_eq!(UTC::from_str(text).unwrap_err(), expected, "{text}");
    }
}

#[test]
fn test_to_iso8601_matches_erfa_d2dtf() {
    let cases = [
        ((2457753.5, 0.000694430620015972), "2016-12-31T00:01:00.000"),
        ((2457753.9999942062, 0.0), "2016-12-31T11:59:59.999"),
        ((2457753.5, 0.4999942072429717), "2016-12-31T12:00:00.000"),
        ((2457753.5, 0.9999884202729136), "2016-12-31T23:59:60.000"),
        ((2457753.5, 0.999999994), "2016-12-31T23:59:60.999"),
        ((2457753.5, 0.999999995), "2017-01-01T00:00:00.000"),
        ((2457752.5, 0.999999994), "2016-12-30T23:59:59.999"),
        ((2457752.5, 0.999999995), "2016-12-31T00:00:00.000"),
        ((2437511.5, 0.5000002893520193), "1961-07-31T12:00:00.025"),
        ((2438333.5, 0.9999988425939321), "1963-10-31T23:59:59.900"),
    ];
    for ((jd1, jd2), expected) in cases {
        let utc = UTC::from_julian_date(JulianDate::new(jd1, jd2));
        assert_eq!(utc.to_iso8601().unwrap(), expected, "({jd1}, {jd2})");
    }
}

#[test]
fn test_to_iso8601_round_trips_a_leap_second() {
    let text = "2016-12-31T23:59:60.500";
    assert_eq!(UTC::from_str(text).unwrap().to_iso8601().unwrap(), text);
}

#[test]
fn test_to_iso8601_rejects_dates_outside_the_calendar() {
    let before = UTC::from_julian_date(JulianDate::new(-68569.5, 0.0));
    assert_eq!(
        before.to_iso8601(),
        Err(TimeError::InvalidDate(
            "-4900-03-01: year is before -4799".to_string()
        ))
    );
    let nan = UTC::from_julian_date(JulianDate::new(f64::NAN, 0.0));
    assert_eq!(
        nan.to_iso8601(),
        Err(TimeError::ConversionError(
            "Julian Date NaN out of valid range [-68569.5, 1000000000]".to_string()
        ))
    );
}

#[test]
fn test_utc_from_calendar_rejects_invalid_input() {
    let cases = [
        (2000, 13, 1, 0, 0, 0.0),
        (2000, 0, 1, 0, 0, 0.0),
        (2001, 2, 29, 0, 0, 0.0),
        (2000, 4, 31, 0, 0, 0.0),
        (-4800, 3, 1, 0, 0, 0.0),
        (2000, 1, 1, 24, 0, 0.0),
        (2000, 1, 1, 0, 60, 0.0),
        (2000, 1, 1, 0, 0, -1e-300),
        (2000, 1, 1, 0, 0, f64::NAN),
        (2000, 1, 1, 0, 0, f64::INFINITY),
        (2016, 12, 31, 23, 59, 61.0),
        (2016, 12, 31, 12, 0, 60.0),
        (2016, 12, 30, 23, 59, 60.0),
        // 1961-07-31 was 0.05 s short.
        (1961, 7, 31, 23, 59, 59.96),
    ];
    for (y, mo, d, h, mi, s) in cases {
        assert!(
            utc_from_calendar(y, mo, d, h, mi, s).is_err(),
            "{y}-{mo:02}-{d:02} {h}:{mi}:{s}"
        );
    }
}

#[test]
fn test_next_calendar_day() {
    assert!(next_calendar_day(2000, 13, 15).is_err());
    assert!(next_calendar_day(2001, 2, 29).is_err());
    assert!(next_calendar_day(i32::MAX, 12, 31).is_err());

    let cases = [
        (2000, 2, 28, (2000, 2, 29)),
        (1999, 2, 28, (1999, 3, 1)),
        (2000, 4, 30, (2000, 5, 1)),
        (2000, 12, 31, (2001, 1, 1)),
    ];

    for (y, m, d, expected) in cases {
        assert_eq!(next_calendar_day(y, m, d).unwrap(), expected);
    }
}
