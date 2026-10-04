use super::*;

#[test]
fn test_julian_date_creation() {
    let jd = JulianDate::new(J2000_JD, 0.5);
    assert_eq!(jd.jd1(), J2000_JD);
    assert_eq!(jd.jd2(), 0.5);
    assert_eq!(jd.to_f64(), 2451545.5);
}

#[test]
fn test_j2000_epoch() {
    let j2000 = JulianDate::j2000();
    assert_eq!(j2000.to_f64(), J2000_JD);
}

#[test]
fn test_unix_epoch() {
    let unix = JulianDate::unix_epoch();
    assert_eq!(unix.to_f64(), crate::constants::UNIX_EPOCH_JD);
}

#[test]
fn test_arithmetic() {
    let jd = JulianDate::new(J2000_JD, 0.0);
    assert_eq!(jd.add_days(1.0).parts(), (2451546.0, 0.0));
    assert_eq!(
        jd.add_seconds(3600.0).parts(),
        (J2000_JD, 0.041666666666666664)
    );
}

#[test]
fn test_add_days_carries_whole_days_into_jd1() {
    let jd = JulianDate::new(J2000_JD, 0.75);
    assert_eq!(jd.add_days(0.5).parts(), (2451546.0, 0.25));
    assert_eq!(jd.add_days(-1.5).parts(), (2451544.0, 0.25));

    let century_later = JulianDate::new(J2000_JD, 0.25).add_days(36525.0);
    assert_eq!(century_later.parts(), (2488070.0, 0.25));
    assert_eq!(
        century_later.add_seconds(1e-7).parts(),
        (2488070.0, 0.2500000000011574)
    );
}

#[test]
fn test_add_seconds_divides_by_the_day_length() {
    // 5 * (1 / 86400) rounds one ULP below 5 / 86400.
    let jd = JulianDate::new(J2000_JD, 0.0).add_seconds(5.0);
    assert_eq!(jd.parts(), (J2000_JD, 5.787037037037037e-5));
}

#[test]
fn test_from_unix_time_keeps_nanoseconds() {
    assert_eq!(
        JulianDate::from_unix_time(1_700_000_000, 0).parts(),
        (2460262.5, 0.9259259259259259)
    );
    assert_eq!(
        JulianDate::from_unix_time(1_700_000_000, 1).parts(),
        (2460262.5, 0.9259259259259376)
    );
    assert_eq!(
        JulianDate::from_unix_time(-1, 500_000_000).parts(),
        (2440587.5, -5.787037037037037e-6)
    );
}

#[test]
fn test_equality_ignores_the_split() {
    let jd = JulianDate::new;
    assert_eq!(jd(J2000_JD, 0.5), jd(2451545.5, 0.0));
    assert_eq!(jd(2400000.5, 51544.75), jd(J2000_JD, 0.25));
    assert_eq!(jd(1e-17, J2000_JD), jd(J2000_JD, 1e-17));
    let tt = crate::scales::tt::TT::from_julian_date;
    assert_eq!(tt(jd(J2000_JD, 0.5)), tt(jd(2451545.5, 0.0)));
}

#[test]
fn test_equality_is_exact_below_an_ulp_of_the_sum() {
    let jd = JulianDate::new;
    // Both sums round to 2451545.0, but the dates are 1e-17 days apart.
    assert_ne!(jd(J2000_JD, 1e-17), jd(J2000_JD, 0.0));
    assert_ne!(jd(J2000_JD, f64::NAN), jd(J2000_JD, f64::NAN));
    assert_eq!(jd(f64::INFINITY, 0.0), jd(f64::INFINITY, 0.0));
    assert_eq!(jd(f64::MAX, f64::MAX), jd(f64::MAX, f64::MAX));
    assert_ne!(jd(f64::MAX, f64::MAX), jd(f64::INFINITY, 0.0));
}

#[test]
fn test_display_shows_both_parts_exactly() {
    for (jd, shown) in [
        (JulianDate::new(J2000_JD, 0.5), "JD 2451545 + 0.5"),
        (JulianDate::new(J2000_JD, 1e-9), "JD 2451545 + 0.000000001"),
        (
            JulianDate::new(2400000.5, -0.123456789012345),
            "JD 2400000.5 + -0.123456789012345",
        ),
    ] {
        assert_eq!(jd.to_string(), shown);
    }
}

#[test]
fn test_julian_year_matches_erfa_epj() {
    let jd = JulianDate::new(2400000.5, 60000.12345678901);
    assert_eq!(jd.to_julian_year(), 2023.1502353368624);
    let jd = JulianDate::new(J2000_JD, -36524.87654321988);
    assert_eq!(jd.to_julian_year(), 1900.0003380062426);
}

#[test]
fn test_from_julian_year_matches_erfa_epj2jd() {
    assert_eq!(
        JulianDate::from_julian_year(2024.123456789).parts(),
        (2400000.5, 60355.592592182285)
    );
    assert_eq!(
        JulianDate::from_julian_year(2000.0).parts(),
        (2400000.5, 51544.5)
    );
}

#[cfg(feature = "serde")]
#[test]
fn test_serde_round_trip() {
    let test_cases = [
        JulianDate::new(J2000_JD, 0.0),
        JulianDate::new(2451545.5, 0.123456789),
        JulianDate::new(2440587.5, 0.0), // Unix epoch
        JulianDate::new(J2000_JD, 0.999999999),
    ];

    for original in test_cases {
        let json = serde_json::to_string(&original).unwrap();
        let deserialized: JulianDate = serde_json::from_str(&json).unwrap();

        assert_eq!(
            original.jd1(),
            deserialized.jd1(),
            "JD1 precision lost in serde round-trip"
        );
        assert_eq!(
            original.jd2(),
            deserialized.jd2(),
            "JD2 precision lost in serde round-trip"
        );
        assert_eq!(
            original, deserialized,
            "JulianDate equality lost in serde round-trip"
        );
    }
}

#[test]
fn test_finite_jd_rejects_either_part_non_finite() {
    for (bad, shown) in [
        (f64::NAN, "NaN"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
    ] {
        let message = |jd1: &str, jd2: &str| {
            TimeError::InvalidEpoch(format!("Julian Date ({}, {}) is not finite", jd1, jd2))
        };
        let jd = JulianDate::new(bad, 0.5);
        assert_eq!(finite_jd(jd), Err(message(shown, "0.5")));
        let jd = JulianDate::new(J2000_JD, bad);
        assert_eq!(finite_jd(jd), Err(message("2451545", shown)));
    }
    let jd = JulianDate::new(J2000_JD, -0.25);
    assert_eq!(finite_jd(jd), Ok(jd));
}

#[test]
fn test_finite_arg_names_the_argument() {
    for (bad, shown) in [
        (f64::NAN, "NaN"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
    ] {
        let expected = format!("dut1_seconds must be finite, got {}", shown);
        assert_eq!(
            finite_arg("dut1_seconds", bad),
            Err(TimeError::ConversionError(expected))
        );
    }
    assert_eq!(finite_arg("dut1_seconds", -0.25), Ok(-0.25));
}

#[test]
fn test_sub_is_the_difference_in_days() {
    let later = JulianDate::new(J2000_JD, 1e-9);
    let earlier = JulianDate::new(J2000_JD, 0.0);
    assert_eq!(later - earlier, 1e-9);
    assert_eq!(
        JulianDate::new(J2000_JD, 0.25) - JulianDate::new(2451544.5, 0.5),
        0.25
    );
}
