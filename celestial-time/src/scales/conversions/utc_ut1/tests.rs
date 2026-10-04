use super::*;
use crate::TimeError;
use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT};

#[test]
fn test_dut1_offset_matches_erfa() {
    // eraUtcut1 and eraUt1utc at J2000.0. UTC to UT1 goes through TAI and picks
    // up its rounding; UT1 to UTC subtracts DUT1 directly.
    let cases = [
        (-0.9, -1.0416666666682758e-5, 1.0416666666666666e-5),
        (0.0, -1.610040226140974e-17, 0.0),
        (0.9, 1.0416666666650557e-5, -1.0416666666666666e-5),
    ];
    for (dut1, ut1_jd2, utc_jd2) in cases {
        let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let ut1 = utc.to_ut1_with_dut1(dut1).unwrap();
        assert_eq!(ut1.to_julian_date().parts(), (J2000_JD, ut1_jd2), "{dut1}");

        let ut1 = UT1::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let utc = ut1.to_utc_with_dut1(dut1).unwrap();
        assert_eq!(utc.to_julian_date().parts(), (J2000_JD, utc_jd2), "{dut1}");
    }

    let jd = JulianDate::new;
    for (input, expected) in [
        (jd(J2000_JD, 0.5), jd(J2000_JD, 0.49999652777777776)),
        (jd(0.1, J2000_JD), jd(0.09999652777777779, J2000_JD)),
    ] {
        let utc = UT1::from_julian_date(input).to_utc_with_dut1(0.3).unwrap();
        assert_eq!(utc.to_julian_date().parts(), expected.parts(), "{input}");
    }
}

#[test]
fn test_utc_ut1_round_trip_matches_erfa() {
    // eraUtcut1 then eraUt1utc, and the reverse, as (input, DUT1, UTC back,
    // UT1 back). Either way can come back an ulp or two off.
    let jd = JulianDate::new;
    let cases = [
        (
            jd(J2000_JD, 0.123456789),
            -0.9,
            jd(J2000_JD, 0.123456789),
            jd(J2000_JD, 0.12345678900000001),
        ),
        (jd(J2000_JD, 0.5), 0.0, jd(J2000_JD, 0.5), jd(J2000_JD, 0.5)),
        (
            jd(J2000_JD, 0.5),
            0.9,
            jd(J2000_JD, 0.5),
            jd(J2000_JD, 0.5000000000000001),
        ),
        (
            jd(J2000_JD, 0.999999999),
            -0.9,
            jd(J2000_JD, 0.9999999990000001),
            jd(J2000_JD, 0.999999999),
        ),
        (
            jd(J2000_JD, 0.25),
            0.9,
            jd(J2000_JD, 0.24999999999999994),
            jd(J2000_JD, 0.24999999999999997),
        ),
        (
            jd(0.1, J2000_JD),
            0.0,
            jd(0.09999999999999996, J2000_JD),
            jd(0.09999999999999996, J2000_JD),
        ),
    ];
    for (input, dut1, utc_back, ut1_back) in cases {
        let utc = UTC::from_julian_date(input);
        let back = utc.to_ut1_with_dut1(dut1).unwrap();
        let back = back.to_utc_with_dut1(dut1).unwrap();
        assert_eq!(
            back.to_julian_date().parts(),
            utc_back.parts(),
            "{input}, {dut1}"
        );

        let ut1 = UT1::from_julian_date(input);
        let back = ut1.to_utc_with_dut1(dut1).unwrap();
        let back = back.to_ut1_with_dut1(dut1).unwrap();
        assert_eq!(
            back.to_julian_date().parts(),
            ut1_back.parts(),
            "{input}, {dut1}"
        );
    }
}

#[test]
fn test_zero_dut1_at_leap_second_is_identity() {
    // eraUt1utc takes the leap second off and adds all of it back here.
    for jd1 in [2441499.5, 2441683.5, 2442048.5] {
        for jd2 in [0.0, 0.001] {
            let input = JulianDate::new(jd1, jd2);
            let utc = UT1::from_julian_date(input).to_utc_with_dut1(0.0).unwrap();
            assert_eq!(utc.to_julian_date().parts(), input.parts(), "{input}");
        }
    }
}

#[test]
fn test_leap_second_day_utc_ut1_matches_erfa() {
    // The UTC inputs of the utc_tai leap-second-day test. UT1 from eraUtcut1 and
    // back from eraUt1utc. DUT1 is held at +0.1 s across the 1972 leap second,
    // so 23:59:60.5 and the next day's 00:00:00.5 share a UT1 and the way back
    // gives the latter. Elsewhere eraUt1utc's iteration can end a bit off.
    let cases = [
        (
            2441498.5,
            0.9999942130299417,
            0.1,
            1.0000069444444444,
            1.000005787037037,
        ),
        (
            2457753.5,
            0.4999942130299418,
            -0.5888,
            0.49999318518518515,
            0.49999421302994174,
        ),
        (
            2457753.5,
            0.9999884260598836,
            -0.5888,
            0.9999931851851852,
            0.9999884260598836,
        ),
        (
            2457753.5,
            0.9999942130299417,
            -0.5888,
            0.9999989722222221,
            0.9999942130299417,
        ),
        (
            2457753.5,
            0.9999999884260599,
            -0.5888,
            1.0000047476851852,
            0.9999999884260599,
        ),
        (
            2457754.5,
            2.8935185185185184e-6,
            0.4112,
            7.652777777777757e-6,
            2.8935185185184985e-6,
        ),
    ];
    for (jd1, utc_jd2, dut1, ut1_jd2, back_jd2) in cases {
        let utc = UTC::from_julian_date(JulianDate::new(jd1, utc_jd2));
        let ut1 = utc.to_ut1_with_dut1(dut1).unwrap();
        assert_eq!(ut1.to_julian_date().jd1(), jd1);
        assert_eq!(
            ut1.to_julian_date().jd2(),
            ut1_jd2,
            "UTC ({jd1}, {utc_jd2})"
        );

        let back = ut1.to_utc_with_dut1(dut1).unwrap();
        assert_eq!(back.to_julian_date().jd1(), jd1);
        assert_eq!(
            back.to_julian_date().jd2(),
            back_jd2,
            "UT1 ({jd1}, {ut1_jd2})"
        );
    }
}

#[test]
fn test_pre_1972_utc_ut1_matches_erfa() {
    // eraUtcut1, then eraUt1utc on its result, with DUT1 = 0.123456789 s.
    // TAI-UTC drifts through these days and 1971-12-31 ends in a 0.108 s
    // step, so the round trip doesn't come back.
    let cases = [
        (36934.25, 36934.25000143264, 36934.25000000375),
        (39130.999999, 39131.0000004589, 39130.999999030006),
        (41316.75, 41316.7500023868, 41316.7500009579),
    ];
    for (utc_mjd, ut1_mjd, back_mjd) in cases {
        let utc = UTC::from_julian_date(JulianDate::new(MJD_ZERO_POINT, utc_mjd));
        let ut1 = utc.to_ut1_with_dut1(0.123456789).unwrap();
        assert_eq!(ut1.to_julian_date().parts(), (MJD_ZERO_POINT, ut1_mjd));

        let back = ut1.to_utc_with_dut1(0.123456789).unwrap();
        assert_eq!(back.to_julian_date().parts(), (MJD_ZERO_POINT, back_mjd));
    }
}

#[test]
fn test_non_finite_julian_date_is_rejected() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let jd = JulianDate::new(bad, 0.0);
        let expected = TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
        let ut1 = UT1::from_julian_date(jd);
        assert_eq!(ut1.to_utc_with_dut1(0.25).unwrap_err(), expected);
    }
}

#[test]
fn test_non_finite_dut1_is_rejected() {
    let jd = JulianDate::new(J2000_JD, 0.0);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let expected =
            TimeError::ConversionError(format!("dut1_seconds must be finite, got {}", bad));
        let utc = UTC::from_julian_date(jd);
        assert_eq!(utc.to_ut1_with_dut1(bad).unwrap_err(), expected);
        let ut1 = UT1::from_julian_date(jd);
        assert_eq!(ut1.to_utc_with_dut1(bad).unwrap_err(), expected);
    }
}
