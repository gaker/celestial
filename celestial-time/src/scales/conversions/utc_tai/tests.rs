use super::*;
use celestial_core::constants::J2000_JD;

#[test]
fn test_utc_tai_leap_second_offset() {
    assert_eq!(get_tai_utc_offset(1972, 1, 1, 0.0), Ok(10.0));
    assert_eq!(get_tai_utc_offset(1980, 1, 1, 0.0), Ok(19.0));
    assert_eq!(get_tai_utc_offset(1999, 1, 1, 0.0), Ok(32.0));
    assert_eq!(get_tai_utc_offset(2017, 1, 1, 0.0), Ok(37.0));
    assert_eq!(get_tai_utc_offset(1970, 6, 15, 0.5), Ok(8.429058000000001));
    // UTC begins in 1960.
    assert_eq!(get_tai_utc_offset(1959, 12, 31, 0.0), Ok(0.0));
}

#[test]
fn test_tai_utc_offset_rejects_invalid_input() {
    let cases = [
        (2000, 1, 1, 1.5),
        (2000, 1, 1, -0.5),
        (2000, 1, 1, f64::NAN),
        (1960, 0, 1, 0.0),
        (2000, 13, 1, 0.0),
        (2001, 2, 29, 0.5),
        (2000, 4, 31, 0.5),
        (-4800, 3, 1, 0.0),
    ];
    for (y, m, d, frac) in cases {
        assert!(
            get_tai_utc_offset(y, m, d, frac).is_err(),
            "{y}-{m:02}-{d:02} + {frac}"
        );
    }
}

#[test]
fn test_pre_1972_drift_matches_erfa_dat() {
    // One date in each drift period, from eraDat.
    let cases: &[(i32, i32, i32, f64, f64)] = &[
        (1960, 6, 15, 0.5, 1.1592660000000001),
        (1961, 3, 1, 0.25, 1.499606),
        (1961, 11, 30, 0.75, 1.805358),
        (1962, 7, 4, 0.5, 2.0530884),
        (1963, 12, 1, 0.5, 2.7315364),
        (1964, 2, 29, 0.5, 2.842906),
        (1964, 6, 1, 0.5, 3.063434),
        (1964, 10, 1, 0.5, 3.321546),
        (1965, 2, 1, 0.5, 3.580954),
        (1965, 5, 1, 0.5, 3.796298),
        (1965, 8, 1, 0.5, 4.01553),
        (1965, 10, 1, 0.5, 4.194586),
        (1966, 1, 1, 0.0, 4.31317),
        (1967, 6, 15, 0.5, 5.688226),
        (1968, 1, 31, 0.75, 6.2850340000000005),
        (1968, 2, 1, 0.0, 6.185682),
        (1971, 12, 31, 0.5, 9.890946),
    ];
    for &(y, m, d, frac, expected) in cases {
        assert_eq!(
            get_tai_utc_offset(y, m, d, frac),
            Ok(expected),
            "{y}-{m:02}-{d:02} + {frac}"
        );
    }
}

#[test]
fn test_pre_1972_utc_tai_matches_erfa() {
    use celestial_core::constants::MJD_ZERO_POINT;
    // From eraUtctai. 1965-12-31 rescales its day length using the 1966 rate.
    let cases: &[(f64, f64)] = &[
        (39125.5, 39125.50004991345),
        (39125.999, 39125.99904992094),
        (39126.0, 39126.00004992095),
        (39500.25, 39500.25006114845),
        (39886.75, 39886.750071875394),
    ];
    for &(utc_mjd, tai_mjd) in cases {
        let tai = utc_to_tai(JulianDate::new(MJD_ZERO_POINT, utc_mjd)).unwrap();
        assert_eq!(tai.to_julian_date().jd1(), MJD_ZERO_POINT);
        assert_eq!(tai.to_julian_date().jd2(), tai_mjd, "UTC MJD {utc_mjd}");

        let utc = tai_to_utc(tai.to_julian_date()).unwrap();
        assert_eq!(utc.to_julian_date().jd1(), MJD_ZERO_POINT);
        assert_eq!(utc.to_julian_date().jd2(), utc_mjd, "TAI MJD {tai_mjd}");
    }
}

#[test]
fn test_leap_second_day_utc_tai_matches_erfa() {
    // UTC from eraDtf2d: 1972-06-30T23:59:60.5; 2016-12-31T12:00:00, 23:59:60,
    // 23:59:60.5 and 23:59:60.999; 2017-01-01T00:00:00.25. TAI from eraUtctai.
    // The way back is eraTaiutc, whose iteration can end a bit off the input.
    let cases = [
        (
            2441498.5,
            0.9999942130299417,
            1.0001215277777777,
            0.9999942130299417,
        ),
        (
            2457753.5,
            0.4999942130299418,
            0.5004166666666666,
            0.49999421302994185,
        ),
        (
            2457753.5,
            0.9999884260598836,
            1.0004166666666667,
            0.9999884260598837,
        ),
        (
            2457753.5,
            0.9999942130299417,
            1.0004224537037036,
            0.9999942130299417,
        ),
        (
            2457753.5,
            0.9999999884260599,
            1.0004282291666666,
            0.9999999884260597,
        ),
        (
            2457754.5,
            2.8935185185185184e-6,
            0.00043113425925925925,
            2.893518518518501e-6,
        ),
    ];
    for (jd1, utc_jd2, tai_jd2, back_jd2) in cases {
        let tai = utc_to_tai(JulianDate::new(jd1, utc_jd2)).unwrap();
        assert_eq!(tai.to_julian_date().jd1(), jd1);
        assert_eq!(
            tai.to_julian_date().jd2(),
            tai_jd2,
            "UTC ({jd1}, {utc_jd2})"
        );

        let utc = tai_to_utc(tai.to_julian_date()).unwrap();
        assert_eq!(utc.to_julian_date().jd1(), jd1);
        assert_eq!(
            utc.to_julian_date().jd2(),
            back_jd2,
            "TAI ({jd1}, {tai_jd2})"
        );
    }
}

#[test]
fn test_utc_tai_round_trip_matches_erfa() {
    // eraUtctai then eraTaiutc, and the reverse. The last two cases take the
    // alternate split in utc_to_tai and tai_to_utc.
    let jd = JulianDate::new;
    let cases = [
        (
            jd(J2000_JD, 0.123456789),
            jd(J2000_JD, 0.12345678900000001),
            jd(J2000_JD, 0.12345678900000001),
        ),
        (jd(J2000_JD, 0.0), jd(J2000_JD, 0.0), jd(J2000_JD, 0.0)),
        (
            jd(J2000_JD, 0.999999),
            jd(J2000_JD, 0.9999990000000001),
            jd(J2000_JD, 0.999999),
        ),
        (jd(0.5, J2000_JD), jd(0.5, J2000_JD), jd(0.5, J2000_JD)),
        (
            jd(0.1, J2000_JD),
            jd(0.09999999999999998, J2000_JD),
            jd(0.09999999999999998, J2000_JD),
        ),
    ];
    for (input, utc_back, tai_back) in cases {
        let back = UTC::from_julian_date(input).to_tai().unwrap();
        let utc = back.to_utc().unwrap();
        assert_eq!(utc.to_julian_date().parts(), utc_back.parts(), "{input}");
        let back = TAI::from_julian_date(input).to_utc().unwrap();
        let tai = back.to_tai().unwrap();
        assert_eq!(tai.to_julian_date().parts(), tai_back.parts(), "{input}");
    }
}

#[test]
fn test_calendar_helper_functions() {
    assert_eq!(next_calendar_day(2000, 1, 31).unwrap(), (2000, 2, 1));
    assert_eq!(next_calendar_day(2000, 12, 31).unwrap(), (2001, 1, 1));
    assert_eq!(next_calendar_day(2000, 2, 29).unwrap(), (2000, 3, 1));

    assert_eq!(next_calendar_day(2000, 2, 28).unwrap(), (2000, 2, 29));
    assert_eq!(next_calendar_day(1900, 2, 28).unwrap(), (1900, 3, 1));

    assert_eq!(next_calendar_day(2000, 4, 30).unwrap(), (2000, 5, 1));
    assert_eq!(next_calendar_day(2000, 6, 30).unwrap(), (2000, 7, 1));
    assert_eq!(next_calendar_day(2000, 9, 30).unwrap(), (2000, 10, 1));
    assert_eq!(next_calendar_day(2000, 11, 30).unwrap(), (2000, 12, 1));

    assert!(next_calendar_day(2000, 0, 1).is_err());
    assert!(next_calendar_day(2000, 13, 1).is_err());
    assert!(next_calendar_day(2000, -1, 1).is_err());

    assert!(julian_to_calendar(1e10, 0.0).is_err());
    assert!(julian_to_calendar(-1e6, 0.0).is_err());

    // eraJd2cal, through the negative fraction correction, the Kahan summation
    // branch where |frac_2| > |sum|, and the near-1.0 fraction correction
    assert_eq!(
        julian_to_calendar(celestial_core::constants::J2000_JD, -0.6),
        Ok((1999, 12, 31, 0.9))
    );
    assert_eq!(
        julian_to_calendar(2451544.6, 0.2),
        Ok((2000, 1, 1, 0.30000000009313227))
    );
    assert_eq!(julian_to_calendar(2451544.75, 0.75), Ok((2000, 1, 2, 0.0)));
}
