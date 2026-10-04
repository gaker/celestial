use super::*;
use crate::julian::JulianDate;
use crate::TimeError;
use celestial_core::constants::J2000_JD;

#[test]
fn test_ut1_tai_offset_applied_correctly() {
    let test_dates = [
        (J2000_JD, "J2000.0"),
        (2455197.5, "2010-01-01"),
        (2459580.5, "2022-01-01"),
    ];
    let ut1_tai_offset = -32.3;
    let delta_t = 69.0;

    for (jd, description) in test_dates {
        // UT1 -> TAI: TAI = UT1 - (UT1-TAI), so TAI should be ahead by 32.3s
        let ut1 = UT1::from_julian_date(JulianDate::new(jd, 0.0));
        let tai = ut1.to_tai_with_offset(ut1_tai_offset).unwrap();

        let ut1_jd = ut1.to_julian_date();
        let tai_jd = tai.to_julian_date();

        let offset_days = (tai_jd.jd1() - ut1_jd.jd1()) + (tai_jd.jd2() - ut1_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        assert_eq!(
            offset_seconds, -ut1_tai_offset,
            "{}: UT1->TAI offset must be exactly {} seconds",
            description, -ut1_tai_offset
        );

        // TAI -> UT1: UT1 = TAI + (UT1-TAI), so UT1 should be behind by 32.3s
        let tai = TAI::from_julian_date(JulianDate::new(jd, 0.0));
        let ut1 = tai.to_ut1_with_offset(ut1_tai_offset).unwrap();

        let tai_jd = tai.to_julian_date();
        let ut1_jd = ut1.to_julian_date();

        let offset_days = (tai_jd.jd1() - ut1_jd.jd1()) + (tai_jd.jd2() - ut1_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        assert_eq!(
            offset_seconds, -ut1_tai_offset,
            "{}: TAI->UT1 means TAI is {} seconds ahead",
            description, -ut1_tai_offset
        );

        // UT1 -> TT: TT = UT1 + Delta-T, so TT should be ahead by 69s
        let ut1 = UT1::from_julian_date(JulianDate::new(jd, 0.0));
        let tt = ut1.to_tt_with_delta_t(delta_t).unwrap();

        let ut1_jd = ut1.to_julian_date();
        let tt_jd = tt.to_julian_date();

        let offset_days = (tt_jd.jd1() - ut1_jd.jd1()) + (tt_jd.jd2() - ut1_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        assert_eq!(
            offset_seconds, delta_t,
            "{}: UT1->TT offset must be exactly {} seconds",
            description, delta_t
        );

        // TT -> UT1: UT1 = TT - Delta-T, so UT1 should be behind by 69s
        let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
        let ut1 = tt.to_ut1_with_delta_t(delta_t).unwrap();

        let tt_jd = tt.to_julian_date();
        let ut1_jd = ut1.to_julian_date();

        let offset_days = (tt_jd.jd1() - ut1_jd.jd1()) + (tt_jd.jd2() - ut1_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        assert_eq!(
            offset_seconds, delta_t,
            "{}: TT->UT1 means TT is {} seconds ahead",
            description, delta_t
        );
    }
}

#[test]
fn test_ut1_tai_round_trip_matches_erfa() {
    // eraUt1tai then eraTaiut1, and the reverse. TAI always comes back exactly;
    // UT1 with a 0.5 day part comes back 1 ulp low for some offsets.
    let jd = JulianDate::new;
    let cases = [
        (jd(J2000_JD, 0.0), -31.8, jd(J2000_JD, 0.0)),
        (jd(J2000_JD, 0.5), -32.0, jd(J2000_JD, 0.5)),
        (jd(J2000_JD, 0.5), -31.8, jd(J2000_JD, 0.49999999999999994)),
        (
            jd(J2000_JD, 0.123456789012345),
            -31.8,
            jd(J2000_JD, 0.123456789012345),
        ),
        (
            jd(J2000_JD, -0.123456789012345),
            -31.8,
            jd(J2000_JD, -0.123456789012345),
        ),
        (jd(J2000_JD, 0.987654321), -31.8, jd(J2000_JD, 0.987654321)),
        (jd(0.5, J2000_JD), -31.8, jd(0.49999999999999994, J2000_JD)),
    ];
    for (input, offset, ut1_back) in cases {
        let ut1 = UT1::from_julian_date(input);
        let back = ut1.to_tai_with_offset(offset).unwrap();
        let back = back.to_ut1_with_offset(offset).unwrap();
        assert_eq!(
            back.to_julian_date().parts(),
            ut1_back.parts(),
            "{input}, {offset}"
        );

        let tai = TAI::from_julian_date(input);
        let back = tai.to_ut1_with_offset(offset).unwrap();
        let back = back.to_tai_with_offset(offset).unwrap();
        assert_eq!(
            back.to_julian_date().parts(),
            input.parts(),
            "{input}, {offset}"
        );
    }
}

#[test]
fn test_ut1_tt_round_trip_matches_erfa() {
    // eraUt1tt then eraTtut1, and the reverse. TT always comes back exactly;
    // UT1 with a 0.5 day part comes back 1 ulp low for some offsets.
    let jd = JulianDate::new;
    let cases = [
        (jd(J2000_JD, 0.0), 63.8, jd(J2000_JD, 0.0)),
        (jd(J2000_JD, 0.5), 69.0, jd(J2000_JD, 0.5)),
        (jd(J2000_JD, 0.5), 63.8, jd(J2000_JD, 0.49999999999999994)),
        (
            jd(J2000_JD, 0.123456789012345),
            63.8,
            jd(J2000_JD, 0.123456789012345),
        ),
        (
            jd(J2000_JD, -0.123456789012345),
            63.8,
            jd(J2000_JD, -0.123456789012345),
        ),
        (jd(J2000_JD, 0.987654321), 63.8, jd(J2000_JD, 0.987654321)),
        (jd(0.5, J2000_JD), 63.8, jd(0.49999999999999994, J2000_JD)),
    ];
    for (input, delta_t, ut1_back) in cases {
        let ut1 = UT1::from_julian_date(input);
        let back = ut1.to_tt_with_delta_t(delta_t).unwrap();
        let back = back.to_ut1_with_delta_t(delta_t).unwrap();
        assert_eq!(
            back.to_julian_date().parts(),
            ut1_back.parts(),
            "{input}, {delta_t}"
        );

        let tt = TT::from_julian_date(input);
        let back = tt.to_ut1_with_delta_t(delta_t).unwrap();
        let back = back.to_tt_with_delta_t(delta_t).unwrap();
        assert_eq!(
            back.to_julian_date().parts(),
            input.parts(),
            "{input}, {delta_t}"
        );
    }
}

#[test]
fn test_non_finite_julian_date_is_rejected() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let jd = JulianDate::new(bad, 0.0);
        let expected = TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
        let ut1 = UT1::from_julian_date(jd);
        assert_eq!(ut1.to_tai_with_offset(-32.0).unwrap_err(), expected);
        assert_eq!(ut1.to_tt_with_delta_t(64.0).unwrap_err(), expected);
        let tai = TAI::from_julian_date(jd);
        assert_eq!(tai.to_ut1_with_offset(-32.0).unwrap_err(), expected);
        let tt = TT::from_julian_date(jd);
        assert_eq!(tt.to_ut1_with_delta_t(64.0).unwrap_err(), expected);
    }
}

#[test]
fn test_non_finite_offset_is_rejected() {
    let jd = JulianDate::new(J2000_JD, 0.0);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let offset = TimeError::ConversionError(format!(
            "ut1_tai_offset_seconds must be finite, got {}",
            bad
        ));
        let delta_t =
            TimeError::ConversionError(format!("delta_t_seconds must be finite, got {}", bad));
        let ut1 = UT1::from_julian_date(jd);
        assert_eq!(ut1.to_tai_with_offset(bad).unwrap_err(), offset);
        assert_eq!(ut1.to_tt_with_delta_t(bad).unwrap_err(), delta_t);
        let tai = TAI::from_julian_date(jd);
        assert_eq!(tai.to_ut1_with_offset(bad).unwrap_err(), offset);
        let tt = TT::from_julian_date(jd);
        assert_eq!(tt.to_ut1_with_delta_t(bad).unwrap_err(), delta_t);
    }
}
