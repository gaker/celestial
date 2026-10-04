use super::*;
use crate::julian::JulianDate;
use crate::TimeError;
use celestial_core::constants::J2000_JD;

#[test]
fn test_tai_tt_offset_32_184_seconds() {
    let test_dates = [
        (J2000_JD, "J2000.0"),
        (2455197.5, "2010-01-01"),
        (2459580.5, "2022-01-01"),
        (2440587.5, "1970-01-01 Unix epoch"),
    ];

    for (jd, description) in test_dates {
        let tai = TAI::from_julian_date(JulianDate::new(jd, 0.0));
        let tt = tai.to_tt().unwrap();

        let tai_jd = tai.to_julian_date();
        let tt_jd = tt.to_julian_date();

        let offset_days = (tt_jd.jd1() - tai_jd.jd1()) + (tt_jd.jd2() - tai_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        assert_eq!(
            offset_seconds, 32.184,
            "{}: TAI->TT offset must be exactly 32.184 seconds",
            description
        );

        let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
        let tai = tt.to_tai().unwrap();

        let tt_jd = tt.to_julian_date();
        let tai_jd = tai.to_julian_date();

        let offset_days = (tt_jd.jd1() - tai_jd.jd1()) + (tt_jd.jd2() - tai_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        assert_eq!(
            offset_seconds, 32.184,
            "{}: TT->TAI means TT is 32.184 seconds ahead",
            description
        );
    }
}

#[test]
fn test_tai_tt_round_trip_precision() {
    let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];

    for jd2 in test_jd2_values {
        let original_tai = TAI::from_julian_date(JulianDate::new(J2000_JD, jd2));
        let tt = original_tai.to_tt().unwrap();
        let round_trip_tai = tt.to_tai().unwrap();

        assert_eq!(
            original_tai.to_julian_date().jd1(),
            round_trip_tai.to_julian_date().jd1(),
            "TAI->TT->TAI JD1 must be exact for jd2={}",
            jd2
        );
        assert_eq!(
            original_tai.to_julian_date().jd2(),
            round_trip_tai.to_julian_date().jd2(),
            "TAI->TT->TAI JD2 must be exact for jd2={}",
            jd2
        );

        let original_tt = TT::from_julian_date(JulianDate::new(J2000_JD, jd2));
        let tai = original_tt.to_tai().unwrap();
        let round_trip_tt = tai.to_tt().unwrap();

        assert_eq!(
            original_tt.to_julian_date().jd1(),
            round_trip_tt.to_julian_date().jd1(),
            "TT->TAI->TT JD1 must be exact for jd2={}",
            jd2
        );
        assert_eq!(
            original_tt.to_julian_date().jd2(),
            round_trip_tt.to_julian_date().jd2(),
            "TT->TAI->TT JD2 must be exact for jd2={}",
            jd2
        );
    }

    let alt_tai = TAI::from_julian_date(JulianDate::new(0.5, J2000_JD));
    let alt_tt = alt_tai.to_tt().unwrap();
    let alt_round_trip = alt_tt.to_tai().unwrap();

    assert_eq!(
        alt_tai.to_julian_date().jd1(),
        alt_round_trip.to_julian_date().jd1(),
        "Alternate split TAI->TT->TAI JD1 must be exact"
    );
    assert_eq!(
        alt_tai.to_julian_date().jd2(),
        alt_round_trip.to_julian_date().jd2(),
        "Alternate split TAI->TT->TAI JD2 must be exact"
    );
}

#[test]
fn test_tai_tcg_chain_round_trip_matches_erfa() {
    // eraTaitt and eraTttcg, then eraTcgtt and eraTttai, and the reverse. TAI
    // always comes back exactly; TCG at two of these splits doesn't.
    let cases = [
        (0.0, 0.0, 2.541098841762901e-21),
        (0.123456789, 0.123456789, 0.123456789),
        (0.5, 0.5, 0.5),
        (-0.25, -0.25, -0.24999999999999997),
    ];
    for (jd2, tai_back, tcg_back) in cases {
        let tai = TAI::from_julian_date(JulianDate::new(J2000_JD, jd2));
        let back = tai.to_tcg().unwrap().to_tai().unwrap();
        assert_eq!(back.to_julian_date().parts(), (J2000_JD, tai_back));

        let tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, jd2));
        let back = tcg.to_tai().unwrap().to_tcg().unwrap();
        assert_eq!(back.to_julian_date().parts(), (J2000_JD, tcg_back));
    }
}

#[test]
fn test_non_finite_julian_date_is_rejected() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let jd = JulianDate::new(bad, 0.0);
        let expected = TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
        assert_eq!(TAI::from_julian_date(jd).to_tt().unwrap_err(), expected);
        assert_eq!(TT::from_julian_date(jd).to_tai().unwrap_err(), expected);
    }
}
