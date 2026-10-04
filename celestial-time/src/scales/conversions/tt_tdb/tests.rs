use super::*;
use crate::TimeError;
use celestial_core::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD};

#[test]
fn test_tt_tdb_with_location_matches_erfa() {
    // eraTttdb and eraTdbtt with UT1 = TT and dtr from eraDtdb, evaluated at
    // the input date. At a (JD, 0) split jd2 is the offset itself.
    let greenwich = greenwich_location().unwrap();
    let tokyo = Location::from_degrees(35.6762, 139.6503, 40.0).unwrap();
    let sydney = Location::from_degrees(-33.8688, 151.2093, 58.0).unwrap();
    let cases = [
        (J2000_JD, greenwich, 1.1510214835963987e-9),
        (
            J2000_JD - DAYS_PER_JULIAN_CENTURY,
            greenwich,
            3.871874123266885e-10,
        ),
        (J2000_JD + 18262.5, greenwich, 9.293054160090208e-10),
        (
            J2000_JD + DAYS_PER_JULIAN_CENTURY,
            greenwich,
            8.769430614331916e-10,
        ),
        (J2000_JD, tokyo, 1.1632537479395206e-9),
        (J2000_JD, sydney, 1.1580498041164867e-9),
    ];
    for (jd, location, offset_days) in cases {
        let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
        let tdb = tt.to_tdb_with_location(&location).unwrap();
        assert_eq!(tdb.to_julian_date().parts(), (jd, -offset_days));

        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let tt = tdb.to_tt_with_location(&location).unwrap();
        assert_eq!(tt.to_julian_date().parts(), (jd, offset_days));
    }
}

#[test]
fn test_tt_tdb_round_trip_matches_erfa() {
    // eraTttdb then eraTdbtt, and the reverse, at Greenwich with UT1 = TT and
    // dtr from eraDtdb at each step's input. Only the (JD, 0) split doesn't
    // come back exactly, by under 3e-19 d.
    let cases = [
        (J2000_JD, 0.0, 2.7330316272451114e-19, 2.733033695196643e-19),
        (J2000_JD, 0.5, 0.5, 0.5),
        (
            J2000_JD,
            0.123456789012345,
            0.123456789012345,
            0.123456789012345,
        ),
        (
            J2000_JD,
            -0.123456789012345,
            -0.123456789012345,
            -0.123456789012345,
        ),
    ];
    for (jd1, jd2, tt_back, tdb_back) in cases {
        let tt = TT::from_julian_date(JulianDate::new(jd1, jd2));
        let back = tt.to_tdb_greenwich().unwrap().to_tt_greenwich().unwrap();
        assert_eq!(back.to_julian_date().parts(), (jd1, tt_back));

        let tdb = TDB::from_julian_date(JulianDate::new(jd1, jd2));
        let back = tdb.to_tt_greenwich().unwrap().to_tdb_greenwich().unwrap();
        assert_eq!(back.to_julian_date().parts(), (jd1, tdb_back));
    }

    let tt = TT::from_julian_date(JulianDate::new(0.5, J2000_JD));
    let tdb = tt.to_tdb_greenwich().unwrap();
    assert_eq!(tdb.to_julian_date().parts(), (0.4999999990169487, J2000_JD));
    let back = tdb.to_tt_greenwich().unwrap();
    assert_eq!(back.to_julian_date().parts(), (0.5, J2000_JD));
}

#[test]
fn test_tt_tdb_with_ut1_offset_matches_erfa() {
    // TDB from eraTttdb, with dtr from eraDtdb given the UT1 day fraction
    // counted from midnight (UT1 from eraTtut1). TT back from eraTdbtt, with
    // dtr evaluated at TDB.
    let greenwich = Location::from_degrees(51.477928, 0.0, 46.0).unwrap();
    let tokyo = Location::from_degrees(35.6762, 139.6503, 40.0).unwrap();
    let sydney = Location::from_degrees(-33.8688, 151.2093, 58.0).unwrap();
    let mauna_kea = Location::from_degrees(19.8207, -155.4681, 4205.0).unwrap();
    let cases = [
        (
            J2000_JD,
            0.0,
            greenwich,
            0.0,
            -1.1510214835963987e-9,
            2.7330316272451114e-19,
        ),
        (J2000_JD, 0.25, tokyo, -63.8285, 0.24999999894881492, 0.25),
        (
            2400000.5,
            60000.75,
            sydney,
            -69.184,
            60000.75000001504,
            60000.75,
        ),
        (2460000.5, 0.375, mauna_kea, 0.3, 0.3750000149283695, 0.375),
        (2415020.5, 0.9, greenwich, -2.5, 0.9000000000856577, 0.9),
    ];
    for (jd1, jd2, location, ut1_minus_tt, tdb_jd2, tt_jd2) in cases {
        let tt = TT::from_julian_date(JulianDate::new(jd1, jd2));

        let tdb = tt
            .to_tdb_with_location_and_ut1_offset(&location, ut1_minus_tt)
            .unwrap();
        assert_eq!(tdb.to_julian_date().jd1(), jd1);
        assert_eq!(tdb.to_julian_date().jd2(), tdb_jd2, "TT ({jd1}, {jd2})");

        let back = tdb
            .to_tt_with_location_and_ut1_offset(&location, ut1_minus_tt)
            .unwrap();
        assert_eq!(back.to_julian_date().jd1(), jd1);
        assert_eq!(
            back.to_julian_date().jd2(),
            tt_jd2,
            "TDB ({jd1}, {tdb_jd2})"
        );
    }
}

#[test]
fn test_greenwich_conversions_match_erfa() {
    // eraTttdb, with dtr from eraDtdb at 51.477928° N, 0°, 46 m and UT1 = TT.
    let cases = [
        (J2000_JD, 0.123456, 0.12345599887953732),
        (J2000_JD, 0.987654, 0.9876539991812064),
        (2400000.5, 60000.75, 60000.750000014996),
    ];
    for (jd1, tt_jd2, tdb_jd2) in cases {
        let tt = TT::from_julian_date(JulianDate::new(jd1, tt_jd2));

        let tdb = tt.to_tdb_greenwich().unwrap();
        assert_eq!(tdb.to_julian_date().jd1(), jd1);
        assert_eq!(tdb.to_julian_date().jd2(), tdb_jd2, "TT ({jd1}, {tt_jd2})");

        let back = tdb.to_tt_greenwich().unwrap();
        assert_eq!(back.to_julian_date().jd1(), jd1);
        assert_eq!(
            back.to_julian_date().jd2(),
            tt_jd2,
            "TDB ({jd1}, {tdb_jd2})"
        );
    }
}

#[test]
fn test_non_finite_julian_date_is_rejected() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let jd = JulianDate::new(bad, 0.0);
        let expected = TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
        let location = Location::greenwich();
        let tt = TT::from_julian_date(jd);
        assert_eq!(tt.to_tdb_greenwich().unwrap_err(), expected);
        assert_eq!(tt.to_tdb_with_offset(0.001).unwrap_err(), expected);
        let tdb = TDB::from_julian_date(jd);
        assert_eq!(tdb.to_tt_greenwich().unwrap_err(), expected);
        assert_eq!(tdb.to_tt_with_offset(0.001).unwrap_err(), expected);
        let offset = compute_tdb_tt_offset(&jd, 0.5, &location);
        assert_eq!(offset.unwrap_err(), expected);
    }
}

#[test]
fn test_non_finite_offset_is_rejected() {
    let location = Location::greenwich();
    let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let dtr = TimeError::ConversionError(format!("dtr_seconds must be finite, got {}", bad));
        assert_eq!(tt.to_tdb_with_offset(bad).unwrap_err(), dtr);
        assert_eq!(tdb.to_tt_with_offset(bad).unwrap_err(), dtr);
        let ut1 =
            TimeError::ConversionError(format!("ut1_minus_tt_seconds must be finite, got {}", bad));
        let to_tdb = tt.to_tdb_with_location_and_ut1_offset(&location, bad);
        assert_eq!(to_tdb.unwrap_err(), ut1);
        let to_tt = tdb.to_tt_with_location_and_ut1_offset(&location, bad);
        assert_eq!(to_tt.unwrap_err(), ut1);
        let fraction =
            TimeError::ConversionError(format!("ut1_fraction must be finite, got {}", bad));
        let offset = compute_tdb_tt_offset(&tt.to_julian_date(), bad, &location);
        assert_eq!(offset.unwrap_err(), fraction);
    }
}
