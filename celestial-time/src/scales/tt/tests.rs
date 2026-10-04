use super::*;
use crate::julian::JulianDate;
use crate::TimeError;
use celestial_core::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD};

#[test]
fn test_tt_julian_year_and_centuries() {
    let tt = TT::j2000();
    assert_eq!(tt.julian_year(), 2000.0);
    assert_eq!(tt.centuries_since_j2000().unwrap(), 0.0);

    let tt_plus_century = tt.add_days(celestial_core::constants::DAYS_PER_JULIAN_CENTURY);
    assert_eq!(tt_plus_century.centuries_since_j2000().unwrap(), 1.0);

    // Values from eraEpj and the (jd1 - J2000) + jd2 form ERFA uses for T.
    let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 9131.987654321));
    assert_eq!(tt.julian_year(), 2025.0020195874633);
    assert_eq!(tt.centuries_since_j2000().unwrap(), 0.2500201958746338);
}

#[test]
fn test_centuries_since_j2000_is_limited_to_twenty() {
    let at = |jd1: f64, jd2: f64| {
        TT::from_julian_date(JulianDate::new(jd1, jd2)).centuries_since_j2000()
    };
    let too_far = |msg: &str| Err(TimeError::InvalidEpoch(msg.into()));

    assert_eq!(at(J2000_JD, 20.0 * DAYS_PER_JULIAN_CENTURY), Ok(20.0));
    assert_eq!(at(J2000_JD, -20.0 * DAYS_PER_JULIAN_CENTURY), Ok(-20.0));
    assert_eq!(
        at(J2000_JD, 21.0 * DAYS_PER_JULIAN_CENTURY),
        too_far("Epoch too far from J2000.0 for the IAU models: 21.0 centuries")
    );
    // Both parts finite, but their sum overflows.
    assert_eq!(
        at(f64::MAX, f64::MAX),
        too_far("Epoch too far from J2000.0 for the IAU models: inf centuries")
    );
    assert_eq!(
        at(J2000_JD, f64::INFINITY),
        too_far("Julian Date (2451545, inf) is not finite")
    );
}
