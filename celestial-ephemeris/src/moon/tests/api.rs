use celestial_core::constants::J2000_JD;
use celestial_core::errors::{AstroError, AstroResult, MathErrorKind};
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

use super::tdb;
use crate::moon::ElpMpp02Moon;

#[test]
fn default_is_the_llr_fit() {
    let t = tdb(J2000_JD + 1000.0);
    assert_eq!(
        ElpMpp02Moon::default().geocentric_state(&t).unwrap(),
        ElpMpp02Moon::new().geocentric_state(&t).unwrap()
    );
}

#[test]
fn positions_are_the_state_positions() {
    let t = tdb(J2000_JD + 1000.0);
    for moon in [ElpMpp02Moon::new(), ElpMpp02Moon::with_de405_fit()] {
        let ecliptic = moon.ecliptic_state(&t).unwrap();
        assert_eq!(moon.ecliptic_position(&t).unwrap(), ecliptic[..3]);
        let (position, _) = moon.geocentric_state(&t).unwrap();
        assert_eq!(moon.geocentric_position(&t).unwrap(), position);
    }
}

fn error_kind<T: std::fmt::Debug>(jd: f64, result: AstroResult<T>) -> MathErrorKind {
    match result {
        Err(AstroError::MathError { kind, .. }) => kind,
        other => panic!("JD {}: expected a MathError, got {:?}", jd, other),
    }
}

fn all_methods(moon: &ElpMpp02Moon, jd1: f64, jd2: f64) -> [MathErrorKind; 4] {
    let tdb = TDB::from_julian_date(JulianDate::new(jd1, jd2));
    let jd = jd1 + jd2;
    [
        error_kind(jd, moon.geocentric_position(&tdb)),
        error_kind(jd, moon.geocentric_state(&tdb)),
        error_kind(jd, moon.ecliptic_position(&tdb)),
        error_kind(jd, moon.ecliptic_state(&tdb)),
    ]
}

#[test]
fn non_finite_epochs_are_errors() {
    let epochs = [
        (f64::NAN, 0.0),
        (J2000_JD, f64::NAN),
        (f64::INFINITY, 0.0),
        (J2000_JD, f64::NEG_INFINITY),
    ];
    for moon in [ElpMpp02Moon::new(), ElpMpp02Moon::with_de405_fit()] {
        for (jd1, jd2) in epochs {
            assert_eq!(all_methods(&moon, jd1, jd2), [MathErrorKind::NotFinite; 4]);
        }
    }
}

#[test]
fn epochs_outside_the_theory_range_are_errors() {
    for moon in [ElpMpp02Moon::new(), ElpMpp02Moon::with_de405_fit()] {
        for jd in [1e8, -1e8, 625294.0, 2816796.0] {
            assert_eq!(all_methods(&moon, jd, 0.0), [MathErrorKind::OutOfRange; 4]);
        }
    }
}
