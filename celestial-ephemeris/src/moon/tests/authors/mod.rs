use celestial_core::constants::AU_KM;
use celestial_core::errors::AstroResult;
use celestial_time::scales::tdb::TDB;

use super::tdb;
use crate::moon::ElpMpp02Moon;

mod vectors;
use vectors::{DE405_ECLIPTIC, DE405_ICRS, EPOCHS, LLR_ECLIPTIC, LLR_ICRS};

type State = fn(&ElpMpp02Moon, &TDB) -> AstroResult<[f64; 6]>;

fn check(moon: &ElpMpp02Moon, state: State, expected: &[[f64; 6]; 8]) {
    for (&jd, want) in EPOCHS.iter().zip(expected) {
        assert_eq!(state(moon, &tdb(jd)).unwrap(), *want, "JD {}", jd);
    }
}

fn icrs_au(moon: &ElpMpp02Moon, tdb: &TDB) -> AstroResult<[f64; 6]> {
    let (p, v) = moon.geocentric_state(tdb)?;
    Ok([p.x, p.y, p.z, v.x, v.y, v.z])
}

fn in_au(rows: &[[f64; 6]; 8]) -> [[f64; 6]; 8] {
    rows.map(|row| row.map(|v| v / AU_KM))
}

#[test]
fn llr_ecliptic_state_matches_the_authors_code() {
    let moon = ElpMpp02Moon::new();
    check(&moon, ElpMpp02Moon::ecliptic_state, &LLR_ECLIPTIC);
}

#[test]
fn de405_ecliptic_state_matches_the_authors_code() {
    let moon = ElpMpp02Moon::with_de405_fit();
    check(&moon, ElpMpp02Moon::ecliptic_state, &DE405_ECLIPTIC);
}

#[test]
fn llr_icrs_state_uses_the_iau_2006_obliquity() {
    let moon = ElpMpp02Moon::new();
    check(&moon, icrs_au, &in_au(&LLR_ICRS));
}

#[test]
fn de405_icrs_state_uses_the_angles_fitted_with_it() {
    let moon = ElpMpp02Moon::with_de405_fit();
    check(&moon, icrs_au, &in_au(&DE405_ICRS));
}
