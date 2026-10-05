use crate::planets::*;
use celestial_core::constants::J2000_JD;
use celestial_core::errors::{AstroError, AstroResult, MathErrorKind};
use celestial_core::matrix::Vector3;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

type Position = fn(&TDB) -> AstroResult<Vector3>;
type State = fn(&TDB) -> AstroResult<(Vector3, Vector3)>;

// A valid Earth state, so that the planet's own epoch check is what fails.
const EARTH: (Vector3, Vector3) = (
    Vector3 {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    },
    Vector3 {
        x: 0.0,
        y: 0.0172,
        z: 0.0,
    },
);

const HELIO: [(&str, Position); 9] = [
    ("mercury", |t| Vsop2013Mercury.heliocentric_position(t)),
    ("venus", |t| Vsop2013Venus.heliocentric_position(t)),
    ("emb", |t| Vsop2013Emb.heliocentric_position(t)),
    ("mars", |t| Vsop2013Mars.heliocentric_position(t)),
    ("jupiter", |t| Vsop2013Jupiter.heliocentric_position(t)),
    ("saturn", |t| Vsop2013Saturn.heliocentric_position(t)),
    ("uranus", |t| Vsop2013Uranus.heliocentric_position(t)),
    ("neptune", |t| Vsop2013Neptune.heliocentric_position(t)),
    ("pluto", |t| Vsop2013Pluto.heliocentric_position(t)),
];

const OTHER_POSITIONS: [(&str, Position); 2] = [
    ("mars geocentric", |t| {
        Vsop2013Mars.geocentric_position(t, &EARTH.0)
    }),
    ("pluto geocentric", |t| {
        Vsop2013Pluto.geocentric_position(t, &EARTH.0)
    }),
];

const STATES: [(&str, State); 4] = [
    ("mars heliocentric", |t| Vsop2013Mars.heliocentric_state(t)),
    ("mars geocentric", |t| {
        Vsop2013Mars.geocentric_state(t, &EARTH)
    }),
    ("pluto heliocentric", |t| {
        Vsop2013Pluto.heliocentric_state(t)
    }),
    ("pluto geocentric", |t| {
        Vsop2013Pluto.geocentric_state(t, &EARTH)
    }),
];

fn tdb(jd1: f64, jd2: f64) -> TDB {
    TDB::from_julian_date(JulianDate::new(jd1, jd2))
}

fn kind<T: std::fmt::Debug>(name: &str, jd: f64, result: AstroResult<T>) -> MathErrorKind {
    match result {
        Err(AstroError::MathError { kind, .. }) => kind,
        other => panic!(
            "{} at JD {}: expected a MathError, got {:?}",
            name, jd, other
        ),
    }
}

#[test]
fn non_finite_epochs_are_errors() {
    let epochs = [
        (f64::NAN, 0.0),
        (J2000_JD, f64::NAN),
        (f64::INFINITY, 0.0),
        (J2000_JD, f64::NEG_INFINITY),
    ];
    for (jd1, jd2) in epochs {
        let t = tdb(jd1, jd2);
        for (name, f) in HELIO.iter().chain(OTHER_POSITIONS.iter()) {
            assert_eq!(kind(name, jd1 + jd2, f(&t)), MathErrorKind::NotFinite);
        }
        for (name, f) in STATES {
            assert_eq!(kind(name, jd1 + jd2, f(&t)), MathErrorKind::NotFinite);
        }
    }
}

#[test]
fn far_epochs_are_errors_not_hangs() {
    for jd in [1e12, -1e12] {
        for (name, f) in HELIO {
            assert_eq!(kind(name, jd, f(&tdb(jd, 0.0))), MathErrorKind::OutOfRange);
        }
    }
    let jd = -4853455.0;
    let pluto = Vsop2013Pluto.heliocentric_position(&tdb(jd, 0.0));
    assert_eq!(kind("pluto", jd, pluto), MathErrorKind::OutOfRange);
}

#[test]
fn epochs_outside_the_theory_range_are_errors() {
    let cases: [(&str, Position, f64); 8] = [
        ("pluto", HELIO[8].1, 7930295.0),
        ("pluto", HELIO[8].1, -3027205.0),
        ("neptune", HELIO[7].1, 1e8),
        ("jupiter", HELIO[4].1, -15810955.0),
        ("mars", HELIO[3].1, 260044.0),
        ("mars", HELIO[3].1, 4643046.0),
        ("pluto", HELIO[8].1, 1721044.0),
        ("pluto", HELIO[8].1, 3182046.0),
    ];
    for (name, f, jd) in cases {
        assert_eq!(kind(name, jd, f(&tdb(jd, 0.0))), MathErrorKind::OutOfRange);
    }
}

#[test]
fn every_planet_is_finite_across_its_range() {
    for (name, f) in HELIO {
        let (first, last) = if name == "pluto" {
            (1721045.0, 3182045.0)
        } else {
            (260045.0, 4643045.0)
        };
        for i in 0..=60 {
            let jd = first + (last - first) * i as f64 / 60.0;
            let p = f(&tdb(jd, 0.0)).unwrap_or_else(|e| panic!("{} at JD {}: {}", name, jd, e));
            let r = libm::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
            assert!(r > 0.3 && r < 60.0, "{} at JD {}: r = {} AU", name, jd, r);
        }
    }
}
