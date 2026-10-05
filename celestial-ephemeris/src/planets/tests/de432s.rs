use crate::earth::Vsop2013Earth;
use crate::jpl::bodies;
use crate::jpl::spk::SpkFile;
use crate::planets::*;
use celestial_core::constants::{AU_KM, J2000_JD, SECONDS_PER_DAY_F64};
use celestial_core::errors::AstroResult;
use celestial_core::matrix::Vector3;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

struct Body {
    name: &'static str,
    id: i32,
    state: fn(&TDB) -> AstroResult<(Vector3, Vector3)>,
    position_bound_km: f64,
    velocity_bound_km_s: f64,
}

// Each bound is the largest miss at EPOCHS, rounded up to two significant
// figures. Truncation and the theory's own error both count.
const BODIES: [Body; 10] = [
    Body {
        name: "mercury",
        id: bodies::MERCURY_BARYCENTER,
        state: |t| Vsop2013Mercury.heliocentric_state(t),
        position_bound_km: 1.2,
        velocity_bound_km_s: 2.6e-6,
    },
    Body {
        name: "venus",
        id: bodies::VENUS_BARYCENTER,
        state: |t| Vsop2013Venus.heliocentric_state(t),
        position_bound_km: 0.56,
        velocity_bound_km_s: 4.4e-7,
    },
    Body {
        name: "emb",
        id: bodies::EARTH_MOON_BARYCENTER,
        state: |t| Vsop2013Emb.heliocentric_state(t),
        position_bound_km: 0.77,
        velocity_bound_km_s: 5.6e-7,
    },
    Body {
        name: "earth",
        id: bodies::EARTH,
        state: |t| Vsop2013Earth::new().heliocentric_state(t),
        position_bound_km: 0.77,
        velocity_bound_km_s: 5.6e-7,
    },
    Body {
        name: "mars",
        id: bodies::MARS_BARYCENTER,
        state: |t| Vsop2013Mars.heliocentric_state(t),
        position_bound_km: 0.72,
        velocity_bound_km_s: 4.5e-7,
    },
    Body {
        name: "jupiter",
        id: bodies::JUPITER_BARYCENTER,
        state: |t| Vsop2013Jupiter.heliocentric_state(t),
        position_bound_km: 67.0,
        velocity_bound_km_s: 9.0e-6,
    },
    Body {
        name: "saturn",
        id: bodies::SATURN_BARYCENTER,
        state: |t| Vsop2013Saturn.heliocentric_state(t),
        position_bound_km: 23.0,
        velocity_bound_km_s: 3.9e-6,
    },
    Body {
        name: "uranus",
        id: bodies::URANUS_BARYCENTER,
        state: |t| Vsop2013Uranus.heliocentric_state(t),
        position_bound_km: 6000.0,
        velocity_bound_km_s: 3.9e-4,
    },
    Body {
        name: "neptune",
        id: bodies::NEPTUNE_BARYCENTER,
        state: |t| Vsop2013Neptune.heliocentric_state(t),
        position_bound_km: 2900.0,
        velocity_bound_km_s: 2.9e-4,
    },
    Body {
        name: "pluto",
        id: bodies::PLUTO_BARYCENTER,
        state: |t| Vsop2013Pluto.heliocentric_state(t),
        position_bound_km: 43000.0,
        velocity_bound_km_s: 1.9e-3,
    },
];

const EPOCHS: [f64; 3] = [J2000_JD, 2440000.5, 2460000.5];

fn misses(spk: &SpkFile, body: &Body, jd: f64) -> (f64, f64) {
    let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
    let (de_position, de_velocity) = spk.compute_state(body.id, bodies::SUN, &tdb).unwrap();
    let (position, velocity) = (body.state)(&tdb).unwrap();
    let to_km_s = AU_KM / SECONDS_PER_DAY_F64;
    (
        (position * AU_KM - de_position).magnitude(),
        (velocity * to_km_s - de_velocity).magnitude(),
    )
}

#[test]
fn states_stay_near_de432s() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/de432s.bsp");
    let spk = SpkFile::open(path).unwrap();
    for body in &BODIES {
        for jd in EPOCHS {
            let (position, velocity) = misses(&spk, body, jd);
            assert!(
                position < body.position_bound_km,
                "{} at JD {jd}: {position} km",
                body.name
            );
            assert!(
                velocity < body.velocity_bound_km_s,
                "{} at JD {jd}: {velocity} km/s",
                body.name
            );
        }
    }
}
