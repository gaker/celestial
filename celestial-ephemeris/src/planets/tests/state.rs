use crate::earth::Vsop2013Earth;
use crate::planets::*;
use celestial_core::errors::AstroResult;
use celestial_core::matrix::Vector3;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

type State = (Vector3, Vector3);

struct Planet {
    name: &'static str,
    position: fn(&TDB) -> AstroResult<Vector3>,
    state: fn(&TDB) -> AstroResult<State>,
    geocentric_position: fn(&TDB, &Vector3) -> AstroResult<Vector3>,
    geocentric_state: fn(&TDB, &State) -> AstroResult<State>,
}

macro_rules! planet {
    ($name:expr, $body:ident) => {
        Planet {
            name: $name,
            position: |t| $body.heliocentric_position(t),
            state: |t| $body.heliocentric_state(t),
            geocentric_position: |t, e| $body.geocentric_position(t, e),
            geocentric_state: |t, e| $body.geocentric_state(t, e),
        }
    };
}

const PLANETS: [Planet; 9] = [
    planet!("mercury", Vsop2013Mercury),
    planet!("venus", Vsop2013Venus),
    planet!("emb", Vsop2013Emb),
    planet!("mars", Vsop2013Mars),
    planet!("jupiter", Vsop2013Jupiter),
    planet!("saturn", Vsop2013Saturn),
    planet!("uranus", Vsop2013Uranus),
    planet!("neptune", Vsop2013Neptune),
    planet!("pluto", Vsop2013Pluto),
];

// A two-part epoch, which the summed Julian date would round.
const EPOCH: (f64, f64) = (2445665.5, 0.45541555912031195);

// From this crate, in the order of PLANETS. An independent evaluation of the
// same truncated series, with the series differentiated term by term and the
// Kepler map by Richardson extrapolation, gives these velocities to 1.5e-10 of
// the speed; a 40-digit evaluation of the series gives these positions to
// 0.03 m.
const STATES: [[f64; 6]; 9] = [
    [
        0.19412436955238996,
        -0.34079073692149586,
        -0.20216997207380485,
        0.019598155153997415,
        0.012988283102935603,
        0.004903948242564971,
    ],
    [
        -0.38900173473034066,
        0.5413065489806661,
        0.26812090403987127,
        -0.017067855074653207,
        -0.010490861559844676,
        -0.0036387266931846574,
    ],
    [
        0.4192444489013819,
        0.8194750798011596,
        0.35532170668282403,
        -0.015855459385801914,
        0.006648057827724448,
        0.0028825199250823445,
    ],
    [
        -1.5075114721669827,
        0.6280303845354325,
        0.32886728147153843,
        -0.005415376706374548,
        -0.010486843879051536,
        -0.004663270178425945,
    ],
    [
        -0.8445865696708112,
        -4.814079745319168,
        -2.0429575818597328,
        0.007357913230046792,
        -0.0007210726910086498,
        -0.0004883817020829518,
    ],
    [
        -7.699108810818928,
        -5.749469929566868,
        -2.0435255490243023,
        0.0031690328475300713,
        -0.004010196365014657,
        -0.0017923546501991397,
    ],
    [
        -6.635187641838182,
        -16.322713735334307,
        -7.054740061428107,
        0.0036638314694881183,
        -0.0014108552326361704,
        -0.0006698855014339356,
    ],
    [
        -0.5066575285276077,
        -28.009923781384988,
        -11.452162303310196,
        0.00312663226565551,
        -5.731058127211438e-06,
        -8.019072946872956e-05,
    ],
    [
        -24.750505712765893,
        -16.536533108890097,
        2.294742704549952,
        0.0018613286346112739,
        -0.002649074098208093,
        -0.0013853598822493765,
    ],
];

fn tdb((jd1, jd2): (f64, f64)) -> TDB {
    TDB::from_julian_date(JulianDate::new(jd1, jd2))
}

fn flat((p, v): State) -> [f64; 6] {
    [p.x, p.y, p.z, v.x, v.y, v.z]
}

#[test]
fn states_match_the_pinned_values() {
    let t = tdb(EPOCH);
    for (planet, want) in PLANETS.iter().zip(STATES) {
        assert_eq!(flat((planet.state)(&t).unwrap()), want, "{}", planet.name);
    }
}

#[test]
fn positions_are_the_state_positions() {
    for epoch in [EPOCH, (2451545.0, 0.0), (2816795.5, -0.125)] {
        let t = tdb(epoch);
        for planet in &PLANETS {
            let (position, _) = (planet.state)(&t).unwrap();
            let got = (planet.position)(&t).unwrap();
            assert_eq!(got, position, "{} at {:?}", planet.name, epoch);
        }
    }
}

#[test]
fn geocentric_is_heliocentric_less_the_earth() {
    let t = tdb(EPOCH);
    let earth = Vsop2013Earth::new().heliocentric_state(&t).unwrap();
    for planet in &PLANETS {
        let (position, velocity) = (planet.state)(&t).unwrap();
        let got = (planet.geocentric_state)(&t, &earth).unwrap();
        assert_eq!(
            got,
            (position - earth.0, velocity - earth.1),
            "{}",
            planet.name
        );
        let got = (planet.geocentric_position)(&t, &earth.0).unwrap();
        assert_eq!(got, position - earth.0, "{}", planet.name);
    }
}

#[test]
fn the_epoch_keeps_both_parts() {
    // 1e-10 day is below the resolution of the summed date.
    let jd1 = 2460000.5;
    assert_eq!(jd1 + 1e-10, jd1);
    for planet in &PLANETS {
        let at = |jd2| (planet.state)(&tdb((jd1, jd2))).unwrap();
        let (moved, still) = (at(1e-10), at(0.0));
        assert_ne!(moved.0, still.0, "{}", planet.name);
        assert_ne!(moved.1, still.1, "{}", planet.name);
    }
}
