use crate::planets::Vsop2013Saturn;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const SATURN_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, -8.5151099046, 3.2825565780, 1.7200893604),
    (2415545.0, 2.3770579174, -8.9996266643, -3.8169867762),
    (2419545.0, 5.2172917580, 6.9815327926, 2.6579637103),
    (2423545.0, -9.1381889728, -3.0134497674, -0.8520810205),
    (2427545.0, 7.7506635038, -5.5025418118, -2.6046062335),
    (2431545.0, -1.7986439511, 8.1545304723, 3.4443131779),
    (2435545.0, -5.2734506611, -7.8866680116, -3.0300154867),
    (2439545.0, 9.5012015992, 0.4527726650, -0.2214084263),
    (2443545.0, -7.7367291471, 4.5131027003, 2.1962088916),
    (2447545.0, 1.0138796256, -9.2197781359, -3.8512664946),
    (2451545.0, 6.4064088704, 6.1746578061, 2.2747707349),
];

#[test]
fn vsop2013_vs_reference() {
    let saturn = Vsop2013Saturn;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in SATURN_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = saturn.heliocentric_position(&tdb).unwrap();

        let dx = pos[0] - x_exp;
        let dy = pos[1] - y_exp;
        let dz = pos[2] - z_exp;
        let error_km = (dx * dx + dy * dy + dz * dz).sqrt() * AU_KM;

        if error_km > max_error_km {
            max_error_km = error_km;
        }
    }

    assert!(
        max_error_km < 50_000.0,
        "Max error {:.0} km exceeds 50,000 km threshold",
        max_error_km
    );
}

#[test]
fn vsop2013_j2000() {
    let saturn = Vsop2013Saturn;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = saturn.heliocentric_position(&tdb).unwrap();

    let expected = (6.4064088704, 6.1746578061, 2.2747707349);
    let dx = pos[0] - expected.0;
    let dy = pos[1] - expected.1;
    let dz = pos[2] - expected.2;
    let error_km = (dx * dx + dy * dy + dz * dz).sqrt() * AU_KM;

    assert!(
        error_km < 20_000.0,
        "Error {:.0} km exceeds threshold",
        error_km
    );
}

#[test]
fn geocentric_distance_range() {
    let saturn = Vsop2013Saturn;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 20).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = saturn.geocentric_position(&tdb).unwrap();
        let dist = (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z).sqrt();
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 7.9,
        "Min distance {:.2} AU below expected 7.9 AU",
        min_dist
    );
    assert!(
        max_dist <= 11.1,
        "Max distance {:.2} AU above expected 11.1 AU",
        max_dist
    );
}
