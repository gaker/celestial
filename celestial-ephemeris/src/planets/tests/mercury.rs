use crate::planets::Vsop2013Mercury;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const MERCURY_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, 0.3493878714, -0.1302077267, -0.1058730361),
    (2415545.0, -0.3953232726, -0.0832703775, -0.0033538163),
    (2419545.0, 0.2950960118, -0.2441772970, -0.1610737357),
    (2423545.0, -0.3676232407, 0.0409400626, 0.0600717273),
    (2427545.0, 0.2077238019, -0.3312635001, -0.1984945320),
    (2431545.0, -0.2846205184, 0.1582302603, 0.1140691451),
    (2435545.0, 0.1004920218, -0.3870987319, -0.2171791680),
    (2439545.0, -0.1477140412, 0.2442561947, 0.1457920872),
    (2443545.0, -0.0153852754, -0.4101850758, -0.2174925982),
    (2447545.0, 0.0231249166, 0.2719576619, 0.1428649538),
    (2451545.0, -0.1300936046, -0.4005937206, -0.2004893069),
];

#[test]
fn vsop2013_vs_reference() {
    let mercury = Vsop2013Mercury;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in MERCURY_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = mercury.heliocentric_position(&tdb).unwrap();

        let dx = pos[0] - x_exp;
        let dy = pos[1] - y_exp;
        let dz = pos[2] - z_exp;
        let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;

        if error_km > max_error_km {
            max_error_km = error_km;
        }
    }

    assert!(
        max_error_km < 5_000.0,
        "Max error {:.0} km exceeds 5,000 km threshold",
        max_error_km
    );
}

#[test]
fn vsop2013_j2000() {
    let mercury = Vsop2013Mercury;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = mercury.heliocentric_position(&tdb).unwrap();

    let expected = (-0.1300936046, -0.4005937206, -0.2004893069);
    let dx = pos[0] - expected.0;
    let dy = pos[1] - expected.1;
    let dz = pos[2] - expected.2;
    let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;

    assert!(
        error_km < 2_000.0,
        "Error {:.0} km exceeds threshold",
        error_km
    );
}

#[test]
fn geocentric_distance_range() {
    let mercury = Vsop2013Mercury;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 4).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = mercury.geocentric_position(&tdb).unwrap();
        let dist = libm::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 0.54,
        "Min distance {:.2} AU below expected 0.54 AU",
        min_dist
    );
    assert!(
        max_dist <= 1.48,
        "Max distance {:.2} AU above expected 1.48 AU",
        max_dist
    );
}
