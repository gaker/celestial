use crate::planets::Vsop2013Jupiter;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const JUPITER_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, 2.9837884053, -3.7723816270, -1.6901903627),
    (2415545.0, 0.7069925496, -4.7445829284, -2.0513981365),
    (2419545.0, -1.7382880384, -4.6471837333, -1.9499553658),
    (2423545.0, -3.8212272643, -3.5628413190, -1.4341490904),
    (2427545.0, -5.1281915292, -1.7498001410, -0.6250207356),
    (2431545.0, -5.4147718783, 0.4194201580, 0.3118141091),
    (2435545.0, -4.6121619014, 2.4953666775, 1.1821569858),
    (2439545.0, -2.8475608552, 4.0524137157, 1.8065920154),
    (2443545.0, -0.4627146537, 4.7108655871, 2.0307146958),
    (2447545.0, 2.0318689801, 4.2537285712, 1.7738331572),
    (2451545.0, 4.0011771819, 2.7365785897, 1.0755125254),
];

#[test]
fn vsop2013_vs_reference() {
    let jupiter = Vsop2013Jupiter;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in JUPITER_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = jupiter.heliocentric_position(&tdb).unwrap();

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
    let jupiter = Vsop2013Jupiter;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = jupiter.heliocentric_position(&tdb).unwrap();

    let expected = (4.0011771819, 2.7365785897, 1.0755125254);
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
    let jupiter = Vsop2013Jupiter;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 12).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = jupiter.geocentric_position(&tdb).unwrap();
        let dist = (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z).sqrt();
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 3.9,
        "Min distance {:.2} AU below expected 3.9 AU",
        min_dist
    );
    assert!(
        max_dist <= 6.5,
        "Max distance {:.2} AU above expected 6.5 AU",
        max_dist
    );
}
