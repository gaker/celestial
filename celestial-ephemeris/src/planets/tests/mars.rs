use crate::planets::Vsop2013Mars;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const MARS_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, -0.1474461434, -1.3278375334, -0.6049474049),
    (2415545.0, -1.4949210753, -0.5677164770, -0.2196035382),
    (2419545.0, -1.4097862365, 0.7887552075, 0.4001552102),
    (2423545.0, -0.0652445211, 1.4309576308, 0.6580631191),
    (2427545.0, 1.2751062367, 0.5929892819, 0.2372980895),
    (2431545.0, 0.9034102302, -0.9540459212, -0.4621103276),
    (2435545.0, -0.7578584474, -1.2100135320, -0.5344072732),
    (2439545.0, -1.6471756551, -0.0694871483, 0.0127846723),
    (2443545.0, -1.0232278535, 1.1621468626, 0.5607395621),
    (2447545.0, 0.5197104541, 1.3038943187, 0.5839816939),
    (2451545.0, 1.3907159214, 0.0014012149, -0.0369601677),
];

#[test]
fn vsop2013_vs_reference() {
    let mars = Vsop2013Mars;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in MARS_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = mars.heliocentric_position(&tdb).unwrap();

        let dx = pos[0] - x_exp;
        let dy = pos[1] - y_exp;
        let dz = pos[2] - z_exp;
        let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;

        if error_km > max_error_km {
            max_error_km = error_km;
        }
    }

    assert!(
        max_error_km < 10_000.0,
        "Max error {:.0} km exceeds 10,000 km threshold",
        max_error_km
    );
}

#[test]
fn vsop2013_j2000() {
    let mars = Vsop2013Mars;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = mars.heliocentric_position(&tdb).unwrap();

    let expected = (1.3907159214, 0.0014012149, -0.0369601677);
    let dx = pos[0] - expected.0;
    let dy = pos[1] - expected.1;
    let dz = pos[2] - expected.2;
    let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;

    assert!(
        error_km < 5_000.0,
        "Error {:.0} km exceeds threshold",
        error_km
    );
}

#[test]
fn geocentric_distance_range() {
    let mars = Vsop2013Mars;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 15).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = mars.geocentric_position(&tdb).unwrap();
        let dist = libm::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 0.37,
        "Min distance {:.2} AU below expected 0.37 AU",
        min_dist
    );
    assert!(
        max_dist <= 2.68,
        "Max distance {:.2} AU above expected 2.68 AU",
        max_dist
    );
}
