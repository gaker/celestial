use crate::planets::Vsop2013Neptune;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const NEPTUNE_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, 12.1323801234, 25.3226220555, 10.0628233249),
    (2415545.0, -0.1360931298, 27.6520562678, 11.3214894525),
    (2419545.0, -12.3856600756, 25.1528348111, 10.6033851946),
    (2423545.0, -22.4955863232, 18.2841471648, 8.0436712472),
    (2427545.0, -28.7387978923, 8.2766903355, 4.1028690883),
    (2431545.0, -30.1064751093, -3.1360332460, -0.5343057384),
    (2435545.0, -26.3968325711, -14.0307173050, -5.0858772277),
    (2439545.0, -18.2523466372, -22.5660704974, -8.7822092401),
    (2443545.0, -7.0565022877, -27.3194094185, -11.0063258261),
    (2447545.0, 5.3316427610, -27.4815799442, -11.3810940223),
    (2451545.0, 16.8120479567, -22.9801038994, -9.8244204429),
];

#[test]
fn vsop2013_vs_reference() {
    let neptune = Vsop2013Neptune;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in NEPTUNE_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = neptune.heliocentric_position(&tdb).unwrap();

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
    let neptune = Vsop2013Neptune;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = neptune.heliocentric_position(&tdb).unwrap();

    let expected = (16.8120479567, -22.9801038994, -9.8244204429);
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
    let neptune = Vsop2013Neptune;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 20).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = neptune.geocentric_position(&tdb).unwrap();
        let dist = (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z).sqrt();
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 28.7,
        "Min distance {:.2} AU below expected 28.7 AU",
        min_dist
    );
    assert!(
        max_dist <= 31.4,
        "Max distance {:.2} AU above expected 31.4 AU",
        max_dist
    );
}
