use crate::planets::Vsop2013Venus;
use astro_core::constants::{AU_KM, J2000_JD};
use astro_time::julian::JulianDate;
use astro_time::TDB;

const VENUS_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, -0.7178452043, 0.0139241146, 0.0517532468),
    (2415545.0, -0.1846601606, 0.6292563055, 0.2945858040),
    (2419545.0, 0.6061342811, 0.3741316134, 0.1298042734),
    (2423545.0, 0.5740437374, -0.3938385063, -0.2134485135),
    (2427545.0, -0.2299142475, -0.6331803015, -0.2701579257),
    (2431545.0, -0.7188539175, -0.0162132431, 0.0382334294),
    (2435545.0, -0.2164924053, 0.6199261557, 0.2925136930),
    (2439545.0, 0.5874139778, 0.3983146995, 0.1419598205),
    (2443545.0, 0.5935137101, -0.3693448084, -0.2037088644),
    (2447545.0, -0.1985961433, -0.6413855038, -0.2759551516),
    (astro_core::constants::J2000_JD, -0.7183022964, -0.0462742464, 0.0246406381),
];

#[test]
fn vsop2013_vs_reference() {
    let venus = Vsop2013Venus;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in VENUS_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = venus.heliocentric_position(&tdb).unwrap();

        let dx = pos[0] - x_exp;
        let dy = pos[1] - y_exp;
        let dz = pos[2] - z_exp;
        let error_km = (dx * dx + dy * dy + dz * dz).sqrt() * AU_KM;

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
    let venus = Vsop2013Venus;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = venus.heliocentric_position(&tdb).unwrap();

    let expected = (-0.7183022964, -0.0462742464, 0.0246406381);
    let dx = pos[0] - expected.0;
    let dy = pos[1] - expected.1;
    let dz = pos[2] - expected.2;
    let error_km = (dx * dx + dy * dy + dz * dz).sqrt() * AU_KM;

    assert!(
        error_km < 2_000.0,
        "Error {:.0} km exceeds threshold",
        error_km
    );
}

#[test]
fn geocentric_distance_range() {
    let venus = Vsop2013Venus;
    let start_jd = astro_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 8).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = venus.geocentric_position(&tdb).unwrap();
        let dist = (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z).sqrt();
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 0.26,
        "Min distance {:.2} AU below expected 0.26 AU",
        min_dist
    );
    assert!(
        max_dist <= 1.74,
        "Max distance {:.2} AU above expected 1.74 AU",
        max_dist
    );
}
