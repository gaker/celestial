use crate::planets::Vsop2013Pluto;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const PLUTO_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, 17.6213054610, 43.9412232438, 8.4004533071),
    (2415545.0, 9.1596791660, 44.4875512788, 11.1198669305),
    (2419545.0, 0.2778860239, 42.9818583697, 13.3261124440),
    (2423545.0, -8.6130696510, 39.2370669120, 14.8346938492),
    (2427545.0, -16.9607735628, 33.1345593015, 15.4466786780),
    (2431545.0, -24.0688624312, 24.6550011074, 14.9422382300),
    (2435545.0, -29.0024768171, 13.9948403475, 13.1021316473),
    (2439545.0, -30.6589088179, 1.8198599397, 9.8043229541),
    (2443545.0, -28.0398775393, -10.5511024785, 5.1547952585),
    (2447545.0, -20.7970730811, -21.1416023203, -0.3309844592),
    (2451545.0, -9.8753695808, -27.9789262247, -5.7537118247),
];

#[test]
fn vsop2013_vs_reference() {
    let pluto = Vsop2013Pluto;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in PLUTO_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = pluto.heliocentric_position(&tdb).unwrap();

        let dx = pos[0] - x_exp;
        let dy = pos[1] - y_exp;
        let dz = pos[2] - z_exp;
        let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;

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
    let pluto = Vsop2013Pluto;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = pluto.heliocentric_position(&tdb).unwrap();

    let expected = (-9.8753695808, -27.9789262247, -5.7537118247);
    let dx = pos[0] - expected.0;
    let dy = pos[1] - expected.1;
    let dz = pos[2] - expected.2;
    let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;

    assert!(
        error_km < 20_000.0,
        "Error {:.0} km exceeds threshold",
        error_km
    );
}

#[test]
fn geocentric_distance_range() {
    let pluto = Vsop2013Pluto;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 20).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = pluto.geocentric_position(&tdb).unwrap();
        let dist = libm::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 28.6,
        "Min distance {:.2} AU below expected 28.6 AU",
        min_dist
    );
    assert!(
        max_dist <= 50.5,
        "Max distance {:.2} AU above expected 50.5 AU",
        max_dist
    );
}
