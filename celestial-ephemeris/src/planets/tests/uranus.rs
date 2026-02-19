use crate::planets::Vsop2013Uranus;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

const URANUS_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, -16.4159097008, -7.7936503250, -3.1804126500),
    (2415545.0, -4.5170400561, -17.0119258599, -7.3868736325),
    (2419545.0, 10.4746921789, -15.3005867330, -6.8500911185),
    (2423545.0, 19.4279722859, -4.5746255506, -2.2790261954),
    (2427545.0, 17.5128800292, 8.7263262152, 3.5736558506),
    (2431545.0, 5.5570811634, 16.9308451272, 7.3365906981),
    (2435545.0, -9.7178358757, 14.4546349234, 6.4686670469),
    (2439545.0, -18.1327079367, 2.0613281354, 1.1596916860),
    (2443545.0, -13.4822148842, -11.8237898401, -4.9874755478),
    (2447545.0, 0.5625964390, -17.6846990306, -7.7533189975),
    (2451545.0, 14.4318565807, -12.5062632452, -5.6816829828),
];

#[test]
fn vsop2013_vs_reference() {
    let uranus = Vsop2013Uranus;
    let mut max_error_km = 0.0;

    for (jd, x_exp, y_exp, z_exp) in URANUS_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = uranus.heliocentric_position(&tdb).unwrap();

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
    let uranus = Vsop2013Uranus;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
    let pos = uranus.heliocentric_position(&tdb).unwrap();

    let expected = (14.4318565807, -12.5062632452, -5.6816829828);
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
    let uranus = Vsop2013Uranus;
    let start_jd = celestial_core::constants::J2000_JD;

    let (min_dist, max_dist) = (0..12 * 20).fold((f64::MAX, 0.0f64), |(min, max), i| {
        let jd = start_jd + (i * 30) as f64;
        let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
        let pos = uranus.geocentric_position(&tdb).unwrap();
        let dist = libm::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
        (min.min(dist), max.max(dist))
    });

    assert!(
        min_dist >= 17.2,
        "Min distance {:.2} AU below expected 17.2 AU",
        min_dist
    );
    assert!(
        max_dist <= 21.1,
        "Max distance {:.2} AU above expected 21.1 AU",
        max_dist
    );
}
