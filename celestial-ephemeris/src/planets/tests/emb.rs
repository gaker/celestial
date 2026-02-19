use crate::planets::Vsop2013Emb;
use celestial_core::constants::{AU_KM, J2000_JD};
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

// VSOP2013.ctl reference values (ICRS, from official Fortran output)
// Line 3 for each epoch: Equatorial Heliocentric Coordinates X,Y,Z (au) - ICRS Frame J2000
const EMB_VSOP2013_REF: &[(f64, f64, f64, f64)] = &[
    (2411545.0, 0.1117527004, -0.9270100498, -0.4021802015),
    (2415545.0, -0.1884496475, -0.9153016306, -0.3970809941),
    (2419545.0, -0.4717744111, -0.8220227738, -0.3565918417),
    (2423545.0, -0.7127406201, -0.6548963219, -0.2840791130),
    (2427545.0, -0.8889625915, -0.4282316398, -0.1857477550),
    (2431545.0, -0.9832170198, -0.1622279919, -0.0703684529),
    (2435545.0, -0.9854629414, 0.1188731342, 0.0515417480),
    (2439545.0, -0.8941873500, 0.3887273239, 0.1685624948),
    (2443545.0, -0.7169612488, 0.6210814225, 0.2693037488),
    (2447545.0, -0.4700594621, 0.7930197220, 0.3438366739),
    (2451545.0, -0.1771587839, 0.8874068590, 0.3847367185),
];

#[test]
fn vsop2013_vs_reference() {
    let emb = Vsop2013Emb;

    for (jd, x_exp, y_exp, z_exp) in EMB_VSOP2013_REF.iter() {
        let tdb = TDB::from_julian_date(JulianDate::new(*jd, 0.0));
        let pos = emb.heliocentric_position(&tdb).unwrap();

        let dx = pos[0] - x_exp;
        let dy = pos[1] - y_exp;
        let dz = pos[2] - z_exp;
        let error_au = libm::sqrt(dx * dx + dy * dy + dz * dz);
        let error_km = error_au * AU_KM;

        assert!(
            error_km < 75.0,
            "JD {}: error {:.0} km exceeds 75 km threshold",
            jd,
            error_km
        );
    }
}

#[test]
fn vsop2013_j2000() {
    let emb = Vsop2013Emb;
    let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

    let pos = emb.heliocentric_position(&tdb).unwrap();

    // VSOP2013.ctl reference for J2000 (ICRS)
    let expected = (-0.1771587839, 0.8874068590, 0.3847367185);

    let dx = pos[0] - expected.0;
    let dy = pos[1] - expected.1;
    let dz = pos[2] - expected.2;
    let error_km = libm::sqrt(dx * dx + dy * dy + dz * dz) * AU_KM;
    assert!(
        error_km < 75.0,
        "Error {:.0} km exceeds 75 km threshold",
        error_km
    );
}
