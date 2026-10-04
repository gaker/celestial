use super::*;

// eraEpv00 at J2000.
#[test]
fn test_earth_velocity_j2000() {
    let state = compute_earth_state(&TT::j2000()).unwrap();
    let velocity = Vector3::new(
        -0.01720224630718366,
        -0.002904925940146081,
        -0.0012594275302390552,
    );
    assert_eq!(state.barycentric_velocity, velocity);
    let position = Vector3::new(-0.17713507281322974, 0.8874285242954301, 0.3847428889988798);
    assert_eq!(state.heliocentric_position, position);
}

#[test]
fn test_aberration_magnitude() {
    let tt = TT::j2000();
    let state = compute_earth_state(&tt).unwrap();

    let p = Vector3::new(1.0, 0.0, 0.0);
    let s = state.heliocentric_position.magnitude();
    // eraAb, with the velocity in units of c and the Lorentz factor from eraApcs.
    assert_eq!(
        apply_aberration(p, state.barycentric_velocity, s),
        Vector3::new(
            0.9999999998327874,
            -1.6778286874271264e-5,
            -7.274208305167048e-6
        )
    );
}

#[test]
fn test_aberration_roundtrip() {
    let tt = TT::j2000();
    let state = compute_earth_state(&tt).unwrap();
    let sun_dist = state.heliocentric_position.magnitude();

    let directions = [
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(1.0, 1.0, 1.0).normalize().unwrap(),
    ];

    for dir in directions {
        let aberrated = apply_aberration(dir, state.barycentric_velocity, sun_dist);
        let recovered = remove_aberration(aberrated, state.barycentric_velocity, sun_dist).unwrap();

        // The inverse stops after a fixed number of iterations, as eraAticq does, which
        // leaves up to about 4e-13.
        let diff = (dir - recovered).magnitude();
        assert!(diff < 1e-12, "Aberration roundtrip error: {:.2e}", diff);
    }
}

#[test]
fn test_remove_aberration_is_inverse() {
    let velocity = Vector3::new(0.01, 0.005, 0.002);
    let sun_dist = 1.0;
    let original = Vector3::new(0.6, 0.7, 0.3).normalize().unwrap();

    let aberrated = apply_aberration(original, velocity, sun_dist);
    let recovered = remove_aberration(aberrated, velocity, sun_dist).unwrap();

    let diff = (original - recovered).magnitude();
    assert!(
        diff < 1e-12,
        "Aberration inverse should be within iterative precision: {:.2e}",
        diff
    );
}

fn spring_equinox_tt() -> TT {
    TT::from_julian_date(celestial_time::julian::JulianDate::from_f64(2460389.5))
}

fn arcsec_diff(a: Angle, b: Angle) -> f64 {
    libm::fabs((a - b).wrapped().unwrap().arcseconds())
}

#[test]
fn annual_aberration_shift_is_on_order_of_twenty_arcsec() {
    let tt = spring_equinox_tt();
    let ra = Angle::from_hours(0.0);
    let dec = Angle::from_degrees(0.0);
    let (ra2, dec2) = apply_annual_aberration(ra, dec, &tt).unwrap();
    let shift = arcsec_diff(ra, ra2) + arcsec_diff(dec, dec2);
    assert!(
        shift > 5.0 && shift < 45.0,
        "expected annual aberration ~20 arcsec, got {}",
        shift
    );
}

#[test]
fn annual_aberration_apply_and_remove_are_inverses() {
    let tt = spring_equinox_tt();
    let ra = Angle::from_hours(12.0);
    let dec = Angle::from_degrees(30.0);
    let (ra2, dec2) = apply_annual_aberration(ra, dec, &tt).unwrap();
    let (ra3, dec3) = remove_annual_aberration(ra2, dec2, &tt).unwrap();
    // The inverse leaves about 3e-14 rad here. The bound is the 1e-12 rad used for the
    // vectors above, 2e-7 arcsec.
    assert!(
        arcsec_diff(ra, ra3) < 2e-7,
        "ra drift {}",
        arcsec_diff(ra, ra3)
    );
    assert!(
        arcsec_diff(dec, dec3) < 2e-7,
        "dec drift {}",
        arcsec_diff(dec, dec3)
    );
}

#[test]
fn annual_aberration_varies_by_sky_position() {
    let tt = spring_equinox_tt();
    let dec = Angle::from_degrees(0.0);
    let (ra_a, dec_a) = apply_annual_aberration(Angle::from_hours(0.0), dec, &tt).unwrap();
    let (ra_b, dec_b) = apply_annual_aberration(Angle::from_hours(12.0), dec, &tt).unwrap();
    let shift_a = arcsec_diff(Angle::from_hours(0.0), ra_a) + arcsec_diff(dec, dec_a);
    let shift_b = arcsec_diff(Angle::from_hours(12.0), ra_b) + arcsec_diff(dec, dec_b);
    assert_ne!(shift_a, shift_b);
}

#[test]
fn compute_earth_state_rejects_non_finite_epoch() {
    for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let tt = TT::from_julian_date(celestial_time::julian::JulianDate::new(jd, 0.0));
        assert!(compute_earth_state(&tt).is_err(), "JD {jd}");
    }
}

#[test]
fn annual_aberration_rejects_non_finite_place() {
    let tt = spring_equinox_tt();
    let places = [(f64::NAN, 0.0), (0.0, f64::NAN), (f64::INFINITY, 0.0)];
    for (ra, dec) in places {
        let (ra, dec) = (Angle::from_radians(ra), Angle::from_radians(dec));
        assert!(
            apply_annual_aberration(ra, dec, &tt).is_err(),
            "{ra:?} {dec:?}"
        );
        assert!(
            remove_annual_aberration(ra, dec, &tt).is_err(),
            "{ra:?} {dec:?}"
        );
    }
}
