use celestial_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};
use celestial_core::Angle;

use crate::common::{check_nonzero_param, native_coord_from_radians};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_bon(native: NativeCoord, theta_1_deg: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_1 = theta_1_deg * DEG_TO_RAD;

    check_nonzero_param(theta_1, "BON projection: theta_1")?;

    let cot_theta_1 = libm::cos(theta_1) / libm::sin(theta_1);
    let y0 = cot_theta_1 + theta_1;
    let r_theta = y0 - theta;

    let a = phi * libm::cos(theta) / r_theta;
    let (a_sin, a_cos) = libm::sincos(a);

    let x = r_theta * a_sin * RAD_TO_DEG;
    let y = (y0 - r_theta * a_cos) * RAD_TO_DEG;

    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_bon(inter: IntermediateCoord, theta_1_deg: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_1 = theta_1_deg * DEG_TO_RAD;

    check_nonzero_param(theta_1, "BON projection: theta_1")?;

    let cot_theta_1 = libm::cos(theta_1) / libm::sin(theta_1);
    let y0 = cot_theta_1 + theta_1;
    let y_offset = y0 - y;

    let r_unsigned = libm::sqrt(x * x + y_offset * y_offset);
    let r = theta_1.signum() * r_unsigned;

    let theta = y0 - r;

    let a = libm::atan2(theta_1.signum() * x, theta_1.signum() * y_offset);

    let cos_theta = libm::cos(theta);
    if cos_theta.abs() < 1e-15 {
        return Err(WcsError::singularity(
            "BON deprojection: cos(theta) ~ 0 (theta near +/-90)",
        ));
    }
    let phi = a * r / cos_theta;

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_pco(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta.abs() < 1e-10 {
        let x = phi * RAD_TO_DEG;
        let y = 0.0;
        return Ok(IntermediateCoord::new(x, y));
    }

    let (sin_theta, cos_theta) = libm::sincos(theta);
    let tan_theta = sin_theta / cos_theta;
    let e = phi * sin_theta;
    let (e_sin, e_cos) = libm::sincos(e);

    let x = e_sin / tan_theta * RAD_TO_DEG;
    let y = (theta + (1.0 - e_cos) / tan_theta) * RAD_TO_DEG;

    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_pco(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;

    if y.abs() < 1e-10 && x.abs() < 1e-10 {
        return Ok(NativeCoord::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
        ));
    }

    if y.abs() < 1e-10 {
        return Ok(native_coord_from_radians(x, 0.0));
    }

    let theta = solve_pco_inverse(x, y)?;

    if theta.abs() < 1e-10 {
        return Ok(native_coord_from_radians(x, 0.0));
    }

    let sin_theta = libm::sin(theta);
    let tan_theta = libm::tan(theta);

    let sin_e = x * tan_theta;
    if sin_e.abs() > 1.0 {
        return Err(WcsError::out_of_bounds("PCO deprojection: |sin(E)| > 1"));
    }

    let e = libm::asin(sin_e);
    let phi = e / sin_theta;

    Ok(native_coord_from_radians(phi, theta))
}

fn solve_pco_inverse(x: f64, y: f64) -> WcsResult<f64> {
    const MAX_ITER: usize = 100;
    const TOL: f64 = 1e-12;

    let mut theta = y;

    for _ in 0..MAX_ITER {
        if theta.abs() < 1e-10 {
            return Ok(0.0);
        }

        let (sin_theta, cos_theta) = libm::sincos(theta);

        if cos_theta.abs() < 1e-15 {
            return Err(WcsError::singularity(
                "PCO inverse: singularity at theta = +/-90",
            ));
        }

        let tan_theta = sin_theta / cos_theta;
        let sin_e = x * tan_theta;

        if sin_e.abs() > 1.0 {
            theta *= 0.9;
            continue;
        }

        let e = libm::asin(sin_e);
        let (e_sin, e_cos) = libm::sincos(e);

        let f = theta + (1.0 - e_cos) / tan_theta - y;

        let de_dtheta = x / (cos_theta * cos_theta * libm::sqrt(1.0 - sin_e * sin_e).max(1e-15));
        let d_cos_e_dtheta = -e_sin * de_dtheta;

        let df_dtheta = 1.0 - d_cos_e_dtheta / tan_theta - (1.0 - e_cos) / (sin_theta * sin_theta);

        if df_dtheta.abs() < 1e-15 {
            theta += 0.01 * y.signum();
            continue;
        }

        let delta = f / df_dtheta;
        theta -= delta;

        theta = theta.clamp(-HALF_PI + 0.01, HALF_PI - 0.01);

        if delta.abs() < TOL {
            return Ok(theta);
        }
    }

    Err(WcsError::convergence_failure(
        "PCO inverse: Newton-Raphson did not converge",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;

    // Per-projection native_reference checks are covered by
    // spherical::tests::test_all_projections_map_reference_to_origin.

    #[test]
    fn test_bon_reference_maps_to_origin() {
        // Spec-anchored: per Paper II Eq. (145), (phi=0, theta=0) must map to
        // (x=0, y=0) for every theta_1.  Forward at the reference is a direct
        // check on Y_0 = (180/pi) cot(theta_1) + theta_1 (Eq. 150).
        for theta_1 in [30.0, 45.0, 60.0, 75.0, -45.0] {
            let proj = Projection::bon(theta_1);
            let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
            let inter = proj.project(native).unwrap();
            assert!(
                inter.x_deg().abs() < 1e-10,
                "x not zero for BON theta_1={}: {}", theta_1, inter.x_deg(),
            );
            assert!(
                inter.y_deg().abs() < 1e-10,
                "y not zero for BON theta_1={}: {}", theta_1, inter.y_deg(),
            );
        }
    }

    #[test]
    fn test_bon_roundtrip_sweep() {
        // Sweep theta_1 (including negative) x phi x theta in one shot.
        for theta_1 in [-45.0, 30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::bon(theta_1);
            // Mirror theta sign for negative theta_1 so we stay on the valid side.
            let thetas: &[f64] = if theta_1 > 0.0 {
                &[0.0, 10.0, 20.0, 30.0, 50.0, 70.0]
            } else {
                &[0.0, -10.0, -20.0, -30.0, -50.0, -70.0]
            };
            for &phi_deg in &[-60.0, -30.0, 0.0, 30.0, 60.0] {
                for &theta_deg in thetas {
                    let original = NativeCoord::new(
                        Angle::from_degrees(phi_deg),
                        Angle::from_degrees(theta_deg),
                    );
                    let inter = proj.project(original).unwrap();
                    let recovered = proj.deproject(inter).unwrap();
                    // The theta=0 equator case has a separate code path that
                    // doesn't recover ULP-tight; assert the absolute floor.
                    assert!(
                        (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-8,
                        "phi (theta_1={}, {}, {})", theta_1, phi_deg, theta_deg,
                    );
                    assert!(
                        (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-8,
                        "theta (theta_1={}, {}, {})", theta_1, phi_deg, theta_deg,
                    );
                }
            }
        }
    }

    #[test]
    fn test_pco_equator_known_value_and_symmetry() {
        let proj = Projection::pco();

        // On the equator PCO collapses to (phi, 0) - a degree-identity in x.
        let native = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert_ulp_lt!(inter.x_deg(), 45.0, 2);
        assert!(inter.y_deg().abs() < 1e-10);

        // Deprojection round trip on the equator.
        let inter = IntermediateCoord::new(30.0, 0.0);
        let result = proj.deproject(inter).unwrap();
        assert_ulp_lt!(result.phi().degrees(), 30.0, 2);
        assert!(result.theta().degrees().abs() < 1e-10);

        // Symmetry about the central meridian: (-phi, theta) -> (-x, y).
        let native_pos = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(45.0));
        let native_neg = NativeCoord::new(Angle::from_degrees(-30.0), Angle::from_degrees(45.0));
        let inter_pos = proj.project(native_pos).unwrap();
        let inter_neg = proj.project(native_neg).unwrap();
        assert_ulp_lt!(inter_pos.x_deg(), -inter_neg.x_deg(), 2);
        assert_ulp_lt!(inter_pos.y_deg(), inter_neg.y_deg(), 2);
    }

    #[test]
    fn test_pco_roundtrip() {
        let proj = Projection::pco();
        for phi_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
            for theta_deg in [-70.0, -45.0, -20.0, 15.0, 20.0, 30.0, 45.0, 60.0, 70.0] {
                let original = NativeCoord::new(
                    Angle::from_degrees(phi_deg),
                    Angle::from_degrees(theta_deg),
                );
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert!(
                    (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-8,
                    "phi mismatch at ({}, {})", phi_deg, theta_deg,
                );
                assert!(
                    (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-8,
                    "theta mismatch at ({}, {})", phi_deg, theta_deg,
                );
            }
        }
    }

    #[test]
    fn test_polyconic_error_cases() {
        // BON with theta_1 = 0 makes cot(theta_1) undefined.
        let bon = Projection::bon(0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        assert!(bon.project(native).is_err());
    }
}
