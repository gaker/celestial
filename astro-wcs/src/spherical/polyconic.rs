use astro_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};
use astro_core::Angle;

use crate::common::{check_nonzero_param, native_coord_from_radians};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_bon(native: NativeCoord, theta_1_deg: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_1 = theta_1_deg * DEG_TO_RAD;

    check_nonzero_param(theta_1, "BON projection: theta_1")?;

    let (theta_1_s, theta_1_c) = theta_1.sin_cos();
    let cot_theta_1 = theta_1_c / theta_1_s;

    let y0 = cot_theta_1;

    if theta.abs() < 1e-10 {
        let r_theta = cot_theta_1 + theta_1;
        let a = phi * theta_1_s / r_theta;
        let (a_sin, a_cos) = a.sin_cos();
        let x = r_theta * a_sin * RAD_TO_DEG;
        let y = (y0 - r_theta * a_cos) * RAD_TO_DEG;
        return Ok(IntermediateCoord::new(x, y));
    }

    let r_theta = cot_theta_1 + theta_1 - theta;

    let a = phi * theta.sin() / r_theta;
    let (a_sin, a_cos) = a.sin_cos();

    let x = r_theta * a_sin * RAD_TO_DEG;
    let y = (y0 - r_theta * a_cos) * RAD_TO_DEG;

    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_bon(inter: IntermediateCoord, theta_1_deg: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_1 = theta_1_deg * DEG_TO_RAD;

    check_nonzero_param(theta_1, "BON projection: theta_1")?;

    let (theta_1_s, theta_1_c) = theta_1.sin_cos();
    let cot_theta_1 = theta_1_c / theta_1_s;
    let y0 = cot_theta_1;
    let y_offset = y0 - y;

    let r_unsigned = (x * x + y_offset * y_offset).sqrt();
    let r = theta_1.signum() * r_unsigned;

    let theta = cot_theta_1 + theta_1 - r;

    let a = (theta_1.signum() * x).atan2(theta_1.signum() * y_offset);

    if theta.abs() < 1e-10 {
        let r_at_equator = cot_theta_1 + theta_1;
        let phi = a * r_at_equator / theta_1_s;
        return Ok(native_coord_from_radians(phi, theta));
    }

    let phi = a * r / theta.sin();

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

    let (sin_theta, cos_theta) = theta.sin_cos();
    let tan_theta = sin_theta / cos_theta;
    let e = phi * sin_theta;
    let (e_sin, e_cos) = e.sin_cos();

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

    let sin_theta = theta.sin();
    let tan_theta = theta.tan();

    let sin_e = x * tan_theta;
    if sin_e.abs() > 1.0 {
        return Err(WcsError::out_of_bounds("PCO deprojection: |sin(E)| > 1"));
    }

    let e = sin_e.asin();
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

        let (sin_theta, cos_theta) = theta.sin_cos();

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

        let e = sin_e.asin();
        let (e_sin, e_cos) = e.sin_cos();

        let f = theta + (1.0 - e_cos) / tan_theta - y;

        let de_dtheta = x / (cos_theta * cos_theta * (1.0 - sin_e * sin_e).sqrt().max(1e-15));
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
    use astro_core::assert_ulp_lt;
    use astro_core::Angle;

    #[test]
    fn test_bon_reference_point() {
        let theta_1 = 45.0;
        let proj = Projection::bon(theta_1);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_1));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_bon_native_reference() {
        let theta_1 = 45.0;
        let proj = Projection::bon(theta_1);
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, theta_1);
    }

    #[test]
    fn test_bon_roundtrip() {
        let proj = Projection::bon(45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
    }

    #[test]
    fn test_bon_roundtrip_various_theta_1() {
        for theta_1 in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::bon(theta_1);
            let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(50.0));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 20);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 20);
        }
    }

    #[test]
    fn test_bon_roundtrip_various_angles() {
        let proj = Projection::bon(45.0);
        for phi_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
            for theta_deg in [20.0, 40.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 20);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 20);
            }
        }
    }

    #[test]
    fn test_bon_theta_1_zero() {
        let proj = Projection::bon(0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_bon_equator_handling() {
        let proj = Projection::bon(45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(0.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-8,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-8,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_bon_wide_latitude_range() {
        let proj = Projection::bon(45.0);
        for theta_deg in [10.0, 20.0, 30.0, 50.0, 70.0] {
            let original =
                NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(theta_deg));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 20);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 20);
        }
    }

    #[test]
    fn test_bon_negative_theta_1() {
        let proj = Projection::bon(-45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(-60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 20);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 20);
    }

    #[test]
    fn test_pco_reference_point() {
        let proj = Projection::pco();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_pco_native_reference() {
        let proj = Projection::pco();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_pco_equator() {
        let proj = Projection::pco();
        let native = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert_ulp_lt!(inter.x_deg(), 45.0, 2);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_pco_roundtrip() {
        let proj = Projection::pco();
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(45.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 50);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 50);
    }

    #[test]
    fn test_pco_roundtrip_various_angles() {
        let proj = Projection::pco();
        for phi_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
            for theta_deg in [-60.0, -30.0, 15.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert!(
                    (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-8,
                    "phi mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.phi().degrees(),
                    recovered.phi().degrees()
                );
                assert!(
                    (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-8,
                    "theta mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.theta().degrees(),
                    recovered.theta().degrees()
                );
            }
        }
    }

    #[test]
    fn test_pco_deproject_origin() {
        let proj = Projection::pco();
        let inter = IntermediateCoord::new(0.0, 0.0);
        let result = proj.deproject(inter).unwrap();
        assert_eq!(result.phi().degrees(), 0.0);
        assert_eq!(result.theta().degrees(), 0.0);
    }

    #[test]
    fn test_pco_deproject_equator() {
        let proj = Projection::pco();
        let inter = IntermediateCoord::new(30.0, 0.0);
        let result = proj.deproject(inter).unwrap();
        assert_ulp_lt!(result.phi().degrees(), 30.0, 2);
        assert!(result.theta().degrees().abs() < 1e-10);
    }

    #[test]
    fn test_pco_symmetric() {
        let proj = Projection::pco();
        let native_pos = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(45.0));
        let native_neg = NativeCoord::new(Angle::from_degrees(-30.0), Angle::from_degrees(45.0));
        let inter_pos = proj.project(native_pos).unwrap();
        let inter_neg = proj.project(native_neg).unwrap();
        assert_ulp_lt!(inter_pos.x_deg(), -inter_neg.x_deg(), 2);
        assert_ulp_lt!(inter_pos.y_deg(), inter_neg.y_deg(), 2);
    }

    #[test]
    fn test_pco_wide_latitude_range() {
        let proj = Projection::pco();
        for theta_deg in [-70.0, -45.0, -20.0, 20.0, 45.0, 70.0] {
            let original =
                NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(theta_deg));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert!(
                (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-8,
                "phi mismatch at theta={}: {} vs {}",
                theta_deg,
                original.phi().degrees(),
                recovered.phi().degrees()
            );
            assert!(
                (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-8,
                "theta mismatch at theta={}: {} vs {}",
                theta_deg,
                original.theta().degrees(),
                recovered.theta().degrees()
            );
        }
    }

    #[test]
    fn test_polyconic_projections_native_reference() {
        let theta_1 = 45.0;
        let bon = Projection::bon(theta_1);
        let pco = Projection::pco();

        let (bon_phi0, bon_theta0) = bon.native_reference();
        assert_eq!(bon_phi0, 0.0);
        assert_eq!(bon_theta0, theta_1);

        let (pco_phi0, pco_theta0) = pco.native_reference();
        assert_eq!(pco_phi0, 0.0);
        assert_eq!(pco_theta0, 0.0);
    }

    #[test]
    fn test_bon_reference_maps_to_origin() {
        for theta_1 in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::bon(theta_1);
            let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_1));
            let inter = proj.project(native).unwrap();
            assert!(
                inter.x_deg().abs() < 1e-10,
                "x not zero for BON with theta_1={}",
                theta_1
            );
            assert!(
                inter.y_deg().abs() < 1e-10,
                "y not zero for BON with theta_1={}",
                theta_1
            );
        }
    }
}
