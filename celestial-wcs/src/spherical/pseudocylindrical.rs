use celestial_core::constants::{DEG_TO_RAD, HALF_PI, PI, RAD_TO_DEG, SQRT2};

use crate::common::{asin_safe, native_coord_from_radians, newton_raphson_1d, NewtonConfig};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_sfl(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let x = phi * libm::cos(theta) * RAD_TO_DEG;
    let y = theta * RAD_TO_DEG;
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_sfl(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;

    let theta = y;

    let cos_theta = libm::cos(theta);
    if cos_theta.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "SFL deprojection: singularity at theta = +/-90",
        ));
    }

    let phi = x / cos_theta;

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_par(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let scale = 2.0 * libm::cos(2.0 * theta / 3.0) - 1.0;
    let x = phi * scale * RAD_TO_DEG;
    let y = 180.0 * libm::sin(theta / 3.0);
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_par(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg();
    let y = inter.y_deg();

    let sin_theta_3 = y / 180.0;
    if sin_theta_3.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "PAR deprojection: |y| > 180 degrees",
        ));
    }

    let theta_3 = libm::asin(sin_theta_3);
    let theta = 3.0 * theta_3;

    let scale = 2.0 * libm::cos(2.0 * theta / 3.0) - 1.0;
    if scale.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "PAR deprojection: singularity at theta = +/-90",
        ));
    }

    let phi = x * DEG_TO_RAD / scale;

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_mol(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let gamma = solve_mollweide_gamma(theta)?;

    let sqrt_8_over_pi = libm::sqrt(8.0_f64) / PI;
    let (gamma_sin, gamma_cos) = libm::sincos(gamma);
    let x = sqrt_8_over_pi * phi * gamma_cos * RAD_TO_DEG;
    let y = SQRT2 * 90.0 * gamma_sin;
    Ok(IntermediateCoord::new(x, y))
}

fn solve_mollweide_gamma(theta: f64) -> WcsResult<f64> {
    if theta.abs() >= HALF_PI - 1e-10 {
        return Ok(theta.signum() * HALF_PI);
    }

    let pi_sin_theta = PI * libm::sin(theta);

    const CONFIG: NewtonConfig = NewtonConfig::new((-HALF_PI, HALF_PI), "MOL forward");
    newton_raphson_1d(
        theta,
        pi_sin_theta,
        |gamma| 2.0 * gamma + libm::sin(2.0 * gamma),
        |gamma| 2.0 + 2.0 * libm::cos(2.0 * gamma),
        &CONFIG,
    )
}

pub(crate) fn deproject_mol(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg();
    let y = inter.y_deg();

    let sqrt_2_times_90 = SQRT2 * 90.0;
    let sin_gamma = y / sqrt_2_times_90;

    if sin_gamma.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "MOL deprojection: point outside projection boundary",
        ));
    }

    let gamma = libm::asin(sin_gamma);
    let cos_gamma = libm::cos(gamma);

    if cos_gamma.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "MOL deprojection: singularity at gamma = +/-90",
        ));
    }

    let sin_theta = (2.0 * gamma + libm::sin(2.0 * gamma)) / PI;
    let theta = asin_safe(sin_theta);

    let sqrt_8_over_pi = libm::sqrt(8.0_f64) / PI;
    let phi = x * DEG_TO_RAD / (sqrt_8_over_pi * cos_gamma);

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_ait(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let (sin_theta, cos_theta) = libm::sincos(theta);
    let half_phi = phi / 2.0;
    let cos_theta_cos_half_phi = cos_theta * libm::cos(half_phi);

    let denom = 1.0 + cos_theta_cos_half_phi;
    if denom < 1e-10 {
        return Err(WcsError::singularity(
            "AIT projection: singularity at antipodal point",
        ));
    }

    let gamma = libm::sqrt(2.0 / denom);
    let x = 2.0 * gamma * cos_theta * libm::sin(half_phi) * RAD_TO_DEG;
    let y = gamma * sin_theta * RAD_TO_DEG;
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_ait(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;

    let x_scaled = x / 4.0;
    let y_scaled = y / 2.0;

    let z_sq = 1.0 - x_scaled * x_scaled - y_scaled * y_scaled;
    if z_sq < 0.0 {
        return Err(WcsError::out_of_bounds(
            "AIT deprojection: point outside projection boundary",
        ));
    }

    let z = libm::sqrt(z_sq);

    let sin_theta = y * z;
    let theta = asin_safe(sin_theta);

    let phi = 2.0 * libm::atan2(x * z / 2.0, 2.0 * z * z - 1.0);

    Ok(native_coord_from_radians(phi, theta))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;

    #[test]
    fn test_sfl_reference_point() {
        let proj = Projection::sfl();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), 0.0);
    }

    #[test]
    fn test_sfl_native_reference() {
        let proj = Projection::sfl();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_sfl_roundtrip() {
        let proj = Projection::sfl();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_sfl_roundtrip_various_angles() {
        let proj = Projection::sfl();
        for phi_deg in [-180.0, -90.0, -45.0, 0.0, 45.0, 90.0, 180.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 5);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 5);
            }
        }
    }

    #[test]
    fn test_sfl_known_value() {
        let proj = Projection::sfl();
        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!((inter.x_deg() - 90.0).abs() < 1e-10);
        assert!((inter.y_deg() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_sfl_singularity_deproject() {
        let proj = Projection::sfl();
        let inter = IntermediateCoord::new(10.0, 90.0);
        let result = proj.deproject(inter);
        assert!(result.is_err());
    }

    #[test]
    fn test_par_reference_point() {
        let proj = Projection::par();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), 0.0);
    }

    #[test]
    fn test_par_native_reference() {
        let proj = Projection::par();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_par_roundtrip() {
        let proj = Projection::par();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_par_roundtrip_various_angles() {
        let proj = Projection::par();
        for phi_deg in [-180.0, -90.0, -45.0, 0.0, 45.0, 90.0, 180.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 5);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 5);
            }
        }
    }

    #[test]
    fn test_par_boundary_y() {
        let proj = Projection::par();
        let inter = IntermediateCoord::new(0.0, 200.0);
        let result = proj.deproject(inter);
        assert!(result.is_err());
    }

    #[test]
    fn test_mol_reference_point() {
        let proj = Projection::mol();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!((inter.x_deg()).abs() < 1e-10);
        assert!((inter.y_deg()).abs() < 1e-10);
    }

    #[test]
    fn test_mol_native_reference() {
        let proj = Projection::mol();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_mol_roundtrip() {
        let proj = Projection::mol();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_mol_roundtrip_various_angles() {
        let proj = Projection::mol();
        for phi_deg in [-180.0, -90.0, -45.0, 0.0, 45.0, 90.0, 180.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
            }
        }
    }

    #[test]
    fn test_mol_poles() {
        let proj = Projection::mol();
        let north_pole = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let inter = proj.project(north_pole).unwrap();
        assert!((inter.x_deg()).abs() < 1e-10);
        let expected_y = std::f64::consts::SQRT_2 * 90.0;
        assert!((inter.y_deg() - expected_y).abs() < 1e-10);
    }

    #[test]
    fn test_mol_boundary() {
        let proj = Projection::mol();
        let sqrt_2_times_90 = std::f64::consts::SQRT_2 * 90.0;
        let inter = IntermediateCoord::new(0.0, sqrt_2_times_90 + 10.0);
        let result = proj.deproject(inter);
        assert!(result.is_err());
    }

    #[test]
    fn test_ait_reference_point() {
        let proj = Projection::ait();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!((inter.x_deg()).abs() < 1e-10);
        assert!((inter.y_deg()).abs() < 1e-10);
    }

    #[test]
    fn test_ait_native_reference() {
        let proj = Projection::ait();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_ait_roundtrip() {
        let proj = Projection::ait();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_ait_roundtrip_various_angles() {
        let proj = Projection::ait();
        for phi_deg in [-150.0, -90.0, -45.0, 0.0, 45.0, 90.0, 150.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
            }
        }
    }

    #[test]
    fn test_ait_poles() {
        let proj = Projection::ait();
        let north_pole = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let inter = proj.project(north_pole).unwrap();
        assert!((inter.x_deg()).abs() < 1e-10);
        let expected_y = std::f64::consts::SQRT_2 * RAD_TO_DEG;
        assert!((inter.y_deg() - expected_y).abs() < 1e-8);
    }

    #[test]
    fn test_ait_boundary() {
        let proj = Projection::ait();
        let inter = IntermediateCoord::new(400.0, 200.0);
        let result = proj.deproject(inter);
        assert!(result.is_err());
    }

    #[test]
    fn test_ait_equator() {
        let proj = Projection::ait();
        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(native.phi().degrees(), recovered.phi().degrees(), 1);
        assert_ulp_lt!(native.theta().degrees(), recovered.theta().degrees(), 1);
    }
}
