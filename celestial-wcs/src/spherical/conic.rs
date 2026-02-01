use celestial_core::constants::{DEG_TO_RAD, HALF_PI};
use celestial_core::Angle;

use crate::common::{
    check_nonzero_param, deproject_conic_polar, native_coord_from_radians, project_conic_xy,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_cop(native: NativeCoord, theta_a_deg: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COP projection: theta_a")?;

    let (sigma_s, sigma_c) = theta_a.sin_cos();
    let c = sigma_s;

    if theta.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "COP projection: singularity at theta = 0",
        ));
    }

    let (theta_s, theta_c) = theta.sin_cos();
    let r_theta = sigma_s * theta_c / theta_s;
    let y0 = sigma_c; // sigma / tan(theta_a) = sin(theta_a) * cos(theta_a) / sin(theta_a) = cos(theta_a)

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_cop(inter: IntermediateCoord, theta_a_deg: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COP projection: theta_a")?;

    let (sigma_s, sigma_c) = theta_a.sin_cos();
    let y0 = sigma_c; // sigma / tan(theta_a) = cos(theta_a)

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a);

    if r_unsigned < 1e-15 {
        return Ok(NativeCoord::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0 * theta_a.signum()),
        ));
    }

    let theta = (sigma_s.abs() / r_unsigned).atan() * theta_a.signum();

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_coe(native: NativeCoord, theta_a_deg: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COE projection: theta_a")?;

    let sin_theta_a = theta_a.sin();
    let sin_theta = theta.sin();

    let gamma = sin_theta_a * (2.0 / (1.0 + sin_theta_a * sin_theta_a)).sqrt();
    let c = gamma;

    let s = 1.0 + sin_theta_a;
    let r_theta_a_sq = 2.0 * s * (1.0 - sin_theta_a) / (gamma * gamma);
    let r_theta_sq = r_theta_a_sq + 2.0 * (sin_theta_a - sin_theta) / (gamma * gamma);

    if r_theta_sq < 0.0 {
        return Err(WcsError::out_of_bounds(
            "COE projection: point outside valid region",
        ));
    }

    let r_theta = r_theta_sq.sqrt();
    let y0 = r_theta_a_sq.sqrt();

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_coe(inter: IntermediateCoord, theta_a_deg: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COE projection: theta_a")?;

    let sin_theta_a = theta_a.sin();
    let gamma = sin_theta_a * (2.0 / (1.0 + sin_theta_a * sin_theta_a)).sqrt();
    let c = gamma;

    let s = 1.0 + sin_theta_a;
    let r_theta_a_sq = 2.0 * s * (1.0 - sin_theta_a) / (gamma * gamma);
    let y0 = r_theta_a_sq.sqrt();

    let y_offset = y0 - y;
    let r_sq = x * x + y_offset * y_offset;

    let phi = x.atan2(y_offset) / c;

    let sin_theta = sin_theta_a - 0.5 * gamma * gamma * (r_sq - r_theta_a_sq);

    if sin_theta.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "COE deprojection: point outside valid region",
        ));
    }

    let theta = sin_theta.asin();

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_cod(native: NativeCoord, theta_a_deg: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COD projection: theta_a")?;

    let sigma = theta_a.sin();
    let c = sigma;

    let r_theta_a = theta_a / sigma;
    let r_theta = r_theta_a + theta_a - theta;
    let y0 = r_theta_a;

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_cod(inter: IntermediateCoord, theta_a_deg: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COD projection: theta_a")?;

    let sigma = theta_a.sin();
    let r_theta_a = theta_a / sigma;
    let y0 = r_theta_a;

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a);
    let r = theta_a.signum() * r_unsigned;

    let theta = r_theta_a + theta_a - r;

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_coo(native: NativeCoord, theta_a_deg: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COO projection: theta_a")?;

    if theta.abs() >= HALF_PI - 1e-10 && theta.signum() != theta_a.signum() {
        return Err(WcsError::singularity(
            "COO projection: singularity at opposite pole",
        ));
    }

    let sigma = theta_a.sin();
    let c = sigma;

    let tan_half_theta_a = ((HALF_PI - theta_a) / 2.0).tan();
    if tan_half_theta_a.abs() < 1e-15 {
        return Err(WcsError::singularity(
            "COO projection: singularity at theta_a = 90",
        ));
    }

    let psi = 1.0 / (sigma * tan_half_theta_a.powf(sigma));

    let tan_half_theta = ((HALF_PI - theta) / 2.0).tan();
    let r_theta = psi * tan_half_theta.powf(sigma);
    let y0 = psi * tan_half_theta_a.powf(sigma);

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_coo(inter: IntermediateCoord, theta_a_deg: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COO projection: theta_a")?;

    let sigma = theta_a.sin();

    let tan_half_theta_a = ((HALF_PI - theta_a) / 2.0).tan();
    if tan_half_theta_a.abs() < 1e-15 {
        return Err(WcsError::singularity(
            "COO projection: singularity at theta_a = 90",
        ));
    }

    let psi = 1.0 / (sigma * tan_half_theta_a.powf(sigma));
    let y0 = psi * tan_half_theta_a.powf(sigma);

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a);

    if r_unsigned < 1e-15 {
        return Ok(NativeCoord::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0 * theta_a.signum()),
        ));
    }

    let tan_half_theta = (r_unsigned / psi.abs()).powf(1.0 / sigma.abs());
    let theta = theta_a.signum() * (HALF_PI - 2.0 * tan_half_theta.atan());

    Ok(native_coord_from_radians(phi, theta))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;

    #[test]
    fn test_cop_reference_point() {
        let theta_a = 45.0;
        let proj = Projection::cop(theta_a);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_a));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_cop_native_reference() {
        let theta_a = 45.0;
        let proj = Projection::cop(theta_a);
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, theta_a);
    }

    #[test]
    fn test_cop_roundtrip() {
        let proj = Projection::cop(45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_cop_roundtrip_various_theta_a() {
        for theta_a in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::cop(theta_a);
            let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(50.0));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_cop_roundtrip_various_angles() {
        let proj = Projection::cop(45.0);
        for phi_deg in [-120.0, -60.0, 0.0, 60.0, 120.0] {
            for theta_deg in [20.0, 40.0, 60.0, 80.0] {
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
    fn test_cop_singularity_theta_zero() {
        let proj = Projection::cop(45.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_cop_theta_a_zero() {
        let proj = Projection::cop(0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_cop_wide_latitude_range() {
        let proj = Projection::cop(45.0);
        for theta_deg in [10.0, 20.0, 30.0, 50.0, 70.0, 85.0] {
            let original =
                NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(theta_deg));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_coe_reference_point() {
        let theta_a = 45.0;
        let proj = Projection::coe(theta_a);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_a));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_coe_native_reference() {
        let theta_a = 45.0;
        let proj = Projection::coe(theta_a);
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, theta_a);
    }

    #[test]
    fn test_coe_roundtrip() {
        let proj = Projection::coe(45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_coe_roundtrip_various_theta_a() {
        for theta_a in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::coe(theta_a);
            let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(50.0));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_coe_roundtrip_various_angles() {
        let proj = Projection::coe(45.0);
        for phi_deg in [-120.0, -60.0, 0.0, 60.0, 120.0] {
            for theta_deg in [20.0, 40.0, 60.0, 80.0] {
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
    fn test_coe_theta_a_zero() {
        let proj = Projection::coe(0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_coe_wide_latitude_range() {
        let proj = Projection::coe(45.0);
        for theta_deg in [10.0, 20.0, 30.0, 50.0, 70.0, 85.0] {
            let original =
                NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(theta_deg));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_cod_reference_point() {
        let theta_a = 45.0;
        let proj = Projection::cod(theta_a);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_a));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_cod_native_reference() {
        let theta_a = 45.0;
        let proj = Projection::cod(theta_a);
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, theta_a);
    }

    #[test]
    fn test_cod_roundtrip() {
        let proj = Projection::cod(45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_cod_roundtrip_various_theta_a() {
        for theta_a in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::cod(theta_a);
            let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(50.0));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_cod_roundtrip_various_angles() {
        let proj = Projection::cod(45.0);
        for phi_deg in [-120.0, -60.0, 0.0, 60.0, 120.0] {
            for theta_deg in [20.0, 40.0, 60.0, 80.0] {
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
    fn test_cod_theta_a_zero() {
        let proj = Projection::cod(0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_cod_wide_latitude_range() {
        let proj = Projection::cod(45.0);
        for theta_deg in [10.0, 20.0, 30.0, 50.0, 70.0, 85.0] {
            let original =
                NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(theta_deg));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_coo_reference_point() {
        let theta_a = 45.0;
        let proj = Projection::coo(theta_a);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_a));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_coo_native_reference() {
        let theta_a = 45.0;
        let proj = Projection::coo(theta_a);
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, theta_a);
    }

    #[test]
    fn test_coo_roundtrip() {
        let proj = Projection::coo(45.0);
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_coo_roundtrip_various_theta_a() {
        for theta_a in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::coo(theta_a);
            let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(50.0));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_coo_roundtrip_various_angles() {
        let proj = Projection::coo(45.0);
        for phi_deg in [-120.0, -60.0, 0.0, 60.0, 120.0] {
            for theta_deg in [20.0, 40.0, 60.0, 80.0] {
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
    fn test_coo_theta_a_zero() {
        let proj = Projection::coo(0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_coo_wide_latitude_range() {
        let proj = Projection::coo(45.0);
        for theta_deg in [10.0, 20.0, 30.0, 50.0, 70.0, 85.0] {
            let original =
                NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(theta_deg));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 15);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 15);
        }
    }

    #[test]
    fn test_coo_singularity_opposite_pole() {
        let proj = Projection::coo(45.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        let result = proj.project(native);
        assert!(result.is_err());
    }

    #[test]
    fn test_conic_projections_native_reference() {
        let theta_a = 45.0;
        let projections = [
            Projection::cop(theta_a),
            Projection::coe(theta_a),
            Projection::cod(theta_a),
            Projection::coo(theta_a),
        ];

        for proj in projections {
            let (phi0, theta0) = proj.native_reference();
            assert_eq!(phi0, 0.0);
            assert_eq!(theta0, theta_a);
        }
    }

    #[test]
    fn test_conic_projections_reference_maps_to_origin() {
        let theta_a = 45.0;
        let projections = [
            Projection::cop(theta_a),
            Projection::coe(theta_a),
            Projection::cod(theta_a),
            Projection::coo(theta_a),
        ];

        for proj in &projections {
            let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_a));
            let inter = proj.project(native).unwrap();
            assert!(inter.x_deg().abs() < 1e-10, "x not zero for {:?}", proj);
            assert!(inter.y_deg().abs() < 1e-10, "y not zero for {:?}", proj);
        }
    }
}
