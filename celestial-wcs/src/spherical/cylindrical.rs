use celestial_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};
use celestial_core::Angle;

use crate::common::native_coord_from_radians;
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_car(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().degrees();
    let theta = native.theta().degrees();
    Ok(IntermediateCoord::new(phi, theta))
}

pub(crate) fn deproject_car(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let phi = inter.x_deg();
    let theta = inter.y_deg();
    Ok(NativeCoord::new(
        Angle::from_degrees(phi),
        Angle::from_degrees(theta),
    ))
}

pub(crate) fn project_mer(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().degrees();
    let theta = native.theta().radians();

    if theta.abs() >= HALF_PI - 1e-10 {
        return Err(WcsError::singularity(
            "MER projection undefined at theta = +/-90",
        ));
    }

    let y = (std::f64::consts::FRAC_PI_4 + theta / 2.0).tan().ln() * RAD_TO_DEG;
    Ok(IntermediateCoord::new(phi, y))
}

pub(crate) fn deproject_mer(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let phi = inter.x_deg();
    let y = inter.y_deg() * DEG_TO_RAD;

    let theta = 2.0 * y.exp().atan() - HALF_PI;
    Ok(NativeCoord::new(
        Angle::from_degrees(phi),
        Angle::from_degrees(theta * RAD_TO_DEG),
    ))
}

pub(crate) fn project_cea(native: NativeCoord, lambda: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().degrees();
    let theta = native.theta().radians();

    let y = theta.sin() / lambda * RAD_TO_DEG;
    Ok(IntermediateCoord::new(phi, y))
}

pub(crate) fn deproject_cea(inter: IntermediateCoord, lambda: f64) -> WcsResult<NativeCoord> {
    let phi = inter.x_deg();
    let y = inter.y_deg() * DEG_TO_RAD;

    let sin_theta = lambda * y;
    if sin_theta.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "CEA deprojection: |lambda * y| > 1",
        ));
    }

    let theta = sin_theta.asin();
    Ok(NativeCoord::new(
        Angle::from_degrees(phi),
        Angle::from_degrees(theta * RAD_TO_DEG),
    ))
}

pub(crate) fn project_cyp(
    native: NativeCoord,
    mu: f64,
    lambda: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let (sin_theta, cos_theta) = theta.sin_cos();
    let denom = mu + cos_theta;

    if denom.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "CYP projection singularity: mu + cos(theta) = 0",
        ));
    }

    let x = lambda * phi * RAD_TO_DEG;
    let y = (mu + lambda) * sin_theta / denom * RAD_TO_DEG;
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_cyp(
    inter: IntermediateCoord,
    mu: f64,
    lambda: f64,
) -> WcsResult<NativeCoord> {
    if lambda.abs() < 1e-15 {
        return Err(WcsError::invalid_parameter(
            "CYP deprojection: lambda cannot be zero",
        ));
    }

    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;

    let phi = x / lambda;

    let eta = y / (mu + lambda);

    let a = eta * (mu - 1.0);
    let c = eta * (mu + 1.0);

    let theta = if a.abs() < 1e-15 {
        2.0 * (c / 2.0).atan()
    } else {
        let discriminant = 4.0 - 4.0 * a * c;
        if discriminant < 0.0 {
            return Err(WcsError::out_of_bounds(
                "CYP deprojection: point outside valid region",
            ));
        }
        let t = (2.0 - discriminant.sqrt()) / (2.0 * a);
        2.0 * t.atan()
    };

    Ok(native_coord_from_radians(phi, theta))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;

    #[test]
    fn test_car_reference_point() {
        let proj = Projection::car();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), 0.0);
    }

    #[test]
    fn test_car_roundtrip() {
        let proj = Projection::car();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();

        assert_eq!(original.phi().degrees(), recovered.phi().degrees());
        assert_eq!(original.theta().degrees(), recovered.theta().degrees());
    }

    #[test]
    fn test_car_known_values() {
        let proj = Projection::car();

        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(45.0));
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 90.0);
        assert_eq!(inter.y_deg(), 45.0);

        let native2 = NativeCoord::new(Angle::from_degrees(-120.0), Angle::from_degrees(-60.0));
        let inter2 = proj.project(native2).unwrap();
        assert_ulp_lt!(inter2.x_deg(), -120.0, 1);
        assert_ulp_lt!(inter2.y_deg(), -60.0, 1);
    }

    #[test]
    fn test_car_roundtrip_various_angles() {
        let proj = Projection::car();
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-85.0, -45.0, 0.0, 45.0, 85.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                assert_eq!(original.phi().degrees(), recovered.phi().degrees());
                assert_eq!(original.theta().degrees(), recovered.theta().degrees());
            }
        }
    }

    #[test]
    fn test_car_native_reference() {
        let proj = Projection::car();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_mer_reference_point() {
        let proj = Projection::mer();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_mer_roundtrip() {
        let proj = Projection::mer();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();

        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 1);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_mer_known_value() {
        let proj = Projection::mer();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        let expected_y = ((std::f64::consts::FRAC_PI_4 + 45.0_f64.to_radians() / 2.0)
            .tan()
            .ln())
            * RAD_TO_DEG;
        assert_ulp_lt!(inter.y_deg(), expected_y, 1);
    }

    #[test]
    fn test_mer_singularity_north_pole() {
        let proj = Projection::mer();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let result = proj.project(native);

        assert!(result.is_err());
    }

    #[test]
    fn test_mer_singularity_south_pole() {
        let proj = Projection::mer();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        let result = proj.project(native);

        assert!(result.is_err());
    }

    #[test]
    fn test_mer_roundtrip_various_angles() {
        let proj = Projection::mer();
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-80.0, -45.0, 0.0, 45.0, 80.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
            }
        }
    }

    #[test]
    fn test_mer_native_reference() {
        let proj = Projection::mer();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_cea_reference_point() {
        let proj = Projection::cea();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), 0.0);
    }

    #[test]
    fn test_cea_roundtrip() {
        let proj = Projection::cea();
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();

        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 1);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_cea_known_value() {
        let proj = Projection::cea();
        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(30.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 90.0);
        let expected_y = 30.0_f64.to_radians().sin() * RAD_TO_DEG;
        assert_ulp_lt!(inter.y_deg(), expected_y, 1);
    }

    #[test]
    fn test_cea_with_lambda() {
        let proj = Projection::cea_with_lambda(0.5);
        let original = NativeCoord::new(Angle::from_degrees(60.0), Angle::from_degrees(45.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();

        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 1);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_cea_roundtrip_various_angles() {
        let proj = Projection::cea();
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-85.0, -45.0, 0.0, 45.0, 85.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
            }
        }
    }

    #[test]
    fn test_cea_out_of_bounds() {
        let proj = Projection::cea();
        let inter = IntermediateCoord::new(0.0, 100.0);
        let result = proj.deproject(inter);

        assert!(result.is_err());
    }

    #[test]
    fn test_cea_native_reference() {
        let proj = Projection::cea();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_cyp_reference_point() {
        let proj = Projection::cyp(1.0, 1.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), 0.0);
    }

    #[test]
    fn test_cyp_roundtrip() {
        let proj = Projection::cyp(1.0, 1.0);
        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();

        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_cyp_various_parameters() {
        for mu in [0.5, 1.0, 2.0, 5.0] {
            for lambda in [0.5, 1.0, 2.0] {
                let proj = Projection::cyp(mu, lambda);
                let original =
                    NativeCoord::new(Angle::from_degrees(60.0), Angle::from_degrees(45.0));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 5);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 5);
            }
        }
    }

    #[test]
    fn test_cyp_roundtrip_various_angles() {
        let proj = Projection::cyp(1.0, 1.0);
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
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
    fn test_cyp_singularity() {
        let mu = 0.5;
        let proj = Projection::cyp(mu, 1.0);
        let theta_singular = (-(mu)).acos() * RAD_TO_DEG;
        let native = NativeCoord::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(theta_singular),
        );
        let result = proj.project(native);

        assert!(result.is_err());
    }

    #[test]
    fn test_cyp_native_reference() {
        let proj = Projection::cyp(1.0, 1.0);
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_cylindrical_projections_native_reference() {
        let projections = [
            Projection::car(),
            Projection::mer(),
            Projection::cea(),
            Projection::cyp(1.0, 1.0),
        ];

        for proj in projections {
            let (phi0, theta0) = proj.native_reference();
            assert_eq!(phi0, 0.0);
            assert_eq!(theta0, 0.0);
        }
    }

    #[test]
    fn test_cyp_lambda_zero() {
        let proj = Projection::cyp(1.0, 0.0);
        let inter = IntermediateCoord::new(10.0, 10.0);
        let result = proj.deproject(inter);
        assert!(result.is_err());
    }

    #[test]
    fn test_cyp_mu_equals_one_linear_case() {
        let proj = Projection::cyp(1.0, 1.0);
        let native = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(30.0));
        let inter = proj.project(native).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(native.phi().degrees(), recovered.phi().degrees(), 2);
        assert_ulp_lt!(native.theta().degrees(), recovered.theta().degrees(), 2);
    }
}
