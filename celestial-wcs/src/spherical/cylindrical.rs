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

    let y = libm::log(libm::tan(std::f64::consts::FRAC_PI_4 + theta / 2.0)) * RAD_TO_DEG;
    Ok(IntermediateCoord::new(phi, y))
}

pub(crate) fn deproject_mer(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let phi = inter.x_deg();
    let y = inter.y_deg() * DEG_TO_RAD;

    let theta = 2.0 * libm::atan(libm::exp(y)) - HALF_PI;
    Ok(NativeCoord::new(
        Angle::from_degrees(phi),
        Angle::from_degrees(theta * RAD_TO_DEG),
    ))
}

pub(crate) fn project_cea(native: NativeCoord, lambda: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().degrees();
    let theta = native.theta().radians();

    let y = libm::sin(theta) / lambda * RAD_TO_DEG;
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

    let theta = libm::asin(sin_theta);
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

    let (sin_theta, cos_theta) = libm::sincos(theta);
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
        2.0 * libm::atan(c / 2.0)
    } else {
        let discriminant = 4.0 - 4.0 * a * c;
        if discriminant < 0.0 {
            return Err(WcsError::out_of_bounds(
                "CYP deprojection: point outside valid region",
            ));
        }
        let t = (2.0 - libm::sqrt(discriminant)) / (2.0 * a);
        2.0 * libm::atan(t)
    };

    Ok(native_coord_from_radians(phi, theta))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::constants::{PI, QUARTER_PI};
    use celestial_core::Angle;

    // Per-projection native_reference checks are covered by
    // spherical::tests::test_all_projections_map_reference_to_origin.

    #[test]
    fn test_car_known_values_and_roundtrip() {
        // CAR is a literal pass-through (phi, theta) -> (x, y) in degrees.
        let proj = Projection::car();

        // Spec-anchored: a couple of known points should be byte-exact.
        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(45.0));
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 90.0);
        assert_eq!(inter.y_deg(), 45.0);

        let native2 = NativeCoord::new(Angle::from_degrees(-120.0), Angle::from_degrees(-60.0));
        let inter2 = proj.project(native2).unwrap();
        assert_ulp_lt!(inter2.x_deg(), -120.0, 1);
        assert_ulp_lt!(inter2.y_deg(), -60.0, 1);

        // CAR is an identity in degrees, so the roundtrip must be byte-exact.
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-85.0, -45.0, 0.0, 45.0, 85.0] {
                let original = NativeCoord::new(
                    Angle::from_degrees(phi_deg),
                    Angle::from_degrees(theta_deg),
                );
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_eq!(original.phi().degrees(), recovered.phi().degrees());
                assert_eq!(original.theta().degrees(), recovered.theta().degrees());
            }
        }
    }

    #[test]
    fn test_mer_known_value_and_roundtrip() {
        let proj = Projection::mer();

        // Spec-anchored Paper II Eq. 17: y = (180/pi) * ln(tan(pi/4 + theta/2)).
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 0.0);
        let expected_y = libm::log(libm::tan(QUARTER_PI + (PI / 4.0) / 2.0)) * RAD_TO_DEG;
        assert_ulp_lt!(inter.y_deg(), expected_y, 1);

        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-80.0, -45.0, 0.0, 45.0, 80.0] {
                let original = NativeCoord::new(
                    Angle::from_degrees(phi_deg),
                    Angle::from_degrees(theta_deg),
                );
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
            }
        }
    }

    #[test]
    fn test_cea_known_value_and_roundtrip() {
        let proj = Projection::cea();

        // Spec-anchored Paper II Eq. 19 with lambda = 1: y = (180/pi) * sin(theta).
        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(30.0));
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 90.0);
        let expected_y = libm::sin(PI / 6.0) * RAD_TO_DEG;
        assert_ulp_lt!(inter.y_deg(), expected_y, 1);

        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-85.0, -45.0, 0.0, 45.0, 85.0] {
                let original = NativeCoord::new(
                    Angle::from_degrees(phi_deg),
                    Angle::from_degrees(theta_deg),
                );
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
            }
        }

        // Custom lambda preserves roundtrip.
        let proj_half = Projection::cea_with_lambda(0.5);
        let pt = NativeCoord::new(Angle::from_degrees(60.0), Angle::from_degrees(45.0));
        let inter = proj_half.project(pt).unwrap();
        let recovered = proj_half.deproject(inter).unwrap();
        assert_ulp_lt!(pt.phi().degrees(), recovered.phi().degrees(), 1);
        assert_ulp_lt!(pt.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_cyp_roundtrip_parameter_sweep() {
        // Sweeps mu x lambda x phi x theta in one shot.  The fixed-parameter
        // grid (mu=1, lambda=1) covers many phi/theta points; the parameter
        // sweep covers many (mu, lambda) at one well-behaved point.
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

        let proj = Projection::cyp(1.0, 1.0);
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original = NativeCoord::new(
                    Angle::from_degrees(phi_deg),
                    Angle::from_degrees(theta_deg),
                );
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 5);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 5);
            }
        }
    }

    #[test]
    fn test_cylindrical_error_cases() {
        // MER blows up at the poles (theta = +/- 90).
        let mer = Projection::mer();
        for theta in [90.0, -90.0] {
            let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta));
            assert!(mer.project(native).is_err(), "MER theta = {} should error", theta);
        }

        // CEA deproject: |lambda * y| > 1 is out of range.
        let cea = Projection::cea();
        assert!(cea.deproject(IntermediateCoord::new(0.0, 100.0)).is_err());

        // CYP projects fail at theta = arccos(-mu) (denom mu + cos(theta) -> 0).
        let mu = 0.5;
        let cyp = Projection::cyp(mu, 1.0);
        let theta_singular = libm::acos(-mu) * RAD_TO_DEG;
        let native = NativeCoord::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(theta_singular),
        );
        assert!(cyp.project(native).is_err());

        // CYP deproject with lambda = 0 is an invalid parameter.
        let cyp_zero = Projection::cyp(1.0, 0.0);
        assert!(cyp_zero.deproject(IntermediateCoord::new(10.0, 10.0)).is_err());
    }
}
