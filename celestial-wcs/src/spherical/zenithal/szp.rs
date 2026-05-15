use celestial_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};

use crate::common::{
    native_coord_from_radians, pole_native_coord,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_szp(
    native: NativeCoord,
    mu: f64,
    phi_c_deg: f64,
    theta_c_deg: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta == HALF_PI {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    let phi_c = phi_c_deg * DEG_TO_RAD;
    let theta_c = theta_c_deg * DEG_TO_RAD;

    let (sin_phi_c, cos_phi_c) = libm::sincos(phi_c);
    let (sin_theta_c, cos_theta_c) = libm::sincos(theta_c);

    let xp = -mu * cos_theta_c * sin_phi_c;
    let yp = mu * cos_theta_c * cos_phi_c;
    let zp = mu * sin_theta_c + 1.0;

    if zp.abs() < 1e-10 {
        return Err(WcsError::singularity("SZP projection singularity: zp = 0"));
    }

    let (sin_theta, cos_theta) = libm::sincos(theta);
    let (sin_phi, cos_phi) = libm::sincos(phi);

    let denom = zp - (1.0 - sin_theta);
    if denom.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "SZP projection singularity: denominator = 0",
        ));
    }

    let x = (zp * cos_theta * sin_phi - xp * (1.0 - sin_theta)) / denom * RAD_TO_DEG;
    let y = -(zp * cos_theta * cos_phi + yp * (1.0 - sin_theta)) / denom * RAD_TO_DEG;

    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_szp(
    inter: IntermediateCoord,
    mu: f64,
    phi_c_deg: f64,
    theta_c_deg: f64,
) -> WcsResult<NativeCoord> {
    let x_big = inter.x_deg() * DEG_TO_RAD;
    let y_big = inter.y_deg() * DEG_TO_RAD;

    if x_big == 0.0 && y_big == 0.0 {
        return Ok(pole_native_coord());
    }

    let phi_c = phi_c_deg * DEG_TO_RAD;
    let theta_c = theta_c_deg * DEG_TO_RAD;

    let (sin_phi_c, cos_phi_c) = libm::sincos(phi_c);
    let (sin_theta_c, cos_theta_c) = libm::sincos(theta_c);

    let xp = -mu * cos_theta_c * sin_phi_c;
    let yp = mu * cos_theta_c * cos_phi_c;
    let zp = mu * sin_theta_c + 1.0;

    if zp.abs() < 1e-10 {
        return Err(WcsError::singularity("SZP projection singularity: zp = 0"));
    }

    let x_prime = (x_big - xp) / zp;
    let y_prime = (y_big - yp) / zp;

    let a = x_prime * x_prime + y_prime * y_prime + 1.0;
    let b = x_prime * (x_big - x_prime) + y_prime * (y_big - y_prime);
    let c = (x_big - x_prime) * (x_big - x_prime) + (y_big - y_prime) * (y_big - y_prime) - 1.0;

    let discriminant = b * b - a * c;
    if discriminant < 0.0 {
        return Err(WcsError::out_of_bounds(
            "Point outside SZP projection boundary",
        ));
    }

    let sin_theta_plus = (-b + libm::sqrt(discriminant)) / a;
    let sin_theta_minus = (-b - libm::sqrt(discriminant)) / a;

    let sin_theta = if sin_theta_plus.abs() <= 1.0 + 1e-10 && sin_theta_minus.abs() <= 1.0 + 1e-10 {
        if sin_theta_plus > sin_theta_minus {
            sin_theta_plus
        } else {
            sin_theta_minus
        }
    } else if sin_theta_plus.abs() <= 1.0 + 1e-10 {
        sin_theta_plus
    } else if sin_theta_minus.abs() <= 1.0 + 1e-10 {
        sin_theta_minus
    } else {
        return Err(WcsError::out_of_bounds("Invalid theta in SZP deprojection"));
    };

    let sin_theta_clamped = sin_theta.clamp(-1.0, 1.0);
    let theta = libm::asin(sin_theta_clamped);

    let one_minus_sin_theta = 1.0 - sin_theta_clamped;
    let arg_x = x_big - x_prime * one_minus_sin_theta;
    let arg_y = -(y_big - y_prime * one_minus_sin_theta);
    let phi = libm::atan2(arg_x, arg_y);

    Ok(native_coord_from_radians(phi, theta))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_szp_roundtrip() {
        // (mu, phi_c, theta_c, phi, theta, ulp)
        // Exercises both ± sign branches of the inverse via varying parameter configs.
        let cases: &[(f64, f64, f64, f64, f64, u64)] = &[
            (0.0, 0.0, 90.0, 45.0, 60.0, 2),
            (2.0, 0.0, 90.0, 30.0, 70.0, 5),
            (2.0, 30.0, 60.0, 20.0, 75.0, 10),
            (0.5, 0.0, 90.0, 30.0, 30.0, 10),
            (3.0, 45.0, 45.0, 60.0, 20.0, 32),
            (2.0, 30.0, 70.0, 45.0, 50.0, 10),
            (1.0, 0.0, 90.0, 60.0, 45.0, 3),
            (5.0, 0.0, 90.0, 60.0, 45.0, 5),
            (10.0, 0.0, 90.0, 60.0, 45.0, 5),
        ];
        for (mu, phi_c, theta_c, phi, theta, ulp) in cases {
            let proj = Projection::szp(*mu, *phi_c, *theta_c);
            let original = NativeCoord::new(Angle::from_degrees(*phi), Angle::from_degrees(*theta));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), *ulp,
                "phi (mu={}, phi_c={}, theta_c={})", mu, phi_c, theta_c);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), *ulp,
                "theta (mu={}, phi_c={}, theta_c={})", mu, phi_c, theta_c);
        }
    }

    #[test]
    fn test_szp_singularity_errors() {
        let native_pole = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        let native_normal = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let inter_normal = IntermediateCoord::new(10.0, 10.0);
        let inter_far = IntermediateCoord::new(500.0, 500.0);

        let err = project_szp(native_normal, 1.0, 0.0, -90.0).unwrap_err();
        assert!(matches!(err, WcsError::Singularity { .. }), "zp_singularity_project: {:?}", err);

        let err = project_szp(native_pole, 1.0, 0.0, 90.0).unwrap_err();
        assert!(matches!(err, WcsError::Singularity { .. }), "denominator_singularity: {:?}", err);

        let err = deproject_szp(inter_normal, 1.0, 0.0, -90.0).unwrap_err();
        assert!(matches!(err, WcsError::Singularity { .. }), "zp_singularity_deproject: {:?}", err);

        let err = deproject_szp(inter_far, 2.0, 45.0, 60.0).unwrap_err();
        assert!(matches!(err, WcsError::OutOfBounds { .. }), "negative_discriminant: {:?}", err);
    }
}
