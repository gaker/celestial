use celestial_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};

use crate::common::native_coord_from_radians;
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_sin(native: NativeCoord, xi: f64, eta: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta == HALF_PI {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    let (sin_theta, cos_theta) = libm::sincos(theta);
    let (sin_phi, cos_phi) = libm::sincos(phi);

    let x = (cos_theta * sin_phi + xi * (1.0 - sin_theta)) * RAD_TO_DEG;
    let y = -(cos_theta * cos_phi - eta * (1.0 - sin_theta)) * RAD_TO_DEG;
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_sin(inter: IntermediateCoord, xi: f64, eta: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;

    let a = xi * xi + eta * eta + 1.0;
    let b = xi * (x - xi) + eta * (y - eta);
    let c = (x - xi) * (x - xi) + (y - eta) * (y - eta) - 1.0;

    let discriminant = b * b - a * c;
    if discriminant < 0.0 {
        return Err(WcsError::out_of_bounds(
            "Point outside SIN projection boundary",
        ));
    }

    let sin_theta = (-b + libm::sqrt(discriminant)) / a;
    if sin_theta.abs() > 1.0 {
        return Err(WcsError::out_of_bounds("Invalid theta in SIN deprojection"));
    }

    let theta = libm::asin(sin_theta);
    let x_adj = x - xi * (1.0 - sin_theta);
    let y_adj = y - eta * (1.0 - sin_theta);
    let phi = libm::atan2(x_adj, -y_adj);

    Ok(native_coord_from_radians(phi, theta))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_sin_roundtrip() {
        // Sweeps both standard SIN (xi=eta=0) and slant SIN (NCP-like) configs.
        let cases: &[(f64, f64, f64, f64, u64)] = &[
            (0.0, 0.0, 30.0, 60.0, 1),
            (0.0, 0.0, -90.0, 30.0, 2),
            (0.0, 0.0, 135.0, 75.0, 2),
            (0.1, -0.2, 25.0, 55.0, 2),
            (0.05, 0.05, 60.0, 45.0, 2),
        ];
        for (xi, eta, phi, theta, ulp) in cases {
            let proj = Projection::sin_with_params(*xi, *eta);
            let original = NativeCoord::new(Angle::from_degrees(*phi), Angle::from_degrees(*theta));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), *ulp,
                "phi (xi={}, eta={})", xi, eta);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), *ulp,
                "theta (xi={}, eta={})", xi, eta);
        }
    }

    #[test]
    fn test_sin_out_of_bounds() {
        // Standard SIN: point outside projection boundary.
        let err = deproject_sin(IntermediateCoord::new(100.0, 100.0), 0.0, 0.0).unwrap_err();
        assert!(matches!(err, WcsError::OutOfBounds { .. }), "standard: {:?}", err);

        // Slant SIN: parameter combo that drives sin_theta out of [-1, 1].
        let err = deproject_sin(IntermediateCoord::new(100.0, 100.0), 0.5, 0.5).unwrap_err();
        assert!(matches!(err, WcsError::OutOfBounds { .. }), "slant: {:?}", err);
    }
}
