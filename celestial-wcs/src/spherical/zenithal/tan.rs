use celestial_core::constants::{DEG_TO_RAD, HALF_PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, pole_native_coord,
    radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_tan(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta == HALF_PI {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }
    if theta <= 0.0 {
        return Err(WcsError::singularity(
            "TAN projection undefined at theta <= 0",
        ));
    }
    let (rt_sin, rt_cos) = libm::sincos(theta);
    let r_theta = rt_cos / rt_sin;
    Ok(radial_to_intermediate(r_theta, phi))
}

pub(crate) fn deproject_tan(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let (phi, r_theta, is_pole) = intermediate_to_polar(x, y);

    if is_pole {
        return Ok(pole_native_coord());
    }

    let theta = libm::atan2(1.0_f64, r_theta);

    Ok(native_coord_from_radians(phi, theta))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_tan_roundtrip() {
        let proj = Projection::tan();
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 135.0, 180.0] {
            // TAN only valid for theta > 0; sweep up to near the pole.
            for theta_deg in [5.0, 30.0, 45.0, 60.0, 80.0, 89.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 4,
                    "phi (phi={}, theta={})", phi_deg, theta_deg);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 4,
                    "theta (phi={}, theta={})", phi_deg, theta_deg);
            }
        }
    }

    #[test]
    fn test_tan_singularity() {
        // TAN diverges as theta -> 0 (cot(theta) -> infinity).
        let proj = Projection::tan();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        assert!(proj.project(native).is_err());
    }
}
