use celestial_core::constants::{DEG_TO_RAD, HALF_PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, pole_native_coord,
    radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_stg(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta == HALF_PI {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }
    if theta == -HALF_PI {
        return Err(WcsError::singularity(
            "STG projection diverges at theta = -90",
        ));
    }
    let (theta_s, theta_c) = libm::sincos(theta);
    let r_theta = 2.0 * theta_c / (1.0 + theta_s);
    Ok(radial_to_intermediate(r_theta, phi))
}

pub(crate) fn deproject_stg(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let (phi, r_theta, is_pole) = intermediate_to_polar(x, y);

    if is_pole {
        return Ok(pole_native_coord());
    }

    let theta = HALF_PI - 2.0 * libm::atan(r_theta / 2.0);

    Ok(native_coord_from_radians(phi, theta))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_stg_roundtrip() {
        let proj = Projection::stg();
        for phi_deg in [-180.0, -60.0, 0.0, 45.0, 120.0, 180.0] {
            for theta_deg in [-45.0, 0.0, 30.0, 60.0, 80.0, 89.0] {
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
    fn test_stg_singularity() {
        // STG diverges at the antipodal point (theta = -90).
        let proj = Projection::stg();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        assert!(proj.project(native).is_err());
    }
}
