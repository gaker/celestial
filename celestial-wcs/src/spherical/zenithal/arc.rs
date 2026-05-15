use celestial_core::constants::{DEG_TO_RAD, HALF_PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, pole_native_coord,
    radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::WcsResult;

pub(crate) fn project_arc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let r_theta = HALF_PI - theta;
    Ok(radial_to_intermediate(r_theta, phi))
}

pub(crate) fn deproject_arc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let (phi, r_theta, is_pole) = intermediate_to_polar(x, y);

    if is_pole {
        return Ok(pole_native_coord());
    }

    let theta = HALF_PI - r_theta;

    Ok(native_coord_from_radians(phi, theta))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_arc_roundtrip() {
        let proj = Projection::arc();
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 120.0, 180.0] {
            for theta_deg in [-60.0, 0.0, 30.0, 45.0, 75.0, 89.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 8,
                    "phi (phi={}, theta={})", phi_deg, theta_deg);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 8,
                    "theta (phi={}, theta={})", phi_deg, theta_deg);
            }
        }
    }

    #[test]
    fn test_arc_known_value() {
        // Anchors absolute output. ARC at (phi=0, theta=45) projects to (0, -45)
        // by direct application of the Paper II formula (r_theta = 90 - theta in degrees).
        let proj = Projection::arc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), -45.0);
    }
}
