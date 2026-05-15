use celestial_core::constants::{DEG_TO_RAD, HALF_PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, pole_native_coord,
    radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_zea(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let r_theta = libm::sqrt(2.0 * (1.0 - libm::sin(theta)));
    Ok(radial_to_intermediate(r_theta, phi))
}

pub(crate) fn deproject_zea(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let (phi, r_theta, is_pole) = intermediate_to_polar(x, y);

    if is_pole {
        return Ok(pole_native_coord());
    }

    let rho = r_theta / 2.0;
    if rho > 1.0 {
        return Err(WcsError::out_of_bounds(
            "Point outside ZEA projection boundary",
        ));
    }

    let theta = HALF_PI - 2.0 * libm::asin(rho);

    Ok(native_coord_from_radians(phi, theta))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_zea_roundtrip() {
        let proj = Projection::zea();
        // ZEA is defined for the whole sphere including theta = -90 (antipode).
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 135.0, 180.0] {
            for theta_deg in [-89.0, -45.0, 30.0, 60.0, 89.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();

                // Near the antipode (theta -> -90), the asin in deproject
                // loses precision; allow more ULP slack there.
                let ulp_bar = if theta_deg <= -85.0 { 64 } else { 8 };
                assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), ulp_bar,
                    "phi (phi={}, theta={})", phi_deg, theta_deg);
                assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), ulp_bar,
                    "theta (phi={}, theta={})", phi_deg, theta_deg);
            }
        }
    }

    #[test]
    fn test_zea_out_of_bounds() {
        // ZEA's projection plane is a disk of radius 2·(180/π); points outside are invalid.
        let proj = Projection::zea();
        assert!(proj.deproject(IntermediateCoord::new(200.0, 0.0)).is_err());
    }
}
