use celestial_core::constants::{DEG_TO_RAD, HALF_PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, pole_native_coord, radial_to_intermediate,
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

    native_coord_from_radians(phi, theta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::phi_on_same_edge;
    use crate::Projection;
    use celestial_core::angle::Angle;
    use celestial_core::assert_ulp_le;

    // ARC computes θ = π/2 − r, so θ carries an absolute error of a few ulps of π/2 (2^-52),
    // which a ULP compare against exactly 0 can't express.
    fn assert_theta_recovered(original: Angle, recovered: Angle, context: &str) {
        if original.radians() == 0.0 {
            let error = recovered.radians().abs();
            assert!(
                error <= 8.0 * f64::EPSILON,
                "theta {context}: {error:e} rad from 0"
            );
        } else {
            assert_ulp_le!(
                original.degrees(),
                recovered.degrees(),
                8,
                "theta {}",
                context
            );
        }
    }

    #[test]
    fn test_arc_roundtrip() {
        let proj = Projection::arc();
        for phi_deg in [-180.0, -90.0, 0.0, 45.0, 120.0, 180.0] {
            for theta_deg in [-60.0, 0.0, 30.0, 45.0, 75.0, 89.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let recovered = proj.deproject(proj.project(original).unwrap()).unwrap();
                let context = format!("(phi={}, theta={})", phi_deg, theta_deg);

                let phi_back = phi_on_same_edge(original.phi(), recovered.phi());
                let (phi, phi_back) = (original.phi().degrees(), phi_back.degrees());
                assert_ulp_le!(phi, phi_back, 8, "phi {}", context);
                assert_theta_recovered(original.theta(), recovered.theta(), &context);
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
