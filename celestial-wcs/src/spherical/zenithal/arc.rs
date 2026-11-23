use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::WcsResult;

// Paper II eq. 67, Rθ = 90° − θ, exact in degrees.
pub(crate) fn project_arc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let r_theta = 90.0 - native.theta_deg();
    Ok(radial_to_intermediate(r_theta, native.phi_deg()))
}

pub(crate) fn deproject_arc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    if is_pole {
        return Ok(pole_native_coord());
    }
    native_coord(phi, 90.0 - r_theta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Against a 60-digit evaluation the projected values are within 0.80 ulp
    // and the deprojected ones within 0.81 ulp, both at (120°, 30°). θ = 90° − Rθ
    // carries the absolute rounding of Rθ, so at (45°, 0°) it is 0 against
    // −1.7e-15, 0.12 ulp of 90°.
    #[test]
    fn test_arc_pinned_both_ways() {
        let rows = [
            ((-180.0, -60.0), (0.0, 150.0), (180.0, -60.0)),
            (
                (45.0, 0.0),
                (63.63961030678927, -63.63961030678928),
                (45.0, 0.0),
            ),
            (
                (120.0, 30.0),
                (51.96152422706631, 30.0),
                (120.00000000000001, 30.000000000000007),
            ),
            ((-90.0, 89.0), (-1.0, -0.0), (-90.0, 89.0)),
        ];
        assert_pinned(&Projection::arc(), &rows);
    }

    #[test]
    fn test_arc_known_value() {
        // Anchors absolute output. ARC at (phi=0, theta=45) projects to (0, -45)
        // by direct application of the Paper II formula (r_theta = 90 - theta in degrees).
        let proj = Projection::arc();
        let native = NativeCoord::new(0.0, 45.0);
        let inter = proj.project(native).unwrap();

        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), -45.0);
    }
}
