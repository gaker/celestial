use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{atand, tand, to_degrees, to_radians};

// Paper II eq. 56, Rθ = (360°/π) tan((90° − θ)/2).
pub(crate) fn project_stg(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    if theta == 90.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }
    if theta == -90.0 {
        return Err(WcsError::singularity(
            "STG projection diverges at theta = -90",
        ));
    }
    let r_theta = to_degrees(2.0 * tand((90.0 - theta) / 2.0));
    Ok(radial_to_intermediate(r_theta, native.phi_deg()))
}

// Eq. 57, θ = 90° − 2 tan⁻¹(πRθ/360°).
pub(crate) fn deproject_stg(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    if is_pole {
        return Ok(pole_native_coord());
    }
    native_coord(phi, 90.0 - 2.0 * atand(to_radians(r_theta) / 2.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Against a 60-digit evaluation the projected values are within 2.10 ulp
    // (x at (45°, 80°)) and the deprojected ones within 0.93 ulp. At (120°, 0°)
    // θ is 90° less a value near 90°, so it carries that value's absolute
    // rounding: 0 against −3.5e-15, 0.25 ulp of 90°.
    #[test]
    fn test_stg_pinned_both_ways() {
        let rows = [
            ((-180.0, -45.0), (0.0, 276.6484959344437), (180.0, -45.0)),
            (
                (45.0, 80.0),
                (7.0890724126051605, -7.089072412605161),
                (45.0, 80.0),
            ),
            (
                (120.0, 0.0),
                (99.23920117592256, 57.29577951308232),
                (120.00000000000001, 0.0),
            ),
            (
                (-60.0, 60.0),
                (-26.59106381259818, -15.352357850242075),
                (-59.99999999999999, 60.0),
            ),
            (
                (120.0, 89.0),
                (0.8660473883227904, 0.5000126927791354),
                (120.00000000000001, 89.0),
            ),
        ];
        assert_pinned(&Projection::stg(), &rows);
    }

    #[test]
    fn test_stg_singularity() {
        // STG diverges at the antipodal point (theta = -90).
        let proj = Projection::stg();
        let native = NativeCoord::new(0.0, -90.0);
        assert!(proj.project(native).is_err());
    }
}
