use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{asind, sind, to_degrees, to_radians};

// Paper II eq. 69, Rθ = (360°/π) sin((90° − θ)/2). In degrees 90° − θ is
// exact near the pole, so the radius keeps full precision there.
pub(crate) fn project_zea(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let r_theta = to_degrees(2.0 * sind((90.0 - native.theta_deg()) / 2.0));
    Ok(radial_to_intermediate(r_theta, native.phi_deg()))
}

// Eq. 70, θ = 90° − 2 sin⁻¹(πRθ/360°).
pub(crate) fn deproject_zea(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    if is_pole {
        return Ok(pole_native_coord());
    }
    let rho = to_radians(r_theta) / 2.0;
    if rho > 1.0 {
        return Err(WcsError::out_of_bounds(
            "Point outside ZEA projection boundary",
        ));
    }
    native_coord(phi, 90.0 - 2.0 * asind(rho))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Against a 60-digit evaluation the projected values are within 1.41 ulp
    // (y at (135°, 60°)) and the deprojected ones within 1.35 ulp, except θ at
    // −89°, 50 and 80 ulp low. Near the antipode dθ/dR = −1/cos ξ, and R and θ
    // share a binade, so one ulp of R moves θ by 1/cos ξ ≈ 115 ulps.
    #[test]
    fn test_zea_pinned_both_ways() {
        let rows = [
            (
                (-180.0, -89.0),
                (0.0, 114.5871957307251),
                (180.0, -89.00000000000173),
            ),
            (
                (45.0, -89.0),
                (81.02538313834592, -81.02538313834593),
                (44.99999999999999, -89.00000000000173),
            ),
            (
                (-90.0, 30.0),
                (-57.29577951308232, -0.0),
                (-90.0, 29.999999999999993),
            ),
            (
                (135.0, 60.0),
                (20.971710831420125, 20.97171083142013),
                (135.0, 60.0),
            ),
            ((0.0, 89.0), (0.0, -0.9999873076558378), (0.0, 89.0)),
        ];
        assert_pinned(&Projection::zea(), &rows);
    }

    #[test]
    fn test_zea_near_pole_keeps_the_radius() {
        // 1e-6° from the pole, 1 − sin θ is below the rounding of sin θ, so
        // eq. 69 takes the sine of the half colatitude. Against a 60-digit
        // evaluation x is 0.58 ulp high and y correctly rounded.
        let native = NativeCoord::new(37.5, 89.999999);
        let inter = Projection::zea().project(native).unwrap();
        assert_eq!(
            (inter.x_deg(), inter.y_deg()),
            (6.087614274717458e-7, -7.933533382882105e-7)
        );
    }

    #[test]
    fn test_zea_out_of_bounds() {
        // ZEA's projection plane is a disk of radius 2·(180/π); points outside are invalid.
        let proj = Projection::zea();
        assert!(proj.deproject(IntermediateCoord::new(200.0, 0.0)).is_err());
    }
}
