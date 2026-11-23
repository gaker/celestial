use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{atan2d, cotd, to_degrees, to_radians};

// Paper II eq. 54, Rθ = (180°/π) cot θ.
pub(crate) fn project_tan(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    if theta == 90.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }
    if theta <= 0.0 {
        return Err(WcsError::singularity(
            "TAN projection undefined at theta <= 0",
        ));
    }
    let r_theta = to_degrees(cotd(theta));
    Ok(radial_to_intermediate(r_theta, native.phi_deg()))
}

// Eq. 55, θ = tan⁻¹(180°/(π Rθ)).
pub(crate) fn deproject_tan(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    if is_pole {
        return Ok(pole_native_coord());
    }
    native_coord(phi, atan2d(1.0, to_radians(r_theta)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Against a 60-digit evaluation the projected values are within 1.92 ulp
    // (y at (45°, 5°)) and the deprojected ones within 0.58 ulp.
    #[test]
    fn test_tan_pinned_both_ways() {
        let rows = [
            (
                (-180.0, 30.0),
                (0.0, 99.23920117592256),
                (180.0, 30.000000000000004),
            ),
            (
                (-90.0, 60.0),
                (-33.07973372530753, -0.0),
                (-90.0, 59.99999999999999),
            ),
            (
                (45.0, 5.0),
                (463.07981622173133, -463.0798162217314),
                (45.0, 4.999999999999999),
            ),
            (
                (135.0, 89.0),
                (0.707178588950518, 0.7071785889505181),
                (135.0, 89.0),
            ),
        ];
        assert_pinned(&Projection::tan(), &rows);
    }

    // Eq. 54 at θ = 45° is Rθ = 180°/π, which rounds to 57.29577951308232.
    #[test]
    fn test_tan_known_value() {
        let inter = Projection::tan()
            .project(NativeCoord::new(0.0, 45.0))
            .unwrap();
        assert_eq!((inter.x_deg(), inter.y_deg()), (0.0, -57.29577951308232));
    }

    #[test]
    fn test_tan_singularity() {
        // TAN diverges as theta -> 0 (cot(theta) -> infinity).
        let proj = Projection::tan();
        let native = NativeCoord::new(0.0, 0.0);
        assert!(proj.project(native).is_err());
    }
}
