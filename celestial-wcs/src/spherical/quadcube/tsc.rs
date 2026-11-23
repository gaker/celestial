use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

use super::{front_face, quadcube_face_from_xy};

pub(crate) fn project_tsc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let face = front_face(native, "TSC")?;
    let (chi, psi) = face.plane();
    Ok(face.intermediate(45.0 * chi, 45.0 * psi))
}

pub(crate) fn deproject_tsc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (face, x_rel, y_rel) = quadcube_face_from_xy(inter.x_deg(), inter.y_deg());

    let chi = x_rel / 45.0;
    let psi = y_rel / 45.0;

    if chi.abs() > 1.0 || psi.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "TSC deprojection: point outside cube face",
        ));
    }

    face.native_from_plane(chi, psi)
}

#[cfg(test)]
mod tests {
    use crate::common::assert_pinned;
    use crate::coordinate::IntermediateCoord;
    use crate::spherical::Projection;

    // One row on each face, then (180°, −30°) on the seam and (−135°, −30°) on
    // the edge between faces 3 and 4, which goes to face 3. Against a 60-digit
    // evaluation the projected values are within 2.61 ulp and the deprojected
    // ones within 2.21 ulp. Each offset is within 1.07 times the summed
    // half-ulp rounding of the steps as evaluated; the largest, y at
    // (−30°, 15°), carries the rounding of sin θ, cos φ cos θ and their
    // quotient.
    #[test]
    fn test_tsc_pinned_both_ways() {
        let rows = [
            (
                (-30.0, 15.0),
                (-25.980762113533164, 13.923048454132642),
                (-30.000000000000004, 15.000000000000005),
            ),
            (
                (120.0, 15.0),
                (115.98076211353316, 13.923048454132642),
                (120.00000000000001, 15.000000000000005),
            ),
            (
                (180.0, -30.0),
                (180.0, -25.980762113533164),
                (180.0, -30.000000000000004),
            ),
            (
                (-90.0, -40.0),
                (-90.0, -37.759483402977594),
                (-90.0, -39.99999999999999),
            ),
            (
                (-135.0, 80.0),
                (-5.61069016942973, 95.61069016942973),
                (-135.0, 80.0),
            ),
            (
                (20.0, -75.0),
                (4.12398095396804, -78.66945545071188),
                (20.000000000000004, -75.0),
            ),
            (
                (-135.0, -30.0),
                (225.0, -36.74234614174767),
                (-135.0, -29.99999999999999),
            ),
        ];
        assert_pinned(&Projection::tsc(), &rows);
    }

    #[test]
    fn test_tsc_inverse_near_the_native_pole() {
        // (37.5, 89.999999) projects to this point. Against eqs. 165 to 169
        // evaluated to 60 digits from it, φ is correctly rounded and θ is
        // 0.69 ulp low.
        let inter = IntermediateCoord::new(4.781201070835181e-7, 89.99999937690174);
        let native = Projection::tsc().deproject(inter).unwrap();
        assert_eq!(
            (native.phi_deg(), native.theta_deg()),
            (37.499999804832704, 89.99999899999999)
        );
    }
}
