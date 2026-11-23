use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{asind, atand, sincosd};

use super::{front_face, quadcube_face_from_xy};

pub(crate) fn project_qsc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let face = front_face(native);

    // Near a face centre 1 − ζ falls below the rounding of ζ. The direction
    // cosines are a unit vector, so 1 − ζ = (ξ² + η²)/(1 + ζ), which keeps
    // full precision there.
    let one_minus_zeta = (face.xi * face.xi + face.eta * face.eta) / (1.0 + face.zeta);
    let (x, y) = if face.xi.abs() >= face.eta.abs() {
        qsc_forward_offsets(face.xi, face.eta, one_minus_zeta)
    } else {
        let (u, v) = qsc_forward_offsets(face.eta, face.xi, one_minus_zeta);
        (v, u)
    };
    Ok(face.intermediate(x, y))
}

// Eqs. 177 and 178, where `major` is whichever of ξ and η is larger in
// magnitude. ω = minor/major is defined everywhere but the face centre, where
// both are zero and so are u and v.
fn qsc_forward_offsets(major: f64, minor: f64, one_minus_zeta: f64) -> (f64, f64) {
    let omega = if major == 0.0 { 0.0 } else { minor / major };
    let omega2 = omega * omega;
    let s = if major >= 0.0 { 1.0 } else { -1.0 };
    let u = s * 45.0 * libm::sqrt(one_minus_zeta / (1.0 - 1.0 / libm::sqrt(2.0 + omega2)));
    let v = (u / 15.0) * (atand(omega) - asind(omega / libm::sqrt(2.0 * (1.0 + omega2))));
    (u, v)
}

pub(crate) fn deproject_qsc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (face, x_rel, y_rel) = quadcube_face_from_xy(inter.x_deg(), inter.y_deg());
    if x_rel.abs() > 45.0 || y_rel.abs() > 45.0 {
        return Err(WcsError::out_of_bounds(
            "QSC deprojection: point outside cube face",
        ));
    }

    let (xi, eta, zeta) = if x_rel.abs() >= y_rel.abs() {
        qsc_inverse_cosines(x_rel, y_rel)
    } else {
        let (eta, xi, zeta) = qsc_inverse_cosines(y_rel, x_rel);
        (xi, eta, zeta)
    };
    face.native(xi, eta, zeta)
}

// Eqs. 180 to 185, with u the larger offset; returns (major, minor, ζ). Near a
// face centre ζ rounds towards 1 and eq. 182's 1 − ζ² cancels. 1 − ζ is the
// small term of eq. 180, and 1 − ζ² = (1 − ζ)(1 + ζ) keeps full precision.
fn qsc_inverse_cosines(u: f64, v: f64) -> (f64, f64, f64) {
    let omega = qsc_inverse_omega(u, v);
    let omega2 = omega * omega;
    let one_minus_zeta = (u / 45.0) * (u / 45.0) * (1.0 - 1.0 / libm::sqrt(2.0 + omega2));
    let zeta = 1.0 - one_minus_zeta;
    let major = libm::sqrt(one_minus_zeta * (1.0 + zeta) / (1.0 + omega2));
    let major = if u >= 0.0 { major } else { -major };
    (major, omega * major, zeta)
}

// At the face centre u = v = 0 and ω is taken as 0. Elsewhere |v| ≤ |u|, so
// the angle is within ±15° and the denominator is at least cos 15° − 1/√2.
fn qsc_inverse_omega(u: f64, v: f64) -> f64 {
    if u == 0.0 {
        return 0.0;
    }
    let (sin_ratio, cos_ratio) = sincosd(15.0 * v / u);
    sin_ratio / (cos_ratio - 1.0 / libm::sqrt(2.0_f64))
}

#[cfg(test)]
mod tests {
    use crate::common::assert_pinned;
    use crate::coordinate::{IntermediateCoord, NativeCoord};
    use crate::spherical::Projection;

    // One row on each face, then (180°, −30°) on the seam and (−135°, −30°) on
    // the edge between faces 3 and 4, which goes to face 3. Against a 60-digit
    // evaluation the projected values are within 3.85 ulp and the deprojected
    // ones within 2.79 ulp, each within the summed half-ulp rounding of the
    // steps as evaluated. The largest, x at (20°, −30°), is v of eq. 178, where
    // the rounding of tan⁻¹ ω, of sin⁻¹ and of its argument each count about an
    // ulp. Coming back, θ at (225°, −38°) carries cos(15°v/u) − 1/√2.
    #[test]
    fn test_qsc_pinned_both_ways() {
        let rows = [
            (
                (20.0, -30.0),
                (20.892234256318257, -32.926030420422414),
                (20.000000000000007, -29.999999999999996),
            ),
            (
                (120.0, 15.0),
                (121.26035846906711, 18.078987735415996),
                (120.00000000000001, 14.999999999999998),
            ),
            (
                (180.0, -30.0),
                (180.0, -30.434704104205988),
                (180.0, -29.999999999999996),
            ),
            ((-90.0, -40.0), (-90.0, -40.218376725992606), (-90.0, -40.0)),
            (
                (-135.0, 80.0),
                (-8.531653748553603, 98.53165374855361),
                (-135.00000000000006, 80.0),
            ),
            (
                (20.0, -75.0),
                (5.922641865204558, -75.20462660630123),
                (19.999999999999996, -75.0),
            ),
            (
                (-135.0, -30.0),
                (225.0, -37.99940791954279),
                (-135.0, -29.999999999999982),
            ),
        ];
        assert_pinned(&Projection::qsc(), &rows);
    }

    #[test]
    fn test_qsc_near_face_centre_keeps_the_offset() {
        // 1e-6° from the centres of faces 0 and 1, 1 − ζ is below the rounding
        // of ζ. 2e-9° out, ξ and η are each below 1e-10 and ω is still their
        // ratio. Against a 60-digit evaluation, the 1e-6° face 1 y is 0.85 ulp
        // high and the other 1e-6° values are correctly rounded. At 2e-9°, x is
        // 0.69 ulp low and y 2.19 ulp low: eq. 178 subtracts 18.4° from 26.6°.
        let cases = [
            ((37.5, 89.999999), (7.223694392195705e-7, 89.99999909728538)),
            ((1e-6, 1e-6), (1.2080908165955517e-6, 1.208090816595552e-6)),
            ((2e-9, 1e-9), (2.1509014535795535e-9, 1.1658032647536427e-9)),
        ];
        for ((phi, theta), expected) in cases {
            let native = NativeCoord::new(phi, theta);
            let inter = Projection::qsc().project(native).unwrap();
            assert_eq!((inter.x_deg(), inter.y_deg()), expected, "({phi}, {theta})");
        }
    }

    #[test]
    fn test_qsc_inverse_near_face_centre_keeps_the_offset() {
        // From 1e-6° down to 4e-11° off the centres of faces 0 and 1, where ζ
        // rounds to 1 or nearly. Against eqs. 180 to 185 evaluated to 60 digits,
        // φ is 0.71, 0.73 and 0.48 ulp low and 0.12 ulp high, and θ is 0.17 ulp
        // low, 0.94 ulp high, and 1.16 and 2.10 ulp low.
        let cases = [
            (
                (7.223694392195705e-7, 89.99999909728538),
                (37.5000000989661, 89.999999),
            ),
            (
                (1.2080908165955517e-6, 1.208090816595552e-6),
                (9.999999999999997e-7, 1.0000000000000004e-6),
            ),
            (
                (2.1509014535795535e-9, 1.1658032647536427e-9),
                (1.9999999999999997e-9, 9.999999999999994e-10),
            ),
            (
                (4.301802907159107e-11, 2.331606529507286e-11),
                (4e-11, 1.9999999999999992e-11),
            ),
        ];
        for ((x, y), expected) in cases {
            let inter = IntermediateCoord::new(x, y);
            let native = Projection::qsc().deproject(inter).unwrap();
            assert_eq!(
                (native.phi_deg(), native.theta_deg()),
                expected,
                "({x}, {y})"
            );
        }
    }
}
