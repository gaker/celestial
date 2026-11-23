use crate::common::{native_coord, pole_native_coord};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{atan2d, sincosd, to_degrees, to_radians};

// Paper II eqs. 59 and 60.
pub(crate) fn project_sin(native: NativeCoord, xi: f64, eta: f64) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    if theta == 90.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    let (sin_theta, cos_theta) = sincosd(theta);
    let (sin_phi, cos_phi) = sincosd(native.phi_deg());

    // The limb, θx = −tan⁻¹(ξ sin φ − η cos φ) (eq. 66). Beyond it the far
    // side of the sphere would be projected over the near side.
    if sin_theta + cos_theta * (xi * sin_phi - eta * cos_phi) < 0.0 {
        return Err(WcsError::out_of_bounds(
            "SIN projection: point beyond the limb",
        ));
    }

    let x = to_degrees(cos_theta * sin_phi + xi * (1.0 - sin_theta));
    let y = -to_degrees(cos_theta * cos_phi - eta * (1.0 - sin_theta));
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_sin(inter: IntermediateCoord, xi: f64, eta: f64) -> WcsResult<NativeCoord> {
    let x = to_radians(inter.x_deg());
    let y = to_radians(inter.y_deg());
    if x == 0.0 && y == 0.0 {
        return Ok(pole_native_coord());
    }
    deproject_slant(x, y, xi, eta, "SIN")
}

// Paper II eqs. 38-44, with (X, Y) in radians (eq. 42), which also serve SIN as the µ → ∞ limit of SZP. They
// are solved for s = 1 − sin θ instead of sin θ, because near the origin sin θ
// rounds to 1 and asin of it keeps no digits of 90° − θ. With a from eq. 39
// and b as eq. 40 plus a, s solves a s² − 2b s + x² + y² = 0.
pub(super) fn deproject_slant(
    x: f64,
    y: f64,
    slant_x: f64,
    slant_y: f64,
    code: &str,
) -> WcsResult<NativeCoord> {
    let a = slant_x * slant_x + slant_y * slant_y + 1.0;
    let b = slant_x * x + slant_y * y + 1.0;
    let s = colatitude_root(a, b, x * x + y * y).ok_or_else(|| {
        WcsError::out_of_bounds(format!("Point outside {code} projection boundary"))
    })?;

    let theta = atan2d(1.0 - s, libm::sqrt(s * (2.0 - s)));
    let phi = atan2d(x - slant_x * s, -(y - slant_y * s));
    native_coord(phi, theta)
}

// The valid root nearer the pole, which Paper II picks when both are valid.
// Taking the smaller root as r² / q avoids cancelling b against the square root.
fn colatitude_root(a: f64, b: f64, r2: f64) -> Option<f64> {
    let discriminant = b * b - a * r2;
    if discriminant < 0.0 {
        return None;
    }
    let q = b + libm::copysign(libm::sqrt(discriminant), b);
    let (s1, s2) = (r2 / q, q / a);
    let (near, far) = if s1 <= s2 { (s1, s2) } else { (s2, s1) };
    [near, far].into_iter().find(|s| (0.0..=2.0).contains(s))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Plain SIN and two slant ones. Against a 60-digit evaluation the
    // projected values are within 1.41 ulp (y at (135°, 75°)) and the
    // deprojected ones within 1.61 ulp (θ at (−90°, 30°)).
    #[test]
    fn test_sin_pinned_both_ways() {
        let plain = [
            (
                (30.0, 60.0),
                (14.32394487827058, -24.80980029398064),
                (30.000000000000004, 60.00000000000001),
            ),
            (
                (-90.0, 30.0),
                (-49.61960058796128, -0.0),
                (-90.0, 30.000000000000004),
            ),
            (
                (135.0, 75.0),
                (10.485855415710063, 10.485855415710065),
                (135.0, 75.0),
            ),
        ];
        assert_pinned(&Projection::sin(), &plain);
        let slant = (
            (25.0, 55.0),
            (14.924901520991131, -31.856819063369052),
            (25.0, 55.00000000000001),
        );
        assert_pinned(&Projection::sin_with_params(0.1, -0.2), &[slant]);
        let slant = (
            (60.0, 45.0),
            (35.92543331981606, -19.418039849234262),
            (60.00000000000001, 45.0),
        );
        assert_pinned(&Projection::sin_with_params(0.05, 0.05), &[slant]);
    }

    #[test]
    fn test_sin_near_origin_keeps_the_latitude() {
        // Near the origin sin θ rounds to within an ulp of 1, which leaves asin
        // no digits for 90 − θ. Against a 60-digit evaluation for the radian
        // offsets, the plain SIN φ is 0.81 ulp low and θ 0.53 ulp low; the slant
        // values are correctly rounded.
        let cases = [
            ((0.0, 0.0), (26.565051177077986, 89.99999776393202)),
            ((0.25, -0.125), (26.56505098957799, 89.99999776393203)),
        ];
        for ((xi, eta), expected) in cases {
            let inter = IntermediateCoord::new(1e-6, -2e-6);
            let native = deproject_sin(inter, xi, eta).unwrap();
            let got = (native.phi_deg(), native.theta_deg());
            assert_eq!(got, expected, "(ξ, η) = ({xi}, {eta})");
        }
    }

    #[test]
    fn test_sin_out_of_bounds() {
        // Standard SIN: point outside projection boundary.
        let err = deproject_sin(IntermediateCoord::new(100.0, 100.0), 0.0, 0.0).unwrap_err();
        assert!(
            matches!(err, WcsError::OutOfBounds { .. }),
            "standard: {:?}",
            err
        );

        // Slant SIN: parameter combo that drives sin_theta out of [-1, 1].
        let err = deproject_sin(IntermediateCoord::new(100.0, 100.0), 0.5, 0.5).unwrap_err();
        assert!(
            matches!(err, WcsError::OutOfBounds { .. }),
            "slant: {:?}",
            err
        );
    }

    #[test]
    fn test_sin_rejects_points_beyond_the_limb() {
        // Eq. 66 puts the limb at θx = −tan⁻¹(ξ sin φ − η cos φ): the equator
        // for plain SIN, and for (ξ, η) = (0.1, −0.2) −5.71° at φ = 90° and
        // −11.31° at φ = 0°.
        let cases = [
            ((0.0, 0.0), 37.5, 0.0, true),
            ((0.0, 0.0), 37.5, -1.0, false),
            ((0.1, -0.2), 90.0, -5.5, true),
            ((0.1, -0.2), 90.0, -6.0, false),
            ((0.1, -0.2), 0.0, -11.0, true),
            ((0.1, -0.2), 0.0, -11.5, false),
        ];
        for ((xi, eta), phi, theta, inside) in cases {
            let native = NativeCoord::new(phi, theta);
            let got = project_sin(native, xi, eta);
            let label = format!("(ξ, η) = ({xi}, {eta}) at ({phi}, {theta}): {got:?}");
            assert_eq!(got.is_ok(), inside, "{label}");
            assert!(
                inside || matches!(got, Err(WcsError::OutOfBounds { .. })),
                "{label}"
            );
        }
    }
}
