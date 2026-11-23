use celestial_core::constants::{PI, SQRT2};

use crate::common::native_coord;
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::newton::bracketed_newton;
use crate::trig::{asind, atan2d, sincosd, sind, to_degrees, to_radians};

// Paper II eqs. 90 and 91.
pub(crate) fn project_sfl(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    let x = native.phi_deg() * sincosd(theta).1;
    Ok(IntermediateCoord::new(x, theta))
}

// Eqs. 92 and 93.
pub(crate) fn deproject_sfl(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let y = inter.y_deg();
    native_coord(parallel_longitude(inter.x_deg(), sincosd(y).1)?, y)
}

// Each projection here draws the parallel at θ as x = φ · scale. At a pole
// that parallel is a point, so φ has no value and 0 stands for it, and any
// other x there is off the map.
pub(super) fn parallel_longitude(x: f64, scale: f64) -> WcsResult<f64> {
    if scale != 0.0 {
        Ok(x / scale)
    } else if x == 0.0 {
        Ok(0.0)
    } else {
        Err(WcsError::out_of_bounds(
            "Point beside the pole, off the map",
        ))
    }
}

// Eqs. 94 and 95.
pub(crate) fn project_par(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    let x = native.phi_deg() * (2.0 * sincosd(2.0 * theta / 3.0).1 - 1.0);
    Ok(IntermediateCoord::new(x, 180.0 * sind(theta / 3.0)))
}

// Eqs. 96 and 97. Eq. 97 as printed carries a factor 180°/π that eq. 94 does
// not. This is eq. 94 inverted, with 2 cos(2θ/3) − 1 = 1 − 4 sin²(θ/3).
pub(crate) fn deproject_par(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let sin_third = inter.y_deg() / 180.0;
    if sin_third.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "PAR deprojection: |y| > 180 degrees",
        ));
    }
    let scale = 1.0 - 4.0 * sin_third * sin_third;
    let phi = parallel_longitude(inter.x_deg(), scale)?;
    native_coord(phi, 3.0 * asind(sin_third))
}

// Eqs. 98 and 99.
pub(crate) fn project_mol(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let (sin_gamma, cos_gamma) = mollweide_gamma(native.theta_deg())?;
    let x = 2.0 * SQRT2 / PI * native.phi_deg() * cos_gamma;
    Ok(IntermediateCoord::new(x, SQRT2 * to_degrees(sin_gamma)))
}

// Eq. 100, returning (sin γ, cos γ). Where |γ| ≤ 45°, that is
// |sin θ| ≤ 1/2 + 1/π, it is solved for γ. Nearer the pole sin θ rounds
// toward ±1 and loses 90° − |θ|, and γ itself loses cos γ, so it is solved for
// g = 90° − |γ| instead.
fn mollweide_gamma(theta: f64) -> WcsResult<(f64, f64)> {
    let sin_theta = sind(theta);
    if libm::fabs(sin_theta) <= 0.5 + 1.0 / PI {
        return Ok(sincosd(mollweide_equatorial_gamma(theta, sin_theta)?));
    }
    let (sin_g, cos_g) = sincosd(mollweide_polar_offset(90.0 - libm::fabs(theta))?);
    Ok((libm::copysign(cos_g, theta), sin_g))
}

// sin θ = γ/90° + sin 2γ/π. The right side rises from −1 to 1 as γ runs from
// −90° to 90°.
fn mollweide_equatorial_gamma(theta: f64, sin_theta: f64) -> WcsResult<f64> {
    let f_and_slope = |gamma: f64| {
        let (sin_2gamma, cos_2gamma) = sincosd(2.0 * gamma);
        let f = gamma / 90.0 + sin_2gamma / PI - sin_theta;
        Ok((f, (1.0 + cos_2gamma) / 90.0))
    };
    bracketed_newton((-90.0, 90.0), theta, f_and_slope, "MOL forward")
}

// Eq. 100 subtracted from 1 at |γ| = 90° − g, for t = 90° − |θ|:
// 2 sin²(t/2) = (u − sin u)/π with u = 2g in radians. The right side rises
// with g, with slope 2 sin²(u/2)/90° per degree, and is u³/6π for small u, which
// gives the seed.
fn mollweide_polar_offset(t: f64) -> WcsResult<f64> {
    let sin_half_t = sind(t / 2.0);
    let target = 2.0 * sin_half_t * sin_half_t;
    if target == 0.0 {
        return Ok(0.0);
    }
    let k_and_slope = |g: f64| {
        let u = to_radians(2.0 * g);
        let sin_half_u = libm::sin(u / 2.0);
        let k = u_minus_sin_u(u) / PI - target;
        Ok((k, 2.0 * sin_half_u * sin_half_u / 90.0))
    };
    let seed = to_degrees(libm::cbrt(6.0 * PI * target)) / 2.0;
    bracketed_newton((0.0, 45.0), seed.min(45.0), k_and_slope, "MOL forward")
}

// u − sin u for 0 ≤ u ≤ π/2, which cancels as u → 0. Its Taylor series is
// u³/3! (1 − u²/(4·5) (1 − u²/(6·7) (…))); at u = π/2 the term in u²³ is below
// half an ulp of the sum, so the nesting stops there.
fn u_minus_sin_u(u: f64) -> f64 {
    let u2 = u * u;
    let mut nested = 1.0;
    for k in (2..=11).rev() {
        let n = f64::from(2 * k);
        nested = 1.0 - u2 * nested / (n * (n + 1.0));
    }
    u * u2 / 6.0 * nested
}

// Eqs. 101 and 102, with u = πy/180° = √2 sin γ (eq. 99).
pub(crate) fn deproject_mol(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let y = inter.y_deg();
    let u = to_radians(y);
    if u.abs() > SQRT2 {
        return Err(WcsError::out_of_bounds(
            "MOL deprojection: point outside projection boundary",
        ));
    }
    let (root, theta) = if u.abs() <= 1.0 {
        mollweide_equatorial_inverse(u, y)
    } else {
        mollweide_polar_inverse(u)
    };
    let phi = parallel_longitude(PI * inter.x_deg(), 2.0 * root)?;
    native_coord(phi, theta)
}

// |γ| ≤ 45°: eqs. 101 and 102 as printed, returning (√(2 − u²), θ).
fn mollweide_equatorial_inverse(u: f64, y: f64) -> (f64, f64) {
    let root = libm::sqrt(2.0 - u * u);
    let sin_theta = asind(u / SQRT2) / 90.0 + y / 180.0 * root;
    (root, asind(sin_theta))
}

// |γ| > 45°. Nearer the pole eq. 102's sin θ rounds toward ±1, and eq. 101's
// 2 − u² cancels unless it is formed with one rounding. 2 − u² is 2 cos² γ, so
// g = 90° − |γ| has sin g = √(2 − u²)/√2. Eq. 102 is then the polar form of eq.
// 100 that mollweide_polar_offset solves: 1 − sin|θ| = (v − sin v)/π for v = 2g
// in radians, and 90° − |θ| = 2 sin⁻¹√((1 − sin|θ|)/2).
fn mollweide_polar_inverse(u: f64) -> (f64, f64) {
    let root = libm::sqrt(libm::fma(-u, u, 2.0).max(0.0));
    let g = asind(root / SQRT2);
    let one_minus_sin_theta = u_minus_sin_u(to_radians(2.0 * g)) / PI;
    let t = 2.0 * asind(libm::sqrt(one_minus_sin_theta / 2.0));
    (root, libm::copysign(90.0 - t, u))
}

// Eqs. 103 to 105. 1 + cos θ cos(φ/2) vanishes only at φ = ±360°, which
// Projection::project wraps to 0.
pub(crate) fn project_ait(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let (sin_theta, cos_theta) = sincosd(native.theta_deg());
    let (sin_half_phi, cos_half_phi) = sincosd(native.phi_deg() / 2.0);
    let denom = 1.0 + cos_theta * cos_half_phi;
    let gamma = to_degrees(libm::sqrt(2.0 / denom));
    let x = 2.0 * gamma * cos_theta * sin_half_phi;
    Ok(IntermediateCoord::new(x, gamma * sin_theta))
}

// Eqs. 106 to 108, with u and v the radian forms of x and y.
pub(crate) fn deproject_ait(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let u = to_radians(inter.x_deg());
    let v = to_radians(inter.y_deg());
    let z_sq = 1.0 - (u / 4.0) * (u / 4.0) - (v / 2.0) * (v / 2.0);
    if z_sq < 0.0 {
        return Err(WcsError::out_of_bounds(
            "AIT deprojection: point outside projection boundary",
        ));
    }
    let z = libm::sqrt(z_sq);
    let phi = 2.0 * atan2d(z / 2.0 * u, 2.0 * z_sq - 1.0);
    // v·Z = sin θ, and v²(1 − v²/4) ≤ 1 for every v, so the clamp only removes
    // rounding at the poles.
    native_coord(phi, asind((v * z).clamp(-1.0, 1.0)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Per-projection native_reference checks are covered by
    // spherical::tests::test_all_projections_map_reference_to_origin_and_back.

    #[test]
    fn test_sfl_known_value_and_pins() {
        let proj = Projection::sfl();

        // SFL on the equator collapses to (phi, 0) - an identity in degrees.
        let inter = proj.project(NativeCoord::new(90.0, 0.0)).unwrap();
        assert_eq!((inter.x_deg(), inter.y_deg()), (90.0, 0.0));

        // Against a 60-digit evaluation x is within 0.15 ulp and φ within 0.83
        // ulp. At (180°, 30°) x/cos θ is 5e-15° past the seam; it rounds to an
        // ulp past, so φ comes back wrapped, 2.3e-14° from the truth.
        let rows = [
            (
                (-90.0, -30.0),
                (-77.94228634059948, -30.0),
                (-90.00000000000001, -30.0),
            ),
            ((45.0, 60.0), (22.5, 60.0), (45.0, 60.0)),
            (
                (180.0, 30.0),
                (155.88457268119896, 30.0),
                (-179.99999999999997, 30.0),
            ),
        ];
        assert_pinned(&proj, &rows);
    }

    // Against a 60-digit evaluation the projected values are within 0.86 ulp
    // and the deprojected ones within 1.14 ulp.
    #[test]
    fn test_par_pinned_both_ways() {
        let rows = [
            (
                (-180.0, -60.0),
                (-95.77599952283208, -61.56362579862037),
                (-179.99999999999994, -60.0),
            ),
            (
                (45.0, 30.0),
                (39.572335870731756, 31.25667198004746),
                (45.0, 30.0),
            ),
            (
                (90.0, -30.0),
                (79.14467174146351, -31.25667198004746),
                (90.0, -30.0),
            ),
        ];
        assert_pinned(&Projection::par(), &rows);
    }

    // Against a 60-digit evaluation the projected values are within 3.02 ulp,
    // x at (−90°, −60°), and the deprojected ones within 1.39 ulp.
    #[test]
    fn test_mol_pinned_both_ways() {
        let rows = [
            (
                (-90.0, -60.0),
                (-52.43534023595313, -61.77497708913637),
                (-90.00000000000001, -60.0),
            ),
            (
                (45.0, 30.0),
                (37.061247270088934, 32.73329349706183),
                (45.0, 30.000000000000004),
            ),
            ((90.0, 0.0), (81.02846845413956, 0.0), (90.0, 0.0)),
        ];
        assert_pinned(&Projection::mol(), &rows);
    }

    // Near the pole eq. 100's sin θ rounds to 1 and loses 90° − θ. Against a
    // 60-digit evaluation x at (37.5°, 89.999999°) is 1.03 ulp high and at
    // (120°, −89.9999999999°) 1.31 ulp high.
    #[test]
    fn test_mol_near_pole_keeps_x() {
        let x = |phi, theta| {
            let native = NativeCoord::new(phi, theta);
            Projection::mol().project(native).unwrap().x_deg()
        };
        assert_eq!(x(37.5, 89.999999), 2.399228875532676e-4);
        assert_eq!(x(120.0, -89.9999999999), 1.654093896154317e-6);
    }

    // Near the pole eq. 102's sin θ rounds toward 1 and eq. 101's 2 − u²
    // cancels. Against a 60-digit evaluation θ is within 30 ulp. φ there is
    // proportional to 1/√(2 − u²), so it is 1.5e9 ulp (1.0e-5°) off at
    // 89.999999° and 217 ulp off at 89.9°; half an ulp of y alone moves it 6
    // and 8 times as far.
    #[test]
    fn test_mol_inverse_near_pole() {
        let native = |x, y| {
            let n = Projection::mol()
                .deproject(IntermediateCoord::new(x, y))
                .unwrap();
            (n.phi_deg(), n.theta_deg())
        };
        assert_eq!(
            native(2.399228875532676e-4, 81.02846845209359),
            (37.50008107874383, 89.99999900000324)
        );
        assert_eq!(
            native(0.45485970864058056, 81.01897181574569),
            (32.999999999998515, 89.9)
        );
    }

    // Against a 60-digit evaluation the projected values are within 1.44 ulp
    // and the deprojected ones within 2.06 ulp.
    #[test]
    fn test_ait_pinned_both_ways() {
        let rows = [
            (
                (-150.0, 60.0),
                (-73.64710601305798, 66.03018884746814),
                (-150.0, 59.99999999999999),
            ),
            (
                (-90.0, -60.0),
                (-49.24758482030464, -60.315726937090375),
                (-89.99999999999997, -59.99999999999999),
            ),
            (
                (-45.0, -60.0),
                (-25.645581509177468, -58.036808503858886),
                (-44.99999999999999, -59.99999999999998),
            ),
            (
                (90.0, -30.0),
                (78.15388603121195, -31.906190365349868),
                (89.99999999999999, -29.999999999999996),
            ),
        ];
        assert_pinned(&Projection::ait(), &rows);
    }

    // Every parallel shrinks to a point at a pole, where x = 0 for any φ.
    // Coming back, φ has no value there and 0 stands for it.
    #[test]
    fn test_pseudocylindrical_poles_are_exact() {
        // SFL and PAR put the pole at y = 90°. MOL and AIT put it at
        // √2·180°/π = 81.028468454139546…, and this is 0.71 ulp high because √2
        // and 180°/π both round up.
        let y_pole = 81.02846845413956;
        let cases = [
            (Projection::sfl(), 90.0),
            (Projection::par(), 90.0),
            (Projection::mol(), y_pole),
            (Projection::ait(), y_pole),
        ];
        for (proj, y) in cases {
            for sign in [1.0, -1.0] {
                let inter = proj.project(NativeCoord::new(37.5, sign * 90.0)).unwrap();
                assert_eq!((inter.x_deg(), inter.y_deg()), (0.0, sign * y), "{proj:?}");
                let native = proj
                    .deproject(IntermediateCoord::new(0.0, sign * y))
                    .unwrap();
                let back = (native.phi_deg(), native.theta_deg());
                assert_eq!(back, (0.0, sign * 90.0), "{proj:?}");
            }
        }
    }

    #[test]
    fn test_pseudocylindrical_error_cases() {
        // SFL deproject is singular at |y| = 90 (cos(theta) = 0).
        let sfl = Projection::sfl();
        assert!(sfl.deproject(IntermediateCoord::new(10.0, 90.0)).is_err());

        // PAR's θ = 3 sin⁻¹(y/180°) (Paper II eq. 96) leaves the sphere past
        // |y| = 90, and sin⁻¹ has no value past |y| = 180.
        let par = Projection::par();
        assert!(par.deproject(IntermediateCoord::new(0.0, 100.0)).is_err());
        assert!(par.deproject(IntermediateCoord::new(0.0, 200.0)).is_err());

        // MOL deproject rejects points outside the projection ellipse.
        let mol = Projection::mol();
        let outside = 91.0;
        assert!(mol.deproject(IntermediateCoord::new(0.0, outside)).is_err());

        // AIT deproject rejects points outside the projection ellipse.
        let ait = Projection::ait();
        assert!(ait.deproject(IntermediateCoord::new(400.0, 200.0)).is_err());
    }
}
