use celestial_core::constants::PI;

use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::newton::bracketed_newton;
use crate::trig::{to_degrees, to_radians};

// Paper II eq. 68, Rθ = (180°/π) Σ Pm ((π/180°)(90° − θ))^m.
pub(crate) fn project_zpn(native: NativeCoord, coeffs: &[f64]) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    let p0 = constant_term(coeffs)?;

    if theta == 90.0 && p0 == 0.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    let r_theta = to_degrees(evaluate_polynomial(to_radians(90.0 - theta), coeffs));

    Ok(radial_to_intermediate(r_theta, native.phi_deg()))
}

// Paper II gives no inverse for eq. 68. Where R takes a radius at more than
// one ζ = 90° − θ, this takes the smallest, the one nearest the pole.
pub(crate) fn deproject_zpn(
    inter: IntermediateCoord,
    coeffs: &[f64],
    piece_ends: &WcsResult<Vec<f64>>,
) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    let p0 = constant_term(coeffs)?;

    if is_pole && p0 == 0.0 {
        return Ok(pole_native_coord());
    }

    let piece_ends = piece_ends.as_ref().map_err(Clone::clone)?;
    let zeta = solve_zpn_inverse(to_radians(r_theta), coeffs, piece_ends)?;

    native_coord(phi, 90.0 - to_degrees(zeta))
}

// The ends of the pieces of [0, π] on which R(ζ) is monotonic: 0, each ζ
// where R′ changes sign, and π, ascending.
pub(crate) fn zpn_piece_ends(coeffs: &[f64]) -> WcsResult<Vec<f64>> {
    let mut ends = vec![0.0];
    ends.extend(sign_changes(&derivative(coeffs), (0.0, PI))?);
    ends.push(PI);
    Ok(ends)
}

// P0 is the radius of the pole. When it is not zero the pole maps to a circle
// and the reference point lies outside the projection (Paper II §5.1.7).
// Without a non-zero Pm for m ≥ 1, R is the same at every θ.
fn constant_term(coeffs: &[f64]) -> WcsResult<f64> {
    match coeffs.split_first() {
        Some((&p0, terms)) if terms.iter().any(|&p| p != 0.0) => Ok(p0),
        _ => Err(WcsError::invalid_parameter(
            "ZPN projection requires a non-zero coefficient beyond P0",
        )),
    }
}

fn evaluate_polynomial(x: f64, coeffs: &[f64]) -> f64 {
    let mut result = 0.0;
    for coeff in coeffs.iter().rev() {
        result = result * x + coeff;
    }
    result
}

fn evaluate_derivative(x: f64, coeffs: &[f64]) -> f64 {
    let mut result = 0.0;
    for (m, coeff) in coeffs.iter().enumerate().skip(1).rev() {
        result = result * x + (m as f64) * coeff;
    }
    result
}

fn derivative(coeffs: &[f64]) -> Vec<f64> {
    coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(m, &p)| m as f64 * p)
        .collect()
}

// The first piece whose values span r holds the smallest root.
fn solve_zpn_inverse(r: f64, coeffs: &[f64], piece_ends: &[f64]) -> WcsResult<f64> {
    for piece in piece_ends.windows(2) {
        let (a, b) = (piece[0], piece[1]);
        let (ra, rb) = (
            evaluate_polynomial(a, coeffs),
            evaluate_polynomial(b, coeffs),
        );
        if r == ra {
            return Ok(a);
        }
        if (ra < r && r <= rb) || (rb <= r && r < ra) {
            let sign = if ra < rb { 1.0 } else { -1.0 };
            let seed = a + (b - a) * ((r - ra) / (rb - ra));
            let f_and_slope = |z| {
                let value = evaluate_polynomial(z, coeffs) - r;
                Ok((sign * value, sign * evaluate_derivative(z, coeffs)))
            };
            return bracketed_newton((a, b), seed, f_and_slope, "ZPN inverse");
        }
    }
    Err(WcsError::out_of_bounds(
        "ZPN projection: radius outside the range of the polynomial",
    ))
}

// The x in (lo, hi) where the polynomial p changes sign, ascending. p is
// monotonic between consecutive sign changes of p′, so each such piece holds
// at most one. A zero where the pieces meet counts only when p has opposite
// signs on the two ends around it.
fn sign_changes(p: &[f64], (lo, hi): (f64, f64)) -> WcsResult<Vec<f64>> {
    if p.len() < 2 {
        return Ok(Vec::new());
    }
    let mut ends = vec![lo];
    ends.extend(sign_changes(&derivative(p), (lo, hi))?);
    ends.push(hi);
    let values: Vec<f64> = ends.iter().map(|&x| evaluate_polynomial(x, p)).collect();
    let mut roots = Vec::new();
    for i in 0..ends.len() - 1 {
        let (va, vb) = (values[i], values[i + 1]);
        if (va < 0.0 && vb > 0.0) || (va > 0.0 && vb < 0.0) {
            roots.push(monotonic_root(p, (ends[i], ends[i + 1]), va < 0.0)?);
        } else if i > 0 && va == 0.0 && values[i - 1] * values[i + 1] < 0.0 {
            roots.push(ends[i]);
        }
    }
    Ok(roots)
}

// The one root of p in a bracket where p is monotonic and changes sign.
fn monotonic_root(p: &[f64], bracket: (f64, f64), rising: bool) -> WcsResult<f64> {
    let sign = if rising { 1.0 } else { -1.0 };
    let f_and_slope = |x| {
        Ok((
            sign * evaluate_polynomial(x, p),
            sign * evaluate_derivative(x, p),
        ))
    };
    let seed = (bracket.0 + bracket.1) / 2.0;
    bracketed_newton(bracket, seed, f_and_slope, "ZPN turning point")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;
    use celestial_core::constants::RAD_TO_DEG;

    fn deproject(inter: IntermediateCoord, coeffs: &[f64]) -> WcsResult<NativeCoord> {
        deproject_zpn(inter, coeffs, &zpn_piece_ends(coeffs))
    }

    // P = (0, 1) is ARC, but eq. 68 goes through radians where ARC takes
    // 90° − θ in degrees, so at (90°, 30°) x is an ulp above ARC's 60. Against
    // a 60-digit evaluation the projected values are within 3.15 ulp (x for
    // the fifth-order P at (45°, 85°)) and the deprojected ones within 3.69
    // ulp (θ at (−30°, 10°), where Rθ = 80 is three binades above θ).
    #[test]
    fn test_zpn_pinned_both_ways() {
        let arc = [
            (
                (90.0, 30.0),
                (60.00000000000001, -0.0),
                (90.0, 29.999999999999993),
            ),
            (
                (-30.0, 10.0),
                (-40.0, -69.28203230275508),
                (-30.000000000000004, 10.0),
            ),
            (
                (120.0, -45.0),
                (116.9134295108992, 67.5),
                (120.00000000000001, -45.0),
            ),
        ];
        assert_pinned(&Projection::zpn(vec![0.0, 1.0]), &arc);
        let quadratic = (
            (45.0, 85.0),
            (3.5663872596699475, -3.566387259669948),
            (45.0, 85.0),
        );
        assert_pinned(&Projection::zpn(vec![0.0, 1.0, 0.1]), &[quadratic]);
        let quadratic = (
            (-180.0, 30.0),
            (0.0, 63.141592653589804),
            (180.0, 29.999999999999993),
        );
        assert_pinned(&Projection::zpn(vec![0.0, 1.0, 0.05]), &[quadratic]);
        let fifth_order = [
            (
                (45.0, 85.0),
                (3.5358033572795207, -3.5358033572795216),
                (44.99999999999999, 85.0),
            ),
            ((-90.0, 55.0), (-35.13547831660654, -0.0), (-90.0, 55.0)),
        ];
        let proj = Projection::zpn(vec![0.0, 1.0, 0.0, 0.01, 0.0, 0.001]);
        assert_pinned(&proj, &fifth_order);
    }

    // R = 2ζ − 0.3ζ³ rises to its maximum of 113.9° at ζ = 85.4° and falls
    // after it, so R = 103° is met at ζ = 62.8° and again at ζ = 106.1°. The
    // inverse takes the first, and no ζ reaches 115°. The 60-digit θ is
    // 27.163331447486943, an ulp above this one: 90° − ζ takes the spacing
    // of ζ near 62.8°, which is two ulp of θ.
    #[test]
    fn test_zpn_inverse_stops_at_the_turning_point() {
        let coeffs = [0.0, 2.0, 0.0, -0.3];
        let inside = deproject(IntermediateCoord::new(0.0, -103.0), &coeffs).unwrap();
        assert_eq!(
            (inside.phi_deg(), inside.theta_deg()),
            (0.0, 27.16333144748694)
        );
        let outside = deproject(IntermediateCoord::new(0.0, -115.0), &coeffs);
        assert!(outside.is_err(), "{outside:?}");
    }

    #[test]
    fn test_zpn_offset_pole_is_a_circle() {
        // With P0 ≠ 0 the native pole maps to a circle of radius P0, and the
        // reference point lies outside the projection (Paper II §5.1.7). The
        // 60-digit values, correctly rounded.
        let coeffs = [0.01, 0.95, 0.02];
        let pole = NativeCoord::new(37.5, 90.0);
        let inter = project_zpn(pole, &coeffs).unwrap();
        assert_eq!(
            (inter.x_deg(), inter.y_deg()),
            (0.34879460612552576, -0.4545579806129398)
        );
        let origin = IntermediateCoord::new(0.0, 0.0);
        assert!(deproject(origin, &coeffs).is_err());
    }

    // Paper II eq. 68 with every Pm (m ≥ 1) zero puts every θ at the same R,
    // which has no inverse.
    // P0 = π/180 puts the pole on the circle of radius 1°, so a point on it
    // deprojects to the pole exactly.
    #[test]
    fn test_zpn_point_on_the_pole_circle() {
        let coeffs = [to_radians(1.0), 1.0];
        let native = deproject(IntermediateCoord::new(1.0, 0.0), &coeffs).unwrap();
        assert_eq!((native.phi_deg(), native.theta_deg()), (90.0, 90.0));
    }

    #[test]
    fn test_degenerate_zpn_is_refused() {
        let native = NativeCoord::new(0.0, 45.0);
        let on_circle = IntermediateCoord::new(0.0, -0.5 * RAD_TO_DEG);
        for coeffs in [vec![0.0], vec![0.5], vec![0.5, 0.0, 0.0]] {
            let proj = Projection::zpn(coeffs.clone());
            let results = [
                proj.project(native).map(|_| ()),
                proj.deproject(on_circle).map(|_| ()),
            ];
            for result in results {
                assert!(
                    matches!(result, Err(WcsError::InvalidParameter { .. })),
                    "{coeffs:?}: {result:?}"
                );
            }
        }
    }

    // Projection::zpn trims trailing zeros, which leaves every value as it was.
    #[test]
    fn test_zpn_trailing_zeros_make_no_difference() {
        let (trimmed, padded) = ([0.0, 1.0, 0.1], [0.0, 1.0, 0.1, 0.0, -0.0]);
        assert_eq!(
            Projection::zpn(padded.to_vec()),
            Projection::zpn(trimmed.to_vec())
        );
        for (phi, theta) in [(45.0, 85.0), (-120.0, 10.0), (180.0, -60.0)] {
            let native = NativeCoord::new(phi, theta);
            let inter = project_zpn(native, &trimmed).unwrap();
            assert_eq!(project_zpn(native, &padded).unwrap(), inter);
            let back = deproject(inter, &trimmed).unwrap();
            assert_eq!(deproject(inter, &padded).unwrap(), back);
        }
    }

    #[test]
    fn test_zpn_error_cases() {
        // Project: empty coefficients
        let proj = Projection::zpn(vec![]);
        let native = NativeCoord::new(0.0, 45.0);
        assert!(proj.project(native).is_err());

        // Deproject: empty coefficients
        let err = deproject(IntermediateCoord::new(10.0, 10.0), &[]).unwrap_err();
        assert!(
            matches!(err, WcsError::InvalidParameter { .. }),
            "empty: {:?}",
            err
        );

        // Deproject: R spans 0.5 rad to 0.5 + 3e-20 rad, so no ζ reaches a
        // radius of 14.1°.
        let err = deproject(IntermediateCoord::new(10.0, 10.0), &[0.5, 1e-20]).unwrap_err();
        assert!(
            matches!(err, WcsError::OutOfBounds { .. }),
            "tiny deriv: {:?}",
            err
        );
    }
}
