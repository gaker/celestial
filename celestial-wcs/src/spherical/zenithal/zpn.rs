use celestial_core::constants::PI;

use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::newton::{newton_raphson_1d, NewtonConfig};
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

pub(crate) fn deproject_zpn(inter: IntermediateCoord, coeffs: &[f64]) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    let p0 = constant_term(coeffs)?;

    if is_pole && p0 == 0.0 {
        return Ok(pole_native_coord());
    }

    let theta = solve_zpn_inverse(to_radians(r_theta), coeffs)?;

    native_coord(phi, theta)
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

fn evaluate_polynomial(theta: f64, coeffs: &[f64]) -> f64 {
    let mut result = 0.0;
    for coeff in coeffs.iter().rev() {
        result = result * theta + coeff;
    }
    result
}

fn evaluate_polynomial_derivative(theta: f64, coeffs: &[f64]) -> f64 {
    let mut result = 0.0;
    for (i, coeff) in coeffs.iter().enumerate().skip(1).rev() {
        result = result * theta + (i as f64) * coeff;
    }
    result
}

fn solve_zpn_inverse(r: f64, coeffs: &[f64]) -> WcsResult<f64> {
    const CONFIG: NewtonConfig = NewtonConfig::new((0.0, PI), "ZPN inverse");
    let s = newton_raphson_1d(
        r.clamp(0.0, PI),
        r,
        |s| evaluate_polynomial(s, coeffs),
        |s| evaluate_polynomial_derivative(s, coeffs),
        &CONFIG,
    )?;
    Ok(90.0 - to_degrees(s))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;
    use celestial_core::constants::RAD_TO_DEG;

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
        assert!(deproject_zpn(origin, &coeffs).is_err());
    }

    // Paper II eq. 68 with every Pm (m ≥ 1) zero puts every θ at the same R,
    // which has no inverse.
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
            let back = deproject_zpn(inter, &trimmed).unwrap();
            assert_eq!(deproject_zpn(inter, &padded).unwrap(), back);
        }
    }

    #[test]
    fn test_zpn_error_cases() {
        // Project: empty coefficients
        let proj = Projection::zpn(vec![]);
        let native = NativeCoord::new(0.0, 45.0);
        assert!(proj.project(native).is_err());

        // Deproject: empty coefficients
        let err = deproject_zpn(IntermediateCoord::new(10.0, 10.0), &[]).unwrap_err();
        assert!(
            matches!(err, WcsError::InvalidParameter { .. }),
            "empty: {:?}",
            err
        );

        // Deproject: derivative too small triggers convergence failure
        let err = deproject_zpn(IntermediateCoord::new(10.0, 10.0), &[0.5, 1e-20]).unwrap_err();
        assert!(
            matches!(err, WcsError::ConvergenceFailure { .. }),
            "tiny deriv: {:?}",
            err
        );
    }
}
