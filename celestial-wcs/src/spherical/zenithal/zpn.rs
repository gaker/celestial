use celestial_core::constants::{DEG_TO_RAD, HALF_PI, PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, newton_raphson_1d, pole_native_coord,
    radial_to_intermediate, NewtonConfig,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_zpn(native: NativeCoord, coeffs: &[f64]) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta == HALF_PI {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    if coeffs.is_empty() {
        return Err(WcsError::invalid_parameter(
            "ZPN projection requires at least one coefficient",
        ));
    }

    let r_theta = evaluate_polynomial(HALF_PI - theta, coeffs);

    Ok(radial_to_intermediate(r_theta, phi))
}

pub(crate) fn deproject_zpn(inter: IntermediateCoord, coeffs: &[f64]) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let (phi, r, is_pole) = intermediate_to_polar(x, y);

    if is_pole {
        return Ok(pole_native_coord());
    }

    if coeffs.is_empty() {
        return Err(WcsError::invalid_parameter(
            "ZPN projection requires at least one coefficient",
        ));
    }

    if coeffs.len() == 1 {
        if (r - coeffs[0]).abs() > 1e-10 {
            return Err(WcsError::out_of_bounds(
                "ZPN with constant coefficient: R does not match",
            ));
        }
        return Ok(pole_native_coord());
    }

    let theta = solve_zpn_inverse(r, coeffs)?;

    Ok(native_coord_from_radians(phi, theta))
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
    Ok(HALF_PI - s)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::constants::RAD_TO_DEG;
    use celestial_core::Angle;
    #[test]
    fn test_zpn_arc_equivalent() {
        let zpn = Projection::zpn(vec![0.0, 1.0]);
        let arc = Projection::arc();

        for &(phi_deg, theta_deg) in &[
            (0.0, 80.0),
            (45.0, 60.0),
            (90.0, 30.0),
            (-30.0, 10.0),
            (120.0, -45.0),
        ] {
            let native = NativeCoord::new(
                Angle::from_degrees(phi_deg),
                Angle::from_degrees(theta_deg),
            );

            let zpn_inter = zpn.project(native).unwrap();
            let arc_inter = arc.project(native).unwrap();

            assert_ulp_lt!(zpn_inter.x_deg(), arc_inter.x_deg(), 2);
            assert_ulp_lt!(zpn_inter.y_deg(), arc_inter.y_deg(), 2);
        }
    }

    #[test]
    fn test_zpn_roundtrip() {
        // Sweep multiple polynomial orders (linear through 5th-order) across
        // a φ × θ grid. Catches both polynomial evaluation and Newton inverse.
        let coeff_sets: &[(Vec<f64>, u64)] = &[
            (vec![0.0, 1.0], 2),
            (vec![0.0, 1.0, 0.1], 5),
            (vec![0.0, 1.0, 0.05], 10),
            (vec![0.0, 1.0, 0.0, 0.01, 0.0, 0.001], 10),
        ];
        for (coeffs, ulp) in coeff_sets {
            let proj = Projection::zpn(coeffs.clone());
            for phi_deg in [-180.0, -90.0, -60.0, 0.0, 30.0, 45.0, 90.0, 135.0, 180.0] {
                for theta_deg in [30.0, 45.0, 55.0, 60.0, 70.0, 75.0, 85.0] {
                    let original = NativeCoord::new(
                        Angle::from_degrees(phi_deg),
                        Angle::from_degrees(theta_deg),
                    );
                    let inter = proj.project(original).unwrap();
                    let recovered = proj.deproject(inter).unwrap();
                    assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), *ulp,
                        "phi (coeffs={:?}, phi={}, theta={})", coeffs, phi_deg, theta_deg);
                    assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), *ulp,
                        "theta (coeffs={:?}, phi={}, theta={})", coeffs, phi_deg, theta_deg);
                }
            }
        }
    }

    #[test]
    fn test_zpn_error_cases() {
        // Project: empty coefficients
        let proj = Projection::zpn(vec![]);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        assert!(proj.project(native).is_err());

        // Deproject: empty coefficients
        let err = deproject_zpn(IntermediateCoord::new(10.0, 10.0), &[]).unwrap_err();
        assert!(matches!(err, WcsError::InvalidParameter { .. }), "empty: {:?}", err);

        // Deproject: single (constant) coefficient — matching r maps to pole
        let r_deg = 0.5 * RAD_TO_DEG;
        let coord = deproject_zpn(IntermediateCoord::new(0.0, -r_deg), &[0.5]).unwrap();
        assert_eq!(coord.theta().degrees(), 90.0);

        // Deproject: single coefficient — non-matching r is out of bounds
        let err = deproject_zpn(IntermediateCoord::new(10.0, 10.0), &[0.5]).unwrap_err();
        assert!(matches!(err, WcsError::OutOfBounds { .. }), "single mismatch: {:?}", err);

        // Deproject: derivative too small triggers convergence failure
        let err = deproject_zpn(IntermediateCoord::new(10.0, 10.0), &[0.5, 1e-20]).unwrap_err();
        assert!(matches!(err, WcsError::ConvergenceFailure { .. }), "tiny deriv: {:?}", err);
    }
}
