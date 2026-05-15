use celestial_core::constants::{DEG_TO_RAD, HALF_PI};

use crate::common::{
    intermediate_to_polar, native_coord_from_radians, pole_native_coord,
    radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_air(native: NativeCoord, theta_b: f64) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    if theta == HALF_PI {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    if theta <= -HALF_PI + 1e-10 {
        return Err(WcsError::singularity(
            "AIR projection diverges at theta = -90",
        ));
    }

    let r_theta = compute_air_r_theta(theta, theta_b)?;
    Ok(radial_to_intermediate(r_theta, phi))
}

fn compute_air_r_theta(theta: f64, theta_b: f64) -> WcsResult<f64> {
    let xi = (HALF_PI - theta) / 2.0;
    let xi_b = (HALF_PI - theta_b * DEG_TO_RAD) / 2.0;

    if xi.abs() < 1e-15 {
        return Ok(0.0);
    }

    let cos_xi = libm::cos(xi);
    let tan_xi = libm::tan(xi);

    let term1 = if cos_xi > 0.0 {
        libm::log(cos_xi) / tan_xi
    } else {
        return Err(WcsError::singularity("AIR projection: cos(xi) <= 0"));
    };

    let term2 = if xi_b.abs() < 1e-10 {
        -0.5 * tan_xi
    } else {
        let cos_xi_b = libm::cos(xi_b);
        let tan_xi_b = libm::tan(xi_b);
        if cos_xi_b > 0.0 && tan_xi_b.abs() > 1e-15 {
            libm::log(cos_xi_b) * tan_xi / (tan_xi_b * tan_xi_b)
        } else {
            return Err(WcsError::singularity(
                "AIR projection: invalid theta_b parameter",
            ));
        }
    };

    Ok(-2.0 * (term1 + term2))
}

pub(crate) fn deproject_air(inter: IntermediateCoord, theta_b: f64) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let (phi, r, is_pole) = intermediate_to_polar(x, y);

    if is_pole {
        return Ok(pole_native_coord());
    }

    let theta = solve_air_inverse(r, theta_b)?;

    Ok(native_coord_from_radians(phi, theta))
}

fn solve_air_inverse(r: f64, theta_b: f64) -> WcsResult<f64> {
    const MAX_ITER: usize = 50;
    const RESIDUAL_TOL: f64 = 5e-16;
    const LOWER_BOUND: f64 = -HALF_PI + 0.01;

    let zea_guess = HALF_PI - 2.0 * libm::asin((r / 2.0).clamp(-1.0, 1.0));
    let mut theta = zea_guess.clamp(LOWER_BOUND, HALF_PI);
    let mut prev_abs_f = f64::INFINITY;

    for _ in 0..MAX_ITER {
        let r_theta = compute_air_r_theta(theta, theta_b)?;
        let f = r_theta - r;
        let abs_f = f.abs();

        if abs_f < RESIDUAL_TOL {
            return Ok(theta);
        }
        if abs_f >= prev_abs_f {
            return Ok(theta);
        }
        prev_abs_f = abs_f;

        let f_prime = compute_air_dr_dtheta(theta, theta_b)?;
        if f_prime.abs() < 1e-15 {
            return Err(WcsError::convergence_failure(
                "AIR inverse: derivative too small",
            ));
        }

        theta = (theta - f / f_prime).clamp(LOWER_BOUND, HALF_PI);
    }

    Err(WcsError::convergence_failure(
        "AIR inverse: Newton-Raphson did not converge",
    ))
}

fn compute_air_dr_dtheta(theta: f64, theta_b: f64) -> WcsResult<f64> {
    let xi = (HALF_PI - theta) / 2.0;
    let xi_b = (HALF_PI - theta_b * DEG_TO_RAD) / 2.0;

    if xi.abs() < 1e-15 {
        return Ok(0.0);
    }

    let sin_xi = libm::sin(xi);
    let cos_xi = libm::cos(xi);
    if cos_xi <= 0.0 {
        return Err(WcsError::singularity("AIR derivative: cos(xi) <= 0"));
    }

    let term1 = -1.0 - libm::log(cos_xi) / (sin_xi * sin_xi);

    let term2 = if xi_b.abs() < 1e-10 {
        -0.5 / (cos_xi * cos_xi)
    } else {
        let cos_xi_b = libm::cos(xi_b);
        let tan_xi_b = libm::tan(xi_b);
        if cos_xi_b <= 0.0 || tan_xi_b.abs() < 1e-15 {
            return Err(WcsError::singularity(
                "AIR derivative: invalid theta_b parameter",
            ));
        }
        libm::log(cos_xi_b) / (tan_xi_b * tan_xi_b * cos_xi * cos_xi)
    };

    Ok(term1 + term2)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_air_roundtrip() {
        for theta_b in [30.0, 45.0, 60.0, 75.0, 90.0] {
            let proj = Projection::air(theta_b);
            for phi_deg in [-180.0, -90.0, 0.0, 30.0, 90.0, 135.0, 180.0] {
                for theta_deg in [10.0, 30.0, 45.0, 60.0, 75.0, 85.0, 89.0] {
                    let original = NativeCoord::new(
                        Angle::from_degrees(phi_deg),
                        Angle::from_degrees(theta_deg),
                    );
                    let inter = proj.project(original).unwrap();
                    let recovered = proj.deproject(inter).unwrap();

                    let ulp_bar = if theta_deg >= 85.0 { 64 } else { 20 };
                    assert_ulp_lt!(
                        original.phi().degrees(),
                        recovered.phi().degrees(),
                        ulp_bar,
                        "phi (theta_b={}, phi={}, theta={})",
                        theta_b,
                        phi_deg,
                        theta_deg
                    );
                    assert_ulp_lt!(
                        original.theta().degrees(),
                        recovered.theta().degrees(),
                        ulp_bar,
                        "theta (theta_b={}, phi={}, theta={})",
                        theta_b,
                        phi_deg,
                        theta_deg
                    );
                }
            }
        }
    }

    #[test]
    fn test_air_edge_cases() {
        // South pole: AIR diverges at theta = -90.
        let proj = Projection::air(90.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        assert!(proj.project(native).is_err());

        // theta near pi/2 takes the xi ≈ 0 early return — result must be tiny.
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(89.9999999));
        let inter = project_air(native, 90.0).unwrap();
        assert!(inter.x_deg().abs() < 1e-6);
        assert!(inter.y_deg().abs() < 1e-6);

        // Invalid theta_b (outside valid range) errors.
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let err = project_air(native, -100.0).unwrap_err();
        assert!(matches!(err, WcsError::Singularity { .. }), "invalid theta_b: {:?}", err);
    }
}
