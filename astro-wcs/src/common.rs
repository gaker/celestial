use astro_core::constants::RAD_TO_DEG;
use astro_core::utils::normalize_longitude;
use astro_core::Angle;

use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

#[inline]
pub fn asin_safe(sin_value: f64) -> f64 {
    sin_value.clamp(-1.0, 1.0).asin()
}

#[inline]
pub fn pole_native_coord() -> NativeCoord {
    NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0))
}

#[inline]
pub fn radial_to_intermediate(r_theta: f64, phi_rad: f64) -> IntermediateCoord {
    let (ps, pc) = phi_rad.sin_cos();
    let x = r_theta * ps * RAD_TO_DEG;
    let y = -r_theta * pc * RAD_TO_DEG;
    IntermediateCoord::new(x, y)
}

#[inline]
pub fn native_coord_from_radians(phi_rad: f64, theta_rad: f64) -> NativeCoord {
    let phi_deg = normalize_longitude(phi_rad * RAD_TO_DEG);
    NativeCoord::new(
        Angle::from_degrees(phi_deg),
        Angle::from_degrees(theta_rad * RAD_TO_DEG),
    )
}

#[inline]
pub fn check_nonzero_param(value: f64, context: &str) -> WcsResult<()> {
    if value.abs() < 1e-10 {
        return Err(WcsError::invalid_parameter(format!(
            "{}: parameter cannot be zero",
            context
        )));
    }
    Ok(())
}

#[inline]
pub fn intermediate_to_polar(x_rad: f64, y_rad: f64) -> (f64, f64, bool) {
    let r_theta = (x_rad * x_rad + y_rad * y_rad).sqrt();
    let is_pole = r_theta == 0.0;
    let phi_rad = if is_pole { 0.0 } else { x_rad.atan2(-y_rad) };
    (phi_rad, r_theta, is_pole)
}

pub fn project_conic_xy(r_theta: f64, y0: f64, c: f64, phi: f64) -> IntermediateCoord {
    let (c_phi_s, c_phi_c) = (c * phi).sin_cos();
    let x = r_theta * c_phi_s * RAD_TO_DEG;
    let y = (y0 - r_theta * c_phi_c) * RAD_TO_DEG;
    IntermediateCoord::new(x, y)
}

pub fn deproject_conic_polar(x_rad: f64, y_rad: f64, y0: f64, theta_a: f64) -> (f64, f64) {
    let y_offset = y0 - y_rad;
    let r_unsigned = (x_rad * x_rad + y_offset * y_offset).sqrt();
    let c = theta_a.sin();
    let phi = (theta_a.signum() * x_rad).atan2(theta_a.signum() * y_offset) / c.abs();
    (phi, r_unsigned)
}

/// Configuration for Newton-Raphson 1D solver
pub struct NewtonConfig {
    pub bounds: (f64, f64),
    pub max_iter: usize,
    pub tol: f64,
    pub context: &'static str,
}

impl NewtonConfig {
    pub const DEFAULT_MAX_ITER: usize = 50;
    pub const DEFAULT_TOL: f64 = 1e-12;

    pub const fn new(bounds: (f64, f64), context: &'static str) -> Self {
        Self {
            bounds,
            max_iter: Self::DEFAULT_MAX_ITER,
            tol: Self::DEFAULT_TOL,
            context,
        }
    }
}

pub fn newton_raphson_1d<F, FP>(
    initial: f64,
    target: f64,
    f: F,
    f_prime: FP,
    config: &NewtonConfig,
) -> WcsResult<f64>
where
    F: Fn(f64) -> f64,
    FP: Fn(f64) -> f64,
{
    let mut x = initial.clamp(config.bounds.0, config.bounds.1);

    for _ in 0..config.max_iter {
        let f_val = f(x) - target;
        let f_prime_val = f_prime(x);

        if f_prime_val.abs() < 1e-15 {
            return Err(WcsError::convergence_failure(format!(
                "{}: derivative too small",
                config.context
            )));
        }

        let delta = f_val / f_prime_val;
        x -= delta;
        x = x.clamp(config.bounds.0, config.bounds.1);

        if delta.abs() < config.tol {
            return Ok(x);
        }
    }

    Err(WcsError::convergence_failure(format!(
        "{}: Newton-Raphson did not converge",
        config.context
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_4;

    #[test]
    fn test_asin_safe_clamping() {
        assert_eq!(asin_safe(1.0000000001), std::f64::consts::FRAC_PI_2);
        assert_eq!(asin_safe(-1.0000000001), -std::f64::consts::FRAC_PI_2);
        assert_eq!(asin_safe(0.5), 0.5_f64.asin());
    }

    #[test]
    fn test_pole_native_coord() {
        let pole = pole_native_coord();
        assert_eq!(pole.phi().degrees(), 0.0);
        assert_eq!(pole.theta().degrees(), 90.0);
    }

    #[test]
    fn test_radial_to_intermediate_at_origin() {
        let inter = radial_to_intermediate(0.0, 0.0);
        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(inter.y_deg(), 0.0);
    }

    #[test]
    fn test_radial_to_intermediate() {
        let inter = radial_to_intermediate(1.0, FRAC_PI_4);
        assert!((inter.x_deg() - FRAC_PI_4.sin() * RAD_TO_DEG).abs() < 1e-10);
        assert!((inter.y_deg() + FRAC_PI_4.cos() * RAD_TO_DEG).abs() < 1e-10);
    }

    #[test]
    fn test_native_coord_from_radians() {
        let native = native_coord_from_radians(FRAC_PI_4, FRAC_PI_4);
        assert!((native.phi().degrees() - 45.0).abs() < 1e-10);
        assert!((native.theta().degrees() - 45.0).abs() < 1e-10);
    }

    #[test]
    fn test_check_nonzero_param_pass() {
        assert!(check_nonzero_param(1.0, "test").is_ok());
        assert!(check_nonzero_param(-0.5, "test").is_ok());
    }

    #[test]
    fn test_check_nonzero_param_fail() {
        assert!(check_nonzero_param(0.0, "test").is_err());
        assert!(check_nonzero_param(1e-11, "test").is_err());
    }

    #[test]
    fn test_intermediate_to_polar_at_origin() {
        let (phi, r, is_pole) = intermediate_to_polar(0.0, 0.0);
        assert_eq!(phi, 0.0);
        assert_eq!(r, 0.0);
        assert!(is_pole);
    }

    #[test]
    fn test_intermediate_to_polar_nonzero() {
        let (phi, r, is_pole) = intermediate_to_polar(1.0, -1.0);
        assert!((phi - FRAC_PI_4).abs() < 1e-10);
        assert!((r - std::f64::consts::SQRT_2).abs() < 1e-10);
        assert!(!is_pole);
    }

    #[test]
    fn test_project_conic_xy() {
        let inter = project_conic_xy(1.0, 0.5, 0.5, 0.0);
        assert!(inter.x_deg().abs() < 1e-10);
        assert!((inter.y_deg() - (-0.5 * RAD_TO_DEG)).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_1d_linear() {
        let config = NewtonConfig::new((-10.0, 10.0), "test");
        let result = newton_raphson_1d(0.0, 5.0, |x| 2.0 * x, |_| 2.0, &config);
        assert!(result.is_ok());
        assert!((result.unwrap() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_1d_quadratic() {
        let config = NewtonConfig::new((0.0, 10.0), "test");
        let result = newton_raphson_1d(1.0, 4.0, |x| x * x, |x| 2.0 * x, &config);
        assert!(result.is_ok());
        assert!((result.unwrap() - 2.0).abs() < 1e-10);
    }
}
