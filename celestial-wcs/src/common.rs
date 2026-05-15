use celestial_core::constants::RAD_TO_DEG;
use celestial_core::utils::normalize_longitude;
use celestial_core::Angle;

use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

#[inline]
pub fn asin_safe(sin_value: f64) -> f64 {
    libm::asin(sin_value.clamp(-1.0, 1.0))
}

#[inline]
pub fn pole_native_coord() -> NativeCoord {
    NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0))
}

#[inline]
pub fn radial_to_intermediate(r_theta: f64, phi_rad: f64) -> IntermediateCoord {
    let r_deg = r_theta * RAD_TO_DEG;
    let (ps, pc) = libm::sincos(phi_rad);
    let x = r_deg * ps;
    let y = -r_deg * pc;
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
    let r_theta = libm::sqrt(x_rad * x_rad + y_rad * y_rad);
    let is_pole = r_theta == 0.0;
    let phi_rad = if is_pole {
        0.0
    } else {
        libm::atan2(x_rad, -y_rad)
    };
    (phi_rad, r_theta, is_pole)
}

pub fn project_conic_xy(r_theta: f64, y0: f64, c: f64, phi: f64) -> IntermediateCoord {
    let (c_phi_s, c_phi_c) = libm::sincos(c * phi);
    let x = r_theta * c_phi_s * RAD_TO_DEG;
    let y = (y0 - r_theta * c_phi_c) * RAD_TO_DEG;
    IntermediateCoord::new(x, y)
}

pub fn deproject_conic_polar(x_rad: f64, y_rad: f64, y0: f64, theta_a: f64, c: f64) -> (f64, f64) {
    let y_offset = y0 - y_rad;
    let r_unsigned = libm::sqrt(x_rad * x_rad + y_offset * y_offset);
    let phi = libm::atan2(theta_a.signum() * x_rad, theta_a.signum() * y_offset) / c.abs();
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

pub fn newton_raphson_2d<F>(
    target: (f64, f64),
    initial_guess: (f64, f64),
    f: F,
    max_iter: usize,
    tol: f64,
    context: &str,
) -> WcsResult<(f64, f64)>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    let (tx, ty) = target;
    let (mut x, mut y) = initial_guess;

    for _ in 0..max_iter {
        let (fx, fy) = f(x, y);
        let (dx, dy) = (fx - tx, fy - ty);

        if dx.abs() < tol && dy.abs() < tol {
            return Ok((x, y));
        }

        let (j11, j12, j21, j22) = jacobian_2d(&f, x, y);
        let (delta_x, delta_y) = solve_2x2(j11, j12, j21, j22, dx, dy, context)?;

        x -= delta_x;
        y -= delta_y;
    }

    Err(WcsError::convergence_failure(format!(
        "{}: Newton-Raphson did not converge",
        context
    )))
}

fn jacobian_2d<F>(f: &F, x: f64, y: f64) -> (f64, f64, f64, f64)
where
    F: Fn(f64, f64) -> (f64, f64),
{
    const H: f64 = 1e-8;
    let (fx, fy) = f(x, y);
    let (fx_px, fy_px) = f(x + H, y);
    let (fx_py, fy_py) = f(x, y + H);

    let j11 = (fx_px - fx) / H;
    let j12 = (fx_py - fx) / H;
    let j21 = (fy_px - fy) / H;
    let j22 = (fy_py - fy) / H;

    (j11, j12, j21, j22)
}

fn solve_2x2(
    j11: f64,
    j12: f64,
    j21: f64,
    j22: f64,
    b1: f64,
    b2: f64,
    context: &str,
) -> WcsResult<(f64, f64)> {
    let det = j11 * j22 - j12 * j21;
    if det.abs() < 1e-15 {
        return Err(WcsError::convergence_failure(format!(
            "{}: singular Jacobian matrix",
            context
        )));
    }
    let inv_det = 1.0 / det;
    let x = inv_det * (j22 * b1 - j12 * b2);
    let y = inv_det * (-j21 * b1 + j11 * b2);
    Ok((x, y))
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::{HALF_PI, QUARTER_PI, SQRT2};

    #[test]
    fn test_asin_safe_clamps_above_and_below_unit() {
        assert_eq!(asin_safe(1.0000000001), HALF_PI);
        assert_eq!(asin_safe(-1.0000000001), -HALF_PI);
    }

    #[test]
    fn test_pole_native_coord() {
        let pole = pole_native_coord();
        assert_eq!(pole.phi().degrees(), 0.0);
        assert_eq!(pole.theta().degrees(), 90.0);
    }

    #[test]
    fn test_radial_to_intermediate() {
        // At the origin the projection is (0, 0); the assertion is byte-exact.
        let zero = radial_to_intermediate(0.0, 0.0);
        assert_eq!(zero.x_deg(), 0.0);
        assert_eq!(zero.y_deg(), 0.0);

        // Off-origin: x = r * sin(phi) * RAD_TO_DEG, y = -r * cos(phi) * RAD_TO_DEG.
        let inter = radial_to_intermediate(1.0, QUARTER_PI);
        assert!((inter.x_deg() - libm::sin(QUARTER_PI) * RAD_TO_DEG).abs() < 1e-10);
        assert!((inter.y_deg() + libm::cos(QUARTER_PI) * RAD_TO_DEG).abs() < 1e-10);
    }

    #[test]
    fn test_native_coord_from_radians() {
        let native = native_coord_from_radians(QUARTER_PI, QUARTER_PI);
        assert!((native.phi().degrees() - 45.0).abs() < 1e-10);
        assert!((native.theta().degrees() - 45.0).abs() < 1e-10);
    }

    #[test]
    fn test_check_nonzero_param() {
        assert!(check_nonzero_param(1.0, "test").is_ok());
        assert!(check_nonzero_param(-0.5, "test").is_ok());
        assert!(check_nonzero_param(0.0, "test").is_err());
        assert!(check_nonzero_param(1e-11, "test").is_err());
    }

    #[test]
    fn test_intermediate_to_polar() {
        // Origin is flagged as the pole.
        let (phi, r, is_pole) = intermediate_to_polar(0.0, 0.0);
        assert_eq!(phi, 0.0);
        assert_eq!(r, 0.0);
        assert!(is_pole);

        // Off-origin: (x=1, y=-1) gives phi = atan2(1, 1) = pi/4 and r = sqrt(2).
        let (phi, r, is_pole) = intermediate_to_polar(1.0, -1.0);
        assert!((phi - QUARTER_PI).abs() < 1e-10);
        assert!((r - SQRT2).abs() < 1e-10);
        assert!(!is_pole);
    }

    #[test]
    fn test_project_conic_xy() {
        let inter = project_conic_xy(1.0, 0.5, 0.5, 0.0);
        assert!(inter.x_deg().abs() < 1e-10);
        assert!((inter.y_deg() - (-0.5 * RAD_TO_DEG)).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_1d_converges() {
        // Linear: f(x) = 2x; solving f(x) = 5 gives x = 2.5.
        let config = NewtonConfig::new((-10.0, 10.0), "linear");
        let x = newton_raphson_1d(0.0, 5.0, |x| 2.0 * x, |_| 2.0, &config).unwrap();
        assert!((x - 2.5).abs() < 1e-10);

        // Quadratic: f(x) = x^2 on [0, 10]; solving f(x) = 4 gives x = 2.
        let config = NewtonConfig::new((0.0, 10.0), "quadratic");
        let x = newton_raphson_1d(1.0, 4.0, |x| x * x, |x| 2.0 * x, &config).unwrap();
        assert!((x - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_2d_recovers_distortion() {
        // Two distortion shapes - separable quadratic and mixed cross-term -
        // exercise both diagonal and off-diagonal Jacobian terms.
        let cases: &[(fn(f64, f64) -> (f64, f64), f64, f64)] = &[
            (|x, y| (x + 0.001 * x * x, y + 0.001 * y * y), 100.0, 200.0),
            (|x, y| (x + 0.0001 * x * y, y + 0.0002 * x * x), 50.0, 75.0),
        ];

        for &(distort, orig_x, orig_y) in cases {
            let (dist_x, dist_y) = distort(orig_x, orig_y);
            let (x, y) = newton_raphson_2d(
                (dist_x, dist_y),
                (dist_x, dist_y),
                distort,
                20,
                1e-12,
                "test",
            )
            .unwrap();
            assert!((x - orig_x).abs() < 1e-10);
            assert!((y - orig_y).abs() < 1e-10);
        }
    }

    #[test]
    fn test_newton_raphson_2d_singular_jacobian_errors() {
        // A constant function has a zero Jacobian and must surface a
        // ConvergenceFailure rather than divide by zero.
        let constant = |_x: f64, _y: f64| (5.0, 5.0);
        let err = newton_raphson_2d((0.0, 0.0), (0.0, 0.0), constant, 20, 1e-12, "test")
            .unwrap_err();
        assert!(matches!(err, WcsError::ConvergenceFailure { .. }), "got: {:?}", err);
    }
}
