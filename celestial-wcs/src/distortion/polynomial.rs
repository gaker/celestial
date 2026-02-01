#[inline]
pub fn horner(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

pub fn chebyshev(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let (mut t_prev, mut t_curr) = (1.0, x);
            for _ in 2..=n {
                let t_next = 2.0 * x * t_curr - t_prev;
                t_prev = t_curr;
                t_curr = t_next;
            }
            t_curr
        }
    }
}

pub fn legendre(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let (mut p_prev, mut p_curr) = (1.0, x);
            for k in 1..n {
                let p_next = ((2 * k + 1) as f64 * x * p_curr - k as f64 * p_prev) / (k + 1) as f64;
                p_prev = p_curr;
                p_curr = p_next;
            }
            p_curr
        }
    }
}

#[inline]
pub fn power_term(x: f64, y: f64, p: u32, q: u32) -> f64 {
    x.powi(p as i32) * y.powi(q as i32)
}

pub fn newton_raphson_2d<F>(
    target: (f64, f64),
    initial_guess: (f64, f64),
    distort_fn: F,
    max_iter: usize,
    tolerance: f64,
) -> Result<(f64, f64), &'static str>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    let (tx, ty) = target;
    let (mut x, mut y) = initial_guess;

    for _ in 0..max_iter {
        let (fx, fy) = distort_fn(x, y);
        let (dx, dy) = (fx - tx, fy - ty);

        if dx.abs() < tolerance && dy.abs() < tolerance {
            return Ok((x, y));
        }

        let (j11, j12, j21, j22) = compute_jacobian(&distort_fn, x, y);
        let (delta_x, delta_y) = solve_2x2(j11, j12, j21, j22, dx, dy)?;

        x -= delta_x;
        y -= delta_y;
    }

    Err("Newton-Raphson failed to converge")
}

fn compute_jacobian<F>(f: &F, x: f64, y: f64) -> (f64, f64, f64, f64)
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
) -> Result<(f64, f64), &'static str> {
    let det = j11 * j22 - j12 * j21;
    if det.abs() < 1e-15 {
        return Err("Singular Jacobian matrix");
    }
    let inv_det = 1.0 / det;
    let x = inv_det * (j22 * b1 - j12 * b2);
    let y = inv_det * (-j21 * b1 + j11 * b2);
    Ok((x, y))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_horner_constant() {
        assert_eq!(horner(&[5.0], 10.0), 5.0);
    }

    #[test]
    fn test_horner_linear() {
        assert_eq!(horner(&[2.0, 3.0], 5.0), 17.0);
    }

    #[test]
    fn test_horner_quadratic() {
        let coeffs = [1.0, 2.0, 3.0];
        let x = 2.0;
        let expected = 1.0 + 2.0 * 2.0 + 3.0 * 4.0;
        assert_eq!(horner(&coeffs, x), expected);
    }

    #[test]
    fn test_horner_cubic() {
        let coeffs = [1.0, -1.0, 2.0, -2.0];
        let x = 3.0;
        let expected = 1.0 - 3.0 + 2.0 * 9.0 - 2.0 * 27.0;
        assert_eq!(horner(&coeffs, x), expected);
    }

    #[test]
    fn test_chebyshev_t0() {
        assert_eq!(chebyshev(0, 0.5), 1.0);
        assert_eq!(chebyshev(0, -0.3), 1.0);
    }

    #[test]
    fn test_chebyshev_t1() {
        assert_eq!(chebyshev(1, 0.5), 0.5);
        assert_eq!(chebyshev(1, -0.7), -0.7);
    }

    #[test]
    fn test_chebyshev_t2() {
        let x = 0.5;
        let expected = 2.0 * x * x - 1.0;
        assert_eq!(chebyshev(2, x), expected);
    }

    #[test]
    fn test_chebyshev_t3() {
        let x: f64 = 0.6;
        let expected = 4.0 * x.powi(3) - 3.0 * x;
        assert!((chebyshev(3, x) - expected).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_t4() {
        let x: f64 = 0.4;
        let expected = 8.0 * x.powi(4) - 8.0 * x.powi(2) + 1.0;
        assert!((chebyshev(4, x) - expected).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_t5() {
        let x: f64 = 0.3;
        let expected = 16.0 * x.powi(5) - 20.0 * x.powi(3) + 5.0 * x;
        assert!((chebyshev(5, x) - expected).abs() < 1e-14);
    }

    #[test]
    fn test_legendre_p0() {
        assert_eq!(legendre(0, 0.5), 1.0);
        assert_eq!(legendre(0, -0.3), 1.0);
    }

    #[test]
    fn test_legendre_p1() {
        assert_eq!(legendre(1, 0.5), 0.5);
        assert_eq!(legendre(1, -0.7), -0.7);
    }

    #[test]
    fn test_legendre_p2() {
        let x = 0.5;
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert_eq!(legendre(2, x), expected);
    }

    #[test]
    fn test_legendre_p3() {
        let x: f64 = 0.6;
        let expected = (5.0 * x.powi(3) - 3.0 * x) / 2.0;
        assert!((legendre(3, x) - expected).abs() < 1e-14);
    }

    #[test]
    fn test_legendre_p4() {
        let x: f64 = 0.4;
        let expected = (35.0 * x.powi(4) - 30.0 * x.powi(2) + 3.0) / 8.0;
        assert!((legendre(4, x) - expected).abs() < 1e-14);
    }

    #[test]
    fn test_legendre_p5() {
        let x: f64 = 0.3;
        let expected = (63.0 * x.powi(5) - 70.0 * x.powi(3) + 15.0 * x) / 8.0;
        assert!((legendre(5, x) - expected).abs() < 1e-14);
    }

    #[test]
    fn test_power_term_basic() {
        assert_eq!(power_term(2.0, 3.0, 0, 0), 1.0);
        assert_eq!(power_term(2.0, 3.0, 1, 0), 2.0);
        assert_eq!(power_term(2.0, 3.0, 0, 1), 3.0);
        assert_eq!(power_term(2.0, 3.0, 1, 1), 6.0);
    }

    #[test]
    fn test_power_term_higher() {
        assert_eq!(power_term(2.0, 3.0, 2, 0), 4.0);
        assert_eq!(power_term(2.0, 3.0, 0, 2), 9.0);
        assert_eq!(power_term(2.0, 3.0, 2, 3), 4.0 * 27.0);
        assert_eq!(power_term(2.0, 3.0, 3, 2), 8.0 * 9.0);
    }

    #[test]
    fn test_newton_raphson_identity() {
        let identity = |x: f64, y: f64| (x, y);
        let result = newton_raphson_2d((3.0, 4.0), (3.0, 4.0), identity, 20, 1e-12);
        let (x, y) = result.unwrap();
        assert!((x - 3.0).abs() < 1e-12);
        assert!((y - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_newton_raphson_quadratic_distortion() {
        let distort = |x: f64, y: f64| (x + 0.001 * x * x, y + 0.001 * y * y);

        let (orig_x, orig_y) = (100.0, 200.0);
        let (dist_x, dist_y) = distort(orig_x, orig_y);

        let result = newton_raphson_2d((dist_x, dist_y), (dist_x, dist_y), distort, 20, 1e-12);
        let (x, y) = result.unwrap();

        assert!((x - orig_x).abs() < 1e-10);
        assert!((y - orig_y).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_mixed_distortion() {
        let distort = |x: f64, y: f64| (x + 0.0001 * x * y, y + 0.0002 * x * x);

        let (orig_x, orig_y) = (50.0, 75.0);
        let (dist_x, dist_y) = distort(orig_x, orig_y);

        let result = newton_raphson_2d((dist_x, dist_y), (dist_x, dist_y), distort, 20, 1e-12);
        let (x, y) = result.unwrap();

        assert!((x - orig_x).abs() < 1e-10);
        assert!((y - orig_y).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_convergence() {
        let distort = |x: f64, y: f64| (x * 1.01 + 0.001 * y, y * 1.01 - 0.001 * x);

        let (orig_x, orig_y) = (500.0, 500.0);
        let (dist_x, dist_y) = distort(orig_x, orig_y);

        let result = newton_raphson_2d((dist_x, dist_y), (dist_x, dist_y), distort, 20, 1e-12);
        assert!(result.is_ok());

        let (x, y) = result.unwrap();
        let (check_x, check_y) = distort(x, y);
        assert!((check_x - dist_x).abs() < 1e-12);
        assert!((check_y - dist_y).abs() < 1e-12);
    }

    #[test]
    fn test_newton_raphson_singular_jacobian() {
        // Constant function has zero gradient (singular Jacobian)
        let constant_fn = |_x: f64, _y: f64| (5.0, 5.0);
        let result = newton_raphson_2d((0.0, 0.0), (0.0, 0.0), constant_fn, 20, 1e-12);
        assert!(result.is_err());
    }
}
