use std::ops::ControlFlow;

use crate::common::is_singular;
use crate::error::{WcsError, WcsResult};

// Newton's method on an f that rises through a bracket holding its one root.
// Every evaluation narrows the bracket and a step that would leave it bisects
// instead, so the iteration can neither stall against a bound nor wander to
// another root. It stops once a step no longer moves x.
pub(crate) fn bracketed_newton(
    mut bracket: (f64, f64),
    seed: f64,
    f_and_slope: impl Fn(f64) -> WcsResult<(f64, f64)>,
    context: &str,
) -> WcsResult<f64> {
    const MAX_ITER: usize = 100;
    let mut x = seed;
    for _ in 0..MAX_ITER {
        let (f, slope) = f_and_slope(x)?;
        if f == 0.0 {
            return Ok(x);
        }
        if f.is_nan() {
            break;
        }
        match next_iterate(x, f, slope, &mut bracket) {
            ControlFlow::Continue(next) => x = next,
            ControlFlow::Break(root) => return Ok(root),
        }
    }
    Err(WcsError::convergence_failure(format!(
        "{context}: Newton-Raphson did not converge"
    )))
}

fn next_iterate(x: f64, f: f64, slope: f64, bracket: &mut (f64, f64)) -> ControlFlow<f64, f64> {
    if f > 0.0 {
        bracket.1 = x;
    } else {
        bracket.0 = x;
    }
    let (lo, hi) = *bracket;
    let step = x - f / slope;
    if step == x {
        return ControlFlow::Break(x);
    }
    let next = if step > lo && step < hi {
        step
    } else {
        (lo + hi) / 2.0
    };
    if next == lo || next == hi {
        ControlFlow::Break(next)
    } else {
        ControlFlow::Continue(next)
    }
}

// Newton's method on f(x, y) = target with a forward-difference Jacobian.
// While it converges each step is smaller than the last; once rounding
// dominates the residual the steps stop shrinking, and the point before that
// step is as close as the arithmetic gets. A fixed residual tolerance would
// stop short of that floor wherever the residual's rounding error is small.
pub(crate) fn newton_raphson_2d(
    target: (f64, f64),
    initial_guess: (f64, f64),
    f: impl Fn(f64, f64) -> (f64, f64),
    max_iter: usize,
    context: &str,
) -> WcsResult<(f64, f64)> {
    let (mut x, mut y) = initial_guess;
    let mut last_step = f64::INFINITY;
    for _ in 0..max_iter {
        let (dx, dy) = newton_step(&f, target, (x, y), context)?;
        let step = libm::hypot(dx, dy);
        if !step.is_finite() {
            break;
        }
        if step == 0.0 || (step >= last_step && is_rounding(step, target, (x, y))) {
            return Ok((x, y));
        }
        (x, y, last_step) = (x - dx, y - dy, step);
    }
    Err(WcsError::convergence_failure(format!(
        "{context}: Newton-Raphson did not converge"
    )))
}

// A step that stops shrinking far from the root means Newton has not settled
// yet, so it must go on. Rounding moves the point by about ε of the values in
// the residual; √ε, midway in exponent between that and a step the size of
// the values themselves, separates the two.
fn is_rounding(step: f64, (tx, ty): (f64, f64), (x, y): (f64, f64)) -> bool {
    const SQRT_EPSILON: f64 = 1.0 / (1u64 << 26) as f64;
    step <= SQRT_EPSILON * (libm::hypot(tx, ty) + libm::hypot(x, y))
}

fn newton_step(
    f: &impl Fn(f64, f64) -> (f64, f64),
    (tx, ty): (f64, f64),
    (x, y): (f64, f64),
    context: &str,
) -> WcsResult<(f64, f64)> {
    let value = f(x, y);
    let jacobian = jacobian_2d(f, x, y, value);
    solve_2x2(jacobian, (value.0 - tx, value.1 - ty), context)
}

fn jacobian_2d(
    f: &impl Fn(f64, f64) -> (f64, f64),
    x: f64,
    y: f64,
    (fx, fy): (f64, f64),
) -> [[f64; 2]; 2] {
    const H: f64 = 1e-8;
    let (fx_px, fy_px) = f(x + H, y);
    let (fx_py, fy_py) = f(x, y + H);

    let j11 = (fx_px - fx) / H;
    let j12 = (fx_py - fx) / H;
    let j21 = (fy_px - fy) / H;
    let j22 = (fy_py - fy) / H;

    [[j11, j12], [j21, j22]]
}
fn solve_2x2(m: [[f64; 2]; 2], (b1, b2): (f64, f64), context: &str) -> WcsResult<(f64, f64)> {
    let [[j11, j12], [j21, j22]] = m;
    let det = j11 * j22 - j12 * j21;
    if is_singular(m, det) {
        return Err(WcsError::convergence_failure(format!(
            "{context}: singular Jacobian matrix"
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

    #[test]
    fn test_newton_raphson_2d_recovers_distortion() {
        // Two distortion shapes - separable quadratic and mixed cross-term -
        // exercise both diagonal and off-diagonal Jacobian terms. The exact
        // root of each distorted point is within 6e-17 of the original, whose
        // double is the nearest.
        type Case = (fn(f64, f64) -> (f64, f64), f64, f64);
        let cases: &[Case] = &[
            (|x, y| (x + 0.001 * x * x, y + 0.001 * y * y), 100.0, 200.0),
            (|x, y| (x + 0.0001 * x * y, y + 0.0002 * x * x), 50.0, 75.0),
        ];

        for &(distort, orig_x, orig_y) in cases {
            let (dist_x, dist_y) = distort(orig_x, orig_y);
            let (x, y) =
                newton_raphson_2d((dist_x, dist_y), (dist_x, dist_y), distort, 20, "test").unwrap();
            assert_eq!((x, y), (orig_x, orig_y));
        }
    }

    #[test]
    fn test_newton_raphson_2d_singular_jacobian_errors() {
        // A constant function has a zero Jacobian and must surface a
        // ConvergenceFailure rather than divide by zero.
        let constant = |_x: f64, _y: f64| (5.0, 5.0);
        let err = newton_raphson_2d((0.0, 0.0), (0.0, 0.0), constant, 20, "test").unwrap_err();
        assert!(
            matches!(err, WcsError::ConvergenceFailure { .. }),
            "got: {:?}",
            err
        );
    }
}
