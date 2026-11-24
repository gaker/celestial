//! Least-squares affine map from tangent-plane positions to pixels, with the leverage
//! terms that give its prediction variance.
//!
//! With `n` fitted points `q_i` (mean `q̄`) and scatter `S = Σ (q_i - q̄)(q_i - q̄)ᵀ`, a
//! prediction at `q` has per-axis variance `σ² h(q)` where
//! `h(q) = 1/n + (q - q̄)ᵀ S⁻¹ (q - q̄)`, for independent pixel errors of per-axis
//! variance `σ²` on the fitted points. The distance between a fresh detection and that
//! prediction therefore has per-axis variance `σ² (1 + h(q))`.

pub(super) type Point = [f64; 2];
type Mat2 = [[f64; 2]; 2];

pub(super) struct AffineFit {
    n: f64,
    mean_q: Point,
    mean_p: Point,
    scatter: Mat2,
    scatter_inv: Mat2,
    // pixel = mean_p + linear * (q - mean_q)
    linear: Mat2,
}

impl AffineFit {
    // `None` for fewer than 3 points or points collinear to working precision.
    pub(super) fn new(points: &[(Point, Point)]) -> Option<Self> {
        if points.len() < 3 {
            return None;
        }
        let n = points.len() as f64;
        let mean_q = mean(points.iter().map(|pt| pt.0), n);
        let mean_p = mean(points.iter().map(|pt| pt.1), n);
        let (scatter, cross) = scatter_and_cross(points, mean_q, mean_p);
        let scatter_inv = invert_scatter(scatter)?;
        let coef = mul(scatter_inv, cross);
        let linear = [[coef[0][0], coef[1][0]], [coef[0][1], coef[1][1]]];
        Some(Self {
            n,
            mean_q,
            mean_p,
            scatter,
            scatter_inv,
            linear,
        })
    }

    pub(super) fn predict(&self, q: Point) -> Point {
        let d = sub(q, self.mean_q);
        let l = self.linear;
        [
            self.mean_p[0] + l[0][0] * d[0] + l[0][1] * d[1],
            self.mean_p[1] + l[1][0] * d[0] + l[1][1] * d[1],
        ]
    }

    pub(super) fn leverage(&self, q: Point) -> f64 {
        self.cross_leverage(q, q)
    }

    pub(super) fn determinant(&self) -> f64 {
        let l = self.linear;
        l[0][0] * l[1][1] - l[0][1] * l[1][0]
    }

    // `x_aᵀ (XᵀX)⁻¹ x_b` for design rows `x = (1, q)`.
    fn cross_leverage(&self, a: Point, b: Point) -> f64 {
        let da = sub(a, self.mean_q);
        let db = sub(b, self.mean_q);
        let s = self.scatter_inv;
        let quad = da[0] * (s[0][0] * db[0] + s[0][1] * db[1])
            + da[1] * (s[1][0] * db[0] + s[1][1] * db[1]);
        1.0 / self.n + quad
    }

    // The fit without the fitted point at `q_t`, whose pixel differs from the full
    // model's prediction by `residual`. `None` when the remaining points do not
    // determine a model.
    pub(super) fn left_out(&self, q_t: Point, residual: Point) -> Option<LeftOut> {
        if self.n < 4.0 {
            return None;
        }
        let d = sub(q_t, self.mean_q);
        let w = self.n / (self.n - 1.0);
        let mut rest = self.scatter;
        for k in 0..2 {
            for c in 0..2 {
                rest[k][c] -= w * d[k] * d[c];
            }
        }
        invert_scatter(rest)?;
        Some(LeftOut {
            q_t,
            residual,
            inv_one_minus_h: 1.0 / (1.0 - self.leverage(q_t)),
        })
    }
}

// Sherman-Morrison form of the fit without one point `t`: the prediction at `q` moves by
// `-c e_t / (1 - h_t)` and its leverage grows by `c² / (1 - h_t)`, where
// `c = x_qᵀ (XᵀX)⁻¹ x_t` and `e_t` is the residual of `t` under the full fit.
pub(super) struct LeftOut {
    q_t: Point,
    residual: Point,
    inv_one_minus_h: f64,
}

impl LeftOut {
    // The move of the prediction at `q` and the growth of its leverage.
    pub(super) fn correction(&self, full: &AffineFit, q: Point) -> (Point, f64) {
        let c = full.cross_leverage(q, self.q_t);
        let shift = c * self.inv_one_minus_h;
        let offset = [-(shift * self.residual[0]), -(shift * self.residual[1])];
        (offset, c * shift)
    }
}

// The scatter of the q points about `mean_q`, sum of d d^T, and their cross scatter with
// the p points about `mean_p`, sum of d e^T.
fn scatter_and_cross(points: &[(Point, Point)], mean_q: Point, mean_p: Point) -> (Mat2, Mat2) {
    let mut scatter = [[0.0; 2]; 2];
    let mut cross = [[0.0; 2]; 2];
    for (q, p) in points {
        let d = sub(*q, mean_q);
        let e = sub(*p, mean_p);
        for k in 0..2 {
            for c in 0..2 {
                scatter[k][c] += d[k] * d[c];
                cross[k][c] += d[k] * e[c];
            }
        }
    }
    (scatter, cross)
}

fn mean(points: impl Iterator<Item = Point>, n: f64) -> Point {
    let sum = points.fold([0.0, 0.0], |acc, p| [acc[0] + p[0], acc[1] + p[1]]);
    [sum[0] / n, sum[1] / n]
}

fn sub(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1]]
}

fn mul(a: Mat2, b: Mat2) -> Mat2 {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

// Inverse of a symmetric positive semidefinite scatter matrix, or `None` when its
// determinant is at the rounding level of its entries (the points are collinear).
fn invert_scatter(s: Mat2) -> Option<Mat2> {
    let det = s[0][0] * s[1][1] - s[0][1] * s[1][0];
    let trace = s[0][0] + s[1][1];
    let full_rank = det.is_finite() && det > f64::EPSILON * trace * trace;
    full_rank.then(|| {
        [
            [s[1][1] / det, -s[0][1] / det],
            [-s[1][0] / det, s[0][0] / det],
        ]
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // pixel = 10 q + (5, 7) on the four points (±1, 0), (0, ±1), whose mean is the
    // origin and whose scatter is diag(2, 2).
    fn cross_points() -> Vec<(Point, Point)> {
        [[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]]
            .into_iter()
            .map(|q| (q, [10.0 * q[0] + 5.0, 10.0 * q[1] + 7.0]))
            .collect()
    }

    #[test]
    fn fit_recovers_the_map_and_its_leverage() {
        let fit = AffineFit::new(&cross_points()).unwrap();
        assert_eq!(fit.predict([0.5, 0.25]), [10.0, 9.5]);
        assert_eq!(fit.determinant(), 100.0);
        // 1/n at the mean; 1/n + 1/2 one unit along an axis.
        assert_eq!(fit.leverage([0.0, 0.0]), 0.25);
        assert_eq!(fit.leverage([1.0, 0.0]), 0.75);
    }

    #[test]
    fn mirrored_map_has_negative_determinant() {
        let points: Vec<(Point, Point)> = cross_points()
            .into_iter()
            .map(|(q, p)| (q, [p[0], -p[1]]))
            .collect();
        assert_eq!(AffineFit::new(&points).unwrap().determinant(), -100.0);
    }

    #[test]
    fn fit_needs_three_points_off_a_line() {
        let two = &cross_points()[..2];
        assert!(AffineFit::new(two).is_none());
        let line: Vec<(Point, Point)> = (0..4)
            .map(|i| ([i as f64, 2.0 * i as f64], [i as f64, i as f64]))
            .collect();
        assert!(AffineFit::new(&line).is_none());
    }

    #[test]
    fn left_out_matches_a_fit_to_the_other_points() {
        // Move the first point's pixel by (4, -2): the full fit absorbs part of it,
        // the fit to the other three points does not.
        let mut points = cross_points();
        points[0].1 = [19.0, 5.0];
        let full = AffineFit::new(&points).unwrap();
        let residual = sub(points[0].1, full.predict(points[0].0));
        let left_out = full.left_out(points[0].0, residual).unwrap();
        let q = [1.0, 0.0];
        let (offset, extra) = left_out.correction(&full, q);
        let p = full.predict(q);
        assert_eq!([p[0] + offset[0], p[1] + offset[1]], [15.0, 7.0]);
        // Against the three remaining points, (1, 0) has leverage 1/3 + (4/3)^2 / (2/3).
        assert_eq!(1.0 + full.leverage(q) + extra, 4.0);
    }

    #[test]
    fn left_out_needs_three_remaining_points_off_a_line() {
        let three = &cross_points()[1..];
        let fit = AffineFit::new(three).unwrap();
        assert!(fit.left_out(three[0].0, [0.0, 0.0]).is_none());
        // Without (1, 0), the points (-1, 0), (2, 0) and (3, 0) are collinear.
        let points: Vec<(Point, Point)> = [[1.0, 1.0], [-1.0, 0.0], [2.0, 0.0], [3.0, 0.0]]
            .into_iter()
            .map(|q| (q, q))
            .collect();
        let fit = AffineFit::new(&points).unwrap();
        assert!(fit.left_out(points[0].0, [0.0, 0.0]).is_none());
    }
}
