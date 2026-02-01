use crate::error::{WcsError, WcsResult};

use super::polynomial::newton_raphson_2d;

#[derive(Debug, Clone)]
pub struct TpvDistortion {
    pv1: [f64; 40],
    pv2: [f64; 40],
}

impl TpvDistortion {
    pub fn new() -> Self {
        Self {
            pv1: [0.0; 40],
            pv2: [0.0; 40],
        }
    }

    pub fn identity() -> Self {
        let mut tpv = Self::new();
        tpv.pv1[1] = 1.0; // x term
        tpv.pv2[2] = 1.0; // y term
        tpv
    }

    pub fn set_pv1(&mut self, index: usize, value: f64) {
        if index < 40 {
            self.pv1[index] = value;
        }
    }

    pub fn set_pv2(&mut self, index: usize, value: f64) {
        if index < 40 {
            self.pv2[index] = value;
        }
    }

    pub fn get_pv1(&self, index: usize) -> Option<f64> {
        self.pv1.get(index).copied()
    }

    pub fn get_pv2(&self, index: usize) -> Option<f64> {
        self.pv2.get(index).copied()
    }

    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let r = (x * x + y * y).sqrt();
        let xi = Self::eval_polynomial(&self.pv1, x, y, r);
        let eta = Self::eval_polynomial(&self.pv2, x, y, r);
        (xi, eta)
    }

    pub fn apply_inverse(&self, xi: f64, eta: f64) -> WcsResult<(f64, f64)> {
        let distort_fn = |x: f64, y: f64| self.apply(x, y);

        newton_raphson_2d((xi, eta), (xi, eta), distort_fn, 20, 1e-12).map_err(|msg| {
            WcsError::convergence_failure(format!("TPV inverse distortion: {}", msg))
        })
    }

    fn eval_polynomial(coeffs: &[f64; 40], x: f64, y: f64, r: f64) -> f64 {
        coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * Self::term(i, x, y, r))
            .sum()
    }

    #[inline]
    fn term(i: usize, x: f64, y: f64, r: f64) -> f64 {
        match i {
            0 => 1.0,
            1 => x,
            2 => y,
            3 => r,
            4 => x * x,
            5 => x * y,
            6 => y * y,
            7 => x * x * x,
            8 => x * x * y,
            9 => x * y * y,
            10 => y * y * y,
            11 => r * r * r,
            12 => x * x * x * x,
            13 => x * x * x * y,
            14 => x * x * y * y,
            15 => x * y * y * y,
            16 => y * y * y * y,
            17 => r * r * r * r,
            18 => x * x * x * x * x,
            19 => x * x * x * x * y,
            20 => x * x * x * y * y,
            21 => x * x * y * y * y,
            22 => x * y * y * y * y,
            23 => y * y * y * y * y,
            24 => r * r * r * r * r,
            25 => x * x * x * x * x * x,
            26 => x * x * x * x * x * y,
            27 => x * x * x * x * y * y,
            28 => x * x * x * y * y * y,
            29 => x * x * y * y * y * y,
            30 => x * y * y * y * y * y,
            31 => y * y * y * y * y * y,
            32 => r * r * r * r * r * r,
            33 => x * x * x * x * x * x * x,
            34 => x * x * x * x * x * x * y,
            35 => x * x * x * x * x * y * y,
            36 => x * x * x * x * y * y * y,
            37 => x * x * x * y * y * y * y,
            38 => x * x * y * y * y * y * y,
            _ => 0.0,
        }
    }
}

impl Default for TpvDistortion {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_distortion() {
        let tpv = TpvDistortion::identity();
        let (xi, eta) = tpv.apply(1.5, 2.5);
        assert_eq!(xi, 1.5);
        assert_eq!(eta, 2.5);
    }

    #[test]
    fn test_zero_distortion_returns_zero() {
        let tpv = TpvDistortion::new();
        let (xi, eta) = tpv.apply(1.5, 2.5);
        assert_eq!(xi, 0.0);
        assert_eq!(eta, 0.0);
    }

    #[test]
    fn test_linear_scale() {
        let mut tpv = TpvDistortion::new();
        tpv.set_pv1(1, 2.0); // 2x
        tpv.set_pv2(2, 2.0); // 2y

        let (xi, eta) = tpv.apply(3.0, 4.0);
        assert_eq!(xi, 6.0);
        assert_eq!(eta, 8.0);
    }

    #[test]
    fn test_radial_term() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(3, 0.01); // r term
        tpv.set_pv2(3, 0.01);

        let (x, y): (f64, f64) = (3.0, 4.0);
        let r = (x * x + y * y).sqrt();

        let (xi, eta) = tpv.apply(x, y);
        assert_eq!(xi, x + 0.01 * r);
        assert_eq!(eta, y + 0.01 * r);
    }

    #[test]
    fn test_radial_roundtrip() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(3, 0.001);
        tpv.set_pv2(3, 0.001);

        let (x_orig, y_orig) = (0.5, 0.7);
        let (xi, eta) = tpv.apply(x_orig, y_orig);
        let (x_back, y_back) = tpv.apply_inverse(xi, eta).unwrap();

        assert!((x_back - x_orig).abs() < 1e-12);
        assert!((y_back - y_orig).abs() < 1e-12);
    }

    #[test]
    fn test_quadratic_x_squared() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(4, 0.1); // x^2 term

        let (x, y) = (2.0, 1.0);
        let (xi, eta) = tpv.apply(x, y);

        assert_eq!(xi, x + 0.1 * x * x);
        assert_eq!(eta, y);
    }

    #[test]
    fn test_quadratic_y_squared() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv2(6, 0.1); // y^2 term

        let (x, y): (f64, f64) = (1.0, 3.0);
        let (xi, eta) = tpv.apply(x, y);
        let expected_eta = y + 0.1 * y * y;

        assert_eq!(xi, x);
        assert!((eta - expected_eta).abs() < 1e-14);
    }

    #[test]
    fn test_quadratic_roundtrip() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(4, 0.01);
        tpv.set_pv2(6, 0.01);

        let (x_orig, y_orig) = (0.3, 0.4);
        let (xi, eta) = tpv.apply(x_orig, y_orig);
        let (x_back, y_back) = tpv.apply_inverse(xi, eta).unwrap();

        assert!((x_back - x_orig).abs() < 1e-12);
        assert!((y_back - y_orig).abs() < 1e-12);
    }

    #[test]
    fn test_cross_term_xy() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(5, 0.05); // xy term

        let (x, y) = (2.0, 3.0);
        let (xi, eta) = tpv.apply(x, y);

        assert_eq!(xi, x + 0.05 * x * y);
        assert_eq!(eta, y);
    }

    #[test]
    fn test_cross_term_roundtrip() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(5, 0.02);
        tpv.set_pv2(5, 0.02);

        let (x_orig, y_orig) = (0.5, 0.6);
        let (xi, eta) = tpv.apply(x_orig, y_orig);
        let (x_back, y_back) = tpv.apply_inverse(xi, eta).unwrap();

        assert!((x_back - x_orig).abs() < 1e-12);
        assert!((y_back - y_orig).abs() < 1e-12);
    }

    #[test]
    fn test_higher_order_terms() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(7, 0.001); // x^3
        tpv.set_pv1(11, 0.002); // r^3
        tpv.set_pv2(10, 0.001); // y^3
        tpv.set_pv2(11, 0.002); // r^3

        let (x_orig, y_orig) = (0.2, 0.3);
        let (xi, eta) = tpv.apply(x_orig, y_orig);
        let (x_back, y_back) = tpv.apply_inverse(xi, eta).unwrap();

        assert!((x_back - x_orig).abs() < 1e-12);
        assert!((y_back - y_orig).abs() < 1e-12);
    }

    #[test]
    fn test_mixed_distortion_roundtrip() {
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(3, 0.001); // r
        tpv.set_pv1(4, 0.002); // x^2
        tpv.set_pv1(5, 0.001); // xy
        tpv.set_pv1(11, 0.0005); // r^3
        tpv.set_pv2(3, 0.001);
        tpv.set_pv2(5, 0.001);
        tpv.set_pv2(6, 0.002); // y^2
        tpv.set_pv2(11, 0.0005);

        let test_points = [
            (0.1, 0.1),
            (0.5, 0.3),
            (-0.2, 0.4),
            (0.3, -0.5),
            (-0.4, -0.4),
        ];

        for (x_orig, y_orig) in test_points {
            let (xi, eta) = tpv.apply(x_orig, y_orig);
            let (x_back, y_back) = tpv.apply_inverse(xi, eta).unwrap();

            assert!(
                (x_back - x_orig).abs() < 1e-12,
                "x roundtrip failed for ({}, {}): expected {}, got {}",
                x_orig,
                y_orig,
                x_orig,
                x_back
            );
            assert!(
                (y_back - y_orig).abs() < 1e-12,
                "y roundtrip failed for ({}, {}): expected {}, got {}",
                x_orig,
                y_orig,
                y_orig,
                y_back
            );
        }
    }

    #[test]
    fn test_getters() {
        let mut tpv = TpvDistortion::new();
        tpv.set_pv1(5, 1.23);
        tpv.set_pv2(10, 4.56);

        assert_eq!(tpv.get_pv1(5), Some(1.23));
        assert_eq!(tpv.get_pv2(10), Some(4.56));
        assert_eq!(tpv.get_pv1(40), None);
        assert_eq!(tpv.get_pv2(100), None);
    }

    #[test]
    fn test_out_of_bounds_set_ignored() {
        let mut tpv = TpvDistortion::new();
        tpv.set_pv1(50, 1.0);
        tpv.set_pv2(100, 2.0);

        assert_eq!(tpv.get_pv1(39), Some(0.0));
    }

    #[test]
    fn test_constant_term() {
        let mut tpv = TpvDistortion::new();
        tpv.set_pv1(0, 1.0);
        tpv.set_pv2(0, 2.0);

        let (xi, eta) = tpv.apply(0.0, 0.0);
        assert_eq!(xi, 1.0);
        assert_eq!(eta, 2.0);
    }

    #[test]
    fn test_all_term_indices() {
        let x: f64 = 0.1;
        let y: f64 = 0.2;
        let r = (x * x + y * y).sqrt();

        let expected_terms: [(usize, f64); 40] = [
            (0, 1.0),
            (1, x),
            (2, y),
            (3, r),
            (4, x * x),
            (5, x * y),
            (6, y * y),
            (7, x.powi(3)),
            (8, x * x * y),
            (9, x * y * y),
            (10, y.powi(3)),
            (11, r.powi(3)),
            (12, x.powi(4)),
            (13, x.powi(3) * y),
            (14, x.powi(2) * y.powi(2)),
            (15, x * y.powi(3)),
            (16, y.powi(4)),
            (17, r.powi(4)),
            (18, x.powi(5)),
            (19, x.powi(4) * y),
            (20, x.powi(3) * y.powi(2)),
            (21, x.powi(2) * y.powi(3)),
            (22, x * y.powi(4)),
            (23, y.powi(5)),
            (24, r.powi(5)),
            (25, x.powi(6)),
            (26, x.powi(5) * y),
            (27, x.powi(4) * y.powi(2)),
            (28, x.powi(3) * y.powi(3)),
            (29, x.powi(2) * y.powi(4)),
            (30, x * y.powi(5)),
            (31, y.powi(6)),
            (32, r.powi(6)),
            (33, x.powi(7)),
            (34, x.powi(6) * y),
            (35, x.powi(5) * y.powi(2)),
            (36, x.powi(4) * y.powi(3)),
            (37, x.powi(3) * y.powi(4)),
            (38, x.powi(2) * y.powi(5)),
            (39, 0.0),
        ];

        for (i, expected) in expected_terms {
            let computed = TpvDistortion::term(i, x, y, r);
            assert!(
                (computed - expected).abs() < 1e-15,
                "term {} mismatch: expected {}, got {}",
                i,
                expected,
                computed
            );
        }
    }
}
