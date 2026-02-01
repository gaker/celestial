use std::collections::HashMap;

use crate::error::{WcsError, WcsResult};

use super::polynomial::{newton_raphson_2d, power_term};

#[derive(Debug, Clone)]
pub struct SipDistortion {
    crpix: [f64; 2],
    a_order: u32,
    b_order: u32,
    a_coeffs: HashMap<(u32, u32), f64>,
    b_coeffs: HashMap<(u32, u32), f64>,
    ap_order: Option<u32>,
    bp_order: Option<u32>,
    ap_coeffs: HashMap<(u32, u32), f64>,
    bp_coeffs: HashMap<(u32, u32), f64>,
}

impl SipDistortion {
    pub fn new(crpix: [f64; 2], a_order: u32, b_order: u32) -> Self {
        Self {
            crpix,
            a_order,
            b_order,
            a_coeffs: HashMap::new(),
            b_coeffs: HashMap::new(),
            ap_order: None,
            bp_order: None,
            ap_coeffs: HashMap::new(),
            bp_coeffs: HashMap::new(),
        }
    }

    pub fn set_a(&mut self, p: u32, q: u32, value: f64) {
        if p + q <= self.a_order && value != 0.0 {
            self.a_coeffs.insert((p, q), value);
        }
    }

    pub fn set_b(&mut self, p: u32, q: u32, value: f64) {
        if p + q <= self.b_order && value != 0.0 {
            self.b_coeffs.insert((p, q), value);
        }
    }

    pub fn set_ap(&mut self, p: u32, q: u32, value: f64) {
        if let Some(order) = self.ap_order {
            if p + q <= order && value != 0.0 {
                self.ap_coeffs.insert((p, q), value);
            }
        }
    }

    pub fn set_bp(&mut self, p: u32, q: u32, value: f64) {
        if let Some(order) = self.bp_order {
            if p + q <= order && value != 0.0 {
                self.bp_coeffs.insert((p, q), value);
            }
        }
    }

    pub fn set_inverse_order(&mut self, ap_order: u32, bp_order: u32) {
        self.ap_order = Some(ap_order);
        self.bp_order = Some(bp_order);
    }

    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let u = x - self.crpix[0];
        let v = y - self.crpix[1];

        let f = Self::eval_poly(&self.a_coeffs, u, v);
        let g = Self::eval_poly(&self.b_coeffs, u, v);

        (x + f, y + g)
    }

    pub fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        if self.has_inverse_coeffs() {
            Ok(self.apply_inverse_analytic(x, y))
        } else {
            self.apply_inverse_iterative(x, y)
        }
    }

    fn has_inverse_coeffs(&self) -> bool {
        self.ap_order.is_some() && self.bp_order.is_some()
    }

    fn apply_inverse_analytic(&self, x: f64, y: f64) -> (f64, f64) {
        let u_prime = x - self.crpix[0];
        let v_prime = y - self.crpix[1];

        let f_prime = Self::eval_poly(&self.ap_coeffs, u_prime, v_prime);
        let g_prime = Self::eval_poly(&self.bp_coeffs, u_prime, v_prime);

        (x + f_prime, y + g_prime)
    }

    fn apply_inverse_iterative(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        let distort_fn = |px: f64, py: f64| self.apply(px, py);

        newton_raphson_2d((x, y), (x, y), distort_fn, 20, 1e-12).map_err(|msg| {
            WcsError::convergence_failure(format!("SIP inverse distortion: {}", msg))
        })
    }

    fn eval_poly(coeffs: &HashMap<(u32, u32), f64>, u: f64, v: f64) -> f64 {
        coeffs
            .iter()
            .map(|(&(p, q), &coeff)| coeff * power_term(u, v, p, q))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_distortion() {
        let sip = SipDistortion::new([512.0, 512.0], 2, 2);
        let (x, y) = sip.apply(100.0, 200.0);
        assert_eq!(x, 100.0);
        assert_eq!(y, 200.0);
    }

    #[test]
    fn test_simple_quadratic() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6);

        let (x, y) = (612.0, 612.0);
        let u = x - 512.0;

        let (x_out, y_out) = sip.apply(x, y);

        let expected_x = x + 1e-6 * u * u;
        assert_eq!(x_out, expected_x);
        assert_eq!(y_out, y);
    }

    #[test]
    fn test_roundtrip_with_inverse_coefficients() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(1, 0, 1e-5);
        sip.set_b(0, 1, 1e-5);

        sip.set_inverse_order(2, 2);
        sip.set_ap(1, 0, -1e-5);
        sip.set_bp(0, 1, -1e-5);

        let (x_orig, y_orig) = (562.0, 562.0);
        let (x_dist, y_dist) = sip.apply(x_orig, y_orig);
        let (x_back, y_back) = sip.apply_inverse(x_dist, y_dist).unwrap();

        assert!((x_back - x_orig).abs() < 1e-8);
        assert!((y_back - y_orig).abs() < 1e-8);
    }

    #[test]
    fn test_roundtrip_newton_raphson() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6);
        sip.set_b(0, 2, 1e-6);

        let (x_orig, y_orig) = (612.0, 612.0);
        let (x_dist, y_dist) = sip.apply(x_orig, y_orig);
        let (x_back, y_back) = sip.apply_inverse(x_dist, y_dist).unwrap();

        assert!((x_back - x_orig).abs() < 1e-10);
        assert!((y_back - y_orig).abs() < 1e-10);
    }

    #[test]
    fn test_cross_terms() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(1, 1, 1e-6);
        sip.set_b(1, 1, 2e-6);

        let (x, y) = (612.0, 712.0);
        let u = x - 512.0;
        let v = y - 512.0;

        let (x_out, y_out) = sip.apply(x, y);

        let expected_x = x + 1e-6 * u * v;
        let expected_y = y + 2e-6 * u * v;

        assert_eq!(x_out, expected_x);
        assert_eq!(y_out, expected_y);
    }

    #[test]
    fn test_large_pixel_offsets() {
        let mut sip = SipDistortion::new([512.0, 512.0], 3, 3);
        sip.set_a(2, 0, 1e-7);
        sip.set_a(0, 2, 1e-7);
        sip.set_a(3, 0, 1e-10);
        sip.set_b(2, 0, 1e-7);
        sip.set_b(0, 2, 1e-7);
        sip.set_b(0, 3, 1e-10);

        let test_points = [
            (512.0 + 1000.0, 512.0 + 1000.0),
            (512.0 - 1000.0, 512.0 - 1000.0),
            (512.0 + 1000.0, 512.0 - 1000.0),
            (512.0 - 1000.0, 512.0 + 1000.0),
        ];

        for (x_orig, y_orig) in test_points {
            let (x_dist, y_dist) = sip.apply(x_orig, y_orig);
            let (x_back, y_back) = sip.apply_inverse(x_dist, y_dist).unwrap();

            assert!(
                (x_back - x_orig).abs() < 1e-10,
                "x roundtrip failed for ({}, {})",
                x_orig,
                y_orig
            );
            assert!(
                (y_back - y_orig).abs() < 1e-10,
                "y roundtrip failed for ({}, {})",
                x_orig,
                y_orig
            );
        }
    }

    #[test]
    fn test_zero_coefficient_omitted() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 0.0);
        sip.set_a(1, 1, 1e-6);

        assert!(!sip.a_coeffs.contains_key(&(2, 0)));
        assert!(sip.a_coeffs.contains_key(&(1, 1)));
    }

    #[test]
    fn test_order_constraint() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(3, 0, 1e-6);
        sip.set_a(2, 1, 1e-6);

        assert!(!sip.a_coeffs.contains_key(&(3, 0)));
        assert!(!sip.a_coeffs.contains_key(&(2, 1)));
    }
}
