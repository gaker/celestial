use std::collections::HashMap;

use crate::common::newton_raphson_2d;
use crate::error::WcsResult;

use super::polynomial::power_term;

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

        newton_raphson_2d((x, y), (x, y), distort_fn, 20, 1e-12, "SIP inverse distortion")
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
    fn test_apply_known_polynomial_terms() {
        // Three forward-only assertions on terms predicted by the SIP polynomial:
        //   (a) identity (no coeffs set),
        //   (b) pure quadratic in u via A_2_0,
        //   (c) cross-term u*v via A_1_1 and B_1_1.

        let id = SipDistortion::new([512.0, 512.0], 2, 2);
        assert_eq!(id.apply(100.0, 200.0), (100.0, 200.0));

        let mut quad = SipDistortion::new([512.0, 512.0], 2, 2);
        quad.set_a(2, 0, 1e-6);
        let (x, y) = (612.0, 612.0);
        let u = x - 512.0;
        let (xo, yo) = quad.apply(x, y);
        assert_eq!(xo, x + 1e-6 * u * u);
        assert_eq!(yo, y);

        let mut cross = SipDistortion::new([512.0, 512.0], 2, 2);
        cross.set_a(1, 1, 1e-6);
        cross.set_b(1, 1, 2e-6);
        let (x, y) = (612.0, 712.0);
        let u = x - 512.0;
        let v = y - 512.0;
        let (xo, yo) = cross.apply(x, y);
        assert_eq!(xo, x + 1e-6 * u * v);
        assert_eq!(yo, y + 2e-6 * u * v);
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
    fn test_set_a_filters_zero_and_over_order_terms() {
        // SIP only stores non-zero coefficients up to A_ORDER (p + q <= A_ORDER).
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);

        sip.set_a(2, 0, 0.0);   // zero -> not stored
        sip.set_a(1, 1, 1e-6);  // valid (1 + 1 = 2) -> stored
        sip.set_a(3, 0, 1e-6);  // 3 + 0 > 2 -> rejected
        sip.set_a(2, 1, 1e-6);  // 2 + 1 > 2 -> rejected

        assert!(!sip.a_coeffs.contains_key(&(2, 0)));
        assert!(sip.a_coeffs.contains_key(&(1, 1)));
        assert!(!sip.a_coeffs.contains_key(&(3, 0)));
        assert!(!sip.a_coeffs.contains_key(&(2, 1)));
    }
}
