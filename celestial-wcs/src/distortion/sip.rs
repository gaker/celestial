use crate::common::require_finite;
use crate::error::{WcsError, WcsResult};
use crate::newton::newton_raphson_2d;

use super::polynomial::power_term;

// Sorted by (p, q), so that every instance sums its terms in the same order.
type Terms = Vec<((u32, u32), f64)>;

// A polynomial's name, its order if it has one, and its terms.
pub(crate) type Polynomial<'a> = (&'static str, Option<u32>, &'a [((u32, u32), f64)]);

#[derive(Debug, Clone)]
pub struct SipDistortion {
    crpix: [f64; 2],
    a_order: u32,
    b_order: u32,
    a_coeffs: Terms,
    b_coeffs: Terms,
    ap_order: Option<u32>,
    bp_order: Option<u32>,
    ap_coeffs: Terms,
    bp_coeffs: Terms,
}

impl SipDistortion {
    pub fn new(crpix: [f64; 2], a_order: u32, b_order: u32) -> Self {
        Self {
            crpix,
            a_order,
            b_order,
            a_coeffs: Terms::new(),
            b_coeffs: Terms::new(),
            ap_order: None,
            bp_order: None,
            ap_coeffs: Terms::new(),
            bp_coeffs: Terms::new(),
        }
    }

    pub fn set_a(&mut self, p: u32, q: u32, value: f64) -> WcsResult<()> {
        check_term("A", Some(self.a_order), (p, q), value)?;
        set_term(&mut self.a_coeffs, (p, q), value);
        Ok(())
    }

    pub fn set_b(&mut self, p: u32, q: u32, value: f64) -> WcsResult<()> {
        check_term("B", Some(self.b_order), (p, q), value)?;
        set_term(&mut self.b_coeffs, (p, q), value);
        Ok(())
    }

    pub fn set_ap(&mut self, p: u32, q: u32, value: f64) -> WcsResult<()> {
        check_term("AP", self.ap_order, (p, q), value)?;
        set_term(&mut self.ap_coeffs, (p, q), value);
        Ok(())
    }

    pub fn set_bp(&mut self, p: u32, q: u32, value: f64) -> WcsResult<()> {
        check_term("BP", self.bp_order, (p, q), value)?;
        set_term(&mut self.bp_coeffs, (p, q), value);
        Ok(())
    }

    pub fn set_inverse_order(&mut self, ap_order: u32, bp_order: u32) -> WcsResult<()> {
        check_order_holds("AP", ap_order, &self.ap_coeffs)?;
        check_order_holds("BP", bp_order, &self.bp_coeffs)?;
        self.ap_order = Some(ap_order);
        self.bp_order = Some(bp_order);
        Ok(())
    }

    pub(crate) fn crpix(&self) -> [f64; 2] {
        self.crpix
    }

    // AP and BP have no order until set_inverse_order gives them one.
    pub(crate) fn polynomials(&self) -> [Polynomial<'_>; 4] {
        [
            ("A", Some(self.a_order), self.a_coeffs.as_slice()),
            ("B", Some(self.b_order), self.b_coeffs.as_slice()),
            ("AP", self.ap_order, self.ap_coeffs.as_slice()),
            ("BP", self.bp_order, self.bp_coeffs.as_slice()),
        ]
    }

    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let u = x - self.crpix[0];
        let v = y - self.crpix[1];

        let f = Self::eval_poly(&self.a_coeffs, u, v);
        let g = Self::eval_poly(&self.b_coeffs, u, v);

        (x + f, y + g)
    }

    // AP and BP are a fit to the inverse of A and B, so they miss the root by the
    // fit error; they only give Newton a closer start.
    pub fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        require_finite(&[x, y], "SIP inverse input")?;
        let seed = if self.has_inverse_coeffs() {
            self.ap_bp_estimate(x, y)
        } else {
            (x, y)
        };
        let distort_fn = |px: f64, py: f64| self.apply(px, py);
        newton_raphson_2d((x, y), seed, distort_fn, 20, "SIP inverse distortion")
    }

    fn has_inverse_coeffs(&self) -> bool {
        self.ap_order.is_some() && self.bp_order.is_some()
    }

    fn ap_bp_estimate(&self, x: f64, y: f64) -> (f64, f64) {
        let u_prime = x - self.crpix[0];
        let v_prime = y - self.crpix[1];

        let f_prime = Self::eval_poly(&self.ap_coeffs, u_prime, v_prime);
        let g_prime = Self::eval_poly(&self.bp_coeffs, u_prime, v_prime);

        (x + f_prime, y + g_prime)
    }

    fn eval_poly(coeffs: &Terms, u: f64, v: f64) -> f64 {
        coeffs
            .iter()
            .map(|&((p, q), coeff)| coeff * power_term(u, v, p, q))
            .sum()
    }
}

fn check_term(name: &str, order: Option<u32>, (p, q): (u32, u32), value: f64) -> WcsResult<()> {
    let order =
        order.ok_or_else(|| WcsError::invalid_parameter(format!("{name}_ORDER is not set")))?;
    if p.saturating_add(q) > order {
        let message = format!("{name}_{p}_{q} is beyond {name}_ORDER = {order}");
        return Err(WcsError::invalid_parameter(message));
    }
    if !value.is_finite() {
        return Err(WcsError::invalid_parameter(format!(
            "{name}_{p}_{q} is {value}"
        )));
    }
    Ok(())
}

fn check_order_holds(name: &str, order: u32, terms: &Terms) -> WcsResult<()> {
    match terms.iter().find(|&&((p, q), _)| p + q > order) {
        Some(((p, q), _)) => Err(WcsError::invalid_parameter(format!(
            "{name}_ORDER = {order} is below the stored {name}_{p}_{q}"
        ))),
        None => Ok(()),
    }
}

// A zero coefficient is the absence of the term.
fn set_term(terms: &mut Terms, (p, q): (u32, u32), value: f64) {
    match terms.binary_search_by_key(&(p, q), |&(pq, _)| pq) {
        Ok(i) if value == 0.0 => {
            terms.remove(i);
        }
        Ok(i) => terms[i].1 = value,
        Err(i) if value != 0.0 => terms.insert(i, ((p, q), value)),
        Err(_) => {}
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
        quad.set_a(2, 0, 1e-6).unwrap();
        let (x, y) = (612.0, 612.0);
        let u = x - 512.0;
        let (xo, yo) = quad.apply(x, y);
        assert_eq!(xo, x + 1e-6 * u * u);
        assert_eq!(yo, y);

        let mut cross = SipDistortion::new([512.0, 512.0], 2, 2);
        cross.set_a(1, 1, 1e-6).unwrap();
        cross.set_b(1, 1, 2e-6).unwrap();
        let (x, y) = (612.0, 712.0);
        let u = x - 512.0;
        let v = y - 512.0;
        let (xo, yo) = cross.apply(x, y);
        assert_eq!(xo, x + 1e-6 * u * v);
        assert_eq!(yo, y + 2e-6 * u * v);
    }

    // In both SIP round-trip tests below, each value either way is the nearest
    // double to the exact result.
    #[test]
    fn test_roundtrip_with_inverse_coefficients() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(1, 0, 1e-5).unwrap();
        sip.set_b(0, 1, 1e-5).unwrap();

        sip.set_inverse_order(2, 2).unwrap();
        sip.set_ap(1, 0, -1e-5).unwrap();
        sip.set_bp(0, 1, -1e-5).unwrap();

        assert_eq!(sip.apply(562.0, 562.0), (562.0005, 562.0005));
        assert_eq!(
            sip.apply_inverse(562.0005, 562.0005).unwrap(),
            (562.0, 562.0)
        );
    }

    #[test]
    fn test_roundtrip_newton_raphson() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6).unwrap();
        sip.set_b(0, 2, 1e-6).unwrap();

        assert_eq!(sip.apply(612.0, 612.0), (612.01, 612.01));
        assert_eq!(sip.apply_inverse(612.01, 612.01).unwrap(), (612.0, 612.0));
    }

    // Each expected point is the nearest double to the exact inverse; the
    // comments give the exact inverse to 21 digits.
    #[test]
    fn test_newton_inverse_solves_to_rounding() {
        let mut both = SipDistortion::new([256.0, 256.0], 2, 2);
        for (p, q, a, b) in [(2, 0, 1e-6, 3e-7), (1, 1, -2e-7, 1e-6), (0, 2, 5e-7, -4e-7)] {
            both.set_a(p, q, a).unwrap();
            both.set_b(p, q, b).unwrap();
        }
        // (355.988690045846128114, 330.991751856665096275)
        assert_eq!(
            both.apply_inverse(356.0, 331.0).unwrap(),
            (355.9886900458461, 330.9917518566651)
        );
        // (137.235966956446519223, 211.741295722281002073)
        assert_eq!(
            both.apply_inverse(137.25, 211.75).unwrap(),
            (137.23596695644653, 211.741295722281)
        );

        let mut x_only = SipDistortion::new([512.0, 512.0], 2, 2);
        x_only.set_a(2, 0, 1e-6).unwrap();
        // (137.109457080829769790, 211.75)
        assert_eq!(
            x_only.apply_inverse(137.25, 211.75).unwrap(),
            (137.10945708082977, 211.75)
        );
    }

    // AP and BP are a fit to the inverse, so they miss the root of A and B by
    // the fit error (1.4e-4 px in y at the far corner); the inverse is that root.
    #[test]
    fn test_inverse_with_ap_bp_is_the_root_of_a_b() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6).unwrap();
        sip.set_a(0, 2, -5e-7).unwrap();
        sip.set_b(1, 1, 8e-7).unwrap();
        sip.set_inverse_order(2, 2).unwrap();
        sip.set_ap(2, 0, -1e-6).unwrap();
        sip.set_ap(0, 2, 5e-7).unwrap();
        sip.set_bp(1, 1, -8e-7).unwrap();
        // (137.154592998134190905, 211.659935124914345393)
        assert_eq!(
            sip.apply_inverse(137.25, 211.75).unwrap(),
            (137.15459299813418, 211.65993512491434)
        );
        // (1022.86946623485993462, 1.20875811927763668424)
        assert_eq!(
            sip.apply_inverse(1023.0, 1.0).unwrap(),
            (1022.8694662348599, 1.2087581192776367)
        );
        // (611.992813487337212381, 586.994000911083383074)
        assert_eq!(
            sip.apply_inverse(612.0, 587.0).unwrap(),
            (611.9928134873372, 586.9940009110834)
        );
    }

    #[test]
    fn test_inverse_rejects_non_finite_points() {
        let mut iterative = SipDistortion::new([512.0, 512.0], 2, 2);
        iterative.set_a(2, 0, 1e-6).unwrap();
        let mut analytic = iterative.clone();
        analytic.set_inverse_order(2, 2).unwrap();
        analytic.set_ap(2, 0, -1e-6).unwrap();
        for (x, y) in [
            (f64::NAN, 37.5),
            (37.5, f64::NAN),
            (f64::INFINITY, 37.5),
            (37.5, f64::NEG_INFINITY),
        ] {
            let result = iterative.apply_inverse(x, y);
            assert!(result.is_err(), "iterative ({x}, {y}): {result:?}");
            let result = analytic.apply_inverse(x, y);
            assert!(result.is_err(), "AP/BP ({x}, {y}): {result:?}");
        }
    }

    #[test]
    fn test_large_pixel_offsets() {
        let mut sip = SipDistortion::new([512.0, 512.0], 3, 3);
        sip.set_a(2, 0, 1e-7).unwrap();
        sip.set_a(0, 2, 1e-7).unwrap();
        sip.set_a(3, 0, 1e-10).unwrap();
        sip.set_b(2, 0, 1e-7).unwrap();
        sip.set_b(0, 2, 1e-7).unwrap();
        sip.set_b(0, 3, 1e-10).unwrap();

        // Each value either way is the nearest double to the exact result.
        for (pixel, distorted) in [
            ((1512.0, 1512.0), (1512.3, 1512.3)),
            ((-488.0, -488.0), (-487.9, -487.9)),
            ((1512.0, -488.0), (1512.3, -487.9)),
            ((-488.0, 1512.0), (-487.9, 1512.3)),
        ] {
            assert_eq!(sip.apply(pixel.0, pixel.1), distorted, "{pixel:?}");
            let back = sip.apply_inverse(distorted.0, distorted.1).unwrap();
            assert_eq!(back, pixel, "{distorted:?}");
        }
    }

    fn cubic_survey() -> SipDistortion {
        let mut sip = SipDistortion::new([1024.5, 768.25], 3, 3);
        let terms = [
            (0, 2, 3.71e-6),
            (0, 3, -2.29e-10),
            (1, 1, -1.13e-6),
            (1, 2, 8.87e-11),
            (2, 0, 2.53e-6),
            (2, 1, -4.41e-10),
            (3, 0, 1.97e-10),
        ];
        for (p, q, c) in terms {
            sip.set_a(p, q, c).unwrap();
            sip.set_b(q, p, -0.7 * c).unwrap();
        }
        sip
    }

    #[test]
    fn test_identical_coefficients_give_identical_results() {
        let points = [(10.0, 20.0), (2011.7, 1533.3), (-300.25, 1700.5)];
        let first: Vec<_> = points
            .iter()
            .map(|&(x, y)| cubic_survey().apply(x, y))
            .collect();
        for _ in 0..64 {
            let sip = cubic_survey();
            let again: Vec<_> = points.iter().map(|&(x, y)| sip.apply(x, y)).collect();
            assert_eq!(again, first);
        }
    }

    #[test]
    fn test_setting_a_coefficient_to_zero_removes_it() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6).unwrap();
        sip.set_a(2, 0, 0.0).unwrap();
        sip.set_b(0, 2, 1e-6).unwrap();
        sip.set_b(0, 2, 0.0).unwrap();
        assert_eq!(sip.apply(612.0, 612.0), (612.0, 612.0));
    }

    #[test]
    fn test_set_a_stores_only_nonzero_terms() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 0.0).unwrap();
        sip.set_a(1, 1, 1e-6).unwrap();
        assert_eq!(sip.a_coeffs, vec![((1, 1), 1e-6)]);
    }

    #[test]
    fn test_setters_reject_terms_they_cannot_hold() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        let rejected = [
            sip.set_a(3, 0, 1e-6),
            sip.set_b(2, 1, 1e-6),
            sip.set_a(1, 1, f64::NAN),
            sip.set_b(1, 1, f64::INFINITY),
            sip.set_ap(1, 0, 1e-6),
            sip.set_bp(0, 1, 1e-6),
        ];
        for (i, result) in rejected.into_iter().enumerate() {
            assert!(
                matches!(result, Err(WcsError::InvalidParameter { .. })),
                "{i}: {result:?}"
            );
        }
        sip.set_inverse_order(2, 2).unwrap();
        assert!(sip.set_ap(3, 0, 1e-6).is_err());
        assert_eq!(sip.apply(612.0, 612.0), (612.0, 612.0));
        assert_eq!(sip.apply_inverse(612.0, 612.0).unwrap(), (612.0, 612.0));
    }

    #[test]
    fn test_inverse_order_cannot_drop_below_a_stored_term() {
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_inverse_order(3, 3).unwrap();
        sip.set_ap(3, 0, 1e-9).unwrap();
        assert!(sip.set_inverse_order(2, 3).is_err());
        assert_eq!(sip.ap_order, Some(3));
        sip.set_inverse_order(3, 2).unwrap();
    }
}
