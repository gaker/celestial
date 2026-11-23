use crate::common::require_finite;
use crate::error::{WcsError, WcsResult};
use crate::newton::newton_raphson_2d;

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
        tpv.pv1[1] = 1.0;
        tpv.pv2[1] = 1.0;
        tpv
    }

    pub fn set_pv1(&mut self, index: usize, value: f64) -> WcsResult<()> {
        Self::set_term(&mut self.pv1, "PV1", index, value)
    }

    pub fn set_pv2(&mut self, index: usize, value: f64) -> WcsResult<()> {
        Self::set_term(&mut self.pv2, "PV2", index, value)
    }

    fn set_term(coeffs: &mut [f64; 40], axis: &str, index: usize, value: f64) -> WcsResult<()> {
        let term = coeffs.get_mut(index).ok_or_else(|| {
            WcsError::invalid_parameter(format!("{axis}_{index} is not a TPV term (0 to 39)"))
        })?;
        if !value.is_finite() {
            return Err(WcsError::invalid_parameter(format!(
                "{axis}_{index} is {value}"
            )));
        }
        *term = value;
        Ok(())
    }

    pub fn get_pv1(&self, index: usize) -> Option<f64> {
        self.pv1.get(index).copied()
    }

    pub fn get_pv2(&self, index: usize) -> Option<f64> {
        self.pv2.get(index).copied()
    }

    pub(crate) fn axes(&self) -> [(&'static str, &[f64; 40]); 2] {
        [("PV1", &self.pv1), ("PV2", &self.pv2)]
    }

    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let r = libm::sqrt(x * x + y * y);
        let xi = Self::eval_polynomial(&self.pv1, x, y, r);
        // Axis 2 uses the same terms with ξ and η swapped.
        let eta = Self::eval_polynomial(&self.pv2, y, x, r);
        (xi, eta)
    }

    pub fn apply_inverse(&self, xi: f64, eta: f64) -> WcsResult<(f64, f64)> {
        require_finite(&[xi, eta], "TPV inverse input")?;
        let distort_fn = |x: f64, y: f64| self.apply(x, y);

        newton_raphson_2d(
            (xi, eta),
            (xi, eta),
            distort_fn,
            20,
            "TPV inverse distortion",
        )
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
            17 => x * x * x * x * x,
            18 => x * x * x * x * y,
            19 => x * x * x * y * y,
            20 => x * x * y * y * y,
            21 => x * y * y * y * y,
            22 => y * y * y * y * y,
            23 => r * r * r * r * r,
            24 => x * x * x * x * x * x,
            25 => x * x * x * x * x * y,
            26 => x * x * x * x * y * y,
            27 => x * x * x * y * y * y,
            28 => x * x * y * y * y * y,
            29 => x * y * y * y * y * y,
            30 => y * y * y * y * y * y,
            31 => x * x * x * x * x * x * x,
            32 => x * x * x * x * x * x * y,
            33 => x * x * x * x * x * y * y,
            34 => x * x * x * x * y * y * y,
            35 => x * x * x * y * y * y * y,
            36 => x * x * y * y * y * y * y,
            37 => x * y * y * y * y * y * y,
            38 => y * y * y * y * y * y * y,
            39 => r * r * r * r * r * r * r,
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
    fn test_identity_and_zero_distortion() {
        // identity() seeds PV1_1 = PV2_1 = 1 (ξ on axis 1, η on axis 2), so
        // apply() is an identity.  new() leaves all coefficients at zero, so apply()
        // produces (0, 0) regardless of input.
        let id = TpvDistortion::identity();
        assert_eq!(id.apply(1.5, 2.5), (1.5, 2.5));
        let mut by_keyword = TpvDistortion::new();
        by_keyword.set_pv1(1, 1.0).unwrap();
        by_keyword.set_pv2(1, 1.0).unwrap();
        assert_eq!(by_keyword.apply(1.5, 2.5), (1.5, 2.5));

        let zero = TpvDistortion::new();
        assert_eq!(zero.apply(1.5, 2.5), (0.0, 0.0));
    }

    #[test]
    fn test_apply_uses_each_term_correctly() {
        // Exercise the apply() path with a handful of representative term
        // indices (linear scale, radial, x^2, y^2, xy).
        // test_each_term_follows_tpv_numbering verifies the full term mapping.
        let mut linear = TpvDistortion::new();
        linear.set_pv1(1, 2.0).unwrap();
        linear.set_pv2(1, 2.0).unwrap();
        assert_eq!(linear.apply(3.0, 4.0), (6.0, 8.0));

        let mut radial = TpvDistortion::identity();
        radial.set_pv1(3, 0.01).unwrap();
        radial.set_pv2(3, 0.01).unwrap();
        let (x, y): (f64, f64) = (3.0, 4.0);
        let r = libm::sqrt(x * x + y * y);
        let (xi, eta) = radial.apply(x, y);
        assert_eq!(xi, x + 0.01 * r);
        assert_eq!(eta, y + 0.01 * r);

        let mut xsq = TpvDistortion::identity();
        xsq.set_pv1(4, 0.1).unwrap();
        assert_eq!(xsq.apply(2.0, 1.0), (2.0 + 0.1 * 2.0 * 2.0, 1.0));

        let mut cross = TpvDistortion::identity();
        cross.set_pv1(5, 0.05).unwrap();
        assert_eq!(cross.apply(2.0, 3.0), (2.0 + 0.05 * 2.0 * 3.0, 3.0));
    }

    #[test]
    fn test_mixed_distortion_roundtrip() {
        // A combined model (linear + radial + quadratic + cubic + r^3) subsumes
        // the per-term roundtrip cases.  Five sign-pattern points cover every
        // quadrant of the projection plane.
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(3, 0.001).unwrap(); // r
        tpv.set_pv1(4, 0.002).unwrap(); // x^2
        tpv.set_pv1(5, 0.001).unwrap(); // xy
        tpv.set_pv1(7, 0.001).unwrap(); // x^3
        tpv.set_pv1(11, 0.0005).unwrap(); // r^3
        tpv.set_pv2(3, 0.001).unwrap();
        tpv.set_pv2(5, 0.001).unwrap();
        tpv.set_pv2(4, 0.002).unwrap(); // y^2
        tpv.set_pv2(7, 0.001).unwrap(); // y^3
        tpv.set_pv2(11, 0.0005).unwrap();

        // Each value either way is the nearest double to the exact result.
        for (point, distorted) in [
            ((0.1, 0.1), (0.10017383556979968, 0.10017383556979968)),
            ((0.5, 0.3), (0.5014572213716969, 0.30103922137169686)),
            ((-0.2, 0.4), (-0.19951606504495006, 0.40079593495504995)),
            ((0.3, -0.5), (0.3007392213716969, -0.4990927786283031)),
            ((-0.4, -0.4), (-0.3989278049070589, -0.3989278049070589)),
        ] {
            assert_eq!(tpv.apply(point.0, point.1), distorted, "{point:?}");
            let back = tpv.apply_inverse(distorted.0, distorted.1).unwrap();
            assert_eq!(back, point, "{distorted:?}");
        }
    }

    // A third-order solution of the size SCAMP writes. The inverse of (0, 0)
    // is (−1.19948258790982265930e-5, 7.50353988488803938233e-6); ξ is
    // 0.57 ulp from it and η −0.26. The nearest ξ is 0.43 ulp away, but the
    // residual cancels PV1_0 = 1.2e-5 against PV1_1 ξ, so its rounding is
    // about one ulp of ξ.
    #[test]
    fn test_inverse_solves_to_rounding() {
        let terms = [
            (0, 1.2e-5, -7.5e-6),
            (1, 1.0003, 0.9998),
            (2, -2.1e-4, 1.7e-4),
            (4, 3.4e-4, -2.6e-4),
            (5, -1.8e-4, 2.9e-4),
            (6, 9.0e-5, -1.1e-4),
            (7, -2.2e-2, 3.8e-4),
            (8, 4.1e-4, -2.1e-2),
            (9, -2.0e-2, 5.2e-4),
            (10, 3.3e-4, -2.3e-2),
        ];
        let mut tpv = TpvDistortion::new();
        for (m, pv1, pv2) in terms {
            tpv.set_pv1(m, pv1).unwrap();
            tpv.set_pv2(m, pv2).unwrap();
        }
        assert_eq!(
            tpv.apply_inverse(0.0, 0.0).unwrap(),
            (-1.1994825879098226e-5, 7.503539884888039e-6)
        );
    }

    #[test]
    fn test_inverse_rejects_non_finite_points() {
        for tpv in [TpvDistortion::new(), TpvDistortion::identity()] {
            for (x, y) in [
                (f64::NAN, 37.5),
                (37.5, f64::NAN),
                (f64::INFINITY, 37.5),
                (37.5, f64::NEG_INFINITY),
            ] {
                let result = tpv.apply_inverse(x, y);
                assert!(result.is_err(), "({x}, {y}): {result:?}");
            }
        }
    }

    #[test]
    fn test_getters() {
        let mut tpv = TpvDistortion::new();
        tpv.set_pv1(5, 1.23).unwrap();
        tpv.set_pv2(10, 4.56).unwrap();

        assert_eq!(tpv.get_pv1(5), Some(1.23));
        assert_eq!(tpv.get_pv2(10), Some(4.56));
        assert_eq!(tpv.get_pv1(40), None);
        assert_eq!(tpv.get_pv2(100), None);
    }

    #[test]
    fn test_setters_reject_terms_they_cannot_hold() {
        let mut tpv = TpvDistortion::identity();
        let rejected = [
            tpv.set_pv1(40, 1.0),
            tpv.set_pv2(100, 2.0),
            tpv.set_pv1(3, f64::NAN),
            tpv.set_pv2(3, f64::NEG_INFINITY),
        ];
        for (i, result) in rejected.into_iter().enumerate() {
            assert!(
                matches!(result, Err(WcsError::InvalidParameter { .. })),
                "{i}: {result:?}"
            );
        }
        assert_eq!(tpv.apply(1.5, 2.5), (1.5, 2.5));
    }

    #[test]
    fn test_constant_term() {
        let mut tpv = TpvDistortion::new();
        tpv.set_pv1(0, 1.0).unwrap();
        tpv.set_pv2(0, 2.0).unwrap();

        let (xi, eta) = tpv.apply(0.0, 0.0);
        assert_eq!(xi, 1.0);
        assert_eq!(eta, 2.0);
    }

    // TPV numbering at (ξ, η) = (3, 4), where r = 5 and every term is an exact
    // integer. Even powers of r are not terms; axis 2 swaps ξ and η.
    const TERMS_AT_3_4: [(usize, f64, f64); 40] = [
        (0, 1.0, 1.0),
        (1, 3.0, 4.0),
        (2, 4.0, 3.0),
        (3, 5.0, 5.0),
        (4, 9.0, 16.0),
        (5, 12.0, 12.0),
        (6, 16.0, 9.0),
        (7, 27.0, 64.0),
        (8, 36.0, 48.0),
        (9, 48.0, 36.0),
        (10, 64.0, 27.0),
        (11, 125.0, 125.0),
        (12, 81.0, 256.0),
        (13, 108.0, 192.0),
        (14, 144.0, 144.0),
        (15, 192.0, 108.0),
        (16, 256.0, 81.0),
        (17, 243.0, 1024.0),
        (18, 324.0, 768.0),
        (19, 432.0, 576.0),
        (20, 576.0, 432.0),
        (21, 768.0, 324.0),
        (22, 1024.0, 243.0),
        (23, 3125.0, 3125.0),
        (24, 729.0, 4096.0),
        (25, 972.0, 3072.0),
        (26, 1296.0, 2304.0),
        (27, 1728.0, 1728.0),
        (28, 2304.0, 1296.0),
        (29, 3072.0, 972.0),
        (30, 4096.0, 729.0),
        (31, 2187.0, 16384.0),
        (32, 2916.0, 12288.0),
        (33, 3888.0, 9216.0),
        (34, 5184.0, 6912.0),
        (35, 6912.0, 5184.0),
        (36, 9216.0, 3888.0),
        (37, 12288.0, 2916.0),
        (38, 16384.0, 2187.0),
        (39, 78125.0, 78125.0),
    ];

    #[test]
    fn test_each_term_follows_tpv_numbering() {
        for (m, xi_term, eta_term) in TERMS_AT_3_4 {
            let mut tpv = TpvDistortion::new();
            tpv.set_pv1(m, 1.0).unwrap();
            tpv.set_pv2(m, 1.0).unwrap();
            assert_eq!(
                tpv.apply(3.0, 4.0),
                (xi_term, eta_term),
                "PV1_{m} / PV2_{m}"
            );
        }
    }
}
