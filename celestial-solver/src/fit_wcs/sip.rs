use std::collections::BTreeMap;

/// SIP (Simple Imaging Polynomial) distortion fitted together with its own TAN
/// parameters.
///
/// Stores the CRVAL and CD matrix of that joint fit plus polynomial coefficient maps
/// `A_p_q` and `B_p_q`. The CRVAL and CD matrix go only with these polynomials: they
/// differ from those of the TAN-only
/// [`WcsSolution`](super::solution::WcsSolution) fitted to the same pairs.
/// Coefficient keys are `(p, q)` monomial powers in the detected pixel offsets
/// `(u, v) = (px - CRPIX1, py - CRPIX2)`. The intermediate pixel coordinates are
/// `x = u + A(u, v)` and `y = v + B(u, v)`, the FITS SIP convention.
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
/// if let Some(sip) = &result.sip {
///     println!("SIP order: {}, RMS: {:.3}px", sip.a_order, sip.rms_px);
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct SipSolution {
    /// Reference pixel x, the same as that of the TAN-only WCS.
    pub crpix1: f64,
    /// Reference pixel y.
    pub crpix2: f64,
    /// RA at CRPIX.
    pub crval1: f64,
    /// Dec at CRPIX.
    pub crval2: f64,
    /// CD matrix element (1,1).
    pub cd1_1: f64,
    /// CD matrix element (1,2).
    pub cd1_2: f64,
    /// CD matrix element (2,1).
    pub cd2_1: f64,
    /// CD matrix element (2,2).
    pub cd2_2: f64,
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Polynomial order for the `A` (x-correction) series.
    pub a_order: u32,
    /// Polynomial order for the `B` (y-correction) series.
    pub b_order: u32,
    /// `A_p_q` coefficients keyed by `(p, q)`.
    pub a_coeffs: BTreeMap<(u32, u32), f64>,
    /// `B_p_q` coefficients keyed by `(p, q)`.
    pub b_coeffs: BTreeMap<(u32, u32), f64>,
    /// Number of stars used in the SIP fit.
    pub n_stars: usize,
    /// RMS residual in pixels for the linear (TAN-only) fit before SIP.
    pub linear_rms_px: f64,
    /// RMS residual in pixels after applying SIP.
    pub rms_px: f64,
    /// RMS residual after SIP, with each pair weighted by the inverse of its centroid
    /// variance.
    pub weighted_rms_px: f64,
}

pub(super) fn sip_monomial_terms(order: u32) -> Vec<(u32, u32)> {
    let mut terms = Vec::new();
    for total in 2..=order {
        for p in (0..=total).rev() {
            terms.push((p, total - p));
        }
    }
    terms
}

// The SIP distortion of `order` that bends the TAN WCS `wcs`, held fixed, onto `pairs`
// by unweighted least squares. A test fixture for a distortion fitted without its TAN.
#[cfg(test)]
pub(crate) fn fit_sip_to_tan(
    pairs: &[crate::match_field::StarPair],
    wcs: &super::solution::WcsSolution,
    order: u32,
) -> SipSolution {
    use super::normal_equations::{monomials, NormalEquations};
    let terms = sip_monomial_terms(order);
    let scale = 1000.0;
    let mut normal = NormalEquations::new(terms.len());
    for p in pairs {
        let (x, y) = wcs
            .sky_to_pixel(p.ra_deg, p.dec_deg)
            .expect("pairs project");
        let row = monomials(
            &terms,
            (p.px_x - wcs.crpix1) / scale,
            (p.px_y - wcs.crpix2) / scale,
            order,
        );
        normal.add(&row, 1.0, x - p.px_x, y - p.px_y);
    }
    let (a, b) = normal
        .solve()
        .expect("the fixture pairs determine the terms");
    let unscale = |c: &[f64]| -> BTreeMap<(u32, u32), f64> {
        let unscaled = terms
            .iter()
            .zip(c)
            .map(|(&(p, q), &c)| ((p, q), c / libm::pow(scale, f64::from(p + q))));
        unscaled.collect()
    };
    SipSolution {
        crpix1: wcs.crpix1,
        crpix2: wcs.crpix2,
        crval1: wcs.crval1,
        crval2: wcs.crval2,
        cd1_1: wcs.cd1_1,
        cd1_2: wcs.cd1_2,
        cd2_1: wcs.cd2_1,
        cd2_2: wcs.cd2_2,
        width: wcs.width,
        height: wcs.height,
        a_order: order,
        b_order: order,
        a_coeffs: unscale(&a),
        b_coeffs: unscale(&b),
        n_stars: pairs.len(),
        linear_rms_px: wcs.rms_px,
        rms_px: 0.0,
        weighted_rms_px: 0.0,
    }
}

pub(super) fn eval_sip_poly(coeffs: &BTreeMap<(u32, u32), f64>, u: f64, v: f64) -> f64 {
    coeffs
        .iter()
        .map(|(&(p, q), &c)| c * libm::pow(u, p as f64) * libm::pow(v, q as f64))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sip_monomial_terms_order_2_has_three_terms() {
        let terms = sip_monomial_terms(2);
        // total = 2: (2,0), (1,1), (0,2)
        assert_eq!(terms.len(), 3);
        assert!(terms.contains(&(2, 0)));
        assert!(terms.contains(&(1, 1)));
        assert!(terms.contains(&(0, 2)));
    }

    #[test]
    fn sip_monomial_terms_order_3_has_seven_terms() {
        // order 2 (3 terms) + order 3 (4 terms) = 7
        let terms = sip_monomial_terms(3);
        assert_eq!(terms.len(), 7);
        for (p, q) in [(3, 0), (2, 1), (1, 2), (0, 3)] {
            assert!(terms.contains(&(p, q)), "missing ({}, {})", p, q);
        }
    }

    #[test]
    fn sip_monomial_terms_order_4_has_twelve_terms() {
        // 3 + 4 + 5 = 12
        let terms = sip_monomial_terms(4);
        assert_eq!(terms.len(), 12);
    }

    #[test]
    fn eval_sip_poly_evaluates_known_polynomial() {
        let mut coeffs = BTreeMap::new();
        coeffs.insert((2, 0), 3.0);
        coeffs.insert((0, 2), 5.0);
        coeffs.insert((1, 1), 7.0);
        // at u=2, v=3: 3*4 + 5*9 + 7*6 = 12 + 45 + 42 = 99
        let v = eval_sip_poly(&coeffs, 2.0, 3.0);
        assert_eq!(v, 99.0);
    }

    #[test]
    fn eval_sip_poly_is_identical_across_coefficient_maps() {
        // Order-5 terms with alternating signs and spread magnitudes, so the sum
        // depends on the order the terms are added.
        let terms: Vec<((u32, u32), f64)> = sip_monomial_terms(5)
            .into_iter()
            .enumerate()
            .map(|(i, t)| (t, libm::pow(-3.7, i as f64) * 1e-9))
            .collect();
        let first = eval_sip_poly(&terms.iter().copied().collect(), 1234.5, -987.25);
        for _ in 0..20 {
            let coeffs = terms.iter().copied().collect();
            assert_eq!(
                eval_sip_poly(&coeffs, 1234.5, -987.25).to_bits(),
                first.to_bits()
            );
        }
    }

    #[test]
    fn eval_sip_poly_empty_is_zero() {
        let coeffs: BTreeMap<(u32, u32), f64> = BTreeMap::new();
        assert_eq!(eval_sip_poly(&coeffs, 5.0, 7.0), 0.0);
    }
}
