//! The inverse SIP polynomials AP and BP, fitted to the forward A and B on a grid.
//!
//! The order, grid and terms follow astrometry.net 0.97 (`sip_ensure_inverse_polynomials`
//! and `sip_compute_inverse_polynomials` in sip-utils.c): the inverse order is one above
//! the forward order, the grid has 10 * (inverse order + 1) points along each axis,
//! evenly spaced over FITS pixel coordinates 0 to the image width and 0 to the image
//! height, and the fit includes every term with p + q up to the inverse order, the
//! constant and linear terms too. At each grid point (u, v), an offset from CRPIX, the
//! forward polynomials give U = u + A(u, v) and V = v + B(u, v), and AP and BP are the
//! least-squares fit of u - U and v - V as polynomials in (U, V).

use std::collections::BTreeMap;

use super::normal_equations::{monomials, terms_up_to, NormalEquations};
use super::polynomial_fit::scale_of;
use super::sip::{eval_sip_poly, SipSolution};
use crate::error::SolveError;

// The SIP convention bounds A_ORDER and B_ORDER to 9, and celestial_wcs reads AP_ORDER
// and BP_ORDER only up to 9 too.
const MAX_INVERSE_ORDER: u32 = 9;

pub(super) type Coefficients = Vec<((u32, u32), f64)>;

pub(super) fn inverse_order(sip: &SipSolution) -> u32 {
    (sip.a_order.max(sip.b_order) + 1).min(MAX_INVERSE_ORDER)
}

// AP and BP of `order` for `sip`, over an image of `sip.width` by `sip.height` pixels
// whose CRPIX is `crpix` in FITS pixel coordinates.
pub(super) fn fit_inverse(
    sip: &SipSolution,
    crpix: (f64, f64),
    order: u32,
) -> Result<(Coefficients, Coefficients), SolveError> {
    let terms = terms_up_to(order);
    let grid = distorted_grid(sip, crpix, 10 * (order as usize + 1));
    // Monomials of U and V divided by the largest |U| and |V| of the grid lie in
    // [-1, 1] along both axes, which keeps the normal equations well conditioned at
    // high orders on images of any shape.
    let extent = |m: f64, x: f64| m.max(libm::fabs(x));
    let scale = grid
        .iter()
        .fold((1.0, 1.0), |(su, sv), &(_, _, big_u, big_v)| {
            (extent(su, big_u), extent(sv, big_v))
        });
    let mut normal = NormalEquations::new(terms.len());
    for &(u, v, big_u, big_v) in &grid {
        let row = monomials(&terms, big_u / scale.0, big_v / scale.1, order);
        normal.add(&row, 1.0, u - big_u, v - big_v);
    }
    let (ap, bp) = normal.solve()?;
    Ok((unscaled(&terms, &ap, scale), unscaled(&terms, &bp, scale)))
}

// The `n` by `n` grid of offsets (u, v) from `crpix` over the image, each with
// U = u + A(u, v) and V = v + B(u, v).
fn distorted_grid(sip: &SipSolution, crpix: (f64, f64), n: usize) -> Vec<(f64, f64, f64, f64)> {
    let lo = (-crpix.0, -crpix.1);
    let hi = (sip.width as f64 - crpix.0, sip.height as f64 - crpix.1);
    let mut grid = Vec::with_capacity(n * n);
    for u in evenly_spaced(lo.0, hi.0, n) {
        for v in evenly_spaced(lo.1, hi.1, n) {
            let big_u = u + eval_sip_poly(&sip.a_coeffs, u, v);
            let big_v = v + eval_sip_poly(&sip.b_coeffs, u, v);
            grid.push((u, v, big_u, big_v));
        }
    }
    grid
}

// `n` values from `lo` to `hi`, both included, a constant step apart.
fn evenly_spaced(lo: f64, hi: f64, n: usize) -> impl Iterator<Item = f64> {
    (0..n).map(move |i| (i as f64 * (hi - lo) / (n - 1) as f64) + lo)
}

// A fitted coefficient multiplies (U / scale.0)^p (V / scale.1)^q.
fn unscaled(terms: &[(u32, u32)], fitted: &[f64], scale: (f64, f64)) -> Coefficients {
    let by_term: BTreeMap<(u32, u32), f64> = terms
        .iter()
        .zip(fitted)
        .map(|(&(p, q), &c)| ((p, q), c / scale_of((p, q), scale)))
        .collect();
    by_term.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // The SIP that scratch probes evaluated with astrometry.net 0.97's
    // sip_ensure_inverse_polynomials: a third-order forward polynomial on a 2048 x 1536
    // image with CRPIX (1024.5, 768.5) in FITS pixels.
    fn sample_sip() -> SipSolution {
        let a = [
            ((2, 0), 2e-7),
            ((1, 1), -1e-7),
            ((0, 2), 5e-8),
            ((3, 0), 1e-10),
            ((1, 2), -3e-11),
        ];
        let b = [
            ((0, 2), -1.5e-7),
            ((2, 0), 4e-8),
            ((1, 1), 9e-8),
            ((0, 3), 2e-10),
            ((2, 1), 5e-11),
        ];
        SipSolution {
            crpix1: 1023.5,
            crpix2: 767.5,
            crval1: 180.0,
            crval2: 30.0,
            cd1_1: -2.8e-4,
            cd1_2: 1e-6,
            cd2_1: 1.2e-6,
            cd2_2: 2.8e-4,
            width: 2048,
            height: 1536,
            a_order: 3,
            b_order: 3,
            a_coeffs: a.into_iter().collect(),
            b_coeffs: b.into_iter().collect(),
            n_stars: 0,
            linear_rms_px: 0.0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
        }
    }

    #[test]
    fn inverse_order_is_one_above_the_forward_order_up_to_nine() {
        let mut sip = sample_sip();
        assert_eq!(inverse_order(&sip), 4);
        sip.a_order = 9;
        assert_eq!(inverse_order(&sip), 9);
    }

    // astrometry.net 0.97 gives AP_2_0 = -1.9999994206649226e-07 and BP_0_2 =
    // 1.5000012521959845e-07 for this SIP. Over a 101 x 101 grid of the image its AP and
    // BP and these differ by at most 4.7e-15 px.
    #[test]
    fn fit_inverse_matches_the_astrometry_net_inverse() {
        let sip = sample_sip();
        let (ap, bp) = fit_inverse(&sip, (1024.5, 768.5), inverse_order(&sip)).unwrap();
        assert_eq!(ap.len(), 15);
        assert_eq!(bp.len(), 15);
        let ap: BTreeMap<(u32, u32), f64> = ap.into_iter().collect();
        let bp: BTreeMap<(u32, u32), f64> = bp.into_iter().collect();
        assert_eq!(ap[&(2, 0)], -1.9999994206649486e-7);
        assert_eq!(bp[&(0, 2)], 1.500001252195855e-7);
        assert_eq!(ap[&(0, 0)], -3.2970341392241087e-9);
    }

    // On a 16000 x 300 image each axis's monomials are scaled by that axis's own extent,
    // so the order 5 inverse of an order 4 forward polynomial is determined. The worst
    // round trip u -> U -> u over a 21 x 21 grid of the image.
    #[test]
    fn fit_inverse_round_trips_on_an_elongated_image() {
        let sip = SipSolution {
            crpix1: 7999.5,
            crpix2: 149.5,
            width: 16000,
            height: 300,
            a_order: 4,
            b_order: 4,
            ..sample_sip()
        };
        let crpix = (8000.5, 150.5);
        let (ap, bp) = fit_inverse(&sip, crpix, inverse_order(&sip)).unwrap();
        let eval = |c: &Coefficients, u: f64, v: f64| -> f64 {
            c.iter()
                .map(|&((p, q), k)| k * libm::pow(u, f64::from(p)) * libm::pow(v, f64::from(q)))
                .sum()
        };
        let mut worst = 0.0_f64;
        for i in 0..=20 {
            for j in 0..=20 {
                let u = -crpix.0 + 16000.0 * f64::from(i) / 20.0;
                let v = -crpix.1 + 300.0 * f64::from(j) / 20.0;
                let big_u = u + eval_sip_poly(&sip.a_coeffs, u, v);
                let big_v = v + eval_sip_poly(&sip.b_coeffs, u, v);
                let du = big_u + eval(&ap, big_u, big_v) - u;
                let dv = big_v + eval(&bp, big_u, big_v) - v;
                worst = worst.max(libm::hypot(du, dv));
            }
        }
        assert_eq!(worst, 0.0014152222368832004);
    }
}
