use std::collections::HashMap;

use anyhow::{ensure, Result};
use celestial_catalog::query::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;

use crate::match_field::StarPair;

use super::cd_matrix::solve_normal_equations;
use super::solution::WcsSolution;

#[derive(Debug, Clone, Copy)]
struct Position {
    px_x: f64,
    px_y: f64,
    x_ideal: f64,
    y_ideal: f64,
}

/// SIP (Simple Imaging Polynomial) distortion fit on top of a TAN WCS.
///
/// Stores the TAN WCS parameters plus polynomial coefficient maps `A_p_q` and
/// `B_p_q`. Coefficient keys are `(p, q)` monomial powers in the undistorted pixel
/// offsets `(u, v) = (px - CRPIX1, py - CRPIX2)`. Order matches the FITS SIP
/// convention.
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
/// if let Some(sip) = &result.sip {
///     println!("SIP order: {}, RMS: {:.3}px", sip.a_order, sip.rms_px);
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct SipSolution {
    /// Reference pixel x (same as the underlying TAN WCS).
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
    pub a_coeffs: HashMap<(u32, u32), f64>,
    /// `B_p_q` coefficients keyed by `(p, q)`.
    pub b_coeffs: HashMap<(u32, u32), f64>,
    /// Number of stars used in the SIP fit.
    pub n_stars: usize,
    /// RMS residual in pixels for the linear (TAN-only) fit before SIP.
    pub linear_rms_px: f64,
    /// RMS residual in pixels after applying SIP.
    pub rms_px: f64,
    /// SNR-weighted RMS residual after SIP.
    pub weighted_rms_px: f64,
}

/// Fits SIP polynomial distortion to `pairs` using `wcs` as the linear prior.
///
/// Least-squares fit of the `A` and `B` monomial coefficients up to `order`. Both `A`
/// and `B` use the same order. Produces a [`SipSolution`] with the underlying TAN
/// parameters copied from `wcs` and the distortion coefficients filled in.
///
/// Usually called internally by [`crate::solve()`] when
/// [`crate::SolveParams::sip_order`] is `Some`. Use directly when you want to refit
/// SIP against a modified pair set.
///
/// # Errors
///
/// Returns an error when `order < 2`, fewer than 6 pairs are supplied, fewer pairs
/// than polynomial terms are available, or the normal-equations system is singular.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::fit_sip;
///
/// # let pairs = vec![];
/// # let wcs = unimplemented!();
/// let sip = fit_sip(&pairs, &wcs, 4)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn fit_sip(
    pairs: &[StarPair],
    wcs: &WcsSolution,
    order: u32,
) -> Result<SipSolution> {
    ensure!(order >= 2, "SIP order must be >= 2, got {}", order);
    ensure!(pairs.len() >= 6, "need at least 6 pairs for SIP fit, got {}", pairs.len());

    let terms = sip_monomial_terms(order);
    let n_terms = terms.len();
    ensure!(
        pairs.len() >= n_terms,
        "need at least {} pairs for order {} SIP, got {}",
        n_terms, order, pairs.len()
    );

    let positions = compute_ideal_positions(pairs, wcs);
    ensure!(positions.len() >= n_terms, "too few projectable stars for SIP fit");

    let sip_weights: Vec<f64> = pairs
        .iter()
        .filter_map(|p| {
            tan_project_star(p.ra_deg, p.dec_deg, wcs.crval1, wcs.crval2)
                .map(|_| p.snr * p.snr)
        })
        .collect();

    let mut mask = vec![true; positions.len()];
    let linear_rms = compute_linear_rms(&positions, &mask);

    let mut a_coeffs;
    let mut b_coeffs;

    let max_iters = 10;
    for iter in 0..max_iters {
        a_coeffs = solve_sip_coefficients(
            &positions, &mask, &sip_weights, wcs.crpix1, wcs.crpix2, &terms, 0,
        )?;
        b_coeffs = solve_sip_coefficients(
            &positions, &mask, &sip_weights, wcs.crpix1, wcs.crpix2, &terms, 1,
        )?;

        let errs = compute_sip_residuals(
            &positions, &mask, wcs.crpix1, wcs.crpix2,
            &a_coeffs, &b_coeffs,
        );

        let mut sum_sr2 = 0.0_f64;
        let mut n_active = 0_usize;
        for (i, &e) in errs.iter().enumerate() {
            if !mask[i] { continue; }
            sum_sr2 += sip_weights[i] * e * e;
            n_active += 1;
        }
        if (n_active as f64) < n_terms as f64 {
            break;
        }
        let sr_rms = libm::sqrt(sum_sr2 / n_active as f64);
        let threshold = 2.5 * sr_rms;

        let mut rejected = 0usize;
        for (i, &e) in errs.iter().enumerate() {
            if mask[i] && e * libm::sqrt(sip_weights[i]) > threshold {
                mask[i] = false;
                rejected += 1;
            }
        }

        let n_active = mask.iter().filter(|&&m| m).count();
        let rms = libm::sqrt(
            errs.iter().enumerate()
                .filter(|(i, _)| mask[*i])
                .map(|(_, &e)| e * e)
                .sum::<f64>() / n_active.max(1) as f64
        );
        log::debug!(
            "fit_sip iter {}: {} stars, rms={:.4} px, rejected {}",
            iter, n_active, rms, rejected,
        );

        if rejected == 0 {
            return build_sip_solution(
                wcs, order, a_coeffs, b_coeffs,
                &positions, &mask, &sip_weights, linear_rms,
            );
        }
    }

    a_coeffs = solve_sip_coefficients(
        &positions, &mask, &sip_weights, wcs.crpix1, wcs.crpix2, &terms, 0,
    )?;
    b_coeffs = solve_sip_coefficients(
        &positions, &mask, &sip_weights, wcs.crpix1, wcs.crpix2, &terms, 1,
    )?;

    build_sip_solution(wcs, order, a_coeffs, b_coeffs, &positions, &mask, &sip_weights, linear_rms)
}

fn sip_monomial_terms(order: u32) -> Vec<(u32, u32)> {
    let mut terms = Vec::new();
    for total in 2..=order {
        for p in (0..=total).rev() {
            terms.push((p, total - p));
        }
    }
    terms
}

fn eval_sip_poly(coeffs: &HashMap<(u32, u32), f64>, u: f64, v: f64) -> f64 {
    coeffs
        .iter()
        .map(|(&(p, q), &c)| c * libm::pow(u, p as f64) * libm::pow(v, q as f64))
        .sum()
}

fn compute_ideal_positions(
    pairs: &[StarPair],
    wcs: &WcsSolution,
) -> Vec<Position> {
    pairs
        .iter()
        .filter_map(|p| {
            let (xi_rad, eta_rad) =
                tan_project_star(p.ra_deg, p.dec_deg, wcs.crval1, wcs.crval2)?;
            let xi = xi_rad * RAD_TO_DEG;
            let eta = eta_rad * RAD_TO_DEG;
            let det = wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1;
            let x_ideal = (wcs.cd2_2 * xi - wcs.cd1_2 * eta) / det + wcs.crpix1;
            let y_ideal = (-wcs.cd2_1 * xi + wcs.cd1_1 * eta) / det + wcs.crpix2;
            Some(Position { px_x: p.px_x, px_y: p.px_y, x_ideal, y_ideal })
        })
        .collect()
}

fn solve_sip_coefficients(
    positions: &[Position],
    mask: &[bool],
    weights: &[f64],
    crpix1: f64,
    crpix2: f64,
    terms: &[(u32, u32)],
    axis: usize,
) -> Result<HashMap<(u32, u32), f64>> {
    let n_terms = terms.len();
    let active: Vec<(usize, &Position)> = positions
        .iter()
        .enumerate()
        .zip(mask.iter())
        .filter(|(_, &m)| m)
        .map(|((i, p), _)| (i, p))
        .collect();
    let n_pts = active.len();
    ensure!(n_pts >= n_terms, "need {} points, have {}", n_terms, n_pts);

    let mut ata = vec![0.0_f64; n_terms * n_terms];
    let mut atb = vec![0.0_f64; n_terms];

    for &(idx, row) in &active {
        let w = weights[idx];
        let u = row.px_x - crpix1;
        let v = row.px_y - crpix2;
        let residual = if axis == 0 { row.px_x - row.x_ideal } else { row.px_y - row.y_ideal };

        let basis: Vec<f64> = terms
            .iter()
            .map(|&(p, q)| libm::pow(u, p as f64) * libm::pow(v, q as f64))
            .collect();

        for i in 0..n_terms {
            for j in 0..n_terms {
                ata[i * n_terms + j] += w * basis[i] * basis[j];
            }
            atb[i] += w * basis[i] * residual;
        }
    }

    let coeffs_vec = solve_normal_equations(&mut ata, &mut atb, n_terms)?;

    let mut coeffs = HashMap::new();
    for (i, &c) in coeffs_vec.iter().enumerate() {
        if libm::fabs(c) > 1e-20 {
            coeffs.insert(terms[i], c);
        }
    }
    Ok(coeffs)
}

fn compute_sip_residuals(
    positions: &[Position],
    mask: &[bool],
    crpix1: f64,
    crpix2: f64,
    a_coeffs: &HashMap<(u32, u32), f64>,
    b_coeffs: &HashMap<(u32, u32), f64>,
) -> Vec<f64> {
    positions
        .iter()
        .zip(mask.iter())
        .map(|(row, &m)| {
            if !m {
                return 0.0;
            }
            let u = row.px_x - crpix1;
            let v = row.px_y - crpix2;
            let pred_du = eval_sip_poly(a_coeffs, u, v);
            let pred_dv = eval_sip_poly(b_coeffs, u, v);
            let actual_du = row.px_x - row.x_ideal;
            let actual_dv = row.px_y - row.y_ideal;
            let ex = actual_du - pred_du;
            let ey = actual_dv - pred_dv;
            libm::sqrt(ex * ex + ey * ey)
        })
        .collect()
}

fn compute_linear_rms(positions: &[Position], mask: &[bool]) -> f64 {
    let mut sum_sq = 0.0_f64;
    let mut n = 0usize;
    for (row, &m) in positions.iter().zip(mask.iter()) {
        if !m {
            continue;
        }
        let du = row.px_x - row.x_ideal;
        let dv = row.px_y - row.y_ideal;
        sum_sq += du * du + dv * dv;
        n += 1;
    }
    if n == 0 { return 0.0; }
    libm::sqrt(sum_sq / n as f64)
}

#[allow(clippy::too_many_arguments)]
fn build_sip_solution(
    wcs: &WcsSolution,
    order: u32,
    a_coeffs: HashMap<(u32, u32), f64>,
    b_coeffs: HashMap<(u32, u32), f64>,
    positions: &[Position],
    mask: &[bool],
    weights: &[f64],
    linear_rms: f64,
) -> Result<SipSolution> {
    let errs = compute_sip_residuals(
        positions, mask, wcs.crpix1, wcs.crpix2,
        &a_coeffs, &b_coeffs,
    );
    let mut n = 0_usize;
    let mut sum_sq = 0.0_f64;
    let mut sum_w = 0.0_f64;
    let mut sum_we2 = 0.0_f64;
    for (i, &e) in errs.iter().enumerate() {
        if !mask[i] { continue; }
        n += 1;
        sum_sq += e * e;
        let w = weights[i];
        sum_w += w;
        sum_we2 += w * e * e;
    }
    let rms = if n > 0 { libm::sqrt(sum_sq / n as f64) } else { 0.0 };
    let wrms = if sum_w > 0.0 { libm::sqrt(sum_we2 / sum_w) } else { 0.0 };

    Ok(SipSolution {
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
        a_coeffs,
        b_coeffs,
        n_stars: n,
        linear_rms_px: linear_rms,
        rms_px: rms,
        weighted_rms_px: wrms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::solution::WcsSolution;

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            crpix1: 0.0,
            crpix2: 0.0,
            crval1: 180.0,
            crval2: 0.0,
            cd1_1: 0.001,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            width: 1000,
            height: 1000,
            focal_mm: None,
            pixel_um: None,
            n_stars: 0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
            residuals: Vec::new(),
        }
    }

    fn pair_at_pixel(wcs: &WcsSolution, px_x: f64, px_y: f64) -> StarPair {
        // Invert the WCS to get the sky position that would project exactly to (px_x, px_y).
        let (ra, dec) = wcs.pixel_to_sky(px_x, px_y);
        StarPair { px_x, px_y, ra_deg: ra, dec_deg: dec, votes: 4, snr: 50.0 }
    }

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
        let mut coeffs = HashMap::new();
        coeffs.insert((2, 0), 3.0);
        coeffs.insert((0, 2), 5.0);
        coeffs.insert((1, 1), 7.0);
        // at u=2, v=3: 3*4 + 5*9 + 7*6 = 12 + 45 + 42 = 99
        let v = eval_sip_poly(&coeffs, 2.0, 3.0);
        assert_eq!(v, 99.0);
    }

    #[test]
    fn eval_sip_poly_empty_is_zero() {
        let coeffs: HashMap<(u32, u32), f64> = HashMap::new();
        assert_eq!(eval_sip_poly(&coeffs, 5.0, 7.0), 0.0);
    }

    #[test]
    fn fit_sip_rejects_order_below_two() {
        let wcs = sample_wcs();
        let pairs = vec![pair_at_pixel(&wcs, 0.0, 0.0); 10];
        let err = fit_sip(&pairs, &wcs, 1).unwrap_err();
        assert!(err.to_string().contains("SIP order must be >= 2"));
    }

    #[test]
    fn fit_sip_rejects_too_few_pairs() {
        let wcs = sample_wcs();
        let pairs = vec![pair_at_pixel(&wcs, 0.0, 0.0); 5];
        let err = fit_sip(&pairs, &wcs, 2).unwrap_err();
        assert!(err.to_string().contains("need at least 6 pairs"));
    }

    #[test]
    fn fit_sip_rejects_fewer_pairs_than_monomial_terms() {
        let wcs = sample_wcs();
        // order 4 needs 12 monomial terms, 11 pairs is short.
        // all pairs at the tangent point project to the same xi/eta; we want
        // the pairs.len() < n_terms check to fire.
        let pairs: Vec<StarPair> = (0..11)
            .map(|i| pair_at_pixel(&wcs, i as f64, 0.0))
            .collect();
        let err = fit_sip(&pairs, &wcs, 4).unwrap_err();
        assert!(err.to_string().contains("need at least 12 pairs"));
    }

    #[test]
    fn fit_sip_happy_path_recovers_injected_distortion() {
        // Inject a known order-2 distortion: px_x += a20 * u^2
        let wcs = sample_wcs();
        let a20 = 1e-6_f64;

        // Scatter sample points across the field
        let grid: Vec<(f64, f64)> = (-4..=4)
            .flat_map(|u| (-4..=4).map(move |v| (u as f64 * 20.0, v as f64 * 20.0)))
            .collect();

        // For each undistorted pixel, compute its catalog sky position, then
        // create the pair with the distorted pixel but the same sky position.
        // That's what a real detection looks like when the optics warp positions.
        let pairs: Vec<StarPair> = grid.iter().map(|&(u, v)| {
            let ideal_x = u;
            let ideal_y = v;
            let (ra, dec) = wcs.pixel_to_sky(ideal_x, ideal_y);
            let distorted_x = ideal_x + a20 * u * u;
            StarPair {
                px_x: distorted_x,
                px_y: ideal_y,
                ra_deg: ra,
                dec_deg: dec,
                votes: 4,
                snr: 50.0,
            }
        }).collect();

        let sip = fit_sip(&pairs, &wcs, 2).unwrap();
        assert_eq!(sip.a_order, 2);
        assert_eq!(sip.b_order, 2);
        assert!(sip.n_stars > 0);

        // A_{2,0} should be close to a20
        let recovered = sip.a_coeffs.get(&(2, 0)).copied().unwrap_or(0.0);
        assert!(
            (recovered - a20).abs() < 1e-10,
            "A_2_0 recovered as {:e}, expected {:e}",
            recovered, a20
        );

        // SIP RMS should be much smaller than linear RMS
        assert!(
            sip.rms_px < sip.linear_rms_px,
            "rms {} should be < linear_rms {}",
            sip.rms_px, sip.linear_rms_px
        );
    }
}
