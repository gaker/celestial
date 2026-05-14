//! WCS fitting: gnomonic projection, CD matrix, SIP distortion, residual refinement.
//!
//! [`fit_wcs`] fits a TAN (gnomonic) WCS from matched star pairs. [`refine_wcs`]
//! re-matches against the catalog and refits with outlier rejection. [`fit_sip`]
//! adds polynomial distortion on top of the TAN fit.

mod cd_matrix;
mod display;
mod refine;
mod sip;
pub mod solution;

pub use display::SolveDisplay;
pub use refine::{refine_wcs, RefineParams, RefineResult};
pub use sip::{fit_sip, SipSolution};
pub use solution::{StarResidual, WcsSolution};

use anyhow::{ensure, Result};
use celestial_catalog::query::{tan_deproject_star, tan_project_star};
use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};

use crate::match_field::StarPair;

use cd_matrix::{
    WcsParams,
    ransac_fit, fit_cd_weighted, sigma_clip_loop_weighted,
    compute_residuals, rms_error, weighted_rms_error,
};

/// Fits a TAN (gnomonic) WCS from matched star pairs.
///
/// RANSAC inlier selection, then weighted least-squares CD-matrix fit, then
/// sigma-clipping loop, then a final weighted fit. The CRVAL is refined from the
/// initial `crval_ra_deg` / `crval_dec_deg` prior using the tangent-plane offsets.
/// Weights are per-pair SNR squared.
///
/// Usually called internally by [`crate::solve()`]. Use directly when you already have
/// matched `(image pixel, catalog RA/Dec)` pairs from another source.
///
/// # Errors
///
/// Returns an error when fewer than 3 pairs are supplied, the tangent-plane projection
/// fails on too many of them, or RANSAC fails to find a consistent inlier set.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::fit_wcs;
///
/// # let pairs = vec![];
/// let wcs = fit_wcs(&pairs, 4656, 3520, 180.0, 30.0, Some(500.0), Some(3.76))?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn fit_wcs(
    pairs: &[StarPair],
    width: usize,
    height: usize,
    crval_ra_deg: f64,
    crval_dec_deg: f64,
    focal_mm: Option<f64>,
    pixel_um: Option<f64>,
) -> Result<WcsSolution> {
    ensure!(pairs.len() >= 3, "need at least 3 matched pairs, got {}", pairs.len());

    let crpix1 = width as f64 / 2.0;
    let crpix2 = height as f64 / 2.0;

    let projected: Vec<(f64, f64, f64, f64)> = pairs
        .iter()
        .filter_map(|p| {
            tan_project_star(p.ra_deg, p.dec_deg, crval_ra_deg, crval_dec_deg)
                .map(|(xi, eta)| {
                    (p.px_x - crpix1, p.px_y - crpix2, xi * RAD_TO_DEG, eta * RAD_TO_DEG)
                })
        })
        .collect();

    ensure!(projected.len() >= 3, "too few stars projected onto tangent plane");

    let snr_weights: Vec<f64> = pairs
        .iter()
        .filter_map(|p| {
            tan_project_star(p.ra_deg, p.dec_deg, crval_ra_deg, crval_dec_deg)
                .map(|_| p.snr * p.snr)
        })
        .collect();

    let inlier_px = 5.0;
    let mut mask = ransac_fit(&projected, pairs, inlier_px)?;

    let n_ransac = mask.iter().filter(|&&m| m).count();
    log::debug!("fit_wcs ransac: {} inliers within {:.1}px", n_ransac, inlier_px);

    let mut cd = fit_cd_weighted(&projected, &mask, &snr_weights)?;

    sigma_clip_loop_weighted(&projected, &mut mask, &mut cd, &snr_weights)?;

    cd = fit_cd_weighted(&projected, &mask, &snr_weights)?;

    let (crval1, crval2) = tan_deproject_star(
        cd.xi0 * DEG_TO_RAD, cd.eta0 * DEG_TO_RAD,
        crval_ra_deg, crval_dec_deg,
    );

    let mut kept_pairs = Vec::new();
    let mut kept_weights = Vec::new();
    for (i, (p, &m)) in pairs.iter().zip(mask.iter()).enumerate() {
        if !m { continue; }
        if i < snr_weights.len() {
            kept_pairs.push(p);
            kept_weights.push(snr_weights[i]);
        }
    }

    let wcs_params = WcsParams {
        crpix1, crpix2, crval_ra: crval1, crval_dec: crval2,
        cd1_1: cd.cd1_1, cd1_2: cd.cd1_2, cd2_1: cd.cd2_1, cd2_2: cd.cd2_2,
    };
    let residuals = compute_residuals(&kept_pairs, wcs_params);
    let rms_px = rms_error(&residuals);
    let weighted_rms_px = weighted_rms_error(&residuals, &kept_weights);

    Ok(WcsSolution {
        crpix1, crpix2,
        crval1, crval2,
        cd1_1: cd.cd1_1, cd1_2: cd.cd1_2, cd2_1: cd.cd2_1, cd2_2: cd.cd2_2,
        width, height,
        focal_mm, pixel_um,
        n_stars: residuals.len(),
        rms_px,
        weighted_rms_px,
        residuals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_catalog::query::tan_deproject_star;

    struct SynthGrid {
        crval_ra: f64,
        crval_dec: f64,
        cd1_1: f64,
        cd2_2: f64,
        crpix1: f64,
        crpix2: f64,
        grid_step: f64,
        grid_half: i32,
    }

    fn synthesize_pairs(g: &SynthGrid) -> Vec<StarPair> {
        let mut pairs = Vec::new();
        for gy in -g.grid_half..=g.grid_half {
            for gx in -g.grid_half..=g.grid_half {
                let px_x = g.crpix1 + gx as f64 * g.grid_step;
                let px_y = g.crpix2 + gy as f64 * g.grid_step;
                let u = px_x - g.crpix1;
                let v = px_y - g.crpix2;
                let xi_deg = g.cd1_1 * u;
                let eta_deg = g.cd2_2 * v;
                let (ra, dec) = tan_deproject_star(
                    xi_deg * DEG_TO_RAD,
                    eta_deg * DEG_TO_RAD,
                    g.crval_ra,
                    g.crval_dec,
                );
                pairs.push(StarPair {
                    px_x,
                    px_y,
                    ra_deg: ra,
                    dec_deg: dec,
                    votes: 4,
                    snr: 50.0,
                });
            }
        }
        pairs
    }

    #[test]
    fn fit_wcs_recovers_known_cd_matrix_and_crval() {
        let crval_ra = 180.0;
        let crval_dec = 0.0;
        let cd1_1 = 0.001;
        let cd2_2 = 0.001;
        let width = 1000;
        let height = 1000;
        let pairs = synthesize_pairs(&SynthGrid {
            crval_ra,
            crval_dec,
            cd1_1,
            cd2_2,
            crpix1: width as f64 / 2.0,
            crpix2: height as f64 / 2.0,
            grid_step: 20.0,
            grid_half: 5,
        });
        let wcs = fit_wcs(&pairs, width, height, crval_ra, crval_dec, None, None).unwrap();

        // crpix stored as width/2, height/2
        assert_eq!(wcs.crpix1, 500.0);
        assert_eq!(wcs.crpix2, 500.0);
        // Recovered crval must round-trip to the input
        assert!((wcs.crval1 - crval_ra).abs() < 1e-9);
        assert!((wcs.crval2 - crval_dec).abs() < 1e-9);
        // Diagonal CD recovered
        assert!((wcs.cd1_1 - cd1_1).abs() < 1e-12);
        assert!((wcs.cd2_2 - cd2_2).abs() < 1e-12);
        assert!(wcs.cd1_2.abs() < 1e-12);
        assert!(wcs.cd2_1.abs() < 1e-12);
        // RMS should be essentially zero for noise-free synthesized pairs
        assert!(wcs.rms_px < 1e-9);
    }

    #[test]
    fn fit_wcs_preserves_focal_and_pixel_metadata() {
        let pairs = synthesize_pairs(&SynthGrid {
            crval_ra: 180.0,
            crval_dec: 0.0,
            cd1_1: 0.001,
            cd2_2: 0.001,
            crpix1: 500.0,
            crpix2: 500.0,
            grid_step: 20.0,
            grid_half: 5,
        });
        let wcs = fit_wcs(&pairs, 1000, 1000, 180.0, 0.0, Some(500.0), Some(3.76)).unwrap();
        assert_eq!(wcs.focal_mm, Some(500.0));
        assert_eq!(wcs.pixel_um, Some(3.76));
    }

    #[test]
    fn fit_wcs_rejects_fewer_than_three_pairs() {
        let pairs = vec![
            StarPair { px_x: 0.0, px_y: 0.0, ra_deg: 0.0, dec_deg: 0.0, votes: 4, snr: 10.0 },
            StarPair { px_x: 1.0, px_y: 0.0, ra_deg: 0.1, dec_deg: 0.0, votes: 4, snr: 10.0 },
        ];
        let err = fit_wcs(&pairs, 100, 100, 0.0, 0.0, None, None).unwrap_err();
        assert!(err.to_string().contains("at least 3 matched pairs"));
    }
}
