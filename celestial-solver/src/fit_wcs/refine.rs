use anyhow::{ensure, Result};
use celestial_catalog::query::{
    tan_deproject_star, tan_project_star,
    Catalog, ConeSearchParams, cone_search,
};
use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};
use celestial_time::JulianDate;

use crate::detect::DetectedStar;
use crate::match_field::StarPair;

use super::cd_matrix::{
    WcsParams,
    fit_cd_weighted, sigma_clip_loop_weighted,
    compute_residuals, rms_error, weighted_rms_error,
};
use super::solution::WcsSolution;

/// Tuning for [`refine_wcs`].
///
/// The refinement loop re-matches all detected stars against the catalog using the
/// current WCS, keeps pairs within `tolerance_px` of their predicted positions, and
/// refits. It stops when the RMS improvement between iterations falls below
/// `rms_delta_threshold` or after `max_iterations`.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::RefineParams;
///
/// let params = RefineParams {
///     tolerance_px: 2.0,
///     max_iterations: 3,
///     ..RefineParams::default()
/// };
/// ```
pub struct RefineParams {
    /// Maximum detected-to-predicted distance (pixels) to accept a pair.
    pub tolerance_px: f64,
    /// Maximum re-match iterations.
    pub max_iterations: usize,
    /// Convergence threshold on RMS change between iterations.
    pub rms_delta_threshold: f64,
    /// Target number of catalog stars to pull in the cone search.
    pub target_stars: usize,
}

impl Default for RefineParams {
    fn default() -> Self {
        Self {
            tolerance_px: 3.0,
            max_iterations: 5,
            rms_delta_threshold: 0.01,
            target_stars: 2000,
        }
    }
}

/// Output of [`refine_wcs`].
///
/// Contains the refined [`WcsSolution`] and the matched pairs used in the final fit
/// (typically many more than the initial quad-matching produced).
pub struct RefineResult {
    /// Refined WCS.
    pub wcs: WcsSolution,
    /// Final matched pairs after refinement.
    pub pairs: Vec<StarPair>,
}

/// Re-matches detections against the catalog and refits the WCS.
///
/// Takes the initial WCS from [`crate::fit_wcs::fit_wcs`] and uses it to predict
/// catalog positions in pixel space. Matches every detected star to its nearest
/// catalog neighbor within `tolerance_px`, refits, and iterates until converged.
///
/// Usually called internally by [`crate::solve()`]. Use directly to re-refine after
/// manually adjusting the initial WCS.
///
/// # Errors
///
/// Returns an error when the cone search fails or too few matches are found to fit.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::{refine_wcs, RefineParams};
///
/// # let initial = unimplemented!();
/// # let stars = vec![];
/// # let catalog: celestial_catalog::query::Catalog = unimplemented!();
/// # let epoch = celestial_time::JulianDate::new(2_400_000.5, 0.0);
/// let refined = refine_wcs(&initial, &stars, &catalog, epoch, &RefineParams::default())?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn refine_wcs(
    initial: &WcsSolution,
    stars: &[DetectedStar],
    catalog: &Catalog,
    epoch: JulianDate,
    params: &RefineParams,
) -> Result<RefineResult> {
    let radius_deg = field_radius_deg(initial);
    let (catalog_stars, mag_limit) = search_optimal_magnitude(
        catalog, initial.crval1, initial.crval2,
        radius_deg, epoch, params.target_stars,
    );
    log::debug!(
        "refine_wcs: {} catalog stars to mag {:.2} within {:.3} deg (target {})",
        catalog_stars.len(), mag_limit, radius_deg, params.target_stars,
    );

    let mut wcs = initial.clone();
    let mut final_pairs = Vec::new();
    log::debug!(
        "refine_wcs: initial rms={:.4}px n_stars={}",
        wcs.rms_px, wcs.n_stars,
    );

    for iter in 0..params.max_iterations {
        let m = match_by_projection(&wcs, stars, &catalog_stars, params.tolerance_px);
        if m.pairs.len() < 3 {
            log::debug!("refine iter {}: only {} pairs, stopping", iter, m.pairs.len());
            if iter == 0 {
                return Err(anyhow::anyhow!(
                    "refine could not extend match: only {} projection pairs within tolerance",
                    m.pairs.len(),
                ));
            }
            break;
        }

        let refined = fit_direct(&m.pairs, &m.weights, &wcs)?;

        let delta = libm::fabs(wcs.rms_px - refined.rms_px);
        log::debug!(
            "refine iter {}: {} stars, rms={:.4} px, wrms={:.4} px (delta={:.4})",
            iter, refined.n_stars, refined.rms_px, refined.weighted_rms_px, delta,
        );

        let prev_n = wcs.n_stars;
        final_pairs = m.pairs;
        wcs = refined;

        if delta < params.rms_delta_threshold && wcs.n_stars == prev_n {
            break;
        }
    }

    if final_pairs.is_empty() {
        return Err(anyhow::anyhow!(
            "refine produced no matched pairs over {} iterations",
            params.max_iterations,
        ));
    }

    log::debug!(
        "refine_wcs: end-of-loop wcs.rms={:.4}px final_pairs={}",
        wcs.rms_px, final_pairs.len(),
    );
    Ok(RefineResult { wcs, pairs: final_pairs })
}

fn fit_direct(
    pairs: &[StarPair],
    weights: &[f64],
    prior: &WcsSolution,
) -> Result<WcsSolution> {
    ensure!(pairs.len() >= 3, "need at least 3 matched pairs, got {}", pairs.len());

    let crpix1 = prior.width as f64 / 2.0;
    let crpix2 = prior.height as f64 / 2.0;
    let crval_ra_deg = prior.crval1;
    let crval_dec_deg = prior.crval2;

    let mut projected = Vec::with_capacity(pairs.len());
    let mut w_vec = Vec::with_capacity(pairs.len());
    for (p, &wt) in pairs.iter().zip(weights.iter()) {
        if let Some((xi, eta)) = tan_project_star(p.ra_deg, p.dec_deg, crval_ra_deg, crval_dec_deg) {
            projected.push((p.px_x - crpix1, p.px_y - crpix2, xi * RAD_TO_DEG, eta * RAD_TO_DEG));
            w_vec.push(wt);
        }
    }

    ensure!(projected.len() >= 3, "too few stars projected onto tangent plane");

    let mut mask = vec![true; projected.len()];
    let mut cd = fit_cd_weighted(&projected, &mask, &w_vec)?;
    sigma_clip_loop_weighted(&projected, &mut mask, &mut cd, &w_vec)?;
    cd = fit_cd_weighted(&projected, &mask, &w_vec)?;

    let (crval1, crval2) = tan_deproject_star(
        cd.xi0 * DEG_TO_RAD, cd.eta0 * DEG_TO_RAD,
        crval_ra_deg, crval_dec_deg,
    );

    let mut kept_pairs = Vec::new();
    let mut kept_weights = Vec::new();
    for (i, (p, &m)) in pairs.iter().zip(mask.iter()).enumerate() {
        if !m { continue; }
        if i < w_vec.len() {
            kept_pairs.push(p);
            kept_weights.push(w_vec[i]);
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
        width: prior.width, height: prior.height,
        focal_mm: prior.focal_mm, pixel_um: prior.pixel_um,
        n_stars: residuals.len(),
        rms_px,
        weighted_rms_px,
        residuals,
    })
}

fn field_radius_deg(wcs: &WcsSolution) -> f64 {
    let (sx, sy) = wcs.scale_arcsec();
    let diag = libm::sqrt(
        (wcs.width as f64 * sx) * (wcs.width as f64 * sx)
            + (wcs.height as f64 * sy) * (wcs.height as f64 * sy),
    );
    diag / 3600.0 / 2.0 * 1.2
}

struct CatalogEntry {
    ra_deg: f64,
    dec_deg: f64,
}

fn cone_search_at_mag(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    max_mag: f64,
    epoch: JulianDate,
) -> Vec<CatalogEntry> {
    let params = ConeSearchParams {
        ra_deg,
        dec_deg,
        radius_deg,
        max_mag: Some(max_mag),
        max_results: None,
        epoch: Some(epoch),
    };
    cone_search(catalog, &params)
        .into_iter()
        .map(|r| CatalogEntry { ra_deg: r.ra_deg, dec_deg: r.dec_deg })
        .collect()
}

fn search_optimal_magnitude(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    epoch: JulianDate,
    target: usize,
) -> (Vec<CatalogEntry>, f64) {
    let mut lo = 8.0_f64;
    let mut hi = 20.0_f64;
    let max_iters = 8;

    for _ in 0..max_iters {
        let mid = (lo + hi) / 2.0;
        let n = cone_search_at_mag(catalog, ra_deg, dec_deg, radius_deg, mid, epoch).len();
        log::debug!("  mag search: m={:.2}, {} stars", mid, n);
        if n < target {
            lo = mid;
        } else {
            hi = mid;
        }
        if libm::fabs(hi - lo) < 0.1 {
            break;
        }
    }

    let mag = hi;
    let stars = cone_search_at_mag(catalog, ra_deg, dec_deg, radius_deg, mag, epoch);
    (stars, mag)
}

struct MatchResult {
    pairs: Vec<StarPair>,
    weights: Vec<f64>,
}

fn match_by_projection(
    wcs: &WcsSolution,
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    tolerance_px: f64,
) -> MatchResult {
    let tol_sq = tolerance_px * tolerance_px;
    let w = wcs.width as f64;
    let h = wcs.height as f64;

    let mut n_in_image = 0usize;
    let mut n_matched = 0usize;
    let mut best_det_for_cat: Vec<Option<(usize, f64)>> = vec![None; catalog_stars.len()];

    for (ci, cat) in catalog_stars.iter().enumerate() {
        let Some((pred_x, pred_y)) = wcs.sky_to_pixel(cat.ra_deg, cat.dec_deg) else {
            continue;
        };
        if pred_x < 0.0 || pred_x >= w || pred_y < 0.0 || pred_y >= h {
            continue;
        }
        n_in_image += 1;

        let mut best_dist_sq = tol_sq;
        let mut best_si = None;
        for (si, s) in stars.iter().enumerate() {
            let dx = s.x - pred_x;
            let dy = s.y - pred_y;
            let d2 = dx * dx + dy * dy;
            if d2 < best_dist_sq {
                best_dist_sq = d2;
                best_si = Some(si);
            }
        }
        if let Some(si) = best_si {
            best_det_for_cat[ci] = Some((si, best_dist_sq));
            n_matched += 1;
        }
    }

    let mut best_cat_for_det: Vec<Option<(usize, f64)>> = vec![None; stars.len()];
    for (ci, entry) in best_det_for_cat.iter().enumerate() {
        if let Some((si, d2)) = entry {
            let slot = &mut best_cat_for_det[*si];
            if slot.is_none() || d2 < &slot.unwrap().1 {
                *slot = Some((ci, *d2));
            }
        }
    }

    let mut pairs = Vec::new();
    let mut weights = Vec::new();
    let mut match_dists = Vec::new();
    for (si, slot) in best_cat_for_det.iter().enumerate() {
        if let Some((ci, d2)) = slot {
            let cat = &catalog_stars[*ci];
            let det = &stars[si];
            pairs.push(StarPair {
                px_x: det.x,
                px_y: det.y,
                ra_deg: cat.ra_deg,
                dec_deg: cat.dec_deg,
                votes: 1,
                snr: det.snr,
            });
            weights.push(det.snr * det.snr);
            match_dists.push(libm::sqrt(*d2));
        }
    }

    match_dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if match_dists.is_empty() {
        0.0
    } else {
        match_dists[match_dists.len() / 2]
    };
    log::debug!(
        "  match_by_projection: {} in image, {} matched, {} final pairs, median dist={:.3} px",
        n_in_image, n_matched, pairs.len(), median,
    );

    MatchResult { pairs, weights }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_wcs(width: usize, height: usize, cd_diag: f64) -> WcsSolution {
        WcsSolution {
            crpix1: width as f64 / 2.0,
            crpix2: height as f64 / 2.0,
            crval1: 180.0,
            crval2: 0.0,
            cd1_1: cd_diag,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: cd_diag,
            width,
            height,
            focal_mm: None,
            pixel_um: None,
            n_stars: 0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
            residuals: Vec::new(),
        }
    }

    #[test]
    fn field_radius_matches_expected_formula() {
        // 1000x1000 image, CD=0.001 deg/px = 3.6 arcsec/px
        let wcs = sample_wcs(1000, 1000, 0.001);
        // diag pixels sqrt(1e6 + 1e6)*3.6" = sqrt(2)*1000*3.6" arcsec
        // formula: diag / 3600 / 2 * 1.2 in degrees
        let expected_diag_arcsec = libm::sqrt(
            (1000.0 * 3.6) * (1000.0 * 3.6) + (1000.0 * 3.6) * (1000.0 * 3.6),
        );
        let expected = expected_diag_arcsec / 3600.0 / 2.0 * 1.2;
        let radius = field_radius_deg(&wcs);
        assert!((radius - expected).abs() < 1e-9, "got {}, expected {}", radius, expected);
    }

    #[test]
    fn field_radius_scales_with_image_size() {
        let small = field_radius_deg(&sample_wcs(500, 500, 0.001));
        let big = field_radius_deg(&sample_wcs(2000, 2000, 0.001));
        assert!(big > small);
        // doubling image size and keeping scale fixed doubles the radius
        let ratio = big / small;
        assert!((ratio - 4.0).abs() < 1e-9);
    }

    #[test]
    fn field_radius_scales_with_plate_scale() {
        let fine = field_radius_deg(&sample_wcs(1000, 1000, 0.001));
        let coarse = field_radius_deg(&sample_wcs(1000, 1000, 0.002));
        assert!((coarse / fine - 2.0).abs() < 1e-9);
    }
}
