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
    let (catalog_stars, _mag) = prep_catalog_for_refine(initial, catalog, epoch, params.target_stars);
    let mut wcs = initial.clone();
    let mut final_pairs = Vec::new();
    log::debug!("refine_wcs: initial rms={:.4}px n_stars={}", wcs.rms_px, wcs.n_stars);

    for iter in 0..params.max_iterations {
        match run_refine_iteration(&wcs, stars, &catalog_stars, params, iter)? {
            RefineIterStep::Continue { new_wcs, new_pairs, converged } => {
                final_pairs = new_pairs;
                wcs = new_wcs;
                if converged { break; }
            }
            RefineIterStep::Break => break,
        }
    }

    ensure_pairs_nonempty(&final_pairs, params.max_iterations)?;
    log::debug!("refine_wcs: end-of-loop wcs.rms={:.4}px final_pairs={}", wcs.rms_px, final_pairs.len());
    Ok(RefineResult { wcs, pairs: final_pairs })
}

fn prep_catalog_for_refine(
    initial: &WcsSolution,
    catalog: &Catalog,
    epoch: JulianDate,
    target_stars: usize,
) -> (Vec<CatalogEntry>, f64) {
    let radius_deg = field_radius_deg(initial);
    let (catalog_stars, mag_limit) = search_optimal_magnitude(
        catalog, initial.crval1, initial.crval2, radius_deg, epoch, target_stars,
    );
    log::debug!(
        "refine_wcs: {} catalog stars to mag {:.2} within {:.3} deg (target {})",
        catalog_stars.len(), mag_limit, radius_deg, target_stars,
    );
    (catalog_stars, mag_limit)
}

enum RefineIterStep {
    Continue { new_wcs: WcsSolution, new_pairs: Vec<StarPair>, converged: bool },
    Break,
}

fn run_refine_iteration(
    wcs: &WcsSolution,
    stars: &[DetectedStar],
    cats: &[CatalogEntry],
    params: &RefineParams,
    iter: usize,
) -> Result<RefineIterStep> {
    let m = match_by_projection(wcs, stars, cats, params.tolerance_px);
    if m.pairs.len() < 3 {
        log::debug!("refine iter {}: only {} pairs, stopping", iter, m.pairs.len());
        if iter == 0 {
            return Err(anyhow::anyhow!(
                "refine could not extend match: only {} projection pairs within tolerance",
                m.pairs.len(),
            ));
        }
        return Ok(RefineIterStep::Break);
    }
    let refined = fit_direct(&m.pairs, &m.weights, wcs)?;
    let delta = libm::fabs(wcs.rms_px - refined.rms_px);
    log::debug!(
        "refine iter {}: {} stars, rms={:.4} px, wrms={:.4} px (delta={:.4})",
        iter, refined.n_stars, refined.rms_px, refined.weighted_rms_px, delta,
    );
    let converged = delta < params.rms_delta_threshold && refined.n_stars == wcs.n_stars;
    Ok(RefineIterStep::Continue { new_wcs: refined, new_pairs: m.pairs, converged })
}

fn ensure_pairs_nonempty(pairs: &[StarPair], max_iterations: usize) -> Result<()> {
    ensure!(
        !pairs.is_empty(),
        "refine produced no matched pairs over {} iterations",
        max_iterations,
    );
    Ok(())
}

pub(crate) fn fit_direct(
    pairs: &[StarPair],
    weights: &[f64],
    prior: &WcsSolution,
) -> Result<WcsSolution> {
    ensure!(pairs.len() >= 3, "need at least 3 matched pairs, got {}", pairs.len());
    let crpix1 = prior.width as f64 / 2.0;
    let crpix2 = prior.height as f64 / 2.0;
    let (projected, w_vec) = project_pairs_to_tangent(
        pairs, weights, prior.crval1, prior.crval2, crpix1, crpix2,
    );
    ensure!(projected.len() >= 3, "too few stars projected onto tangent plane");

    let mut mask = vec![true; projected.len()];
    let mut cd = fit_cd_weighted(&projected, &mask, &w_vec)?;
    sigma_clip_loop_weighted(&projected, &mut mask, &mut cd, &w_vec)?;
    cd = fit_cd_weighted(&projected, &mask, &w_vec)?;

    let (crval1, crval2) = tan_deproject_star(
        cd.xi0 * DEG_TO_RAD, cd.eta0 * DEG_TO_RAD, prior.crval1, prior.crval2,
    );
    let (kept_pairs, kept_weights) = extract_kept_pairs(pairs, &mask, &w_vec);
    Ok(assemble_wcs_solution(WcsAssembly {
        prior, crval1, crval2, cd: &cd,
        kept_pairs: &kept_pairs, kept_weights: &kept_weights,
        crpix1, crpix2,
    }))
}

pub(crate) type ProjectedPair = (f64, f64, f64, f64);

pub(crate) fn project_pairs_to_tangent(
    pairs: &[StarPair],
    weights: &[f64],
    crval_ra: f64,
    crval_dec: f64,
    crpix1: f64,
    crpix2: f64,
) -> (Vec<ProjectedPair>, Vec<f64>) {
    let mut projected = Vec::with_capacity(pairs.len());
    let mut w_vec = Vec::with_capacity(pairs.len());
    for (p, &wt) in pairs.iter().zip(weights.iter()) {
        if let Some((xi, eta)) = tan_project_star(p.ra_deg, p.dec_deg, crval_ra, crval_dec) {
            projected.push((p.px_x - crpix1, p.px_y - crpix2, xi * RAD_TO_DEG, eta * RAD_TO_DEG));
            w_vec.push(wt);
        }
    }
    (projected, w_vec)
}

pub(crate) fn extract_kept_pairs<'a>(
    pairs: &'a [StarPair],
    mask: &[bool],
    weights: &[f64],
) -> (Vec<&'a StarPair>, Vec<f64>) {
    let mut kept_pairs = Vec::new();
    let mut kept_weights = Vec::new();
    for (i, (p, &m)) in pairs.iter().zip(mask.iter()).enumerate() {
        if !m { continue; }
        if i < weights.len() {
            kept_pairs.push(p);
            kept_weights.push(weights[i]);
        }
    }
    (kept_pairs, kept_weights)
}

struct WcsAssembly<'a> {
    prior: &'a WcsSolution,
    crval1: f64,
    crval2: f64,
    cd: &'a super::cd_matrix::CdParams,
    kept_pairs: &'a [&'a StarPair],
    kept_weights: &'a [f64],
    crpix1: f64,
    crpix2: f64,
}

fn assemble_wcs_solution(a: WcsAssembly<'_>) -> WcsSolution {
    let wcs_params = WcsParams {
        crpix1: a.crpix1, crpix2: a.crpix2,
        crval_ra: a.crval1, crval_dec: a.crval2,
        cd1_1: a.cd.cd1_1, cd1_2: a.cd.cd1_2, cd2_1: a.cd.cd2_1, cd2_2: a.cd.cd2_2,
    };
    let residuals = compute_residuals(a.kept_pairs, wcs_params);
    let rms_px = rms_error(&residuals);
    let weighted_rms_px = weighted_rms_error(&residuals, a.kept_weights);
    WcsSolution {
        crpix1: a.crpix1, crpix2: a.crpix2,
        crval1: a.crval1, crval2: a.crval2,
        cd1_1: a.cd.cd1_1, cd1_2: a.cd.cd1_2, cd2_1: a.cd.cd2_1, cd2_2: a.cd.cd2_2,
        width: a.prior.width, height: a.prior.height,
        focal_mm: a.prior.focal_mm, pixel_um: a.prior.pixel_um,
        n_stars: residuals.len(),
        rms_px, weighted_rms_px, residuals,
    }
}

fn field_radius_deg(wcs: &WcsSolution) -> f64 {
    let (sx, sy) = wcs.scale_arcsec();
    let diag = libm::sqrt(
        (wcs.width as f64 * sx) * (wcs.width as f64 * sx)
            + (wcs.height as f64 * sy) * (wcs.height as f64 * sy),
    );
    diag / 3600.0 / 2.0 * 1.2
}

pub(crate) struct CatalogEntry {
    pub ra_deg: f64,
    pub dec_deg: f64,
}

pub(crate) fn cone_search_at_mag(
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

pub(crate) fn search_optimal_magnitude(
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

pub(crate) struct MatchResult {
    pub pairs: Vec<StarPair>,
    pub weights: Vec<f64>,
}

pub(crate) fn match_by_projection(
    wcs: &WcsSolution,
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    tolerance_px: f64,
) -> MatchResult {
    let tol_sq = tolerance_px * tolerance_px;
    let best_det_for_cat = assign_cats_to_dets(wcs, stars, catalog_stars, tol_sq);
    let best_cat_for_det = dedupe_best_cat_per_det(&best_det_for_cat, stars.len());
    build_match_pairs(stars, catalog_stars, &best_cat_for_det)
}

pub(crate) fn predict_cat_in_image(wcs: &WcsSolution, cat: &CatalogEntry) -> Option<(f64, f64)> {
    let (px, py) = wcs.sky_to_pixel(cat.ra_deg, cat.dec_deg)?;
    let w = wcs.width as f64;
    let h = wcs.height as f64;
    if px < 0.0 || px >= w || py < 0.0 || py >= h {
        return None;
    }
    Some((px, py))
}

pub(crate) fn nearest_detection(
    stars: &[DetectedStar],
    pred_x: f64,
    pred_y: f64,
    tol_sq: f64,
) -> Option<(usize, f64)> {
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
    best_si.map(|si| (si, best_dist_sq))
}

pub(crate) fn assign_cats_to_dets(
    wcs: &WcsSolution,
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    tol_sq: f64,
) -> Vec<Option<(usize, f64)>> {
    let mut best_det_for_cat = vec![None; catalog_stars.len()];
    for (ci, cat) in catalog_stars.iter().enumerate() {
        let Some((pred_x, pred_y)) = predict_cat_in_image(wcs, cat) else { continue };
        best_det_for_cat[ci] = nearest_detection(stars, pred_x, pred_y, tol_sq);
    }
    best_det_for_cat
}

pub(crate) fn dedupe_best_cat_per_det(
    best_det_for_cat: &[Option<(usize, f64)>],
    n_dets: usize,
) -> Vec<Option<(usize, f64)>> {
    let mut best_cat_for_det: Vec<Option<(usize, f64)>> = vec![None; n_dets];
    for (ci, entry) in best_det_for_cat.iter().enumerate() {
        let Some((si, d2)) = entry else { continue };
        let slot = &mut best_cat_for_det[*si];
        if slot.is_none() || d2 < &slot.unwrap().1 {
            *slot = Some((ci, *d2));
        }
    }
    best_cat_for_det
}

fn build_match_pairs(
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    best_cat_for_det: &[Option<(usize, f64)>],
) -> MatchResult {
    let mut pairs = Vec::new();
    let mut weights = Vec::new();
    let mut match_dists = Vec::new();
    for (si, slot) in best_cat_for_det.iter().enumerate() {
        let Some((ci, d2)) = slot else { continue };
        let cat = &catalog_stars[*ci];
        let det = &stars[si];
        pairs.push(StarPair {
            px_x: det.x, px_y: det.y,
            ra_deg: cat.ra_deg, dec_deg: cat.dec_deg,
            votes: 1, snr: det.snr,
        });
        weights.push(det.snr * det.snr);
        match_dists.push(libm::sqrt(*d2));
    }
    log_match_stats(&match_dists, pairs.len());
    MatchResult { pairs, weights }
}

fn log_match_stats(match_dists: &[f64], n_pairs: usize) {
    let mut sorted = match_dists.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if sorted.is_empty() { 0.0 } else { sorted[sorted.len() / 2] };
    log::debug!(
        "  match_by_projection: {} final pairs, median dist={:.3} px",
        n_pairs, median,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    use celestial_catalog::query::Catalog;

    use crate::detect::DetectedStar;
    use crate::match_field::test_catalog::{build, SynthStar};

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

    fn wcs_centered_at(ra: f64, dec: f64, scale_arcsec: f64, w: usize, h: usize) -> WcsSolution {
        let s = scale_arcsec / 3600.0;
        WcsSolution {
            crpix1: w as f64 / 2.0,
            crpix2: h as f64 / 2.0,
            crval1: ra,
            crval2: dec,
            cd1_1: -s,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: s,
            width: w,
            height: h,
            focal_mm: None,
            pixel_um: None,
            n_stars: 0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
            residuals: Vec::new(),
        }
    }

    fn det(x: f64, y: f64, snr: f64) -> DetectedStar {
        DetectedStar {
            x, y,
            flux: 1000.0,
            snr,
            saturated: false,
            saturated_count: 0,
            background: 0.0,
        }
    }

    fn cat(ra: f64, dec: f64) -> CatalogEntry {
        CatalogEntry { ra_deg: ra, dec_deg: dec }
    }

    fn synth_catalog(stars: Vec<SynthStar>) -> (Catalog, tempfile::NamedTempFile) {
        let file = build(4, &stars);
        let c = Catalog::open(file.path()).unwrap();
        (c, file)
    }

    #[test]
    fn field_radius_matches_expected_formula() {
        let wcs = sample_wcs(1000, 1000, 0.001);
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
        let ratio = big / small;
        assert!((ratio - 4.0).abs() < 1e-9);
    }

    #[test]
    fn field_radius_scales_with_plate_scale() {
        let fine = field_radius_deg(&sample_wcs(1000, 1000, 0.001));
        let coarse = field_radius_deg(&sample_wcs(1000, 1000, 0.002));
        assert!((coarse / fine - 2.0).abs() < 1e-9);
    }

    #[test]
    fn refine_params_defaults_match_documented_values() {
        let p = RefineParams::default();
        assert_eq!(p.tolerance_px, 3.0);
        assert_eq!(p.max_iterations, 5);
        assert_eq!(p.rms_delta_threshold, 0.01);
        assert_eq!(p.target_stars, 2000);
    }

    #[test]
    fn cone_search_at_mag_filters_by_magnitude() {
        // 3 stars: mag 8, 12, 18. Search with max_mag=10 should return only the
        // bright one; max_mag=20 should return all three.
        let stars = vec![
            SynthStar { source_id: 1, ra: 0.0, dec: 0.0, mag: 8.0 },
            SynthStar { source_id: 2, ra: 0.01, dec: 0.0, mag: 12.0 },
            SynthStar { source_id: 3, ra: 0.0, dec: 0.01, mag: 18.0 },
        ];
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let bright_only = cone_search_at_mag(&c, 0.0, 0.0, 1.0, 10.0, epoch);
        assert_eq!(bright_only.len(), 1);
        assert!((bright_only[0].ra_deg - 0.0).abs() < 1e-6);

        let all = cone_search_at_mag(&c, 0.0, 0.0, 1.0, 20.0, epoch);
        assert_eq!(all.len(), 3);
    }

    #[test]
    fn cone_search_at_mag_returns_empty_when_no_stars_in_cone() {
        let stars = vec![SynthStar { source_id: 1, ra: 10.0, dec: 0.0, mag: 8.0 }];
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let out = cone_search_at_mag(&c, 200.0, 0.0, 0.5, 20.0, epoch);
        assert!(out.is_empty());
    }

    #[test]
    fn search_optimal_magnitude_converges_when_target_reachable() {
        // 50 stars spread across magnitudes 8-18. Target 10 → search should
        // settle on a magnitude that yields ≥10 stars.
        let stars: Vec<_> = (0..50)
            .map(|i| SynthStar {
                source_id: i,
                ra: (i as f64) * 0.005,
                dec: 0.0,
                mag: 8.0 + (i as f32) * 0.2,
            })
            .collect();
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let (results, mag) = search_optimal_magnitude(&c, 0.0, 0.0, 1.0, epoch, 10);
        assert!(results.len() >= 10, "got {} stars at mag {mag}", results.len());
        assert!((8.0..=20.0).contains(&mag), "mag {mag} outside expected range");
    }

    #[test]
    fn search_optimal_magnitude_returns_high_cap_when_target_unreachable() {
        // 2 stars in a target-of-100 search → never reaches target; mag should
        // hit the hi cap (20.0).
        let stars = vec![
            SynthStar { source_id: 1, ra: 0.0, dec: 0.0, mag: 8.0 },
            SynthStar { source_id: 2, ra: 0.01, dec: 0.0, mag: 10.0 },
        ];
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let (results, mag) = search_optimal_magnitude(&c, 0.0, 0.0, 1.0, epoch, 100);
        assert_eq!(results.len(), 2);
        assert!(mag > 19.0, "expected mag near 20.0 cap, got {mag}");
    }

    #[test]
    fn match_by_projection_skips_catalog_stars_outside_image() {
        // Wcs centered at (180, 0); place one catalog star far away from the
        // tangent center so its predicted pixel lands outside the 100×100 image,
        // and one near the center that maps inside. Detection is exactly under
        // the near one. Only the in-bounds star should produce a pair.
        let wcs = wcs_centered_at(180.0, 0.0, 3.6, 100, 100);
        let stars = vec![det(50.0, 50.0, 20.0)];
        let cats = vec![
            cat(180.0, 0.0),  // near image center
            cat(0.0, 80.0),    // anti-meridian + high dec → far outside
        ];

        let result = match_by_projection(&wcs, &stars, &cats, 5.0);
        assert_eq!(result.pairs.len(), 1, "only the in-bounds star should match");
        assert_eq!(result.pairs[0].ra_deg, 180.0);
    }

    #[test]
    fn match_by_projection_no_detection_within_tolerance_yields_no_pair() {
        // Catalog star projects to image center, but the only detection is
        // 50 px away — outside tolerance of 5.
        let wcs = wcs_centered_at(180.0, 0.0, 3.6, 100, 100);
        let stars = vec![det(0.0, 0.0, 20.0)];
        let cats = vec![cat(180.0, 0.0)];

        let result = match_by_projection(&wcs, &stars, &cats, 5.0);
        assert!(result.pairs.is_empty());
        assert!(result.weights.is_empty());
    }

    #[test]
    fn match_by_projection_dedupes_multiple_cats_to_closest_per_det() {
        // Two catalog stars both project very near (50, 50). One detection at
        // exactly (50, 50). Each cat-loop assigns the detection to itself, but
        // the second pass keeps only the closest cat. Result: 1 pair.
        let wcs = wcs_centered_at(180.0, 0.0, 3.6, 100, 100);
        let stars = vec![det(50.0, 50.0, 20.0)];
        // Both catalog stars project very close to the image center.
        // Pick two slightly different RA values so they're distinguishable.
        let cats = vec![
            cat(180.0, 0.0),       // projects to exactly (50, 50)
            cat(180.0005, 0.0),    // projects to (50 - small_offset, 50)
        ];

        let result = match_by_projection(&wcs, &stars, &cats, 5.0);
        assert_eq!(result.pairs.len(), 1, "single det shouldn't pair to two cats");
        // The closer cat is the one at RA=180 (zero offset from center).
        assert!((result.pairs[0].ra_deg - 180.0).abs() < 1e-6);
    }

    #[test]
    fn match_by_projection_weight_is_snr_squared() {
        let wcs = wcs_centered_at(180.0, 0.0, 3.6, 100, 100);
        let stars = vec![det(50.0, 50.0, 25.0)];
        let cats = vec![cat(180.0, 0.0)];

        let result = match_by_projection(&wcs, &stars, &cats, 5.0);
        assert_eq!(result.pairs.len(), 1);
        assert_eq!(result.weights[0], 25.0 * 25.0);
        assert_eq!(result.pairs[0].snr, 25.0);
    }

    fn pairs_on_tangent_plane(
        center_ra: f64,
        center_dec: f64,
        scale_arcsec: f64,
        crpix: f64,
    ) -> Vec<StarPair> {
        let s_deg = scale_arcsec / 3600.0;
        let mut pairs = Vec::new();
        for &(dx, dy) in &[
            (-200.0, -200.0), (200.0, -200.0), (-200.0, 200.0), (200.0, 200.0),
            (-100.0, 50.0), (100.0, -50.0), (0.0, 0.0), (150.0, 150.0),
        ] {
            let px_x = crpix + dx;
            let px_y = crpix + dy;
            let xi_deg = -s_deg * dx;
            let eta_deg = s_deg * dy;
            let (ra, dec) = tan_deproject_star(
                xi_deg * DEG_TO_RAD, eta_deg * DEG_TO_RAD, center_ra, center_dec,
            );
            pairs.push(StarPair {
                px_x, px_y, ra_deg: ra, dec_deg: dec, votes: 4, snr: 20.0,
            });
        }
        pairs
    }

    #[test]
    fn fit_direct_recovers_crval_from_clean_pairs() {
        let center_ra = 180.0;
        let center_dec = 30.0;
        let scale = 1.5;
        let prior = wcs_centered_at(center_ra, center_dec, scale, 1024, 1024);
        let pairs = pairs_on_tangent_plane(center_ra, center_dec, scale, 512.0);
        let weights: Vec<f64> = pairs.iter().map(|p| p.snr * p.snr).collect();

        let wcs = fit_direct(&pairs, &weights, &prior)
            .expect("clean pairs must produce a WCS");
        assert!((wcs.crval1 - center_ra).abs() < 0.1);
        assert!((wcs.crval2 - center_dec).abs() < 0.1);
        assert!(wcs.n_stars >= 3, "expected ≥3 inliers, got {}", wcs.n_stars);
        assert!(wcs.rms_px >= 0.0);
    }

    #[test]
    fn fit_direct_errors_on_too_few_pairs() {
        let prior = wcs_centered_at(180.0, 30.0, 1.5, 1024, 1024);
        let pairs = vec![
            StarPair { px_x: 0.0, px_y: 0.0, ra_deg: 180.0, dec_deg: 30.0, votes: 1, snr: 10.0 },
            StarPair { px_x: 10.0, px_y: 10.0, ra_deg: 180.01, dec_deg: 30.0, votes: 1, snr: 10.0 },
        ];
        let weights = vec![100.0, 100.0];
        let err = fit_direct(&pairs, &weights, &prior).expect_err("2 pairs must error");
        assert!(err.to_string().contains("at least 3"));
    }

    #[test]
    fn predict_cat_in_image_returns_none_outside_bounds() {
        let wcs = wcs_centered_at(180.0, 0.0, 3.6, 100, 100);
        assert!(predict_cat_in_image(&wcs, &cat(180.0, 0.0)).is_some());
        // Far from tangent center → projection succeeds but pixel lies outside.
        assert!(predict_cat_in_image(&wcs, &cat(0.0, 80.0)).is_none());
    }

    #[test]
    fn nearest_detection_picks_closest_within_tolerance() {
        let stars = vec![det(50.0, 50.0, 10.0), det(52.0, 50.0, 10.0), det(60.0, 60.0, 10.0)];
        // pred at (51, 50): closest is star 0 (d=1), then star 1 (d=1) — first wins on tie.
        let (idx, d2) = nearest_detection(&stars, 51.0, 50.0, 5.0 * 5.0).expect("must match");
        assert_eq!(idx, 0);
        assert!((d2 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn nearest_detection_returns_none_when_all_outside_tolerance() {
        let stars = vec![det(0.0, 0.0, 10.0), det(100.0, 100.0, 10.0)];
        assert!(nearest_detection(&stars, 50.0, 50.0, 5.0 * 5.0).is_none());
    }

    #[test]
    fn dedupe_best_cat_per_det_keeps_closest_cat() {
        // Cat 0 claims det 0 at d²=4. Cat 1 also claims det 0 at d²=1.
        // After dedupe, det 0's slot should point at cat 1.
        let input = vec![
            Some((0usize, 4.0)),
            Some((0usize, 1.0)),
            None,
        ];
        let out = dedupe_best_cat_per_det(&input, 1);
        assert_eq!(out.len(), 1);
        let (ci, d2) = out[0].expect("det 0 must be claimed");
        assert_eq!(ci, 1);
        assert!((d2 - 1.0).abs() < 1e-9);
    }
}
