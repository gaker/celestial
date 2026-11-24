// Pairs each catalog star with the nearest detection to where a WCS places it.

use celestial_wcs::builder::Wcs;

use crate::detect::DetectedStar;
use crate::match_field::StarPair;
use crate::pixel_frame::sky_to_detection;

use super::super::solution::WcsSolution;
use super::catalog::CatalogEntry;

pub(super) struct MatchResult {
    pub(super) pairs: Vec<StarPair>,
    pub(super) weights: Vec<f64>,
    // The detection and catalog index of each pair.
    pub(super) links: Vec<(usize, usize)>,
}

// Pairs each catalog star that `predictor` places inside the image of `wcs` with the
// nearest detection within `tolerance_px` of it, keeping the nearest catalog star for
// each detection.
pub(super) fn match_by_projection(
    predictor: &Wcs,
    wcs: &WcsSolution,
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    tolerance_px: f64,
) -> MatchResult {
    let tol_sq = tolerance_px * tolerance_px;
    let image = (wcs.width as f64, wcs.height as f64);
    let best_det_for_cat = assign_cats_to_dets(predictor, image, stars, catalog_stars, tol_sq);
    let best_cat_for_det = dedupe_best_cat_per_det(&best_det_for_cat, stars.len());
    build_match_pairs(stars, catalog_stars, &best_cat_for_det)
}

// The detection pixel of `cat` under `predictor`, when it is inside an image of
// `image` (width, height) pixels.
fn predict_cat_in_image(
    predictor: &Wcs,
    (w, h): (f64, f64),
    cat: &CatalogEntry,
) -> Option<(f64, f64)> {
    let (px, py) = sky_to_detection(predictor, cat.ra_deg, cat.dec_deg).ok()?;
    if px < 0.0 || px >= w || py < 0.0 || py >= h {
        return None;
    }
    Some((px, py))
}

fn nearest_detection(
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

fn assign_cats_to_dets(
    predictor: &Wcs,
    image: (f64, f64),
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    tol_sq: f64,
) -> Vec<Option<(usize, f64)>> {
    let mut best_det_for_cat = vec![None; catalog_stars.len()];
    for (ci, cat) in catalog_stars.iter().enumerate() {
        let Some((pred_x, pred_y)) = predict_cat_in_image(predictor, image, cat) else {
            continue;
        };
        best_det_for_cat[ci] = nearest_detection(stars, pred_x, pred_y, tol_sq);
    }
    best_det_for_cat
}

fn dedupe_best_cat_per_det(
    best_det_for_cat: &[Option<(usize, f64)>],
    n_dets: usize,
) -> Vec<Option<(usize, f64)>> {
    let mut best_cat_for_det: Vec<Option<(usize, f64)>> = vec![None; n_dets];
    for (ci, entry) in best_det_for_cat.iter().enumerate() {
        let Some((si, d2)) = entry else { continue };
        let slot = &mut best_cat_for_det[*si];
        let closer = match slot {
            None => true,
            Some((_, best)) => d2 < best,
        };
        if closer {
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
    let mut links = Vec::new();
    let mut match_dists = Vec::new();
    for (si, slot) in best_cat_for_det.iter().enumerate() {
        let Some((ci, d2)) = slot else { continue };
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
        weights.push(super::super::pair_weight(det.snr));
        links.push((si, *ci));
        match_dists.push(libm::sqrt(*d2));
    }
    log_match_stats(&match_dists, pairs.len());
    MatchResult {
        pairs,
        weights,
        links,
    }
}

fn log_match_stats(match_dists: &[f64], n_pairs: usize) {
    let mut sorted = match_dists.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = if sorted.is_empty() {
        0.0
    } else {
        sorted[sorted.len() / 2]
    };
    log::debug!(
        "  match_by_projection: {} final pairs, median dist={:.3} px",
        n_pairs,
        median,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::super::fits_wcs::build_fits_wcs;
    use super::super::test_support::det;

    fn tan(wcs: &WcsSolution) -> Wcs {
        build_fits_wcs(wcs, None).unwrap()
    }

    fn cat(ra: f64, dec: f64) -> CatalogEntry {
        CatalogEntry {
            ra_deg: ra,
            dec_deg: dec,
        }
    }

    #[test]
    fn match_by_projection_skips_catalog_stars_outside_image() {
        // Wcs centered at (180, 0); place one catalog star far away from the
        // tangent center so its predicted pixel lands outside the 100×100 image,
        // and one near the center that maps inside. Detection is exactly under
        // the near one. Only the in-bounds star should produce a pair.
        let wcs = WcsSolution::north_up((180.0, 0.0), 3.6, (100, 100));
        let stars = vec![det(50.0, 50.0, 20.0)];
        let cats = vec![
            cat(180.0, 0.0), // near image center
            cat(0.0, 80.0),  // anti-meridian + high dec → far outside
        ];

        let result = match_by_projection(&tan(&wcs), &wcs, &stars, &cats, 5.0);
        assert_eq!(
            result.pairs.len(),
            1,
            "only the in-bounds star should match"
        );
        assert_eq!(result.pairs[0].ra_deg, 180.0);
    }

    #[test]
    fn match_by_projection_no_detection_within_tolerance_yields_no_pair() {
        // Catalog star projects to image center, but the only detection is
        // 50 px away, outside the tolerance of 5.
        let wcs = WcsSolution::north_up((180.0, 0.0), 3.6, (100, 100));
        let stars = vec![det(0.0, 0.0, 20.0)];
        let cats = vec![cat(180.0, 0.0)];

        let result = match_by_projection(&tan(&wcs), &wcs, &stars, &cats, 5.0);
        assert!(result.pairs.is_empty());
        assert!(result.weights.is_empty());
    }

    #[test]
    fn match_by_projection_dedupes_multiple_cats_to_closest_per_det() {
        // Two catalog stars both project very near (50, 50). One detection at
        // exactly (50, 50). Each cat-loop assigns the detection to itself, but
        // the second pass keeps only the closest cat. Result: 1 pair.
        let wcs = WcsSolution::north_up((180.0, 0.0), 3.6, (100, 100));
        let stars = vec![det(50.0, 50.0, 20.0)];
        // Both catalog stars project very close to the image center.
        // Pick two slightly different RA values so they're distinguishable.
        let cats = vec![
            cat(180.0, 0.0),    // projects to exactly (50, 50)
            cat(180.0005, 0.0), // projects to (50 - small_offset, 50)
        ];

        let result = match_by_projection(&tan(&wcs), &wcs, &stars, &cats, 5.0);
        assert_eq!(
            result.pairs.len(),
            1,
            "single det shouldn't pair to two cats"
        );
        // The closer cat is the one at RA=180 (zero offset from center).
        assert_eq!(result.pairs[0].ra_deg, 180.0);
    }

    #[test]
    fn match_by_projection_weight_is_the_inverse_centroid_variance() {
        let wcs = WcsSolution::north_up((180.0, 0.0), 3.6, (100, 100));
        let stars = vec![det(50.0, 50.0, 25.0)];
        let cats = vec![cat(180.0, 0.0)];

        let result = match_by_projection(&tan(&wcs), &wcs, &stars, &cats, 5.0);
        assert_eq!(result.pairs.len(), 1);
        assert_eq!(result.weights[0], 0.9989824448373895);
        assert_eq!(result.pairs[0].snr, 25.0);
    }

    #[test]
    fn predict_cat_in_image_returns_none_outside_bounds() {
        let wcs = WcsSolution::north_up((180.0, 0.0), 3.6, (100, 100));
        assert!(predict_cat_in_image(&tan(&wcs), (100.0, 100.0), &cat(180.0, 0.0)).is_some());
        // Far from tangent center → projection succeeds but pixel lies outside.
        assert!(predict_cat_in_image(&tan(&wcs), (100.0, 100.0), &cat(0.0, 80.0)).is_none());
    }

    #[test]
    fn nearest_detection_picks_closest_within_tolerance() {
        let stars = vec![
            det(50.0, 50.0, 10.0),
            det(52.0, 50.0, 10.0),
            det(60.0, 60.0, 10.0),
        ];
        // pred at (51, 50): closest is star 0 (d=1), then star 1 (d=1); the first wins the tie.
        let (idx, d2) = nearest_detection(&stars, 51.0, 50.0, 5.0 * 5.0).expect("must match");
        assert_eq!((idx, d2), (0, 1.0));
    }

    #[test]
    fn nearest_detection_returns_none_when_all_outside_tolerance() {
        let stars = vec![det(0.0, 0.0, 10.0), det(100.0, 100.0, 10.0)];
        assert!(nearest_detection(&stars, 50.0, 50.0, 5.0 * 5.0).is_none());
    }

    #[test]
    fn log_match_stats_accepts_a_nan_distance() {
        log_match_stats(&[1.0, f64::NAN, 0.5], 3);
    }

    #[test]
    fn dedupe_best_cat_per_det_keeps_closest_cat() {
        // Cat 0 claims det 0 at d²=4. Cat 1 also claims det 0 at d²=1.
        // After dedupe, det 0's slot should point at cat 1.
        let input = vec![Some((0usize, 4.0)), Some((0usize, 1.0)), None];
        let out = dedupe_best_cat_per_det(&input, 1);
        assert_eq!(out.len(), 1);
        let (ci, d2) = out[0].expect("det 0 must be claimed");
        assert_eq!((ci, d2), (1, 1.0));
    }
}
