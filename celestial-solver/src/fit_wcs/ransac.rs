// RANSAC selection of the pairs that one affine CD fit agrees with.

use crate::error::SolveError;
use crate::match_field::StarPair;
use crate::triple_search::{index_triples, MAX_SAMPLE_POOL, MIN_CONFIRMED_INLIERS};

use super::cd_matrix::{fit_cd_weighted, InverseAffine};

struct InlierSet {
    count: usize,
    mask: Vec<bool>,
}

// `pairs[i]` must describe `projected[i]`.
pub(super) fn ransac_fit(
    projected: &[(f64, f64, f64, f64)],
    pairs: &[&StarPair],
    inlier_threshold: f64,
) -> Result<Vec<bool>, SolveError> {
    let n = projected.len();
    if pairs.len() != n {
        return Err(SolveError::fit_failed(format!(
            "RANSAC input mismatch: {} pairs, {} projected",
            pairs.len(),
            n
        )));
    }
    let best = best_inlier_set(projected, &sample_pool(pairs), inlier_threshold);
    if best.count < MIN_CONFIRMED_INLIERS {
        return Err(SolveError::fit_failed(format!(
            "RANSAC found no consistent model ({n} pairs)"
        )));
    }
    log::debug!("RANSAC: {} inliers from {} pairs", best.count, n);
    Ok(best.mask)
}

// Tries the affine model of every triple from `pool`, in pool order, and keeps the first
// one with strictly the most inliers.
fn best_inlier_set(
    projected: &[(f64, f64, f64, f64)],
    pool: &[usize],
    threshold: f64,
) -> InlierSet {
    let mut best = InlierSet {
        count: 0,
        mask: vec![false; projected.len()],
    };
    for [a, b, c] in index_triples(pool.len()) {
        let sample = [projected[pool[a]], projected[pool[b]], projected[pool[c]]];
        let Ok(cd) = fit_cd_weighted(&sample, &[true; 3], &[1.0; 3]) else {
            continue;
        };
        let model = InverseAffine::new(cd);
        let count = count_inliers(projected, &model, threshold);
        if count > best.count {
            let found = InlierSet {
                count,
                mask: inlier_mask(projected, &model, threshold),
            };
            best = refit_on_inliers(projected, found, threshold);
            if best.count == projected.len() {
                break;
            }
        }
    }
    best
}

// Indices of the pairs that samples are drawn from, in `StarPair::cmp_rank` order so the
// search does not depend on input order. Pairs with at least a third of the top vote
// count qualify; when fewer than 3 do, every pair does.
fn sample_pool(pairs: &[&StarPair]) -> Vec<usize> {
    let max_v = pairs.iter().map(|p| p.votes).max().unwrap_or(0);
    let mut pool: Vec<usize> = (0..pairs.len())
        .filter(|&i| pairs[i].votes >= max_v / 3)
        .collect();
    if pool.len() < 3 {
        pool = (0..pairs.len()).collect();
    }
    pool.sort_by(|&a, &b| pairs[a].cmp_rank(pairs[b]));
    pool.truncate(MAX_SAMPLE_POOL);
    pool
}

// Counts without allocating: this runs once per sampled triple.
fn count_inliers(
    projected: &[(f64, f64, f64, f64)],
    model: &InverseAffine,
    threshold: f64,
) -> usize {
    projected
        .iter()
        .filter(|row| model.pixel_error(row) <= threshold)
        .count()
}

fn inlier_mask(
    projected: &[(f64, f64, f64, f64)],
    model: &InverseAffine,
    threshold: f64,
) -> Vec<bool> {
    projected
        .iter()
        .map(|row| model.pixel_error(row) <= threshold)
        .collect()
}

// A least-squares fit to all inliers can pull in pairs the 3-point model missed; it is
// kept only when it does not lose inliers.
fn refit_on_inliers(
    projected: &[(f64, f64, f64, f64)],
    found: InlierSet,
    threshold: f64,
) -> InlierSet {
    let Ok(cd) = fit_cd_weighted(projected, &found.mask, &vec![1.0; projected.len()]) else {
        return found;
    };
    let model = InverseAffine::new(cd);
    let count = count_inliers(projected, &model, threshold);
    if count < found.count {
        return found;
    }
    InlierSet {
        count,
        mask: inlier_mask(projected, &model, threshold),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use celestial_catalog::query::quad::tan_project_star;
    use celestial_core::constants::RAD_TO_DEG;

    fn refs(pairs: &[StarPair]) -> Vec<&StarPair> {
        pairs.iter().collect()
    }

    // The shared minority scene and its tangent-plane rows as `fit_wcs` builds them.
    type Scene = (Vec<StarPair>, Vec<(f64, f64, f64, f64)>);

    fn minority_scene(n_good: usize, n_junk: usize) -> Scene {
        let pairs = crate::triple_search::test_scene::minority_scene(n_good, n_junk);
        let projected = pairs
            .iter()
            .filter_map(|p| {
                tan_project_star(p.ra_deg, p.dec_deg, 180.0, 0.0)
                    .map(|(xi, eta)| (p.px_x, p.px_y, xi * RAD_TO_DEG, eta * RAD_TO_DEG))
            })
            .collect();
        (pairs, projected)
    }

    #[test]
    fn ransac_fit_finds_minority_model_among_wrong_pairs() {
        let (pairs, projected) = minority_scene(5, 55);
        let mask = ransac_fit(&projected, &refs(&pairs), 5.0).unwrap();
        let expected: Vec<bool> = (0..60).map(|i| i < 5).collect();
        assert_eq!(mask, expected);
    }

    #[test]
    fn ransac_fit_rejects_model_supported_only_by_its_sample() {
        let (pairs, projected) = minority_scene(3, 5);
        let err = ransac_fit(&projected, &refs(&pairs), 5.0).unwrap_err();
        assert_eq!(
            err.to_string(),
            "RANSAC found no consistent model (8 pairs)"
        );
    }

    #[test]
    fn ransac_fit_rejects_exactly_three_pairs() {
        let (pairs, projected) = minority_scene(3, 0);
        let err = ransac_fit(&projected, &refs(&pairs), 5.0).unwrap_err();
        assert_eq!(
            err.to_string(),
            "RANSAC found no consistent model (3 pairs)"
        );
    }

    #[test]
    fn ransac_fit_inlier_set_is_independent_of_pair_order() {
        let (pairs, projected) = minority_scene(5, 45);
        let reference = ransac_fit(&projected, &refs(&pairs), 5.0).unwrap();
        for k in [1, 2, 3, 17, 31, 49] {
            let mut rot_pairs = pairs.clone();
            let mut rot_projected = projected.clone();
            rot_pairs.rotate_left(k);
            rot_projected.rotate_left(k);
            let mut mask = ransac_fit(&rot_projected, &refs(&rot_pairs), 5.0).unwrap();
            mask.rotate_right(k);
            assert_eq!(mask, reference, "rotation {k}");
        }
    }
}
