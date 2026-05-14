//! Geometric verification of candidate pairs via RANSAC affine fit.
//!
//! A correct quad match set will lie on a single affine transform (pixel → tangent
//! plane); a scrambled match set will not. Counting inliers under a random-sampled
//! affine gives a scalar score used to pick the better parity orientation.

use celestial_catalog::query::tan_project_star;

use super::StarPair;

pub(super) fn verify_pairs(pairs: &[StarPair], center_ra: f64, center_dec: f64) -> usize {
    if pairs.len() < 3 {
        return 0;
    }

    let mut projected: Vec<(f64, f64, f64, f64)> = pairs
        .iter()
        .filter_map(|p| {
            tan_project_star(p.ra_deg, p.dec_deg, center_ra, center_dec)
                .map(|(xi, eta)| (p.px_x, p.px_y, xi, eta))
        })
        .collect();
    if projected.len() < 3 {
        return 0;
    }
    projected.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });

    let n = projected.len();
    let inlier_px = 10.0;
    let trials = 1000;
    let mut best_inliers = 0usize;

    let mut rng: u64 = 0xCAFE_BABE_DEAD_BEEF;
    let mut next = || -> usize {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        (rng as usize) % n
    };

    for _ in 0..trials {
        let a = next();
        let b = next();
        let c = next();
        if a == b || b == c || a == c {
            continue;
        }

        let Some(affine) = fit_affine_3pt(&projected, a, b, c) else {
            continue;
        };
        let count = count_affine_inliers(&projected, &affine, inlier_px);
        if count > best_inliers {
            best_inliers = count;
        }
    }

    best_inliers
}

pub(super) fn fit_affine_3pt(
    pts: &[(f64, f64, f64, f64)],
    a: usize,
    b: usize,
    c: usize,
) -> Option<[f64; 6]> {
    let (x0, y0, u0, v0) = pts[a];
    let (x1, y1, u1, v1) = pts[b];
    let (x2, y2, u2, v2) = pts[c];

    let det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    if libm::fabs(det) < 1e-20 {
        return None;
    }
    let inv = 1.0 / det;

    let a11 = ((u1 - u0) * (y2 - y0) - (u2 - u0) * (y1 - y0)) * inv;
    let a12 = ((u2 - u0) * (x1 - x0) - (u1 - u0) * (x2 - x0)) * inv;
    let a21 = ((v1 - v0) * (y2 - y0) - (v2 - v0) * (y1 - y0)) * inv;
    let a22 = ((v2 - v0) * (x1 - x0) - (v1 - v0) * (x2 - x0)) * inv;
    let tx = u0 - a11 * x0 - a12 * y0;
    let ty = v0 - a21 * x0 - a22 * y0;

    Some([a11, a12, tx, a21, a22, ty])
}

pub(super) fn count_affine_inliers(
    pts: &[(f64, f64, f64, f64)],
    affine: &[f64; 6],
    threshold_px: f64,
) -> usize {
    let det = affine[0] * affine[4] - affine[1] * affine[3];
    if libm::fabs(det) < 1e-30 {
        return 0;
    }
    let inv00 = affine[4] / det;
    let inv01 = -affine[1] / det;
    let inv10 = -affine[3] / det;
    let inv11 = affine[0] / det;
    let thr_sq = threshold_px * threshold_px;

    pts.iter()
        .filter(|(px, py, xi, eta)| {
            let pred_xi = affine[0] * px + affine[1] * py + affine[2];
            let pred_eta = affine[3] * px + affine[4] * py + affine[5];
            let du = xi - pred_xi;
            let dv = eta - pred_eta;
            let dpx = inv00 * du + inv01 * dv;
            let dpy = inv10 * du + inv11 * dv;
            dpx * dpx + dpy * dpy < thr_sq
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_catalog::query::tan_deproject_star;

    fn pair(px: f64, py: f64, ra: f64, dec: f64) -> StarPair {
        StarPair { px_x: px, px_y: py, ra_deg: ra, dec_deg: dec, votes: 2, snr: 10.0 }
    }

    #[test]
    fn fit_affine_identity_recovers_identity() {
        // Build a set of points where xi=px, eta=py, pick three non-collinear.
        let pts = vec![
            (0.0, 0.0, 0.0, 0.0),
            (1.0, 0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0, 1.0),
        ];
        let affine = fit_affine_3pt(&pts, 0, 1, 2).unwrap();
        assert_eq!(affine[0], 1.0);
        assert_eq!(affine[1], 0.0);
        assert_eq!(affine[2], 0.0);
        assert_eq!(affine[3], 0.0);
        assert_eq!(affine[4], 1.0);
        assert_eq!(affine[5], 0.0);
    }

    #[test]
    fn fit_affine_recovers_scale_and_translation() {
        // xi = 2*px + 5, eta = 3*py - 4 (pure scale + translation, no rotation).
        let pts = vec![
            (0.0, 0.0, 5.0, -4.0),
            (1.0, 0.0, 7.0, -4.0),
            (0.0, 1.0, 5.0, -1.0),
        ];
        let a = fit_affine_3pt(&pts, 0, 1, 2).unwrap();
        assert_eq!(a[0], 2.0);
        assert_eq!(a[1], 0.0);
        assert_eq!(a[2], 5.0);
        assert_eq!(a[3], 0.0);
        assert_eq!(a[4], 3.0);
        assert_eq!(a[5], -4.0);
    }

    #[test]
    fn fit_affine_collinear_points_return_none() {
        // Three points on the x-axis → zero-area triangle.
        let pts = vec![
            (0.0, 0.0, 0.0, 0.0),
            (1.0, 0.0, 1.0, 0.0),
            (2.0, 0.0, 2.0, 0.0),
        ];
        assert!(fit_affine_3pt(&pts, 0, 1, 2).is_none());
    }

    #[test]
    fn count_affine_inliers_identity_counts_exact_matches() {
        let pts = vec![
            (0.0, 0.0, 0.0, 0.0),
            (10.0, 0.0, 10.0, 0.0),
            (0.0, 10.0, 0.0, 10.0),
            (100.0, 100.0, 100.0, 100.0),
        ];
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        assert_eq!(count_affine_inliers(&pts, &identity, 0.5), 4);
    }

    #[test]
    fn count_affine_inliers_excludes_far_points() {
        // Identity maps px→xi, py→eta. Any offset in xi/eta shows up equally in pixel units.
        let pts = vec![
            (0.0, 0.0, 0.0, 0.0),                // on-model → in
            (10.0, 0.0, 10.0, 0.0),               // on-model → in
            (20.0, 0.0, 20.0 + 5.0, 0.0),         // 5-px offset → out at threshold 1.0
        ];
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        assert_eq!(count_affine_inliers(&pts, &identity, 1.0), 2);
    }

    #[test]
    fn count_affine_inliers_singular_affine_returns_zero() {
        // det(A) = 1*0 - 0*1 = 0 → degenerate model.
        let pts = vec![(0.0, 0.0, 0.0, 0.0)];
        let singular = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        assert_eq!(count_affine_inliers(&pts, &singular, 10.0), 0);
    }

    #[test]
    fn verify_pairs_returns_zero_below_three_pairs() {
        let pairs = vec![pair(0.0, 0.0, 180.0, 0.0), pair(1.0, 0.0, 180.001, 0.0)];
        assert_eq!(verify_pairs(&pairs, 180.0, 0.0), 0);
    }

    #[test]
    fn verify_pairs_scores_consistent_affine_set() {
        // 25 pairs, each constructed so (ra, dec) is the tangent-plane deprojection
        // of (xi, eta) around a known center. All lie on a single affine in pixel↔
        // tangent space — RANSAC must score them as nearly-all inliers.
        let center_ra = 180.0;
        let center_dec = 30.0;
        let mut pairs = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                // xi, eta in radians; verify_pairs re-projects via tan_project_star.
                let xi = (i as f64 - 2.0) * 0.01;
                let eta = (j as f64 - 2.0) * 0.01;
                let (ra, dec) = tan_deproject_star(xi, eta, center_ra, center_dec);
                pairs.push(pair(xi * 1000.0, eta * 1000.0, ra, dec));
            }
        }
        let score = verify_pairs(&pairs, center_ra, center_dec);
        assert!(score >= 20, "expected at least 20 inliers, got {score}");
    }

    #[test]
    fn verify_pairs_scrambled_pairs_score_low() {
        // Random-looking pixel↔sky assignments — no affine fits them all.
        let pairs = vec![
            pair(0.0, 0.0, 180.0, 30.0),
            pair(100.0, 0.0, 45.0, -20.0),
            pair(0.0, 100.0, 270.0, 60.0),
            pair(50.0, 50.0, 10.0, 0.0),
            pair(200.0, 100.0, 120.0, -40.0),
            pair(30.0, 200.0, 300.0, 15.0),
        ];
        let score = verify_pairs(&pairs, 180.0, 30.0);
        // Fixed RNG, so this is deterministic. Scrambled → score below inlier ceiling.
        assert!(score < pairs.len(), "scrambled pairs should not all be inliers; got {score}");
    }
}
