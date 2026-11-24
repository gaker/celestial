//! Exhaustive enumeration of 3-pair samples for affine model searches.
//!
//! The RANSAC inlier selection in `fit_wcs` fits an affine model to every triple drawn
//! from a ranked pool of pairs.

// Samples come from at most this many pairs, so a search fits at most
// C(100, 3) = 161,700 affine models.
pub(crate) const MAX_SAMPLE_POOL: usize = 100;

// Three pairs fix the six affine parameters exactly, so every non-degenerate sample has
// 3 inliers. A model counts only when at least one pair outside its sample agrees. One
// agreeing pair does not separate a real match from a chance one. `solve` verifies its
// pairs by log-odds before and after the fit; a direct caller of `fit_wcs` gets no such
// check.
pub(crate) const MIN_CONFIRMED_INLIERS: usize = 4;

// Every index triple `i < j < k` below `len`, in lexicographic order.
pub(crate) fn index_triples(len: usize) -> impl Iterator<Item = [usize; 3]> {
    (0..len)
        .flat_map(move |i| ((i + 1)..len).flat_map(move |j| ((j + 1)..len).map(move |k| [i, j, k])))
}

#[cfg(test)]
pub(crate) mod test_scene {
    use celestial_catalog::query::quad::tan_deproject_star;
    use celestial_core::constants::DEG_TO_RAD;

    use crate::match_field::StarPair;

    // Deterministic coordinates in [-1000, 1000) from a 64-bit LCG.
    fn lcg_coords(seed: u64, n: usize) -> Vec<(f64, f64)> {
        let mut s = seed;
        let mut next = || {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (s >> 11) as f64 / (1u64 << 53) as f64 * 2000.0 - 1000.0
        };
        (0..n).map(|_| (next(), next())).collect()
    }

    // Sky position of pixel offset (u, v) under a 1.44"/px TAN WCS centred on (180, 0).
    fn sky_of(u: f64, v: f64) -> (f64, f64) {
        tan_deproject_star(0.0004 * u * DEG_TO_RAD, 0.0004 * v * DEG_TO_RAD, 180.0, 0.0)
    }

    fn pair(px_x: f64, px_y: f64, (ra_deg, dec_deg): (f64, f64)) -> StarPair {
        StarPair {
            px_x,
            px_y,
            ra_deg,
            dec_deg,
            votes: 4,
            snr: 20.0,
        }
    }

    // `n_good` pairs that follow the WCS above, then `n_junk` pairs whose pixel and sky
    // positions are independent. All pairs have equal votes and snr.
    pub(crate) fn minority_scene(n_good: usize, n_junk: usize) -> Vec<StarPair> {
        let good = lcg_coords(1, n_good)
            .into_iter()
            .map(|(u, v)| pair(u, v, sky_of(u, v)));
        let junk = lcg_coords(2, n_junk)
            .into_iter()
            .zip(lcg_coords(3, n_junk))
            .map(|((u, v), (su, sv))| pair(u, v, sky_of(su, sv)));
        good.chain(junk).collect()
    }
}
