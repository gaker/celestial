//! Geometric hash index over catalog quads and nearest-neighbor lookup.

use std::collections::HashMap;

use celestial_catalog::query::quad::{
    discrete_hash_key, discrete_hash_key_neighbors, relabelled_quads, Quad,
};

use super::QuadMatch;

pub(super) fn build_hash_index(quads: &[Quad]) -> HashMap<u64, Vec<usize>> {
    let mut map = HashMap::with_capacity(quads.len());
    for (i, q) in quads.iter().enumerate() {
        map.entry(discrete_hash_key(&q.hash))
            .or_insert_with(Vec::new)
            .push(i);
    }
    map
}

// Each image quad is looked up under its own labels and under every relabelling that a
// small position error could have given the same stars, so a catalog quad labelled the
// other way across a boundary is still found. A match carries the image quad with the
// star order of the code that matched.
pub(super) fn find_matches(
    image_quads: &[Quad],
    catalog_quads: &[Quad],
    index: &HashMap<u64, Vec<usize>>,
) -> Vec<QuadMatch> {
    let mut matches = Vec::new();
    for iq in image_quads {
        for quad in std::iter::once(iq.clone()).chain(relabelled_quads(iq)) {
            push_bin_matches(&quad, catalog_quads, index, &mut matches);
        }
    }
    matches
}

fn push_bin_matches(
    image_quad: &Quad,
    catalog_quads: &[Quad],
    index: &HashMap<u64, Vec<usize>>,
    matches: &mut Vec<QuadMatch>,
) {
    for key in &discrete_hash_key_neighbors(&image_quad.hash) {
        for &ci in index.get(key).into_iter().flatten() {
            matches.push(QuadMatch {
                image_quad: image_quad.clone(),
                catalog_quad: catalog_quads[ci].clone(),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::quad_stars::{stars_to_quad_stars, CdSign};
    use super::super::MatchParams;
    use super::*;
    use crate::detect::DetectedStar;
    use celestial_catalog::query::quad::{neighbor_quads, QuadStar};
    use celestial_core::constants::{ARCSEC_TO_RAD, DEG_TO_RAD};

    fn quad(hash: [f64; 4]) -> Quad {
        Quad {
            hash,
            star_indices: [0, 1, 2, 3],
        }
    }

    #[test]
    fn build_hash_index_empty_is_empty() {
        let idx = build_hash_index(&[]);
        assert!(idx.is_empty());
    }

    #[test]
    fn build_hash_index_one_bucket_per_unique_hash() {
        let quads = vec![quad([0.1, 0.2, 0.3, 0.4]), quad([0.5, 0.6, 0.7, 0.8])];
        let idx = build_hash_index(&quads);
        assert_eq!(idx.len(), 2);
    }

    #[test]
    fn build_hash_index_groups_same_hash() {
        let quads = vec![
            quad([0.1, 0.2, 0.3, 0.4]),
            quad([0.1, 0.2, 0.3, 0.4]),
            quad([0.1, 0.2, 0.3, 0.4]),
        ];
        let idx = build_hash_index(&quads);
        assert_eq!(idx.len(), 1);
        let bucket = idx.values().next().unwrap();
        assert_eq!(bucket.len(), 3);
        assert!(bucket.contains(&0));
        assert!(bucket.contains(&1));
        assert!(bucket.contains(&2));
    }

    #[test]
    fn find_matches_empty_image_quads_yields_none() {
        let cat = vec![quad([0.1, 0.2, 0.3, 0.4])];
        let idx = build_hash_index(&cat);
        let m = find_matches(&[], &cat, &idx);
        assert!(m.is_empty());
    }

    #[test]
    fn find_matches_empty_index_yields_none() {
        let img = vec![quad([0.1, 0.2, 0.3, 0.4])];
        let idx = HashMap::new();
        let m = find_matches(&img, &[], &idx);
        assert!(m.is_empty());
    }

    #[test]
    fn find_matches_exact_hash_matches() {
        let cat = vec![quad([0.30, 0.40, 0.70, 0.10])];
        let idx = build_hash_index(&cat);
        let img = vec![quad([0.30, 0.40, 0.70, 0.10])];
        let m = find_matches(&img, &cat, &idx);
        assert!(!m.is_empty());
        assert_eq!(m[0].image_quad.hash, img[0].hash);
        assert_eq!(m[0].catalog_quad.hash, cat[0].hash);
    }

    #[test]
    fn find_matches_near_hash_matches_via_neighbors() {
        // Within ±0.01 of the catalog bucket — neighbors() covers ±1 bucket shift.
        let cat = vec![quad([0.300, 0.400, 0.700, 0.100])];
        let idx = build_hash_index(&cat);
        let img = vec![quad([0.304, 0.397, 0.702, 0.103])];
        let m = find_matches(&img, &cat, &idx);
        assert!(!m.is_empty(), "nearby hash must match via neighbor lookup");
    }

    #[test]
    fn find_matches_far_hash_does_not_match() {
        let cat = vec![quad([0.30, 0.40, 0.70, 0.10])];
        let idx = build_hash_index(&cat);
        let img = vec![quad([0.70, 0.10, 0.30, 0.50])];
        let m = find_matches(&img, &cat, &idx);
        assert!(m.is_empty());
    }

    #[test]
    fn find_matches_emits_one_per_catalog_hit() {
        let cat = vec![
            quad([0.30, 0.40, 0.70, 0.10]),
            quad([0.30, 0.40, 0.70, 0.10]),
            quad([0.30, 0.40, 0.70, 0.10]),
        ];
        let idx = build_hash_index(&cat);
        let img = vec![quad([0.30, 0.40, 0.70, 0.10])];
        let m = find_matches(&img, &cat, &idx);
        // 3 catalog quads sharing the same hash → 3 matches from one image quad
        // (plus any neighbor-bucket duplicates would also appear, but neighbors
        // of an exact hash don't collide here).
        assert!(m.len() >= 3);
    }

    // A 60-star field in tangent-plane radians, about 0.6 degrees across.
    fn field_points() -> Vec<(f64, f64)> {
        let mut state = 5u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..60)
            .map(|_| (0.01 * uniform() - 0.005, 0.01 * uniform() - 0.005))
            .collect()
    }

    // Detections of `points` rotated by `deg`, at 1 arcsec per pixel. A mirrored image
    // has its y axis reversed, so its det(CD) is negative, and stars_to_quad_stars with
    // CdSign::Negative undoes the mirror.
    fn detections(points: &[(f64, f64)], deg: f64, mirrored: bool) -> Vec<DetectedStar> {
        let rad_per_px = ARCSEC_TO_RAD;
        let (s, c) = libm::sincos(deg * DEG_TO_RAD);
        let y_sign = if mirrored { 1.0 } else { -1.0 };
        points
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| DetectedStar {
                x: 2048.0 - (x * c - y * s) / rad_per_px,
                y: 2048.0 + y_sign * (x * s + y * c) / rad_per_px,
                flux: 10_000.0 - i as f64,
                snr: 20.0,
                saturated: false,
                saturated_count: 0,
                background: 0.0,
            })
            .collect()
    }

    // Counts image quads that find the catalog quad built from the same four stars with
    // the same A, B, C, D labels, which extract_pairs relies on. Detection i is catalog
    // star i, and the flux order keeps image star i at index i.
    fn own_quad_matches(image_quads: &[Quad], catalog_quads: &[Quad]) -> usize {
        let index = build_hash_index(catalog_quads);
        let mut found = std::collections::HashSet::new();
        for m in find_matches(image_quads, catalog_quads, &index) {
            if m.image_quad.star_indices == m.catalog_quad.star_indices {
                found.insert(m.image_quad.star_indices);
            }
        }
        found.len()
    }

    fn spine_quad(c: (f64, f64), d: (f64, f64)) -> Vec<Quad> {
        let points = [(0.0, 0.0), (1.0, 0.0), c, d];
        let stars: Vec<QuadStar> = points
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| QuadStar {
                source_id: i as i64,
                x,
                y,
                mag: 0.0,
            })
            .collect();
        neighbor_quads(&stars, 3, 0.0)
    }

    // In each case one star moves by 0.01 between the catalog and the image, which takes
    // the quad across a labelling boundary: cx + dx crosses 1, which swaps A and B, or
    // cx crosses dx, which swaps C and D. The image quad's own code is then far from the
    // catalog's, and only the relabelled code finds it.
    #[test]
    fn find_matches_finds_quads_labelled_across_a_boundary() {
        let cases = [
            (((0.3, 0.2), (0.695, -0.3)), ((0.3, 0.2), (0.705, -0.3))),
            (((0.4, 0.2), (0.405, -0.3)), ((0.41, 0.2), (0.405, -0.3))),
        ];
        for ((cat_c, cat_d), (img_c, img_d)) in cases {
            let catalog_quads = spine_quad(cat_c, cat_d);
            let image_quads = spine_quad(img_c, img_d);
            assert_ne!(catalog_quads[0].star_indices, image_quads[0].star_indices);
            assert_eq!(
                own_quad_matches(&image_quads, &catalog_quads),
                1,
                "{cat_c:?} {cat_d:?}"
            );
        }
    }

    #[test]
    fn mirrored_rotated_field_matches_with_a_negative_cd_sign() {
        let points = field_points();
        let catalog: Vec<QuadStar> = points
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| QuadStar {
                source_id: i as i64,
                x,
                y,
                mag: 0.0,
            })
            .collect();
        let catalog_quads = neighbor_quads(&catalog, 10, 0.0);
        let params = MatchParams::default();
        let mut counts = Vec::new();
        for deg in [0.0, 90.0, 150.0, 180.0] {
            let stars = detections(&points, deg, true);
            let found: Vec<usize> = [CdSign::Negative, CdSign::Positive]
                .iter()
                .map(|&cd_sign| {
                    let img = stars_to_quad_stars(&stars, (4096, 4096), 1.0, &params, cd_sign);
                    own_quad_matches(&neighbor_quads(&img, 10, 0.0), &catalog_quads)
                })
                .collect();
            counts.push((deg, catalog_quads.len(), found));
        }
        assert_eq!(
            counts,
            [
                (0.0, 4663, vec![4663, 0]),
                (90.0, 4663, vec![4663, 0]),
                (150.0, 4663, vec![4663, 0]),
                (180.0, 4663, vec![4663, 0]),
            ]
        );
    }
}
