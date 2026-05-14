//! Geometric hash index over catalog quads and nearest-neighbor lookup.

use std::collections::HashMap;

use celestial_catalog::query::{discrete_hash_key, discrete_hash_key_neighbors, Quad};

use super::QuadMatch;

pub(super) fn build_hash_index(quads: &[Quad]) -> HashMap<u64, Vec<usize>> {
    let mut map = HashMap::with_capacity(quads.len());
    for (i, q) in quads.iter().enumerate() {
        map.entry(discrete_hash_key(&q.hash)).or_insert_with(Vec::new).push(i);
    }
    map
}

pub(super) fn find_matches(
    image_quads: &[Quad],
    catalog_quads: &[Quad],
    index: &HashMap<u64, Vec<usize>>,
) -> Vec<QuadMatch> {
    let mut matches = Vec::new();
    for iq in image_quads {
        let keys = discrete_hash_key_neighbors(&iq.hash);
        for key in &keys {
            if let Some(indices) = index.get(key) {
                for &ci in indices {
                    matches.push(QuadMatch {
                        image_quad: iq.clone(),
                        catalog_quad: catalog_quads[ci].clone(),
                    });
                }
            }
        }
    }
    matches
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad(hash: [f64; 4]) -> Quad {
        Quad { hash, star_indices: [0, 1, 2, 3] }
    }

    #[test]
    fn build_hash_index_empty_is_empty() {
        let idx = build_hash_index(&[]);
        assert!(idx.is_empty());
    }

    #[test]
    fn build_hash_index_one_bucket_per_unique_hash() {
        let quads = vec![
            quad([0.1, 0.2, 0.3, 0.4]),
            quad([0.5, 0.6, 0.7, 0.8]),
        ];
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
}
