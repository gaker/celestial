//! Vote aggregation: reduce per-quad matches to unique image↔catalog pairs.
//!
//! Each quad match casts four votes (one per star in the quad). The same
//! (image_star, catalog_star) pair can be voted for by many quads; this module
//! picks the best catalog for each image star, resolves catalog collisions, and
//! drops singletons.

use std::collections::HashMap;

use celestial_catalog::query::{ConeSearchResult, QuadStar};

use super::{QuadMatch, StarPair};
use crate::detect::DetectedStar;

#[derive(Debug, Clone, Copy)]
pub(super) struct VoteEntry {
    pub px_x: f64,
    pub px_y: f64,
    pub ra: f64,
    pub dec: f64,
    pub snr: f64,
    pub votes: usize,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct BestMatch {
    pub cat_id: i64,
    pub entry: VoteEntry,
}

pub(super) fn extract_pairs(
    matches: &[QuadMatch],
    detected: &[DetectedStar],
    image_stars: &[QuadStar],
    catalog_stars: &[QuadStar],
    catalog_results: &[ConeSearchResult],
) -> Vec<StarPair> {
    let ra_dec: HashMap<i64, (f64, f64)> = catalog_results
        .iter()
        .map(|r| (r.star.source_id, (r.ra_deg, r.dec_deg)))
        .collect();

    let mut votes: HashMap<(i64, i64), VoteEntry> = HashMap::new();
    for m in matches {
        for i in 0..4 {
            let img_idx = m.image_quad.star_indices[i];
            let cat_idx = m.catalog_quad.star_indices[i];
            if img_idx >= image_stars.len() || cat_idx >= catalog_stars.len() {
                continue;
            }
            let img_star = &image_stars[img_idx];
            let cat_id = catalog_stars[cat_idx].source_id;
            let Some(&(ra, dec)) = ra_dec.get(&cat_id) else { continue };
            let det_idx = img_star.source_id as usize;
            let det = &detected[det_idx];
            let key = (img_star.source_id, cat_id);
            let entry = votes.entry(key).or_insert(VoteEntry {
                px_x: det.x, px_y: det.y, ra, dec, snr: det.snr, votes: 0,
            });
            entry.votes += 1;
        }
    }

    let mut by_image: HashMap<i64, BestMatch> = HashMap::new();
    for (&(img_id, cat_id), entry) in &votes {
        by_image.entry(img_id)
            .and_modify(|best| {
                if entry.votes > best.entry.votes
                    || (entry.votes == best.entry.votes && cat_id < best.cat_id)
                {
                    *best = BestMatch { cat_id, entry: *entry };
                }
            })
            .or_insert(BestMatch { cat_id, entry: *entry });
    }

    let mut used_cat: HashMap<i64, (i64, usize)> = HashMap::new();
    for (&img_id, best) in &by_image {
        used_cat.entry(best.cat_id)
            .and_modify(|e| {
                if best.entry.votes > e.1 || (best.entry.votes == e.1 && img_id < e.0) {
                    *e = (img_id, best.entry.votes);
                }
            })
            .or_insert((img_id, best.entry.votes));
    }

    let mut pairs: Vec<StarPair> = by_image
        .iter()
        .filter(|(&img_id, best)| {
            best.entry.votes >= 2 && used_cat.get(&best.cat_id).is_some_and(|e| e.0 == img_id)
        })
        .map(|(_, best)| StarPair {
            px_x: best.entry.px_x,
            px_y: best.entry.px_y,
            ra_deg: best.entry.ra,
            dec_deg: best.entry.dec,
            votes: best.entry.votes,
            snr: best.entry.snr,
        })
        .collect();
    pairs.sort_by_key(|p| std::cmp::Reverse(p.votes));
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_catalog::query::{ConeSearchResult, Quad, QuadStar};

    use crate::detect::DetectedStar;

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

    fn img_star(id: i64) -> QuadStar {
        QuadStar { source_id: id, x: 0.0, y: 0.0, mag: 0.0 }
    }

    fn cat_star(id: i64) -> QuadStar {
        QuadStar { source_id: id, x: 0.0, y: 0.0, mag: 0.0 }
    }

    fn cone_result(source_id: i64, ra: f64, dec: f64) -> ConeSearchResult {
        let mut star: celestial_catalog::query::StarRecord = unsafe { std::mem::zeroed() };
        star.source_id = source_id;
        star.ra = ra;
        star.dec = dec;
        star.mag = 10.0;
        ConeSearchResult { star, ra_deg: ra, dec_deg: dec, distance_deg: 0.0 }
    }

    fn quad_match(img_idx: [usize; 4], cat_idx: [usize; 4]) -> QuadMatch {
        QuadMatch {
            image_quad: Quad { hash: [0.0; 4], star_indices: img_idx },
            catalog_quad: Quad { hash: [0.0; 4], star_indices: cat_idx },
        }
    }

    #[test]
    fn empty_matches_produce_empty_pairs() {
        let detected = vec![det(0.0, 0.0, 10.0)];
        let img = vec![img_star(0)];
        let cat = vec![cat_star(1000)];
        let res = vec![cone_result(1000, 10.0, 20.0)];
        let pairs = extract_pairs(&[], &detected, &img, &cat, &res);
        assert!(pairs.is_empty());
    }

    #[test]
    fn single_quad_match_drops_singletons() {
        // Each star in a single quad gets exactly 1 vote — below the >=2 threshold.
        let detected = vec![det(10.0, 20.0, 30.0); 4];
        let img: Vec<_> = (0..4).map(|i| img_star(i as i64)).collect();
        let cat: Vec<_> = (0..4).map(|i| cat_star(1000 + i as i64)).collect();
        let res: Vec<_> = (0..4)
            .map(|i| cone_result(1000 + i as i64, i as f64, i as f64))
            .collect();
        let m = vec![quad_match([0, 1, 2, 3], [0, 1, 2, 3])];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        assert!(pairs.is_empty(), "singletons must be filtered");
    }

    #[test]
    fn two_overlapping_quads_promote_pair() {
        // Two quad matches that share star (img 0 ↔ cat 1000) → 2 votes → kept.
        let detected = vec![
            det(1.0, 2.0, 50.0),
            det(3.0, 4.0, 40.0),
            det(5.0, 6.0, 30.0),
            det(7.0, 8.0, 20.0),
            det(9.0, 10.0, 10.0),
        ];
        let img: Vec<_> = (0..5).map(|i| img_star(i as i64)).collect();
        let cat: Vec<_> = (0..5).map(|i| cat_star(1000 + i as i64)).collect();
        let res: Vec<_> = (0..5)
            .map(|i| cone_result(1000 + i as i64, 100.0 + i as f64, i as f64))
            .collect();
        let m = vec![
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
            quad_match([0, 1, 2, 4], [0, 1, 2, 4]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        // Stars 0, 1, 2 got 2 votes each; 3 and 4 only 1.
        assert_eq!(pairs.len(), 3);
        for p in &pairs {
            assert_eq!(p.votes, 2);
        }
        // Top pair carries the correct detection-side and catalog-side data.
        let p0 = pairs.iter().find(|p| p.px_x == 1.0).unwrap();
        assert_eq!(p0.px_y, 2.0);
        assert_eq!(p0.snr, 50.0);
        assert_eq!(p0.ra_deg, 100.0);
    }

    #[test]
    fn pairs_sorted_by_descending_votes() {
        let detected = vec![det(1.0, 1.0, 10.0); 6];
        let img: Vec<_> = (0..6).map(|i| img_star(i as i64)).collect();
        let cat: Vec<_> = (0..6).map(|i| cat_star(1000 + i as i64)).collect();
        let res: Vec<_> = (0..6)
            .map(|i| cone_result(1000 + i as i64, i as f64, 0.0))
            .collect();
        // (0↔1000) appears in 3 quads, (1↔1001) in 2, (2↔1002) in 2.
        let m = vec![
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
            quad_match([0, 1, 2, 4], [0, 1, 2, 4]),
            quad_match([0, 3, 4, 5], [0, 3, 4, 5]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        for w in pairs.windows(2) {
            assert!(w[0].votes >= w[1].votes);
        }
        assert_eq!(pairs[0].votes, 3);
    }

    #[test]
    fn out_of_range_indices_are_skipped() {
        let detected = vec![det(0.0, 0.0, 10.0)];
        let img = vec![img_star(0)];
        let cat = vec![cat_star(1000)];
        let res = vec![cone_result(1000, 1.0, 2.0)];
        // Both quads reference indices past the end of img/cat — must not panic.
        let m = vec![
            quad_match([99, 99, 99, 99], [99, 99, 99, 99]),
            quad_match([99, 99, 99, 99], [99, 99, 99, 99]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        assert!(pairs.is_empty());
    }

    #[test]
    fn missing_catalog_ra_dec_is_skipped() {
        // Catalog star has source_id 1000 but ra_dec map only has 2000 → entry skipped.
        let detected = vec![det(0.0, 0.0, 10.0); 4];
        let img: Vec<_> = (0..4).map(|i| img_star(i as i64)).collect();
        let cat: Vec<_> = (0..4).map(|_| cat_star(1000)).collect();
        let res = vec![cone_result(2000, 5.0, 6.0)];
        let m = vec![
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        assert!(pairs.is_empty());
    }

    #[test]
    fn catalog_collision_keeps_higher_vote_image() {
        // Two image stars (0 and 1) both vote for catalog source 7000, with
        // distinct unique catalog source ids for the remaining quad slots so
        // by_image entries for indices 2..=4 don't also collide on 7000.
        let detected = vec![det(0.0, 0.0, 10.0); 5];
        let img: Vec<_> = (0..5).map(|i| img_star(i as i64)).collect();
        // cat[0] collides on source 7000; cat[1..] are distinct.
        let cat = vec![
            cat_star(7000),
            cat_star(7001),
            cat_star(7002),
            cat_star(7003),
            cat_star(7004),
        ];
        let res = vec![
            cone_result(7000, 5.0, 6.0),
            cone_result(7001, 1.0, 1.0),
            cone_result(7002, 2.0, 2.0),
            cone_result(7003, 3.0, 3.0),
            cone_result(7004, 4.0, 4.0),
        ];
        // Img 0 ↔ cat 7000 votes: 3 quads. Img 1 ↔ cat 7000 votes: 2 quads.
        // Other slots vote for distinct cat ids so they don't fight for 7000.
        let m = vec![
            quad_match([0, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([0, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([0, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([1, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([1, 2, 3, 4], [0, 2, 3, 4]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        let cat7000 = pairs.iter().find(|p| p.ra_deg == 5.0);
        assert!(cat7000.is_some(), "img 0 should claim cat 7000");
        assert_eq!(cat7000.unwrap().votes, 3);
        // Img 1's entry for cat 7000 must have been suppressed by used_cat.
        let losers = pairs.iter().filter(|p| p.ra_deg == 5.0).count();
        assert_eq!(losers, 1);
    }

    #[test]
    fn by_image_keeps_highest_vote_catalog_when_image_appears_twice() {
        // Image star 0 votes for two different catalog source ids across two
        // sets of quads. The candidate with more votes must win by_image's slot.
        let detected = vec![det(0.0, 0.0, 10.0); 4];
        let img: Vec<_> = (0..4).map(|i| img_star(i as i64)).collect();
        let cat: Vec<_> = vec![
            cat_star(8000),
            cat_star(8001),
            cat_star(8002),
            cat_star(8003),
            cat_star(9000),
        ];
        let res = vec![
            cone_result(8000, 1.0, 1.0),
            cone_result(8001, 2.0, 2.0),
            cone_result(8002, 3.0, 3.0),
            cone_result(8003, 4.0, 4.0),
            cone_result(9000, 5.0, 5.0),
        ];
        // 3 quads pair img 0 with cat slot 0 (source 8000). 2 quads pair img 0 with
        // cat slot 4 (source 9000). Other slots use disjoint image ids.
        let m = vec![
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
            quad_match([0, 1, 2, 3], [0, 1, 2, 3]),
            quad_match([0, 1, 2, 3], [4, 1, 2, 3]),
            quad_match([0, 1, 2, 3], [4, 1, 2, 3]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        // Img 0 had two catalog candidates; by_image's and_modify branch picks 8000 (3 votes).
        let img0_pair = pairs.iter().find(|p| p.ra_deg == 1.0);
        assert!(img0_pair.is_some(), "img 0 must claim catalog 8000 (3 votes)");
        assert_eq!(img0_pair.unwrap().votes, 3);
        // The losing candidate (9000) must not appear as img 0's pair.
        assert!(pairs.iter().find(|p| p.ra_deg == 5.0).is_none());
    }

    #[test]
    fn catalog_collision_ties_break_on_lowest_image_id() {
        // Same vote count → smaller image source_id wins per the tiebreak rule.
        let detected = vec![det(0.0, 0.0, 10.0); 5];
        let img: Vec<_> = (0..5).map(|i| img_star(i as i64)).collect();
        let cat = vec![
            cat_star(7000),
            cat_star(7001),
            cat_star(7002),
            cat_star(7003),
            cat_star(7004),
        ];
        let res = vec![
            cone_result(7000, 5.0, 6.0),
            cone_result(7001, 1.0, 1.0),
            cone_result(7002, 2.0, 2.0),
            cone_result(7003, 3.0, 3.0),
            cone_result(7004, 4.0, 4.0),
        ];
        let m = vec![
            quad_match([0, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([0, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([1, 2, 3, 4], [0, 2, 3, 4]),
            quad_match([1, 2, 3, 4], [0, 2, 3, 4]),
        ];
        let pairs = extract_pairs(&m, &detected, &img, &cat, &res);
        // Img 0 and img 1 tied at 2 votes for cat 7000; lower id wins → only one pair at ra=5.0.
        let cat7000 = pairs.iter().filter(|p| p.ra_deg == 5.0).count();
        assert_eq!(cat7000, 1);
    }
}
