//! Vote aggregation: reduce per-quad matches to unique image↔catalog pairs.
//!
//! Each quad match casts four votes (one per star in the quad). The same
//! (image_star, catalog_star) pair can be voted for by many quads; this module
//! picks the best catalog for each image star, resolves catalog collisions, and
//! drops singletons.

use std::collections::{BTreeMap, HashMap};

use celestial_catalog::query::cone::ConeSearchResult;
use celestial_catalog::query::quad::QuadStar;

use super::{QuadMatch, StarPair};
use crate::detect::DetectedStar;

#[derive(Debug, Clone, Copy)]
pub(super) struct VoteEntry {
    pub(super) px_x: f64,
    pub(super) px_y: f64,
    pub(super) ra: f64,
    pub(super) dec: f64,
    pub(super) snr: f64,
    votes: usize,
}

#[derive(Debug, Clone, Copy)]
struct BestMatch {
    cat_id: i64,
    entry: VoteEntry,
}

pub(super) fn extract_pairs(matches: &[QuadMatch], sources: &VoteSources<'_>) -> Vec<StarPair> {
    let votes = count_votes(matches, sources);
    let by_image = best_catalog_per_image(&votes);
    let used_cat = claim_catalog_stars(&by_image);

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
    pairs.sort_by(StarPair::cmp_rank);
    pairs
}

// Resolves quad star indices to detections and catalog positions.
pub(super) struct VoteSources<'a> {
    detected: &'a [DetectedStar],
    image_stars: &'a [QuadStar],
    catalog_stars: &'a [QuadStar],
    ra_dec: HashMap<i64, (f64, f64)>,
}

impl<'a> VoteSources<'a> {
    pub(super) fn new(
        detected: &'a [DetectedStar],
        image_stars: &'a [QuadStar],
        catalog_stars: &'a [QuadStar],
        catalog_results: &[ConeSearchResult],
    ) -> Self {
        Self {
            detected,
            image_stars,
            catalog_stars,
            ra_dec: catalog_results
                .iter()
                .map(|r| (r.star.source_id, (r.ra_deg, r.dec_deg)))
                .collect(),
        }
    }

    // The (image id, catalog id) key and a zero-vote entry for one quad slot, or None
    // when an index is out of range or the catalog star has no cone-search position.
    pub(super) fn slot_entry(
        &self,
        img_idx: usize,
        cat_idx: usize,
    ) -> Option<((i64, i64), VoteEntry)> {
        let img_star = self.image_stars.get(img_idx)?;
        let cat_id = self.catalog_stars.get(cat_idx)?.source_id;
        let &(ra, dec) = self.ra_dec.get(&cat_id)?;
        let det = self
            .detected
            .get(usize::try_from(img_star.source_id).ok()?)?;
        let entry = VoteEntry {
            px_x: det.x,
            px_y: det.y,
            ra,
            dec,
            snr: det.snr,
            votes: 0,
        };
        Some(((img_star.source_id, cat_id), entry))
    }
}

// BTreeMap keeps every later step independent of hasher state: iteration runs in
// (image id, catalog id) order.
fn count_votes(
    matches: &[QuadMatch],
    sources: &VoteSources<'_>,
) -> BTreeMap<(i64, i64), VoteEntry> {
    let mut votes = BTreeMap::new();
    for m in matches {
        let slots = m
            .image_quad
            .star_indices
            .iter()
            .zip(&m.catalog_quad.star_indices);
        for (&img_idx, &cat_idx) in slots {
            if let Some((key, fresh)) = sources.slot_entry(img_idx, cat_idx) {
                votes.entry(key).or_insert(fresh).votes += 1;
            }
        }
    }
    votes
}

// For each image star, the catalog star with the most votes; ties go to the lower
// catalog id.
fn best_catalog_per_image(votes: &BTreeMap<(i64, i64), VoteEntry>) -> BTreeMap<i64, BestMatch> {
    let mut by_image: BTreeMap<i64, BestMatch> = BTreeMap::new();
    for (&(img_id, cat_id), entry) in votes {
        let candidate = BestMatch {
            cat_id,
            entry: *entry,
        };
        by_image
            .entry(img_id)
            .and_modify(|best| {
                if entry.votes > best.entry.votes
                    || (entry.votes == best.entry.votes && cat_id < best.cat_id)
                {
                    *best = candidate;
                }
            })
            .or_insert(candidate);
    }
    by_image
}

// For each catalog star, the image star that claims it with the most votes; ties go
// to the lower image id. Values are (image id, votes).
fn claim_catalog_stars(by_image: &BTreeMap<i64, BestMatch>) -> BTreeMap<i64, (i64, usize)> {
    let mut used_cat: BTreeMap<i64, (i64, usize)> = BTreeMap::new();
    for (&img_id, best) in by_image {
        used_cat
            .entry(best.cat_id)
            .and_modify(|e| {
                if best.entry.votes > e.1 || (best.entry.votes == e.1 && img_id < e.0) {
                    *e = (img_id, best.entry.votes);
                }
            })
            .or_insert((img_id, best.entry.votes));
    }
    used_cat
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_catalog::query::cone::ConeSearchResult;
    use celestial_catalog::query::quad::{Quad, QuadStar};

    use crate::detect::DetectedStar;

    fn det(x: f64, y: f64, snr: f64) -> DetectedStar {
        DetectedStar {
            x,
            y,
            flux: 1000.0,
            snr,
            saturated: false,
            saturated_count: 0,
            background: 0.0,
        }
    }

    fn img_star(id: i64) -> QuadStar {
        QuadStar {
            source_id: id,
            x: 0.0,
            y: 0.0,
            mag: 0.0,
        }
    }

    fn cat_star(id: i64) -> QuadStar {
        QuadStar {
            source_id: id,
            x: 0.0,
            y: 0.0,
            mag: 0.0,
        }
    }

    fn cone_result(source_id: i64, ra: f64, dec: f64) -> ConeSearchResult {
        let mut star = celestial_catalog::query::catalog::StarRecord::from_le_bytes(&[0; 56]);
        star.source_id = source_id;
        star.ra = ra;
        star.dec = dec;
        star.mag = 10.0;
        ConeSearchResult {
            star,
            ra_deg: ra,
            dec_deg: dec,
            distance_deg: 0.0,
        }
    }

    fn quad_match(img_idx: [usize; 4], cat_idx: [usize; 4]) -> QuadMatch {
        QuadMatch {
            image_quad: Quad {
                hash: [0.0; 4],
                star_indices: img_idx,
            },
            catalog_quad: Quad {
                hash: [0.0; 4],
                star_indices: cat_idx,
            },
        }
    }

    #[test]
    fn empty_matches_produce_empty_pairs() {
        let detected = vec![det(0.0, 0.0, 10.0)];
        let img = vec![img_star(0)];
        let cat = vec![cat_star(1000)];
        let res = vec![cone_result(1000, 10.0, 20.0)];
        let pairs = extract_pairs(&[], &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
        // Img 0 had two catalog candidates; by_image's and_modify branch picks 8000 (3 votes).
        let img0_pair = pairs.iter().find(|p| p.ra_deg == 1.0);
        assert!(
            img0_pair.is_some(),
            "img 0 must claim catalog 8000 (3 votes)"
        );
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
        let pairs = extract_pairs(&m, &VoteSources::new(&detected, &img, &cat, &res));
        // Img 0 and img 1 tied at 2 votes for cat 7000; lower id wins → only one pair at ra=5.0.
        let cat7000 = pairs.iter().filter(|p| p.ra_deg == 5.0).count();
        assert_eq!(cat7000, 1);
    }

    type PairBits = (u64, u64, u64, u64, usize, u64);

    fn pair_bits(pairs: &[StarPair]) -> Vec<PairBits> {
        pairs
            .iter()
            .map(|p| {
                (
                    p.px_x.to_bits(),
                    p.px_y.to_bits(),
                    p.ra_deg.to_bits(),
                    p.dec_deg.to_bits(),
                    p.votes,
                    p.snr.to_bits(),
                )
            })
            .collect()
    }

    // 16 image stars in four disjoint quads, each quad matched twice, so every pair
    // has 2 votes. Stars 0..8 share snr 50; stars 8..16 have distinct snr.
    struct TiedVotes {
        detected: Vec<DetectedStar>,
        img: Vec<QuadStar>,
        cat: Vec<QuadStar>,
        res: Vec<ConeSearchResult>,
        matches: Vec<QuadMatch>,
    }

    fn tied_votes() -> TiedVotes {
        let detected = (0..16)
            .map(|i| {
                let snr = if i < 8 { 50.0 } else { 10.0 + i as f64 };
                det(100.0 - i as f64, 3.0 * i as f64, snr)
            })
            .collect();
        let matches = (0..4)
            .flat_map(|q| {
                let idx = [4 * q, 4 * q + 1, 4 * q + 2, 4 * q + 3];
                [quad_match(idx, idx), quad_match(idx, idx)]
            })
            .collect();
        TiedVotes {
            detected,
            img: (0..16).map(|i| img_star(i as i64)).collect(),
            cat: (0..16).map(|i| cat_star(1000 + i as i64)).collect(),
            res: (0..16)
                .map(|i| cone_result(1000 + i as i64, i as f64, -(i as f64)))
                .collect(),
            matches,
        }
    }

    #[test]
    fn extract_pairs_is_identical_across_calls() {
        let t = tied_votes();
        let first = extract_pairs(
            &t.matches,
            &VoteSources::new(&t.detected, &t.img, &t.cat, &t.res),
        );
        for _ in 0..20 {
            let again = extract_pairs(
                &t.matches,
                &VoteSources::new(&t.detected, &t.img, &t.cat, &t.res),
            );
            assert_eq!(pair_bits(&again), pair_bits(&first));
        }
    }

    #[test]
    fn extract_pairs_breaks_vote_ties_by_snr_then_pixel_position() {
        let t = tied_votes();
        let pairs = extract_pairs(
            &t.matches,
            &VoteSources::new(&t.detected, &t.img, &t.cat, &t.res),
        );
        // Stars 7..0 (snr 50, ascending px_x), then stars 15..8 (snr 25..18).
        let expected: Vec<f64> = (0..8)
            .rev()
            .chain((8..16).rev())
            .map(|i| i as f64)
            .collect();
        let ra: Vec<f64> = pairs.iter().map(|p| p.ra_deg).collect();
        assert_eq!(ra, expected);
    }
}
