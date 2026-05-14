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
    pairs.sort_by(|a, b| b.votes.cmp(&a.votes));
    pairs
}
