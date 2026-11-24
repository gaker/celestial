//! Verification of quad-match hypotheses for one parity orientation.
//!
//! Each quad match proposes that four detections are four catalog stars. The affine
//! model of those four pairs is verified with `crate::log_odds`, leaving the four
//! detections out of the score, as Lang et al. 2010 (sec. 2.2 and 2.4) verify every
//! matched quad. Matches are tried in the order the hash lookup found them, and the
//! first one whose log-odds passes is accepted.

use std::collections::BTreeSet;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::aggregate::VoteSources;
use super::{QuadMatch, StarPair};
use crate::log_odds::field::Field;
use crate::log_odds::{FitStars, Hypothesis, ACCEPT_LOG_ODDS};

pub(super) struct Accepted {
    pub(super) log_odds: f64,
    // The position of the hypothesis in `FirstPass` order.
    pub(super) position: usize,
    // `det(CD) < 0` for the accepted model, which can differ from the sign the image quads
    // were built for when a quad is nearly symmetric.
    pub(super) flipped: bool,
    // The aggregated pairs that support the accepted model, plus its quad's own pairs.
    pub(super) pairs: Vec<StarPair>,
}

pub(super) struct ParityVerdict {
    // Quad matches with a determined affine model, each verified once.
    pub(super) n_hypotheses: usize,
    pub(super) best_log_odds: f64,
    pub(super) accepted: Option<Accepted>,
}

// The first hypothesis to pass among attempts that verify at the same time, in the order
// that takes hypothesis k of each of the `attempts` in turn: hypothesis k of attempt `a`
// is at position k * attempts + a. An attempt stops before a position later than one
// that has passed. Every hypothesis before the first to pass is therefore verified
// however the attempts are scheduled, and the first to pass does not depend on timing.
pub(super) struct FirstPass {
    attempts: usize,
    // Only the order of this atomic's own writes matters, so its loads and stores are
    // relaxed.
    first: AtomicUsize,
}

impl FirstPass {
    pub(super) fn new(attempts: usize) -> Self {
        Self {
            attempts,
            first: AtomicUsize::new(usize::MAX),
        }
    }

    fn position(&self, attempt: usize, k: usize) -> usize {
        k.saturating_mul(self.attempts).saturating_add(attempt)
    }

    fn follows_a_pass(&self, position: usize) -> bool {
        position > self.first.load(Ordering::Relaxed)
    }

    fn record(&self, position: usize) {
        self.first.fetch_min(position, Ordering::Relaxed);
    }
}

// Verifies `matches` for attempt `attempt` of `first`, until one passes or a hypothesis
// of another attempt passed at an earlier position.
pub(super) fn verify_quad_matches(
    field: &Field<'_>,
    candidates: (&[QuadMatch], &VoteSources<'_>, &[StarPair]),
    first: &FirstPass,
    attempt: usize,
) -> ParityVerdict {
    let (matches, sources, pairs) = candidates;
    let mut verdict = ParityVerdict {
        n_hypotheses: 0,
        best_log_odds: f64::NEG_INFINITY,
        accepted: None,
    };
    let mut seen = BTreeSet::new();
    let fresh = matches
        .iter()
        .filter(|m| seen.insert((m.image_quad.star_indices, m.catalog_quad.star_indices)));
    for m in fresh {
        let position = first.position(attempt, verdict.n_hypotheses);
        if first.follows_a_pass(position) {
            break;
        }
        let quad = quad_pairs(m, sources);
        let Some(hypothesis) = Hypothesis::new(field, &quad) else {
            continue;
        };
        verdict.n_hypotheses += 1;
        let log_odds = hypothesis.log_odds(FitStars::Excluded, true);
        verdict.best_log_odds = f64::max(verdict.best_log_odds, log_odds);
        if log_odds > ACCEPT_LOG_ODDS {
            first.record(position);
            verdict.accepted = Some(Accepted {
                log_odds,
                position,
                flipped: hypothesis.is_flipped(),
                pairs: supporting_pairs(&hypothesis, pairs, quad),
            });
            break;
        }
    }
    verdict
}

fn quad_pairs(m: &QuadMatch, sources: &VoteSources<'_>) -> Vec<StarPair> {
    let slots = m
        .image_quad
        .star_indices
        .iter()
        .zip(&m.catalog_quad.star_indices);
    slots
        .filter_map(|(&img_idx, &cat_idx)| sources.slot_entry(img_idx, cat_idx))
        .map(|(_, e)| StarPair {
            px_x: e.px_x,
            px_y: e.px_y,
            ra_deg: e.ra,
            dec_deg: e.dec,
            votes: 1,
            snr: e.snr,
        })
        .collect()
}

// The pairs the accepted model supports, in `StarPair::cmp_rank` order: the aggregated
// pairs it supports, which keep their votes, then the scored detections it matches to
// its references, then the quad's own pairs. A later pair is kept only when its
// detection and its catalog star are both unused.
fn supporting_pairs(
    hypothesis: &Hypothesis<'_>,
    pairs: &[StarPair],
    quad: Vec<StarPair>,
) -> Vec<StarPair> {
    let mut kept: Vec<StarPair> = pairs
        .iter()
        .filter(|p| hypothesis.supports(p))
        .cloned()
        .collect();
    for p in hypothesis.matched_pairs().into_iter().chain(quad) {
        let taken = kept.iter().any(|k| {
            (k.px_x, k.px_y) == (p.px_x, p.px_y) || (k.ra_deg, k.dec_deg) == (p.ra_deg, p.dec_deg)
        });
        if !taken {
            kept.push(p);
        }
    }
    kept.sort_by(StarPair::cmp_rank);
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_pass_stops_positions_after_the_earliest_pass() {
        let first = FirstPass::new(2);
        assert_eq!(first.position(1, 3), 7);
        assert!(!first.follows_a_pass(7));
        first.record(5);
        assert_eq!(
            (first.follows_a_pass(5), first.follows_a_pass(6)),
            (false, true)
        );
        // A later pass does not move the first one.
        first.record(9);
        assert_eq!(
            (first.follows_a_pass(5), first.follows_a_pass(6)),
            (false, true)
        );
    }
}
