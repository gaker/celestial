//! The search for a verified quad match under one sign of det(CD), and the choice
//! between the two signs.

use celestial_catalog::query::quad::{neighbor_quads, Quad, QuadStar};
use rayon::prelude::*;

use super::aggregate::{extract_pairs, VoteSources};
use super::catalog_side::CatalogSide;
use super::hash_index::find_matches;
use super::quad_stars::{stars_to_quad_stars, CdSign};
use super::verify::{verify_quad_matches, Accepted, FirstPass, ParityVerdict};
use super::{MatchParams, QuadMatch};
use crate::detect::DetectedStar;
use crate::error::SolveError;
use crate::log_odds::field::Field;
use crate::log_odds::ACCEPT_LOG_ODDS;

// Everything one attempt reads; shared by the attempts for both signs of det(CD).
pub(super) struct ParitySearch<'a> {
    pub(super) stars: &'a [DetectedStar],
    pub(super) image_size: (usize, usize),
    pub(super) scale_arcsec: f64,
    pub(super) params: &'a MatchParams,
    pub(super) catalog: &'a CatalogSide,
    pub(super) field: &'a Field<'a>,
}

impl ParitySearch<'_> {
    // The attempts for both signs of det(CD), run in parallel, and the hypothesis that
    // `choose_parity` accepts between them.
    pub(super) fn run(&self) -> Result<(ParityAttempt, Accepted), SolveError> {
        // The two attempts verify at the same time; the first hypothesis to pass in the
        // order that alternates between them is accepted.
        let first_pass = FirstPass::new(2);
        let attempts: Vec<ParityAttempt> = [CdSign::Positive, CdSign::Negative]
            .par_iter()
            .enumerate()
            .map(|(attempt, &cd_sign)| self.attempt(&first_pass, attempt, cd_sign))
            .collect();
        let (attempt, accepted) = choose_parity(attempts)?;
        log::debug!(
            "match_field: accepted under {} at position {}, model det(CD) {} 0, log-odds {:.2}, {} pairs",
            attempt.cd_sign.label(),
            accepted.position,
            if accepted.flipped { "<" } else { ">" },
            accepted.log_odds,
            accepted.pairs.len()
        );
        Ok((attempt, accepted))
    }

    // Attempt number `attempt` of `first_pass`, for image quads built for `cd_sign`.
    fn attempt(&self, first_pass: &FirstPass, attempt: usize, cd_sign: CdSign) -> ParityAttempt {
        let (img_stars, img_quads) = self.image_quads(cd_sign);
        let catalog = self.catalog;
        let matches = find_matches(&img_quads, &catalog.quads, &catalog.index);
        let sources = VoteSources::new(
            self.stars,
            &img_stars,
            &catalog.quad_stars,
            &catalog.results,
        );
        let pairs = extract_pairs(&matches, &sources);
        let candidates = (matches.as_slice(), &sources, pairs.as_slice());
        let verdict = verify_quad_matches(self.field, candidates, first_pass, attempt);
        log::debug!(
            "match_field [{}]: {} img quads, {} raw matches, {} pairs, {}",
            cd_sign.label(),
            img_quads.len(),
            matches.len(),
            pairs.len(),
            describe(&verdict)
        );
        ParityAttempt {
            cd_sign,
            img_stars,
            matches,
            verdict,
        }
    }

    // The image quad stars for `cd_sign` and the quads built from them.
    fn image_quads(&self, cd_sign: CdSign) -> (Vec<QuadStar>, Vec<Quad>) {
        let img_stars = stars_to_quad_stars(
            self.stars,
            self.image_size,
            self.scale_arcsec,
            self.params,
            cd_sign,
        );
        let min_spine_sq = self.catalog.min_spine_sq;
        let img_quads = neighbor_quads(&img_stars, self.params.k_neighbors, min_spine_sq);
        (img_stars, img_quads)
    }
}

pub(super) struct ParityAttempt {
    pub(super) cd_sign: CdSign,
    pub(super) img_stars: Vec<QuadStar>,
    pub(super) matches: Vec<QuadMatch>,
    verdict: ParityVerdict,
}

fn describe(verdict: &ParityVerdict) -> String {
    if verdict.n_hypotheses == 0 {
        return "0 hypotheses".to_string();
    }
    format!(
        "best log-odds {:.2} over {} hypotheses",
        verdict.best_log_odds, verdict.n_hypotheses
    )
}

// Accepts the hypothesis that passed first in `FirstPass` order. As in astrometry.net,
// which stops at the first hypothesis to pass in either parity, the other attempt's
// later hypotheses are not verified. A second attempt can still pass before it learns of
// the first; its hypothesis comes later in that order and is not used.
fn choose_parity(attempts: Vec<ParityAttempt>) -> Result<(ParityAttempt, Accepted), SolveError> {
    let summary = summarize(&attempts);
    let first = attempts
        .into_iter()
        .filter_map(|mut a| a.verdict.accepted.take().map(|acc| (a, acc)))
        .min_by_key(|(_, acc)| acc.position);
    first.ok_or_else(|| {
        SolveError::no_match(format!(
            "no quad match passed verification: {summary} (acceptance needs log-odds above {ACCEPT_LOG_ODDS:.2}, odds of 1e9)",
        ))
    })
}

fn summarize(attempts: &[ParityAttempt]) -> String {
    let parts: Vec<String> = attempts
        .iter()
        .map(|a| format!("{}: {}", a.cd_sign.label(), describe(&a.verdict)))
        .collect();
    parts.join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn passed(cd_sign: CdSign, position: usize, n_pairs: usize) -> ParityAttempt {
        let pair = crate::match_field::StarPair {
            px_x: 0.0,
            px_y: 0.0,
            ra_deg: 0.0,
            dec_deg: 0.0,
            votes: 1,
            snr: 10.0,
        };
        ParityAttempt {
            cd_sign,
            img_stars: Vec::new(),
            matches: Vec::new(),
            verdict: ParityVerdict {
                n_hypotheses: position / 2 + 1,
                best_log_odds: 30.0,
                accepted: Some(Accepted {
                    log_odds: 30.0,
                    position,
                    flipped: cd_sign == CdSign::Negative,
                    pairs: vec![pair; n_pairs],
                }),
            },
        }
    }

    // The attempt with fewer pairs passed first, so it is the one accepted.
    #[test]
    fn choose_parity_takes_the_first_hypothesis_to_pass() {
        let attempts = vec![
            passed(CdSign::Positive, 6, 3),
            passed(CdSign::Negative, 3, 1),
        ];
        let (attempt, accepted) = choose_parity(attempts).unwrap();
        assert_eq!((attempt.cd_sign, accepted.position), (CdSign::Negative, 3));
    }
}
