//! Bayesian verification of a WCS hypothesis (Lang et al. 2010, AJ 139, 1782, sec. 2.4).
//!
//! A hypothesis is an affine map from the tangent plane to pixels, least-squares fitted
//! to a few matched pairs. It predicts where the catalog reference stars fall in the
//! image. A fitted WCS can place them instead (`Hypothesis::of_wcs`); the affine fit to
//! the WCS's pairs then gives only the variance of each placement. Detected stars are
//! then scored in descending flux under two models:
//!
//! - foreground (the hypothesis is true): with probability `DISTRACTOR_FRACTION` a
//!   detection has no counterpart and falls uniformly over the image; otherwise it is
//!   drawn from a Gaussian around one of the `N_R` reference stars in the image, chosen
//!   with equal weight;
//! - background (the hypothesis is false): every detection falls uniformly over the
//!   image of area `A`.
//!
//! Each detection adds `ln(d + (1 - d) * A * f(x) / N_R)` to the log-odds, where `f` is
//! the sum of the reference Gaussians at the detection. The result is the highest
//! running log-odds over the flux-ordered detections, as the paper stops at the first
//! detection that takes the odds over the threshold.
//!
//! The foreground is a normalised density over the image whatever `d` and the Gaussian
//! widths are, so under the background model the running odds is a nonnegative
//! supermartingale with mean at most 1. By Ville's inequality a false hypothesis reaches
//! odds `K` with probability at most `1 / K`. The constants below set how readily a true
//! hypothesis is accepted; they do not change that bound. Cutting the Gaussians off at
//! `TRUNCATION_SIGMAS` only lowers the foreground density, and stopping a scan early only
//! lowers the result, so neither can raise the false-accept rate either.
//!
//! Detections used to fit the hypothesis say nothing about whether it is chance, so
//! they are either left out of the score (`FitStars::Excluded`) or scored against the
//! model refitted without them (`FitStars::LeftOut`).

mod affine;
pub(crate) mod field;
mod fitted;
mod grid;

use celestial_catalog::query::quad::tan_project_star;
use celestial_core::constants::TWOPI;

use crate::detect::DetectedStar;
use crate::fit_wcs::solution::WcsSolution;
use crate::match_field::StarPair;

use affine::{AffineFit, Point};
use field::{place_references, Field, Placed, Placement};
use fitted::FittedPairs;
use grid::RefGrid;

// ln(1e9). Lang et al. 2010 eqs. 11-13: prior odds p(F)/p(B) = 1e-6 and utilities
// u(TP) = u(TN) = +1, u(FN) = -1, u(FP) = -1999 accept when the Bayes factor exceeds 1e9.
pub(crate) const ACCEPT_LOG_ODDS: f64 = 20.723_265_836_946_41;

// ln(1e-100). A scan whose running log-odds falls below this stops: the hypothesis is
// rejected. DEFAULT_BAIL_THRESHOLD in astrometry.net 0.97 (solver.h).
const BAIL_LOG_ODDS: f64 = -230.258_509_299_404_57;

// Fraction of detections with no catalog counterpart under the foreground model. The
// paper describes the mixture without giving a value; 0.25 is DEFAULT_DISTRACTOR_RATIO
// in astrometry.net 0.97 (solver.h), the authors' implementation, which documents it
// as always kept at 0.25 (tweak2.h).
const DISTRACTOR_FRACTION: f64 = 0.25;

// Per-axis error, in pixels, of a detected position relative to the true position of
// its catalog star. 1 px is DEFAULT_VERIFY_PIX in astrometry.net 0.97 (solver.h). The
// fit residuals of the five reference images after refinement are 0.9 to 1.1 px RMS.
pub(crate) const POSITION_SIGMA_PX: f64 = 1.0;

// Lang et al. 2010 sec. 3.1 keep the brightest 300 detections of each image and use all
// of them to check hypotheses.
const MAX_SCORED_STARS: usize = 300;

// The brightest catalog stars inside the image under a hypothesis are its references,
// as many as the detections scored, so that each scored detection can have its own
// counterpart.
pub(crate) const MAX_REFERENCES: usize = MAX_SCORED_STARS;

// A reference Gaussian is not evaluated beyond this many standard deviations: the mass
// it leaves out is exp(-12.5), about 4e-6.
const TRUNCATION_SIGMAS: f64 = 5.0;

// Which detections the score leaves to the fitted pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FitStars {
    // Fitted detections are not scored. Use when the fitted pairs were chosen for their
    // mutual agreement, as the four stars of a hash-matched quad are.
    Excluded,
    // Each fitted detection is scored against the model fitted to the other pairs, and
    // is not scored when the other pairs do not determine a model.
    LeftOut,
}

// A hypothesis fitted to `fit` pairs on the field's tangent plane, with the reference
// stars it places inside the image.
pub(crate) struct Hypothesis<'a> {
    // The affine fit to the pairs. It places the references unless a WCS does; its
    // leverage always gives their prediction variance.
    fit: AffineFit,
    placement: Placement<'a>,
    fitted: FittedPairs,
    // The `MAX_REFERENCES` brightest references whose predicted pixel is in the image.
    placed: Vec<Placed>,
    grid: RefGrid,
    field: &'a Field<'a>,
}

impl<'a> Hypothesis<'a> {
    // `None` when the pairs with a tangent-plane projection do not determine an affine
    // model: fewer than 3, or collinear.
    pub(crate) fn new(field: &'a Field<'a>, fit: &[StarPair]) -> Option<Self> {
        Self::placed_by(field, fit, Placement::Affine)
    }

    // The hypothesis that `wcs`, fitted to `fit`, is the solution: references go where
    // `wcs` projects them, with the prediction variance of the affine fit to `fit`.
    // `None` as for `new`.
    pub(crate) fn of_wcs(
        field: &'a Field<'a>,
        fit: &[StarPair],
        wcs: &'a WcsSolution,
    ) -> Option<Self> {
        Self::placed_by(field, fit, Placement::Wcs(wcs))
    }

    fn placed_by(field: &'a Field<'a>, fit: &[StarPair], placement: Placement<'a>) -> Option<Self> {
        let fitted = FittedPairs::new(fit, field.tangent);
        let model = fitted.affine_fit()?;
        let placed = place_references(field, &model, placement);
        let disks: Vec<(Point, f64)> = placed
            .iter()
            .map(|r| (r.pixel, truncation_radius(r.variance)))
            .collect();
        Some(Self {
            fit: model,
            placement,
            fitted,
            grid: RefGrid::new(field.width, field.height, &disks),
            placed,
            field,
        })
    }

    // `det(CD) < 0`, the convention of `WcsSolution::is_flipped`.
    pub(crate) fn is_flipped(&self) -> bool {
        self.fit.determinant() < 0.0
    }

    // The highest running log-odds over the scored detections in descending flux. The
    // scan ends when the running log-odds falls below `BAIL_LOG_ODDS` and, with
    // `stop_at_accept`, at the first detection that takes it above `ACCEPT_LOG_ODDS`.
    // 0 when no detection is scored.
    pub(crate) fn log_odds(&self, fit_stars: FitStars, stop_at_accept: bool) -> f64 {
        let mut sum = 0.0;
        let mut best = f64::NEG_INFINITY;
        for &si in &self.field.scored {
            let Some(term) = self.score(&self.field.stars[si], fit_stars) else {
                continue;
            };
            sum += term;
            best = best.max(sum);
            if (stop_at_accept && best > ACCEPT_LOG_ODDS) || sum < BAIL_LOG_ODDS {
                break;
            }
        }
        if best == f64::NEG_INFINITY {
            0.0
        } else {
            best
        }
    }

    // True when the pair's detection, scored against the pair's own catalog star alone,
    // raises the odds that the hypothesis is true.
    pub(crate) fn supports(&self, pair: &StarPair) -> bool {
        let sky = (pair.ra_deg, pair.dec_deg);
        let Some(q) = project(sky.0, sky.1, self.field.tangent) else {
            return false;
        };
        let Some(pixel) = self.placement.pixel(&self.fit, sky, q) else {
            return false;
        };
        let x = [pair.px_x, pair.px_y];
        let density = gaussian(x, pixel, 1.0 + self.fit.leverage(q));
        self.log_ratio(density) > 0.0
    }

    // Pairs each placed reference with the scored detection it supports most, when
    // that detection supports the hypothesis by `supports`'s rule. Each detection and
    // each reference is used once; a detection goes to the reference of highest density.
    pub(crate) fn matched_pairs(&self) -> Vec<StarPair> {
        let mut best_for_ref: Vec<Option<(usize, f64)>> = vec![None; self.placed.len()];
        for &si in &self.field.scored {
            let star = &self.field.stars[si];
            let Some((k, density)) = self.nearest_reference([star.x, star.y]) else {
                continue;
            };
            let slot = &mut best_for_ref[k];
            if self.log_ratio(density) > 0.0 && slot.is_none_or(|(_, d)| density > d) {
                *slot = Some((si, density));
            }
        }
        let pairs = best_for_ref.iter().zip(&self.placed);
        pairs
            .filter_map(|(slot, r)| slot.map(|(si, _)| star_pair(&self.field.stars[si], r.sky)))
            .collect()
    }

    // The placed reference of highest density at `x`, and that density.
    fn nearest_reference(&self, x: Point) -> Option<(usize, f64)> {
        let densities = self.grid.candidates(x).iter().map(|&k| {
            let r = &self.placed[k];
            (k, gaussian(x, r.pixel, r.variance))
        });
        densities
            .filter(|&(_, d)| d > 0.0)
            .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))
    }

    fn score(&self, star: &DetectedStar, fit_stars: FitStars) -> Option<f64> {
        let x = [star.x, star.y];
        let Some(fit_index) = self.fitted.find(x) else {
            let density = self
                .grid
                .candidates(x)
                .iter()
                .map(|&k| gaussian(x, self.placed[k].pixel, self.placed[k].variance))
                .sum();
            return Some(self.log_ratio(density));
        };
        match fit_stars {
            FitStars::Excluded => None,
            FitStars::LeftOut => self.score_left_out(x, fit_index),
        }
    }

    // Scores a fitted detection against the model refitted without its pair, using the
    // rank-one downdate of the affine least-squares solution. With a WCS placement the
    // downdate moves the WCS's predictions; it is exact when the WCS equals the affine
    // fit and first order in the difference otherwise.
    fn score_left_out(&self, x: Point, fit_index: usize) -> Option<f64> {
        let left_out = self.fitted.left_out(fit_index, &self.fit, self.placement)?;
        let density = self
            .placed
            .iter()
            .map(|r| {
                let (offset, extra) = left_out.correction(&self.fit, r.q);
                let pred = [r.pixel[0] + offset[0], r.pixel[1] + offset[1]];
                gaussian(x, pred, r.variance + extra)
            })
            .sum();
        Some(self.log_ratio(density))
    }

    // ln(p(x | foreground) / p(x | background)) for a detection whose summed reference
    // Gaussian density is `density`.
    fn log_ratio(&self, density: f64) -> f64 {
        let n_refs = self.placed.len();
        if n_refs == 0 {
            return 0.0;
        }
        let d = DISTRACTOR_FRACTION;
        libm::log(d + (1.0 - d) * self.field.area() * density / n_refs as f64)
    }
}

fn star_pair(star: &DetectedStar, (ra_deg, dec_deg): (f64, f64)) -> StarPair {
    StarPair {
        px_x: star.x,
        px_y: star.y,
        ra_deg,
        dec_deg,
        votes: 1,
        snr: star.snr,
    }
}

fn project(ra: f64, dec: f64, tangent: (f64, f64)) -> Option<Point> {
    tan_project_star(ra, dec, tangent.0, tangent.1).map(|(xi, eta)| [xi, eta])
}

fn truncation_radius(variance_factor: f64) -> f64 {
    TRUNCATION_SIGMAS * POSITION_SIGMA_PX * libm::sqrt(variance_factor)
}

// Isotropic Gaussian density at `x` around `mean`, with per-axis variance
// `POSITION_SIGMA_PX^2 * variance_factor`; 0 beyond `TRUNCATION_SIGMAS` deviations and
// for a non-finite distance, so a NaN position adds no evidence.
fn gaussian(x: Point, mean: Point, variance_factor: f64) -> f64 {
    let var = POSITION_SIGMA_PX * POSITION_SIGMA_PX * variance_factor;
    let dx = x[0] - mean[0];
    let dy = x[1] - mean[1];
    let d2 = dx * dx + dy * dy;
    if !d2.is_finite() || d2 > TRUNCATION_SIGMAS * TRUNCATION_SIGMAS * var {
        return 0.0;
    }
    libm::exp(-d2 / (2.0 * var)) / (TWOPI * var)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::triple_search::index_triples;
    use crate::triple_search::test_scene::minority_scene;

    const TANGENT: (f64, f64) = (180.0, 0.0);
    const SIZE: usize = 2000;

    // `minority_scene` pixels lie in [-1000, 1000); this moves them into a 2000 x 2000
    // image.
    fn scene(n_good: usize, n_junk: usize) -> Vec<StarPair> {
        minority_scene(n_good, n_junk)
            .into_iter()
            .map(|p| StarPair {
                px_x: p.px_x + 1000.0,
                px_y: p.px_y + 1000.0,
                ..p
            })
            .collect()
    }

    // One detection per pair, in descending flux in pair order.
    fn detections(pairs: &[StarPair]) -> Vec<DetectedStar> {
        pairs
            .iter()
            .enumerate()
            .map(|(i, p)| DetectedStar {
                x: p.px_x,
                y: p.px_y,
                flux: 1000.0 - i as f64,
                snr: p.snr,
                saturated: false,
                saturated_count: 0,
                background: 0.0,
            })
            .collect()
    }

    fn references(pairs: &[StarPair]) -> Vec<(f64, f64)> {
        pairs.iter().map(|p| (p.ra_deg, p.dec_deg)).collect()
    }

    #[test]
    fn accept_threshold_is_ln_of_1e9() {
        assert_eq!(ACCEPT_LOG_ODDS, libm::log(1.0e9));
    }

    #[test]
    fn bail_threshold_is_ln_of_1e_minus_100() {
        assert_eq!(BAIL_LOG_ODDS, libm::log(1.0e-100));
    }

    #[test]
    fn only_the_brightest_detections_are_scored() {
        let pairs = scene(0, 400);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let expected: Vec<usize> = (0..MAX_SCORED_STARS).collect();
        assert_eq!(field.scored, expected);
    }

    #[test]
    fn references_are_the_brightest_the_hypothesis_places_in_the_image() {
        // The first five pairs fix the true model; every pair's sky position then lands
        // in the image, and the first `MAX_REFERENCES` in the given order are kept.
        let pairs = scene(5, 400);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[..5]).unwrap();
        let placed: Vec<(f64, f64)> = h.placed.iter().map(|r| r.sky).collect();
        assert_eq!(placed, references(&pairs)[..MAX_REFERENCES]);
    }

    #[test]
    fn matched_pairs_join_each_confirming_detection_to_its_reference() {
        let pairs = scene(5, 45);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[..3]).unwrap();
        let key = |p: &StarPair| (p.px_x, p.px_y, p.ra_deg, p.dec_deg, p.votes);
        let matched: Vec<_> = h.matched_pairs().iter().map(key).collect();
        let expected: Vec<_> = pairs[..5]
            .iter()
            .map(|p| StarPair {
                votes: 1,
                ..p.clone()
            })
            .map(|p| key(&p))
            .collect();
        assert_eq!(matched, expected);
    }

    #[test]
    fn a_detection_at_nan_adds_a_distractor_term() {
        let pairs = scene(5, 0);
        let mut stars = detections(&pairs);
        stars[3].x = f64::NAN;
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[..3]).unwrap();
        // Star 3 scores ln(0.25); star 4 then confirms the model.
        assert_eq!(h.log_odds(FitStars::Excluded, false), 7.936865369176321);
    }

    #[test]
    fn hypothesis_needs_three_pairs() {
        let pairs = scene(2, 0);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        assert!(Hypothesis::new(&field, &pairs).is_none());
    }

    #[test]
    fn no_triple_of_unrelated_pairs_passes() {
        // Every one of the C(40, 3) = 9880 models, each verified on the other 37 pairs.
        let pairs = scene(0, 40);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let best = index_triples(pairs.len())
            .filter_map(|t| {
                let fit = t.map(|i| pairs[i].clone());
                Hypothesis::new(&field, &fit)
            })
            .map(|h| h.log_odds(FitStars::Excluded, false))
            .fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(best, 7.178196018187171);
    }

    #[test]
    fn exactly_fitted_pairs_carry_no_evidence() {
        // Three pairs fix the affine model, so a solution resting on three pairs has no
        // support whatever they are.
        let pairs = scene(0, 3);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs).unwrap();
        assert_eq!(h.log_odds(FitStars::LeftOut, false), 0.0);
        assert_eq!(h.log_odds(FitStars::Excluded, false), 0.0);
    }

    #[test]
    fn five_true_pairs_pass_left_out_among_unrelated_ones() {
        let pairs = scene(5, 45);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[..5]).unwrap();
        assert_eq!(h.log_odds(FitStars::LeftOut, false), 40.58049496497802);
        assert!(!h.is_flipped());
    }

    #[test]
    fn a_true_triple_needs_more_than_two_confirming_stars() {
        // Only pairs 3 and 4 confirm the model of pairs 0 to 2.
        let pairs = scene(5, 45);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[..3]).unwrap();
        assert_eq!(h.log_odds(FitStars::Excluded, false), 13.27950577689512);
        let supported: Vec<bool> = pairs.iter().map(|p| h.supports(p)).collect();
        assert_eq!(supported[..5], [true; 5]);
        assert_eq!(supported[5..], [false; 45]);
    }

    #[test]
    fn an_unrelated_triple_peaks_at_its_first_distractor() {
        let pairs = scene(5, 45);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[5..8]).unwrap();
        assert_eq!(h.log_odds(FitStars::Excluded, false), libm::log(0.25));
    }

    #[test]
    fn scan_stops_at_the_first_detection_past_the_threshold() {
        let pairs = scene(8, 0);
        let stars = detections(&pairs);
        let field = Field::new(&stars, (SIZE, SIZE), TANGENT, &references(&pairs));
        let h = Hypothesis::new(&field, &pairs[..3]).unwrap();
        // Five stars confirm the model; the third takes it past ln(1e9).
        assert_eq!(h.log_odds(FitStars::Excluded, true), 27.146421757590627);
        assert_eq!(h.log_odds(FitStars::Excluded, false), 45.51419889497292);
    }
}
