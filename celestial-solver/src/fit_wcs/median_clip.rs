// Clipping of the pairs whose residual is an outlier against the median residual.

use crate::detect::stats::median_in_place;
use crate::error::SolveError;
use crate::log_odds::POSITION_SIGMA_PX;

// For a residual size e = sqrt(dx^2 + dy^2) with independent Gaussian errors of standard
// deviation σ on each axis, e follows a Rayleigh distribution, whose median is
// σ·sqrt(2 ln 2). Dividing the median residual by this factor estimates σ.
fn rayleigh_median_per_sigma() -> f64 {
    libm::sqrt(2.0 * libm::log(2.0))
}

// The residual above which a pair is clipped: 3 times the per-axis σ estimated from the
// median residual of the active pairs, and never below `POSITION_SIGMA_PX`, the
// centroid error that no residual below it can tell apart from a correct pair. `None`
// with fewer than 3 active pairs.
fn clip_threshold(errs: &[f64], mask: &[bool]) -> Option<f64> {
    let mut active: Vec<f64> = errs
        .iter()
        .zip(mask)
        .filter(|(_, &on)| on)
        .map(|(&e, _)| e)
        .collect();
    if active.len() < 3 {
        return None;
    }
    let sigma = median_in_place(&mut active)? / rayleigh_median_per_sigma();
    Some((3.0 * sigma).max(POSITION_SIGMA_PX))
}

// The clip of `median_clip_loop` for any model: `errors` gives the residual in pixels of
// each active pair under a model, and `refit` fits a model to the active pairs.
pub(super) fn median_clip<M>(
    mask: &mut [bool],
    model: &mut M,
    errors: impl Fn(&M, &[bool]) -> Vec<f64>,
    mut refit: impl FnMut(&[bool]) -> Result<M, SolveError>,
) -> Result<(), SolveError> {
    for iter in 0..10 {
        let errs = errors(model, mask);
        let Some(threshold) = clip_threshold(&errs, mask) else {
            break;
        };
        let mut rejected = 0usize;
        for (on, &e) in mask.iter_mut().zip(&errs) {
            if *on && e > threshold {
                *on = false;
                rejected += 1;
            }
        }
        let n_active = mask.iter().filter(|&&m| m).count();
        log::debug!(
            "fit_wcs iter {iter}: {n_active} stars, clip above {threshold:.3}px, rejected {rejected}"
        );
        if rejected == 0 || n_active < 3 {
            break;
        }
        *model = refit(mask)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clip_threshold_is_three_rayleigh_sigmas_of_the_active_residuals() {
        let errs = [1.0, 2.0, 3.0, 4.0, 5.0, 100.0];
        let mask = [true, true, true, true, true, false];
        assert_eq!(clip_threshold(&errs, &mask), Some(7.643896202592171));
    }

    #[test]
    fn clip_threshold_is_at_least_the_position_sigma() {
        let errs = [0.1, 0.2, 0.3];
        assert_eq!(clip_threshold(&errs, &[true; 3]), Some(1.0));
    }

    #[test]
    fn clip_threshold_needs_three_active_residuals() {
        assert_eq!(clip_threshold(&[1.0, 2.0, 3.0], &[true, false, true]), None);
    }
}
