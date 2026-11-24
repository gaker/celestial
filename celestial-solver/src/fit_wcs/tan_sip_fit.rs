//! The TAN WCS and its SIP distortion fitted to matched pairs.
//!
//! In the SIP convention (Shupe et al. 2005, eq. 1) the pixel offset (u, v) from CRPIX
//! has the tangent-plane position CD (u + A(u, v), v + B(u, v)), where A and B have
//! terms of order 2 and above. That position is a polynomial in (u, v): its constant
//! terms are the offset of CRVAL, its linear terms are the CD matrix, and its higher
//! terms are CD times those of A and B. The polynomial is linear in its coefficients, so
//! the CD matrix, the CRVAL offset and the SIP terms are one least-squares fit. A SIP
//! fitted with the CD matrix held fixed cannot return the part of the distortion that
//! the CD fit absorbed, since that part is linear in (u, v).
//!
//! The pairs are weighted by `pair_weight` and clipped by `median_clip` on the residuals
//! of the whole model, so a distortion that a lower order cannot represent is not
//! clipped as an outlier. The order is chosen by `polynomial_fit::fit_best_order`.
//!
//! The SIP solution keeps the CRVAL and CD matrix of that fit, its value and derivative
//! at CRPIX. The TAN WCS is a separate first-order fit to the pairs the clip keeps, the linear
//! map that best fits the whole field, for the scale, rotation and log-odds.

use std::collections::BTreeMap;

use celestial_catalog::query::quad::tan_deproject_star;
use celestial_core::constants::DEG_TO_RAD;

use crate::error::SolveError;
use crate::match_field::StarPair;

use super::cd_matrix::{compute_residuals, rms_error, weighted_rms_error, WcsParams};
use super::median_clip::median_clip;
use super::polynomial_fit::{fit_best_order, fit_order, pixel_errors, Polynomial};
use super::sip::{sip_monomial_terms, SipSolution};
use super::solution::WcsSolution;
use super::tangent_projection::{extract_kept_pairs, project_pairs_to_tangent, TangentProjection};
use super::{check_pair_count, check_projected_count};

// A TAN WCS and, when the chosen order is 2 or above, its SIP distortion.
pub(crate) struct TanSip {
    pub(crate) wcs: WcsSolution,
    pub(crate) sip: Option<SipSolution>,
}

// The largest distance in pixels from CRPIX to an image edge along each axis, at least
// 1. Monomials of (u, v) divided by them lie in [-1, 1] along both axes, which keeps the
// normal equations well conditioned at high orders on images of any shape.
fn monomial_scale(prior: &WcsSolution) -> (f64, f64) {
    let extent = |crpix: f64, size: usize| {
        let far = libm::fabs(size as f64 - crpix);
        libm::fabs(crpix).max(far).max(1.0)
    };
    (
        extent(prior.crpix1, prior.width),
        extent(prior.crpix2, prior.height),
    )
}

// The SIP distortion with its own TAN parameters fitted to `pairs` with `weights`, the
// order chosen up to `max_order` (`None` for TAN only) and the pairs clipped by
// `median_clip` on the residuals of the chosen model, and the TAN WCS fitted alone to
// the pairs that clip keeps. Both are fitted about the CRPIX of `prior`, with CRVAL at
// their fitted sky position of CRPIX.
pub(crate) fn fit_tan_sip(
    pairs: &[StarPair],
    weights: &[f64],
    prior: &WcsSolution,
    max_order: Option<u32>,
) -> Result<TanSip, SolveError> {
    check_pair_count(pairs.len())?;
    let scale = monomial_scale(prior);
    let tangent = (prior.crval1, prior.crval2);
    let (crpix1, crpix2) = (prior.crpix1, prior.crpix2);
    let tp = project_pairs_to_tangent(pairs, weights, tangent.0, tangent.1, crpix1, crpix2);
    check_projected_count(tp.rows.len())?;
    let (order, mask) = clipped_order(&tp, max_order.unwrap_or(1), scale)?;
    let (kept, kept_weights) = extract_kept_pairs(&tp, &mask);
    let kept: Vec<StarPair> = kept.into_iter().cloned().collect();
    let model = Kept {
        pairs: &kept,
        weights: &kept_weights,
        crpix: (crpix1, crpix2),
        scale,
    };
    let wcs = tan_wcs(prior, &model.fit_about(tangent, 1)?);
    let sip = match order {
        1 => None,
        order => Some(sip_solution(&wcs, &model.fit_about(tangent, order)?)),
    };
    Ok(TanSip { wcs, sip })
}

// The order `fit_best_order` chooses up to `max_order` and the pairs `median_clip` keeps
// under the fit of that order.
fn clipped_order(
    tp: &TangentProjection<'_>,
    max_order: u32,
    scale: (f64, f64),
) -> Result<(u32, Vec<bool>), SolveError> {
    let mut mask = vec![true; tp.rows.len()];
    let mut fit = fit_best_order(tp, &mask, max_order, scale)?;
    median_clip(
        &mut mask,
        &mut fit,
        |fit, mask| pixel_errors(fit, tp, mask),
        |mask| fit_best_order(tp, mask, max_order, scale),
    )?;
    Ok((fit.order, mask))
}

// The pairs a fit keeps, with their weights, the CRPIX and the monomial scale.
struct Kept<'a> {
    pairs: &'a [StarPair],
    weights: &'a [f64],
    crpix: (f64, f64),
    scale: (f64, f64),
}

// A polynomial fitted about a tangent point, its CRVAL, and the pairs projected about
// that tangent point.
struct FittedAbout<'a> {
    fit: Polynomial,
    crval: (f64, f64),
    tp: TangentProjection<'a>,
}

impl<'a> Kept<'a> {
    // The polynomial of `order` fitted about `tangent`, then fitted again about its
    // sky position of CRPIX, which becomes its CRVAL. A tangent point away from that
    // position changes the projection by an amount that grows with the separation, so
    // the second fit leaves the CRVAL a second-order distance from its tangent point.
    fn fit_about(&self, tangent: (f64, f64), order: u32) -> Result<FittedAbout<'a>, SolveError> {
        let (first, _) = self.fit_at(tangent, order)?;
        let centre = offset_sky(&first, tangent);
        let (fit, tp) = self.fit_at(centre, order)?;
        let crval = offset_sky(&fit, centre);
        Ok(FittedAbout { fit, crval, tp })
    }

    fn fit_at(
        &self,
        tangent: (f64, f64),
        order: u32,
    ) -> Result<(Polynomial, TangentProjection<'a>), SolveError> {
        let (crpix1, crpix2) = self.crpix;
        let tp = project_pairs_to_tangent(
            self.pairs,
            self.weights,
            tangent.0,
            tangent.1,
            crpix1,
            crpix2,
        );
        let all = vec![true; tp.rows.len()];
        let (fit, _) = fit_order(&tp, &all, order, self.scale).map_err(|e| {
            let n = tp.rows.len();
            SolveError::fit_failed(format!("order {order} fit to the {n} kept pairs: {e}"))
        })?;
        Ok((fit, tp))
    }
}

// The sky position of the polynomial's value at CRPIX, its constant terms, about
// `tangent`.
fn offset_sky(fit: &Polynomial, tangent: (f64, f64)) -> (f64, f64) {
    let (xi0, eta0) = fit.coefficients((0, 0));
    tan_deproject_star(xi0 * DEG_TO_RAD, eta0 * DEG_TO_RAD, tangent.0, tangent.1)
}

// The TAN WCS of a first-order fit, with its residuals over the pairs it was fitted on.
fn tan_wcs(prior: &WcsSolution, tan: &FittedAbout<'_>) -> WcsSolution {
    let [cd1_1, cd1_2, cd2_1, cd2_2] = tan.fit.cd();
    let params = WcsParams {
        crpix1: prior.crpix1,
        crpix2: prior.crpix2,
        crval_ra: tan.crval.0,
        crval_dec: tan.crval.1,
        cd1_1,
        cd1_2,
        cd2_1,
        cd2_2,
    };
    let residuals = compute_residuals(&tan.tp.pairs, params);
    let sizes = || residuals.iter().map(|r| r.err_px);
    WcsSolution {
        crval1: tan.crval.0,
        crval2: tan.crval.1,
        cd1_1,
        cd1_2,
        cd2_1,
        cd2_2,
        n_stars: residuals.len(),
        rms_px: rms_error(sizes()),
        weighted_rms_px: weighted_rms_error(sizes(), &tan.tp.weights),
        residuals,
        ..prior.clone()
    }
}

// The SIP distortion of `model`, CD⁻¹ times its terms of order 2 and above, with its
// own CRVAL and CD matrix and its residuals over the pairs it was fitted on. `tan` is
// the TAN WCS fitted alone to the same pairs.
fn sip_solution(tan: &WcsSolution, model: &FittedAbout<'_>) -> SipSolution {
    let fit = &model.fit;
    let mut a_coeffs = BTreeMap::new();
    let mut b_coeffs = BTreeMap::new();
    for term in sip_monomial_terms(fit.order) {
        let (a, b) = fit.pixel_coefficients(term);
        a_coeffs.insert(term, a);
        b_coeffs.insert(term, b);
    }
    let errors = || model.tp.rows.iter().map(|row| fit.pixel_error(row));
    let [cd1_1, cd1_2, cd2_1, cd2_2] = fit.cd();
    SipSolution {
        crpix1: tan.crpix1,
        crpix2: tan.crpix2,
        crval1: model.crval.0,
        crval2: model.crval.1,
        cd1_1,
        cd1_2,
        cd2_1,
        cd2_2,
        width: tan.width,
        height: tan.height,
        a_order: fit.order,
        b_order: fit.order,
        a_coeffs,
        b_coeffs,
        n_stars: model.tp.rows.len(),
        linear_rms_px: tan.rms_px,
        rms_px: rms_error(errors()),
        weighted_rms_px: weighted_rms_error(errors(), &model.tp.weights),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::pair_weight;
    use super::super::sip::eval_sip_poly;

    fn pairs_on_tangent_plane(
        center_ra: f64,
        center_dec: f64,
        scale_arcsec: f64,
        crpix: f64,
    ) -> Vec<StarPair> {
        let s_deg = scale_arcsec / 3600.0;
        let mut pairs = Vec::new();
        for &(dx, dy) in &[
            (-200.0, -200.0),
            (200.0, -200.0),
            (-200.0, 200.0),
            (200.0, 200.0),
            (-100.0, 50.0),
            (100.0, -50.0),
            (0.0, 0.0),
            (150.0, 150.0),
        ] {
            let px_x = crpix + dx;
            let px_y = crpix + dy;
            let xi_deg = -s_deg * dx;
            let eta_deg = s_deg * dy;
            let (ra, dec) = tan_deproject_star(
                xi_deg * DEG_TO_RAD,
                eta_deg * DEG_TO_RAD,
                center_ra,
                center_dec,
            );
            pairs.push(StarPair {
                px_x,
                px_y,
                ra_deg: ra,
                dec_deg: dec,
                votes: 4,
                snr: 20.0,
            });
        }
        pairs
    }

    #[test]
    fn fit_tan_sip_recovers_crval_from_clean_pairs() {
        let center_ra = 180.0;
        let center_dec = 30.0;
        let scale = 1.5;
        let prior = WcsSolution::north_up((center_ra, center_dec), scale, (1024, 1024));
        let pairs = pairs_on_tangent_plane(center_ra, center_dec, scale, 512.0);
        let weights: Vec<f64> = pairs.iter().map(|p| pair_weight(p.snr)).collect();

        let wcs = fit_tan_sip(&pairs, &weights, &prior, None)
            .expect("clean pairs must produce a WCS")
            .wcs;
        // The 8 pairs were placed about CRVAL (180, 30) with no noise. The fitted Dec
        // is 3.6e-15 deg from it and the RMS is 2.0e-11 px.
        assert_eq!((wcs.crval1, wcs.crval2), (center_ra, 29.999999999999996));
        assert_eq!((wcs.n_stars, wcs.rms_px), (8, 2.0032315610125392e-11));
    }

    #[test]
    fn fit_tan_sip_ignores_pairs_behind_the_tangent_plane() {
        let prior = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let clean = pairs_on_tangent_plane(180.0, 30.0, 1.5, 512.0);
        // More than 90 degrees from the tangent point, so it has no gnomonic projection.
        let mut pairs = vec![StarPair {
            px_x: 1.0,
            px_y: 1.0,
            ra_deg: 0.0,
            dec_deg: -30.0,
            votes: 4,
            snr: 20.0,
        }];
        pairs.extend(clean.iter().cloned());
        let weights: Vec<f64> = pairs.iter().map(|p| pair_weight(p.snr)).collect();

        let wcs = fit_tan_sip(&pairs, &weights, &prior, None).unwrap().wcs;
        let fitted: Vec<(f64, f64)> = wcs.residuals.iter().map(|r| (r.px_x, r.px_y)).collect();
        let expected: Vec<(f64, f64)> = clean.iter().map(|p| (p.px_x, p.px_y)).collect();
        assert_eq!(fitted, expected);
    }

    #[test]
    fn fit_tan_sip_errors_on_too_few_pairs() {
        let prior = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let pairs = vec![
            StarPair {
                px_x: 0.0,
                px_y: 0.0,
                ra_deg: 180.0,
                dec_deg: 30.0,
                votes: 1,
                snr: 10.0,
            },
            StarPair {
                px_x: 10.0,
                px_y: 10.0,
                ra_deg: 180.01,
                dec_deg: 30.0,
                votes: 1,
                snr: 10.0,
            },
        ];
        let weights = vec![100.0, 100.0];
        let err = fit_tan_sip(&pairs, &weights, &prior, None)
            .err()
            .expect("2 pairs must error");
        assert!(err.to_string().contains("at least 3"));
    }

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            crpix1: 0.0,
            crpix2: 0.0,
            ..WcsSolution::diagonal((180.0, 0.0), (0.001, 0.001), (1000, 1000))
        }
    }

    fn pair_at_pixel(wcs: &WcsSolution, px_x: f64, px_y: f64) -> StarPair {
        // Invert the WCS to get the sky position that would project exactly to (px_x, px_y).
        let (ra, dec) = wcs.pixel_to_sky(px_x, px_y);
        StarPair {
            px_x,
            px_y,
            ra_deg: ra,
            dec_deg: dec,
            votes: 4,
            snr: 50.0,
        }
    }

    // The SIP distortion `fit_tan_sip` returns for `pairs` weighted by `pair_weight`, its
    // order chosen up to `max_order`.
    fn sip_up_to(pairs: &[StarPair], wcs: &WcsSolution, max_order: u32) -> Option<SipSolution> {
        let weights: Vec<f64> = pairs.iter().map(|p| pair_weight(p.snr)).collect();
        fit_tan_sip(pairs, &weights, wcs, Some(max_order))
            .unwrap()
            .sip
    }

    // Five pairs give no SIP term a leave-one-out prediction: order 2 has six
    // coefficients per axis, so its fit interpolates them.
    #[test]
    fn fit_tan_sip_with_fewer_pairs_than_order_2_terms_is_tan_only() {
        let wcs = sample_wcs();
        let pixels = [
            (-300.0, -200.0),
            (250.0, -150.0),
            (-100.0, 300.0),
            (350.0, 280.0),
            (20.0, 40.0),
        ];
        let pairs: Vec<StarPair> = pixels
            .iter()
            .map(|&(x, y)| pair_at_pixel(&wcs, x, y))
            .collect();
        assert!(sip_up_to(&pairs, &wcs, 2).is_none());
    }

    // Nine pairs, the order 2 distortion of the noisy-pairs test below with offsets of
    // up to 0.15 px. Order 3 has 10 terms per axis, so its fit interpolates the nine
    // pairs and has no leave-one-out residual to score.
    #[test]
    fn fit_tan_sip_does_not_choose_an_order_with_more_terms_than_pairs() {
        let wcs = sample_wcs();
        let rows = [
            (
                -434.3923094787462,
                -263.4569402268622,
                179.56643677099555,
                -0.26377890337197196,
            ),
            (
                -237.64040026926313,
                -388.62988735252065,
                179.76245609831219,
                -0.38916647241648644,
            ),
            (
                -413.0014901291966,
                -314.7388255773179,
                179.5875764816143,
                -0.3150859762030186,
            ),
            (
                -187.7466594795252,
                45.622753318515606,
                179.8124946276227,
                0.04556549801473213,
            ),
            (
                -296.98849789073023,
                170.29307188914072,
                179.70330373190043,
                0.17017504404868547,
            ),
            (
                -276.28081128579197,
                266.79835248586176,
                179.7240040735587,
                0.26664932933733465,
            ),
            (
                396.06507781039136,
                -57.60806874931899,
                180.39661147305665,
                -0.0575684468817908,
            ),
            (
                -473.1142360594911,
                160.21285272924138,
                179.52778824669133,
                0.16023806864970228,
            ),
            (
                -227.07744916626314,
                -286.4738854599543,
                179.77308315421993,
                -0.28662046640754996,
            ),
        ];
        let pairs: Vec<StarPair> = rows
            .iter()
            .map(|&(px_x, px_y, ra_deg, dec_deg)| StarPair {
                px_x,
                px_y,
                ra_deg,
                dec_deg,
                votes: 4,
                snr: 50.0,
            })
            .collect();
        for max_order in [3, 4] {
            let sip = sip_up_to(&pairs, &wcs, max_order);
            let chosen = sip.as_ref().map(|s| {
                let a = eval_sip_poly(&s.a_coeffs, 500.0, 500.0);
                let b = eval_sip_poly(&s.b_coeffs, 500.0, 500.0);
                (s.a_order, s.n_stars, a, b)
            });
            assert_eq!(chosen, None, "max order {max_order}");
        }
    }

    // Six pairs from an order 2 truth with no noise: order 2 has six terms per axis, as
    // many as the pairs, so it is not chosen even though it fits them exactly.
    #[test]
    fn fit_tan_sip_does_not_choose_an_order_with_as_many_terms_as_pairs() {
        let wcs = sample_wcs();
        let pixels = [
            (-300.0, -200.0),
            (250.0, -150.0),
            (-100.0, 300.0),
            (350.0, 280.0),
            (20.0, 40.0),
            (-420.0, 90.0),
        ];
        let pairs: Vec<StarPair> = pixels
            .iter()
            .map(|&(u, v)| {
                let (ra, dec) = wcs.pixel_to_sky(u + 4.0e-6 * u * u, v - 3.0e-6 * v * v);
                StarPair {
                    px_x: u,
                    px_y: v,
                    ra_deg: ra,
                    dec_deg: dec,
                    votes: 4,
                    snr: 50.0,
                }
            })
            .collect();
        let sip = sip_up_to(&pairs, &wcs, 2);
        assert_eq!(sip.map(|s| (s.a_order, s.n_stars)), None);
    }

    // Pairs on a 5 x 5 grid from a TAN-SIP truth of order 2, with offsets of up to
    // 0.3 px standing in for centroid noise. Order 4 fits the offsets as well as the
    // distortion; its leave-one-out residuals are larger than those of order 2.
    #[test]
    fn fit_tan_sip_chooses_a_lower_order_than_the_highest_for_few_noisy_pairs() {
        let wcs = sample_wcs();
        let (a20, b02, a11) = (4.0e-6_f64, -3.0e-6_f64, 2.0e-6_f64);
        let pairs: Vec<StarPair> = (0..25)
            .map(|k| {
                let (u, v) = (
                    f64::from(k % 5) * 200.0 - 400.0,
                    f64::from(k / 5) * 200.0 - 400.0,
                );
                let noise = |phase: f64| 0.3 * libm::sin(f64::from(k) * 2.399 + phase);
                let x = u + a20 * u * u + a11 * u * v + noise(0.0);
                let y = v + b02 * v * v + noise(1.3);
                let (ra, dec) = wcs.pixel_to_sky(x, y);
                StarPair {
                    px_x: u,
                    px_y: v,
                    ra_deg: ra,
                    dec_deg: dec,
                    votes: 4,
                    snr: 50.0,
                }
            })
            .collect();
        let sip = sip_up_to(&pairs, &wcs, 4).expect("an order 2 distortion");
        assert_eq!((sip.a_order, sip.b_order, sip.n_stars), (2, 2, 25));
    }

    // Pairs from a TAN-SIP truth in the FITS SIP convention (Shupe et al. 2005, eq. 1):
    // the detected pixel offset (u, v) maps to the intermediate pixel coordinate
    // x = u + A(u, v), y = v + B(u, v), and the catalog position is the TAN
    // deprojection of CD * (x, y).
    #[test]
    fn fit_tan_sip_happy_path_recovers_injected_distortion() {
        let wcs = sample_wcs();
        let (a20, b02) = (2.0e-7_f64, -1.5e-7_f64);
        let pairs: Vec<StarPair> = (-4..=4)
            .flat_map(|i| (-4..=4).map(move |j| (i as f64 * 100.0, j as f64 * 100.0)))
            .map(|(u, v)| {
                let (x, y) = (u + a20 * u * u, v + b02 * v * v);
                let (ra, dec) = wcs.pixel_to_sky(x, y);
                StarPair {
                    px_x: u,
                    px_y: v,
                    ra_deg: ra,
                    dec_deg: dec,
                    votes: 4,
                    snr: 50.0,
                }
            })
            .collect();

        let sip = sip_up_to(&pairs, &wcs, 2).expect("an order 2 distortion");
        assert_eq!((sip.a_order, sip.b_order, sip.n_stars), (2, 2, 81));
        // The injected values, to the rounding of the TAN projection.
        assert_eq!(
            sip.a_coeffs.get(&(2, 0)).copied(),
            Some(1.999999999884922e-7)
        );
        assert_eq!(
            sip.b_coeffs.get(&(0, 2)).copied(),
            Some(-1.5000000000003845e-7)
        );
        assert!(
            sip.rms_px < sip.linear_rms_px,
            "rms {} should be < linear_rms {}",
            sip.rms_px,
            sip.linear_rms_px
        );
    }
}
