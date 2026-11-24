// The stages of a solve, from detection to the verified result.

use celestial_catalog::query::catalog::Catalog;
use celestial_time::julian::JulianDate;

use celestial_images::formats::image::Image;
use celestial_images::formats::pixel_data::PixelData;

use crate::detect::{find_bright_stars, DetectedStar, DetectionParams};
use crate::error::SolveError;
use crate::fit_wcs::fit_wcs;
use crate::fit_wcs::refine::{refine_wcs, RefineParams};
use crate::fit_wcs::solution::WcsSolution;
use crate::match_field::{match_field, StarPair};
use crate::metadata::image_metadata::ImageMetadata;
use crate::verify_solution::{verify_outcome, Evidence};

use super::sanity::{check_center_offset_sanity, check_params, check_scale_sanity};
use super::{RefineOutcome, Refinement, SolveParams, SolveResult};

macro_rules! time_stage {
    ($label:literal, $value:ident => $detail:expr, $body:expr) => {{
        let __t = std::time::Instant::now();
        let $value = $body;
        log::info!("[solve] {}: {:?} {}", $label, __t.elapsed(), $detail);
        $value
    }};
}

pub(super) fn run_pipeline(
    image: &Image,
    catalog: &Catalog,
    meta: &ImageMetadata,
    params: &SolveParams,
) -> Result<SolveResult, SolveError> {
    check_params(params)?;
    let w = image.width();
    let h = image.height();

    let stars = time_stage!("detect", v => format!("({} stars)", v.len()),
        detect_stars(image, &params.detection)?);
    let field_match = time_stage!("match", fm => format!("({} quads, {} pairs)", fm.matches.len(), fm.pairs.len()),
        match_field(&stars, image, catalog, &meta.hint, meta.scale_arcsec, meta.epoch, &params.matching)?);
    let initial_wcs = time_stage!("fit_wcs", _w => "",
        fit_initial_wcs(&field_match.pairs, w, h, meta, params)?);
    let outcome = time_stage!(
        "refine", o => format!("(refined={}, {} pairs)", o.refinement.is_refined(), o.pairs.len()),
        refine_or_fallback(&initial_wcs, &stars, catalog, &field_match.pairs, meta.epoch, &params.refine, params.sip_order)?
    );
    let evidence = Evidence {
        catalog,
        stars: &stars,
        epoch: meta.epoch,
        min_pairs: params.min_refined_pairs,
    };
    let verified = time_stage!("verify", v => format!("(log-odds {:.2}, sip={})", v.log_odds, v.outcome.sip.is_some()),
        verify_outcome(&evidence, outcome, &initial_wcs, &field_match.pairs)?);

    Ok(SolveResult {
        stars,
        field_match,
        initial_wcs,
        wcs: verified.outcome.wcs,
        sip: verified.outcome.sip,
        refinement: verified.outcome.refinement,
        pairs: verified.outcome.pairs,
        fits_wcs: verified.fits_wcs,
        log_odds: verified.log_odds,
    })
}

fn detect_stars(image: &Image, params: &DetectionParams) -> Result<Vec<DetectedStar>, SolveError> {
    let w = image.width();
    let h = image.height();
    let stars = match &image.pixels {
        PixelData::U8(v) => find_bright_stars(v, w, h, params)?,
        PixelData::I8(v) => find_bright_stars(v, w, h, params)?,
        PixelData::U16(v) => find_bright_stars(v, w, h, params)?,
        PixelData::I16(v) => find_bright_stars(v, w, h, params)?,
        PixelData::U32(v) => find_bright_stars(v, w, h, params)?,
        PixelData::I32(v) => find_bright_stars(v, w, h, params)?,
        PixelData::U64(v) => find_bright_stars(v, w, h, params)?,
        PixelData::I64(v) => find_bright_stars(v, w, h, params)?,
        PixelData::F32(v) => find_bright_stars(v, w, h, params)?,
        PixelData::F64(v) => find_bright_stars(v, w, h, params)?,
    };
    if stars.is_empty() {
        return Err(SolveError::NoStars);
    }
    Ok(stars)
}

pub(super) fn fit_initial_wcs(
    pairs: &[StarPair],
    width: usize,
    height: usize,
    meta: &ImageMetadata,
    params: &SolveParams,
) -> Result<WcsSolution, SolveError> {
    let wcs = fit_wcs(
        pairs,
        width,
        height,
        meta.hint.ra().degrees(),
        meta.hint.dec().degrees(),
        meta.focal_mm,
        meta.pixel_um,
    )?;
    check_scale_sanity(&wcs, meta, params, "initial")?;
    check_center_offset_sanity(&wcs, meta, width, height, params, "initial")?;
    Ok(wcs)
}

// The refined solution, or the initial WCS with `fallback_pairs` when the refinement
// cannot fit a model. A catalog error, or any other error that is not about the fit, is
// returned.
fn refine_or_fallback(
    initial_wcs: &WcsSolution,
    stars: &[DetectedStar],
    catalog: &Catalog,
    fallback_pairs: &[StarPair],
    epoch: JulianDate,
    refine_params: &RefineParams,
    sip_order: Option<u32>,
) -> Result<RefineOutcome, SolveError> {
    match refine_wcs(initial_wcs, stars, catalog, epoch, refine_params, sip_order) {
        Ok(r) => Ok(RefineOutcome {
            wcs: r.wcs,
            sip: r.sip,
            pairs: r.pairs,
            refinement: Refinement::Refined,
        }),
        Err(e @ (SolveError::FitFailed { .. } | SolveError::Wcs(_))) => {
            log::warn!("refine_wcs failed: {}", e.chain());
            Ok(RefineOutcome {
                wcs: initial_wcs.clone(),
                sip: None,
                pairs: fallback_pairs.to_vec(),
                refinement: Refinement::Failed(e),
            })
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::solve::test_support::{
        build_catalog_keep_alive, meta_at, pairs_on_tangent_plane, zero_image,
    };

    fn image_with_one_star(pixels: PixelData, w: usize, h: usize) -> Image {
        Image::new(pixels, [w, h])
    }

    #[test]
    fn detect_stars_zero_image_errors_with_no_stars_detected() {
        let img = zero_image(128, 128);
        let err =
            detect_stars(&img, &DetectionParams::default()).expect_err("zero image must error");
        assert!(err.to_string().contains("no stars detected"));
    }

    #[test]
    fn detect_stars_finds_star_in_each_pixel_type() {
        // Build a synthetic star in f32 once, cast/scale into each PixelData
        // variant, and verify the detector finds it through all six match arms.
        let w = 128;
        let h = 128;
        let base =
            crate::detect::test_images::gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);

        // (variant_name, pixel_data), built so each variant carries a detectable star.
        let cases: Vec<(&str, PixelData)> = vec![
            (
                "U8",
                PixelData::from(
                    base.iter()
                        .map(|&v| (v / 25.0).clamp(0.0, 255.0) as u8)
                        .collect::<Vec<_>>(),
                ),
            ),
            (
                "U16",
                PixelData::from(
                    base.iter()
                        .map(|&v| v.clamp(0.0, 65535.0) as u16)
                        .collect::<Vec<_>>(),
                ),
            ),
            (
                "I16",
                PixelData::from(
                    base.iter()
                        .map(|&v| v.clamp(0.0, 32767.0) as i16)
                        .collect::<Vec<_>>(),
                ),
            ),
            (
                "I32",
                PixelData::from(base.iter().map(|&v| v as i32).collect::<Vec<_>>()),
            ),
            ("F32", PixelData::from(base.clone())),
            (
                "F64",
                PixelData::from(base.iter().map(|&v| v as f64).collect::<Vec<_>>()),
            ),
        ];

        for (name, pixels) in cases {
            let img = image_with_one_star(pixels, w, h);
            let stars = detect_stars(&img, &DetectionParams::default())
                .unwrap_or_else(|e| panic!("{name}: detect_stars failed: {e}"));
            assert!(!stars.is_empty(), "{name}: expected ≥1 star");
            // The painted star sits at (64, 64).
            assert_eq!((stars[0].x, stars[0].y), (64.0, 64.0), "{name}");
        }
    }

    #[test]
    fn fit_initial_wcs_recovers_hint_from_clean_pairs() {
        let center_ra = 180.0;
        let center_dec = 30.0;
        let scale = 1.5;
        let pairs = pairs_on_tangent_plane(center_ra, center_dec, scale);
        let meta = meta_at(center_ra, center_dec, scale);
        let params = SolveParams::default();

        let wcs = fit_initial_wcs(&pairs, 1024, 1024, &meta, &params)
            .expect("clean pairs must produce a valid WCS");

        // The pairs were placed about CRVAL (180, 30) at 1.5"/px. The fitted Dec is
        // 3.6e-15 deg from it and the scales are 3.1e-14 and 1.1e-14 arcsec/px from it.
        assert_eq!((wcs.crval1, wcs.crval2), (center_ra, 29.999999999999996));
        assert_eq!(wcs.scale_arcsec(), (1.5000000000000306, 1.5000000000000109));
    }

    #[test]
    fn fit_initial_wcs_propagates_fit_wcs_error_on_too_few_pairs() {
        // fit_wcs requires >=3 pairs; pass 2 and check the error surfaces from
        // fit_initial_wcs rather than panicking or returning a bogus WCS.
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
                dec_deg: 30.01,
                votes: 1,
                snr: 10.0,
            },
        ];
        let meta = meta_at(180.0, 30.0, 1.5);
        let params = SolveParams::default();
        let err =
            fit_initial_wcs(&pairs, 1024, 1024, &meta, &params).expect_err("2 pairs must error");
        assert!(err.to_string().contains("at least 3"), "got: {err}");
    }

    #[test]
    fn refine_or_fallback_returns_initial_wcs_and_pairs_on_refine_err() {
        // Empty stars slice → refine_wcs can find no projection matches and
        // returns Err on iteration 0. The fallback path must produce a clone
        // of initial_wcs, a clone of fallback_pairs, and refined=false.
        let (cat, _file) = build_catalog_keep_alive();
        let initial = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let fallback = vec![
            StarPair {
                px_x: 1.0,
                px_y: 2.0,
                ra_deg: 180.0,
                dec_deg: 30.0,
                votes: 3,
                snr: 50.0,
            },
            StarPair {
                px_x: 5.0,
                px_y: 6.0,
                ra_deg: 180.01,
                dec_deg: 30.01,
                votes: 2,
                snr: 30.0,
            },
        ];
        let epoch = JulianDate::new(2451545.0, 0.0);
        let refine_params = RefineParams::default();

        let outcome = refine_or_fallback(
            &initial,
            &[],
            &cat,
            &fallback,
            epoch,
            &refine_params,
            Some(4),
        )
        .unwrap();

        assert!(
            matches!(
                outcome.refinement,
                Refinement::Failed(SolveError::FitFailed { .. })
            ),
            "empty stars must take fallback branch"
        );
        assert_eq!(outcome.wcs.crval1, initial.crval1);
        assert_eq!(outcome.wcs.crval2, initial.crval2);
        assert_eq!(outcome.pairs.len(), 2);
        assert_eq!(outcome.pairs[0].px_x, 1.0);
        assert_eq!(outcome.pairs[1].px_x, 5.0);
        assert!(outcome.sip.is_none());
    }

    // Only a failure to fit makes the solve keep the initial WCS. A bad parameter is
    // returned, where it was logged and the initial WCS kept.
    #[test]
    fn refine_or_fallback_returns_errors_that_are_not_about_the_fit() {
        let (cat, _file) = build_catalog_keep_alive();
        let initial = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let epoch = JulianDate::new(2451545.0, 0.0);
        let params = RefineParams {
            tolerance_px: f64::NAN,
            ..RefineParams::default()
        };
        let outcome = refine_or_fallback(&initial, &[], &cat, &[], epoch, &params, Some(4));
        assert_eq!(
            outcome.err().map(|e| e.describe()).as_deref(),
            Some("InvalidParameter: invalid parameter tolerance_px: must be finite and positive, got NaN")
        );
    }
}
