// Checks of the solve parameters, and of the initial WCS against the hint.

use celestial_core::{angle::Angle, math::angular_separation};

use crate::error::SolveError;
use crate::fit_wcs::solution::WcsSolution;
use crate::input;
use crate::match_field::quad_stars::field_diagonal_deg;
use crate::metadata::image_metadata::ImageMetadata;

use super::SolveParams;

// The parameters `run` would otherwise first read after detection and matching. The
// largest factor between the fitted and hinted scales is at least 1, so a limit below 1
// rejects every solution.
pub(super) fn check_params(params: &SolveParams) -> Result<(), SolveError> {
    input::sip_order("sip_order", params.sip_order)?;
    // Infinity is allowed in both limits and turns its check off.
    let factor = params.max_scale_factor;
    if factor.is_nan() || factor < 1.0 {
        return Err(SolveError::InvalidParameter {
            name: "max_scale_factor",
            reason: format!("must be at least 1, got {factor}"),
        });
    }
    let offset = params.max_center_offset_fov;
    if offset.is_nan() || offset <= 0.0 {
        return Err(SolveError::InvalidParameter {
            name: "max_center_offset_fov",
            reason: format!("must be positive, got {offset}"),
        });
    }
    input::match_params(&params.matching)?;
    input::positive("tolerance_px", params.refine.tolerance_px)
}

pub(super) fn check_scale_sanity(
    wcs: &WcsSolution,
    meta: &ImageMetadata,
    params: &SolveParams,
    label: &str,
) -> Result<(), SolveError> {
    let fitted_scale = wcs.mean_scale_arcsec();
    let hint_scale = meta.scale_arcsec;
    let factor = (fitted_scale / hint_scale).max(hint_scale / fitted_scale);
    if factor <= params.max_scale_factor {
        return Ok(());
    }
    Err(SolveError::no_match(format!(
        "{label} WCS plate scale {fitted_scale:.3}\"/px diverges from hint {hint_scale:.3}\"/px by factor {factor:.2} (limit {:.2})",
        params.max_scale_factor,
    )))
}

pub(super) fn check_center_offset_sanity(
    wcs: &WcsSolution,
    meta: &ImageMetadata,
    width: usize,
    height: usize,
    params: &SolveParams,
    label: &str,
) -> Result<(), SolveError> {
    let sep_deg = wcs_to_hint_separation_deg(wcs, meta);
    let fov_diag_deg = field_diagonal_deg(width, height, wcs.mean_scale_arcsec());
    let max_offset_deg = params.max_center_offset_fov * fov_diag_deg;
    if sep_deg <= max_offset_deg {
        return Ok(());
    }
    Err(SolveError::no_match(format!(
        "{label} WCS center {sep_deg:.3}° from hint exceeds {:.3}° ({}× FOV diagonal {fov_diag_deg:.3}°)",
        max_offset_deg,
        params.max_center_offset_fov,
    )))
}

fn wcs_to_hint_separation_deg(wcs: &WcsSolution, meta: &ImageMetadata) -> f64 {
    let separation = angular_separation(
        meta.hint.ra().radians(),
        meta.hint.dec().radians(),
        Angle::from_degrees(wcs.crval1).radians(),
        Angle::from_degrees(wcs.crval2).radians(),
    );
    Angle::from_radians(separation).degrees()
}

#[cfg(test)]
mod tests {
    use super::*;

    use celestial_catalog::query::catalog::Catalog;
    use celestial_time::julian::JulianDate;

    use crate::fit_wcs::refine::RefineParams;
    use crate::match_field::MatchParams;
    use crate::solve::solve;
    use crate::solve::test_support::meta_at;

    #[test]
    fn check_scale_sanity_accepts_matching_scale() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let params = SolveParams::default();
        assert!(check_scale_sanity(&wcs, &meta, &params, "test").is_ok());
    }

    #[test]
    fn check_center_offset_sanity_accepts_matching_center() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let params = SolveParams::default();
        assert!(check_center_offset_sanity(&wcs, &meta, 1024, 1024, &params, "test").is_ok());
    }

    #[test]
    fn check_scale_sanity_rejects_scale_too_large() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = WcsSolution::north_up((180.0, 30.0), 6.0, (1024, 1024));
        let params = SolveParams::default();
        let err = check_scale_sanity(&wcs, &meta, &params, "label-A").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("label-A"), "label propagated: {msg}");
        assert!(msg.contains("plate scale"), "scale error wording: {msg}");
    }

    #[test]
    fn check_scale_sanity_rejects_scale_too_small() {
        // Exercises the other half of the .max() in the factor computation.
        let meta = meta_at(180.0, 30.0, 6.0);
        let wcs = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let params = SolveParams::default();
        let err = check_scale_sanity(&wcs, &meta, &params, "test").unwrap_err();
        assert!(err.to_string().contains("plate scale"));
    }

    #[test]
    fn check_center_offset_sanity_rejects_center_too_far() {
        // 30° offset on a ~0.43° FOV (1024px × 1.5"/px) far exceeds 1.0 FOV.
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = WcsSolution::north_up((180.0 + 30.0, 30.0), 1.5, (1024, 1024));
        let params = SolveParams::default();
        let err = check_center_offset_sanity(&wcs, &meta, 1024, 1024, &params, "centerlabel")
            .unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("centerlabel"));
        assert!(msg.contains("center"), "center error wording: {msg}");
    }

    #[test]
    fn check_center_offset_sanity_respects_relaxed_thresholds() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = WcsSolution::north_up((180.0 + 30.0, 30.0), 1.5, (1024, 1024));
        let params = SolveParams {
            max_center_offset_fov: 1000.0,
            ..SolveParams::default()
        };
        assert!(check_center_offset_sanity(&wcs, &meta, 1024, 1024, &params, "x").is_ok());
    }

    // The FITS reader of celestial_wcs reads SIP orders up to 9, and a SIP fit needs
    // order 2 or above.
    #[test]
    fn solve_rejects_a_sip_order_outside_2_to_9() {
        let scene = crate::sky_scene::star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        for order in [10, 1] {
            let params = SolveParams {
                sip_order: Some(order),
                ..SolveParams::default()
            };
            let outcome = solve(&scene.image(), &catalog)
                .with_metadata(scene.metadata())
                .params(params)
                .run();
            let Err(err) = outcome else {
                panic!("SIP order {order} was accepted");
            };
            let expected = format!("must be from 2 to 9, got {order}");
            assert!(matches!(
                err,
                SolveError::InvalidParameter { name: "sip_order", ref reason } if *reason == expected
            ));
        }
    }

    #[test]
    fn solve_rejects_non_finite_metadata_and_parameters() {
        let scene = crate::sky_scene::star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let image = scene.image();
        let run = |meta: ImageMetadata, params: SolveParams| {
            let out = solve(&image, &catalog)
                .with_metadata(meta)
                .params(params)
                .run();
            out.err().map(|e| e.describe())
        };
        let meta = || scene.metadata();
        let defaults = SolveParams::default;
        let cases = [
            (
                run(
                    ImageMetadata {
                        scale_arcsec: f64::NAN,
                        ..meta()
                    },
                    defaults(),
                ),
                "Metadata: plate scale must be finite, got NaN",
            ),
            (
                run(
                    ImageMetadata {
                        epoch: JulianDate::new(f64::NAN, 0.0),
                        ..meta()
                    },
                    defaults(),
                ),
                "Metadata: epoch must be finite, got NaN",
            ),
            (
                run(
                    ImageMetadata {
                        focal_mm: Some(f64::INFINITY),
                        ..meta()
                    },
                    defaults(),
                ),
                "Metadata: focal length must be finite, got inf",
            ),
            (
                run(
                    meta(),
                    SolveParams {
                        max_scale_factor: f64::NAN,
                        ..defaults()
                    },
                ),
                "InvalidParameter: invalid parameter max_scale_factor: must be at least 1, got NaN",
            ),
            (
                run(
                    meta(),
                    SolveParams {
                        max_center_offset_fov: -1.0,
                        ..defaults()
                    },
                ),
                "InvalidParameter: invalid parameter max_center_offset_fov: must be positive, got -1",
            ),
            (
                run(
                    meta(),
                    SolveParams {
                        matching: MatchParams {
                            search_radius_deg: Some(f64::NAN),
                            ..MatchParams::default()
                        },
                        ..defaults()
                    },
                ),
                "InvalidParameter: invalid parameter search_radius_deg: must be finite and positive, got NaN",
            ),
            (
                run(
                    meta(),
                    SolveParams {
                        refine: RefineParams {
                            tolerance_px: 0.0,
                            ..RefineParams::default()
                        },
                        ..defaults()
                    },
                ),
                "InvalidParameter: invalid parameter tolerance_px: must be finite and positive, got 0",
            ),
        ];
        let (got, want): (Vec<_>, Vec<_>) = cases
            .iter()
            .map(|(err, expected)| (err.as_deref(), Some(*expected)))
            .unzip();
        assert_eq!(got, want);
    }

    #[test]
    fn infinite_sanity_limits_switch_the_checks_off() {
        let scene = crate::sky_scene::star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let params = SolveParams {
            max_scale_factor: f64::INFINITY,
            max_center_offset_fov: f64::INFINITY,
            ..SolveParams::default()
        };
        let out = solve(&scene.image(), &catalog)
            .with_metadata(scene.metadata())
            .params(params)
            .run();
        assert_eq!(out.err().map(|e| e.describe()), None);
    }

    #[test]
    fn solve_rejects_parameters_that_form_no_quad() {
        let scene = crate::sky_scene::star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let image = scene.image();
        let run = |matching: MatchParams| {
            let params = SolveParams {
                matching,
                ..SolveParams::default()
            };
            let out = solve(&image, &catalog)
                .with_metadata(scene.metadata())
                .params(params)
                .run();
            out.err().map(|e| e.describe())
        };
        let defaults = MatchParams::default;
        let got = [
            run(MatchParams {
                max_stars: 3,
                ..defaults()
            }),
            run(MatchParams {
                k_neighbors: 2,
                ..defaults()
            }),
            run(MatchParams {
                search_radius_deg: Some(180.5),
                ..defaults()
            }),
        ];
        let want = [
            "InvalidParameter: invalid parameter max_stars: must be at least 4, the stars of one quad, got 3",
            "InvalidParameter: invalid parameter k_neighbors: must be at least 3, the neighbours a quad adds to its first star, got 2",
            "InvalidParameter: invalid parameter search_radius_deg: must be at most 180, got 180.5",
        ];
        assert_eq!(got.each_ref().map(|g| g.as_deref()), want.map(Some));
    }
}
