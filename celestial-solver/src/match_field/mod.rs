//! Star-quad hash matching against a catalog.
//!
//! Builds four-star quads from detections and from a cone-searched catalog region,
//! then matches them by geometric hash. Each quad match is a hypothesis, verified by
//! the Bayesian log-odds that it is real rather than chance. The output is a set of
//! individual [`StarPair`]s (image pixel ↔ catalog RA/Dec), aggregated across all
//! matched quads with vote counts, that support the accepted hypothesis.

mod aggregate;
mod catalog_side;
mod hash_index;
mod parity;
pub(crate) mod quad_stars;
mod uniform;
mod verify;

use celestial_catalog::query::quad::Quad;
use celestial_coords::frames::icrs::ICRSPosition;
use celestial_images::formats::image::Image;
use celestial_time::julian::JulianDate;

use crate::detect::DetectedStar;
use crate::error::SolveError;
use crate::input;
use crate::log_odds::field::Field;

use catalog_side::search_catalog;
use parity::ParitySearch;

/// Tuning for the quad matching stage of the solve.
///
/// `max_stars` caps the number of detections used for quads, the brightest of each cell
/// of a grid over the image. More cost matching time; fewer risk missing the right
/// solution in sparse fields.
/// `search_radius_deg` overrides the cone search radius (default: auto-computed from
/// image dimensions and plate scale).
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::match_field::MatchParams;
///
/// let params = MatchParams {
///     max_stars: 200,
///     ..MatchParams::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct MatchParams {
    /// Maximum number of detections to build image quads from: the brightest of each cell
    /// of a grid over the image, in cells of about `k_neighbors` of them. At least 4.
    pub max_stars: usize,
    /// How many nearest-neighbor stars to consider when building quads per anchor. At
    /// least 3.
    pub k_neighbors: usize,
    /// Catalog cone search radius override, in degrees, above 0 and at most 180. `None`
    /// auto-computes.
    pub search_radius_deg: Option<f64>,
}

impl Default for MatchParams {
    fn default() -> Self {
        Self {
            max_stars: 100,
            k_neighbors: 10,
            search_radius_deg: None,
        }
    }
}

/// A geometric quad match: one image-side quad paired with a catalog-side quad.
///
/// Produced by the quad matching stage. The two quads have similar four-star topology (edge
/// ratios and interior angles) under the current hashing scheme.
#[derive(Debug, Clone)]
pub struct QuadMatch {
    /// Quad in image pixel space.
    pub image_quad: Quad,
    /// Quad in catalog sky space.
    pub catalog_quad: Quad,
}

/// A single matched star pair: detected pixel ↔ catalog sky position.
///
/// Aggregated across all [`QuadMatch`]es. Multiple quads touching the same pair bump
/// its `votes` count; the solver weights pairs by vote count and SNR.
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
/// let top = &result.pairs[0];
/// println!("pixel ({:.1}, {:.1}) = sky ({:.4}°, {:.4}°)", top.px_x, top.px_y, top.ra_deg, top.dec_deg);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct StarPair {
    /// Detected pixel x.
    pub px_x: f64,
    /// Detected pixel y.
    pub px_y: f64,
    /// Catalog RA in degrees.
    pub ra_deg: f64,
    /// Catalog Dec in degrees.
    pub dec_deg: f64,
    /// Number of quad matches that contributed this pair.
    pub votes: usize,
    /// SNR from the detection.
    pub snr: f64,
}

impl StarPair {
    // Total order used wherever pair order affects a result: votes descending, then
    // snr descending, then pixel position, then sky position.
    pub(crate) fn cmp_rank(&self, other: &Self) -> std::cmp::Ordering {
        other
            .votes
            .cmp(&self.votes)
            .then(other.snr.total_cmp(&self.snr))
            .then(self.px_x.total_cmp(&other.px_x))
            .then(self.px_y.total_cmp(&other.px_y))
            .then(self.ra_deg.total_cmp(&other.ra_deg))
            .then(self.dec_deg.total_cmp(&other.dec_deg))
    }
}

/// Output of the quad matching stage, in [`SolveResult::field_match`].
///
/// [`SolveResult::field_match`]: crate::solve::SolveResult::field_match
///
/// Exposes the raw quad-star sets and quad matches alongside the aggregated
/// [`StarPair`]s. Pairs are what the WCS fitter consumes; the other fields are
/// exposed for debug rendering and analysis.
#[derive(Debug)]
pub struct FieldMatch {
    /// Image-side quad-star set (brightest detections with quad-relevant bookkeeping).
    pub image_stars: Vec<celestial_catalog::query::quad::QuadStar>,
    /// Catalog-side quad-star set (cone-searched stars with quad-relevant bookkeeping).
    pub catalog_stars: Vec<celestial_catalog::query::quad::QuadStar>,
    /// Every quad match found.
    pub matches: Vec<QuadMatch>,
    /// Pixel↔sky pairs that support the verified hypothesis: the aggregated quad-match
    /// pairs it supports, the detections it places on catalog stars, and its own quad
    /// pairs.
    pub pairs: Vec<StarPair>,
}

/// Matches detected stars against a catalog via four-star quad hashing.
///
/// Builds quads from up to `params.max_stars` detections, the brightest of each cell of a
/// grid over the image, cone-searches the catalog around `hint` with a radius derived
/// from the image geometry, builds catalog quads from the brightest stars of cells of
/// the same area at the same density, and matches by geometric hash. Each quad match in
/// both parity orientations is a hypothesis, accepted only when the Bayesian log-odds
/// that it is real rather than chance exceeds ln(1e9) (Lang et al. 2010, sec. 2.4).
/// Returns the pairs that support the first accepted hypothesis.
///
/// # Errors
///
/// Returns [`SolveError::InvalidParameter`] when `scale_arcsec` is not finite and
/// positive, a set `params.search_radius_deg` is not above 0 and at most 180,
/// `params.max_stars` is below 4, `params.k_neighbors` is below 3, `epoch` is not
/// finite, or a star's position or SNR is not finite; [`SolveError::NoStars`] when
/// `stars` is empty; [`SolveError::Catalog`] when the cone search fails; and
/// [`SolveError::NoMatch`] when no quad match passes verification.
pub(crate) fn match_field(
    stars: &[DetectedStar],
    image: &Image,
    catalog: &celestial_catalog::query::catalog::Catalog,
    hint: &ICRSPosition,
    scale_arcsec: f64,
    epoch: JulianDate,
    params: &MatchParams,
) -> Result<FieldMatch, SolveError> {
    check_inputs(stars, scale_arcsec, epoch, params)?;
    let size = (image.width(), image.height());
    let centre = (hint.ra().degrees(), hint.dec().degrees());
    let side = search_catalog(
        catalog,
        centre,
        size,
        scale_arcsec,
        epoch,
        stars.len(),
        params,
    )?;
    let field = Field::new(stars, size, centre, &side.references);
    let search = ParitySearch {
        stars,
        image_size: size,
        scale_arcsec,
        params,
        catalog: &side,
        field: &field,
    };
    let (attempt, accepted) = search.run()?;
    Ok(FieldMatch {
        image_stars: attempt.img_stars,
        catalog_stars: side.quad_stars,
        matches: attempt.matches,
        pairs: accepted.pairs,
    })
}

fn check_inputs(
    stars: &[DetectedStar],
    scale_arcsec: f64,
    epoch: JulianDate,
    params: &MatchParams,
) -> Result<(), SolveError> {
    input::match_params(params)?;
    input::positive("scale_arcsec", scale_arcsec)?;
    input::epoch(epoch)?;
    input::stars(stars)?;
    if stars.is_empty() {
        return Err(SolveError::NoStars);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use celestial_catalog::query::catalog::Catalog;
    use celestial_images::formats::pixel_data::PixelData;
    use celestial_time::julian::JulianDate;

    use crate::detect::DetectedStar;
    use crate::sky_scene::{build, SynthStar};

    use super::quad_stars::{catalog_cone_search, search_radius, QuadSelection};

    fn det(x: f64, y: f64, flux: f64) -> DetectedStar {
        DetectedStar {
            x,
            y,
            flux,
            snr: 20.0,
            saturated: false,
            saturated_count: 0,
            background: 0.0,
        }
    }

    fn empty_image(w: usize, h: usize) -> Image {
        Image::new(PixelData::from(vec![0u16; w * h]), [w, h])
    }

    #[test]
    fn match_params_default_values() {
        let p = MatchParams::default();
        assert_eq!(p.max_stars, 100);
        assert_eq!(p.k_neighbors, 10);
        assert!(p.search_radius_deg.is_none());
    }

    #[test]
    fn match_params_clone() {
        let p = MatchParams {
            max_stars: 250,
            k_neighbors: 20,
            search_radius_deg: Some(3.5),
        };
        let q = p.clone();
        assert_eq!(q.max_stars, 250);
        assert_eq!(q.k_neighbors, 20);
        assert_eq!(q.search_radius_deg, Some(3.5));
    }

    #[test]
    fn star_pair_clone_preserves_fields() {
        let p = StarPair {
            px_x: 1.0,
            px_y: 2.0,
            ra_deg: 10.0,
            dec_deg: 20.0,
            votes: 5,
            snr: 30.0,
        };
        let q = p.clone();
        assert_eq!(q.px_x, 1.0);
        assert_eq!(q.px_y, 2.0);
        assert_eq!(q.ra_deg, 10.0);
        assert_eq!(q.dec_deg, 20.0);
        assert_eq!(q.votes, 5);
        assert_eq!(q.snr, 30.0);
    }

    #[test]
    fn quad_match_clone_preserves_quads() {
        use celestial_catalog::query::quad::Quad;
        let qm = QuadMatch {
            image_quad: Quad {
                hash: [0.1, 0.2, 0.3, 0.4],
                star_indices: [1, 2, 3, 4],
            },
            catalog_quad: Quad {
                hash: [0.5, 0.6, 0.7, 0.8],
                star_indices: [5, 6, 7, 8],
            },
        };
        let c = qm.clone();
        assert_eq!(c.image_quad.hash, [0.1, 0.2, 0.3, 0.4]);
        assert_eq!(c.catalog_quad.star_indices, [5, 6, 7, 8]);
    }

    #[test]
    fn match_field_no_detections_is_an_error() {
        // Catalog with a handful of stars, but zero detections: no image quads can be
        // built, which is reported before the catalog search.
        let center_ra = 180.0;
        let center_dec = 0.0;
        let stars: Vec<SynthStar> = (0..30)
            .map(|i| {
                let i = i as f64;
                SynthStar {
                    source_id: 1000 + i as i64,
                    ra: center_ra + (i * 0.01),
                    dec: center_dec + (i * 0.01 - 0.15),
                    mag: 8.0 + (i * 0.05) as f32,
                }
            })
            .collect();
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();

        let img = empty_image(1024, 1024);
        let hint = ICRSPosition::from_degrees(center_ra, center_dec).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let err = match_field(
            &[],
            &img,
            &catalog,
            &hint,
            2.0,
            epoch,
            &MatchParams::default(),
        )
        .unwrap_err();
        assert_eq!(err.describe(), "NoStars: no stars detected");
    }

    #[test]
    fn match_field_rejects_unrelated_star_fields() {
        for seed in 1..=6 {
            let scene = crate::sky_scene::unrelated_field(seed, 100);
            let catalog = Catalog::open(scene.catalog.path()).unwrap();
            let outcome = match_field(
                &scene.detections(),
                &empty_image(1024, 1024),
                &catalog,
                &crate::sky_scene::hint(),
                crate::sky_scene::SCALE_ARCSEC,
                crate::sky_scene::epoch(),
                &MatchParams::default(),
            );
            let Err(err) = outcome else {
                panic!("seed {seed}: stars unrelated to the catalog matched");
            };
            assert!(
                err.to_string()
                    .starts_with("no quad match passed verification: det(CD) > 0: "),
                "seed {seed}: {err}"
            );
        }
    }

    #[test]
    fn match_field_search_radius_override_is_used() {
        // Tiny override radius → very few catalog stars in the cone → catalog_stars short.
        let center_ra = 0.0;
        let center_dec = 0.0;
        let stars: Vec<SynthStar> = (0..50)
            .map(|i| SynthStar {
                source_id: 2000 + i,
                ra: center_ra + (i as f64) * 0.2,
                dec: center_dec + (i as f64) * 0.01,
                mag: 9.0,
            })
            .collect();
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let cone = |radius_deg: f64| {
            let radius = search_radius(1024, 1024, 1.0, Some(radius_deg));
            let selection = QuadSelection {
                count: MatchParams::default().max_stars,
                cell_rad: 1.0,
            };
            catalog_cone_search(
                &catalog,
                (center_ra, center_dec),
                radius,
                &selection,
                0,
                epoch,
            )
            .unwrap()
            .quad_stars
        };

        // Wider radius must reach at least as many catalog stars.
        assert!(cone(5.0).len() >= cone(0.05).len());
    }

    #[test]
    fn match_field_errors_when_no_parity_verifies() {
        // Detections on a line, unrelated to a ring of catalog stars: both parities are
        // tried and neither verifies.
        let center_ra = 90.0;
        let center_dec = 45.0;
        let mut stars: Vec<SynthStar> = Vec::new();
        for i in 0..40 {
            let theta = (i as f64) * 0.4;
            stars.push(SynthStar {
                source_id: 3000 + i as i64,
                ra: center_ra + 0.1 * libm::cos(theta),
                dec: center_dec + 0.1 * libm::sin(theta),
                mag: 9.0 + (i as f32) * 0.05,
            });
        }
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();

        let detections: Vec<DetectedStar> = (0..20)
            .map(|i| {
                let i = i as f64;
                det(100.0 + i * 30.0, 200.0 + i * 25.0, 5000.0 - i * 100.0)
            })
            .collect();

        let img = empty_image(1024, 1024);
        let hint = ICRSPosition::from_degrees(center_ra, center_dec).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let err = match_field(
            &detections,
            &img,
            &catalog,
            &hint,
            2.0,
            epoch,
            &MatchParams::default(),
        )
        .unwrap_err();
        assert!(
            err.to_string()
                .starts_with("no quad match passed verification: det(CD) > 0: "),
            "{err}"
        );
    }

    #[test]
    fn match_field_finds_sparse_field_among_unrelated_stars() {
        for flip in [false, true] {
            let scene = crate::sky_scene::star_field(3, 40, 30, flip);
            let catalog = Catalog::open(scene.catalog.path()).unwrap();
            let out = match_field(
                &scene.detections(),
                &empty_image(1024, 1024),
                &catalog,
                &crate::sky_scene::hint(),
                crate::sky_scene::SCALE_ARCSEC,
                crate::sky_scene::epoch(),
                &MatchParams::default(),
            )
            .unwrap_or_else(|e| panic!("flip={flip}: {e:#}"));
            scene.assert_true_pairs(&out.pairs);
        }
    }

    fn match_sparse_scene(
        flip: bool,
        max_stars: usize,
    ) -> (crate::sky_scene::Scene, Result<FieldMatch, SolveError>) {
        let scene = crate::sky_scene::star_field(3, 40, 30, flip);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let params = MatchParams {
            max_stars,
            ..MatchParams::default()
        };
        let out = match_field(
            &scene.detections(),
            &empty_image(1024, 1024),
            &catalog,
            &crate::sky_scene::hint(),
            crate::sky_scene::SCALE_ARCSEC,
            crate::sky_scene::epoch(),
            &params,
        );
        (scene, out)
    }

    #[test]
    fn match_field_pairs_cover_every_catalog_star_in_the_cone() {
        for flip in [false, true] {
            let (scene, out) = match_sparse_scene(flip, MatchParams::default().max_stars);
            let out = out.unwrap_or_else(|e| panic!("flip={flip}: {e:#}"));
            scene.assert_true_pairs(&out.pairs);
            // 13 of the 18 scene stars lie inside the cone about the hint; the other 5
            // are not catalog references.
            assert_eq!(scene.truth.len(), 18);
            assert_eq!(out.pairs.len(), 13, "flip={flip}");
        }
    }

    fn match_wide_scene(
        scene: &crate::sky_scene::Scene,
        search_radius_deg: Option<f64>,
    ) -> Result<FieldMatch, SolveError> {
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let params = MatchParams {
            search_radius_deg,
            ..MatchParams::default()
        };
        match_field(
            &scene.detections(),
            &empty_image(scene.size.0, scene.size.1),
            &catalog,
            &crate::sky_scene::hint(),
            crate::sky_scene::SCALE_ARCSEC,
            crate::sky_scene::epoch(),
            &params,
        )
    }

    // The hint does not give the rotation, so the catalog stars for quads cover the circle
    // the strip can turn in: with the search margin, 34 times the strip's area.
    #[test]
    fn match_field_finds_a_strip() {
        let scene = crate::sky_scene::wide_field(5, 5600, 2.1, (6000, 200));
        assert_eq!(scene.truth.len(), 97);
        let out = match_wide_scene(&scene, None).unwrap_or_else(|e| panic!("{e:#}"));
        scene.assert_true_pairs(&out.pairs);
    }

    // A 1.5 degree search radius is a cone 22 times the area of the 0.57 degree square
    // image.
    #[test]
    fn match_field_finds_the_image_in_a_wide_search_cone() {
        let scene = crate::sky_scene::wide_field(7, 1365, 1.6, (1024, 1024));
        assert_eq!(scene.truth.len(), 45);
        let out = match_wide_scene(&scene, Some(1.5)).unwrap_or_else(|e| panic!("{e:#}"));
        scene.assert_true_pairs(&out.pairs);
    }

    #[test]
    fn match_field_verifies_against_stars_beyond_its_quad_stars() {
        // 12 quad stars on each side; the field holds 18 catalog stars.
        let (scene, out) = match_sparse_scene(false, 12);
        let out = out.unwrap_or_else(|e| panic!("{e:#}"));
        scene.assert_true_pairs(&out.pairs);
    }

    // Each input that is not finite, or not positive where it must be, is rejected before
    // the catalog is read.
    #[test]
    fn match_field_rejects_non_finite_inputs() {
        let scene = crate::sky_scene::star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let mut bad_star = stars.clone();
        bad_star[0] = det(f64::NAN, 10.0, 5000.0);
        let scale = crate::sky_scene::SCALE_ARCSEC;
        let epoch = crate::sky_scene::epoch();
        let run = |stars: &[DetectedStar], scale: f64, epoch: JulianDate, radius: Option<f64>| {
            let params = MatchParams {
                search_radius_deg: radius,
                ..MatchParams::default()
            };
            let image = empty_image(1024, 1024);
            let hint = crate::sky_scene::hint();
            let out = match_field(stars, &image, &catalog, &hint, scale, epoch, &params);
            out.err().map(|e| e.describe())
        };
        let cases = [
            (
                run(&stars, f64::NAN, epoch, None),
                "InvalidParameter: invalid parameter scale_arcsec: must be finite and positive, got NaN",
            ),
            (
                run(&stars, 0.0, epoch, None),
                "InvalidParameter: invalid parameter scale_arcsec: must be finite and positive, got 0",
            ),
            (
                run(&stars, scale, JulianDate::new(f64::NAN, 0.0), None),
                "InvalidParameter: invalid parameter epoch: must be finite, got (NaN, 0)",
            ),
            (
                run(&stars, scale, epoch, Some(f64::INFINITY)),
                "InvalidParameter: invalid parameter search_radius_deg: must be finite and positive, got inf",
            ),
            (
                run(&bad_star, scale, epoch, None),
                "InvalidParameter: invalid parameter stars: star 0 is not finite: pixel (NaN, 10), snr 20",
            ),
        ];
        let (got, want): (Vec<_>, Vec<_>) = cases
            .iter()
            .map(|(err, expected)| (err.as_deref(), Some(*expected)))
            .unzip();
        assert_eq!(got, want);
    }

    #[test]
    fn match_field_rejects_parameters_that_form_no_quad() {
        let scene = crate::sky_scene::star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let run = |stars: &[DetectedStar], params: MatchParams| {
            let image = empty_image(1024, 1024);
            let hint = crate::sky_scene::hint();
            let scale = crate::sky_scene::SCALE_ARCSEC;
            let epoch = crate::sky_scene::epoch();
            let out = match_field(stars, &image, &catalog, &hint, scale, epoch, &params);
            out.err().map(|e| e.describe())
        };
        let defaults = MatchParams::default;
        let got = [
            run(
                &stars,
                MatchParams {
                    max_stars: 3,
                    ..defaults()
                },
            ),
            run(
                &stars,
                MatchParams {
                    k_neighbors: 2,
                    ..defaults()
                },
            ),
            run(
                &stars,
                MatchParams {
                    search_radius_deg: Some(180.5),
                    ..defaults()
                },
            ),
            run(&[], defaults()),
        ];
        let want = [
            "InvalidParameter: invalid parameter max_stars: must be at least 4, the stars of one quad, got 3",
            "InvalidParameter: invalid parameter k_neighbors: must be at least 3, the neighbours a quad adds to its first star, got 2",
            "InvalidParameter: invalid parameter search_radius_deg: must be at most 180, got 180.5",
            "NoStars: no stars detected",
        ];
        assert_eq!(got.each_ref().map(|g| g.as_deref()), want.map(Some));
    }
}
