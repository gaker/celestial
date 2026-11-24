//! Quad-star construction for image detections and catalog cone searches.

use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::cone::{cone_search, ConeSearchParams, ConeSearchResult};
use celestial_catalog::query::quad::{tan_project_star, QuadStar};
use celestial_core::constants::{ARCSEC_TO_RAD, DEG_TO_RAD, HALF_PI, PI};
use celestial_time::julian::JulianDate;

use super::uniform::{brightest_in_cells, disc_cell_share, image_cell_px};
use super::MatchParams;
use crate::detect::{image_centre, DetectedStar};
use crate::error::SolveError;
use crate::log_odds::MAX_REFERENCES;

pub(crate) fn search_radius(
    w: usize,
    h: usize,
    scale_arcsec: f64,
    override_deg: Option<f64>,
) -> f64 {
    if let Some(r) = override_deg {
        return r;
    }
    field_diagonal_deg(w, h, scale_arcsec) / 2.0 * CONE_MARGIN
}

// How much larger than the image's half-diagonal a catalog cone about the field centre
// is, so that the cone still covers the image when the centre or the scale is off.
const CONE_MARGIN: f64 = 1.2;

// The image diagonal in degrees at `scale_arcsec` per pixel.
pub(crate) fn field_diagonal_deg(w: usize, h: usize, scale_arcsec: f64) -> f64 {
    libm::sqrt((w * w + h * h) as f64) * scale_arcsec / 3600.0
}

// A cone's stars, brightest first.
pub(super) struct ConeStars {
    // The `QuadSelection::count` stars brightest within their cells, with tangent-plane
    // positions about the cone centre.
    pub(super) quad_stars: Vec<QuadStar>,
    // The cone-search results of `quad_stars`, in the same order.
    pub(super) results: Vec<ConeSearchResult>,
    // (ra, dec) in degrees of the `n_references` brightest.
    pub(super) references: Vec<(f64, f64)>,
}

// Which catalog stars to build quads from: `count` stars, the brightest within square
// cells `cell_rad` wide on the tangent plane.
pub(super) struct QuadSelection {
    pub(super) count: usize,
    pub(super) cell_rad: f64,
}

// The query returns every star down to a magnitude limit, so it holds at least twice
// `selection.count` stars: a cell then has about twice the stars its sweeps take, and
// few cells run out of stars before the faint limit of the query.
pub(super) fn catalog_cone_search(
    catalog: &Catalog,
    (ra_deg, dec_deg): (f64, f64),
    radius_deg: f64,
    selection: &QuadSelection,
    n_references: usize,
    epoch: JulianDate,
) -> Result<ConeStars, SolveError> {
    let need = n_references.max(selection.count.saturating_mul(2));
    let results = brightest_in_cone(catalog, (ra_deg, dec_deg), radius_deg, need, epoch)?;
    let references = sky_positions(&results, n_references);
    let projected = project_about(results, (ra_deg, dec_deg));
    let points: Vec<(f64, f64)> = projected.iter().map(|(_, q)| (q.x, q.y)).collect();
    let cell = (selection.cell_rad, selection.cell_rad);
    let share = disc_cell_share(tangent_radius(radius_deg), selection.cell_rad);
    let (results, quad_stars) = brightest_in_cells(&points, cell, selection.count, share)
        .into_iter()
        .map(|i| projected[i].clone())
        .unzip();
    Ok(ConeStars {
        quad_stars,
        results,
        references,
    })
}

// Each result with a gnomonic projection about `(ra_deg, dec_deg)`, paired with its
// quad star at that tangent-plane position.
fn project_about(
    results: Vec<ConeSearchResult>,
    (ra_deg, dec_deg): (f64, f64),
) -> Vec<(ConeSearchResult, QuadStar)> {
    results
        .into_iter()
        .filter_map(|r| {
            let (x, y) = tan_project_star(r.ra_deg, r.dec_deg, ra_deg, dec_deg)?;
            let star = QuadStar {
                source_id: r.star.source_id,
                x,
                y,
                mag: r.star.mag,
            };
            Some((r, star))
        })
        .collect()
}

// The radius on the tangent plane of a cone of `radius_deg` about the tangent point:
// gnomonic projection puts a star `r` from the centre at tan(r). A cone of 90 degrees or
// more covers the whole plane.
fn tangent_radius(radius_deg: f64) -> f64 {
    let radius = radius_deg * DEG_TO_RAD;
    if radius < HALF_PI {
        libm::tan(radius)
    } else {
        f64::INFINITY
    }
}

// The (ra, dec) in degrees of the `n_references` brightest stars in a cone about
// (`ra_deg`, `dec_deg`) of `radius_deg`.
pub(crate) fn reference_stars(
    catalog: &Catalog,
    (ra_deg, dec_deg): (f64, f64),
    radius_deg: f64,
    n_references: usize,
    epoch: JulianDate,
) -> Result<Vec<(f64, f64)>, SolveError> {
    let results = brightest_in_cone(catalog, (ra_deg, dec_deg), radius_deg, n_references, epoch)?;
    Ok(sky_positions(&results, n_references))
}

// How many of a cone's brightest stars to take as references: as many as put
// `MAX_REFERENCES` inside an area the size of the image at uniform star density.
// `scale_arcsec` is the plate scale in arcseconds per pixel.
pub(crate) fn reference_count(w: usize, h: usize, scale_arcsec: f64, radius_deg: f64) -> usize {
    cone_share(MAX_REFERENCES, w, h, scale_arcsec, radius_deg)
}

// The number of stars a cone of `radius_deg` holds at the density of `n_image` stars over
// the image, and at least `n_image`.
pub(super) fn cone_share(
    n_image: usize,
    w: usize,
    h: usize,
    scale_arcsec: f64,
    radius_deg: f64,
) -> usize {
    let image_deg2 = (w * h) as f64 * (scale_arcsec / 3600.0) * (scale_arcsec / 3600.0);
    let cone_deg2 = PI * radius_deg * radius_deg;
    let ratio = (cone_deg2 / image_deg2).max(1.0);
    let count = libm::ceil(n_image as f64 * ratio);
    if count.is_finite() {
        count as usize
    } else {
        usize::MAX
    }
}

// The first magnitude limit `brightest_in_cone` queries with: the naked-eye limit, with
// about 5000 stars over the whole sky.
const FIRST_MAG_LIMIT: f64 = 6.0;

// Each step of the limit about doubles the number of stars a query returns.
const MAG_LIMIT_STEP: f64 = 1.0;

// No limit is tried fainter than this, whatever a catalog header says, so the steps end;
// Gaia reaches about magnitude 21.
const FAINTEST_MAG_LIMIT: f64 = 30.0;

// Every star of the cone down to the brightest magnitude limit, in steps of
// `MAG_LIMIT_STEP`, that gives at least `need` stars, or the whole cone when no limit
// brighter than the catalog's own does; brightest first. The `need` brightest are the
// same as those of the whole cone.
fn brightest_in_cone(
    catalog: &Catalog,
    (ra_deg, dec_deg): (f64, f64),
    radius_deg: f64,
    need: usize,
    epoch: JulianDate,
) -> Result<Vec<ConeSearchResult>, SolveError> {
    let catalog_limit = f64::from(catalog.header().mag_limit).min(FAINTEST_MAG_LIMIT);
    let mut max_mag = FIRST_MAG_LIMIT;
    loop {
        let limit = (max_mag < catalog_limit).then_some(max_mag);
        let params = ConeSearchParams {
            ra_deg,
            dec_deg,
            radius_deg,
            max_mag: limit,
            max_results: None,
            epoch: Some(epoch),
        };
        let mut results = cone_search(catalog, &params)?;
        if results.len() >= need || limit.is_none() {
            results.sort_by(|a, b| a.star.mag.total_cmp(&b.star.mag));
            return Ok(results);
        }
        max_mag += MAG_LIMIT_STEP;
    }
}

fn sky_positions(results: &[ConeSearchResult], n: usize) -> Vec<(f64, f64)> {
    results
        .iter()
        .take(n)
        .map(|r| (r.ra_deg, r.dec_deg))
        .collect()
}

// The sign of det(CD), the determinant of the map from pixel offsets to the tangent
// plane, that one search for a quad match assumes. `WcsSolution::is_flipped` is true for
// `Negative`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CdSign {
    Positive,
    Negative,
}

impl CdSign {
    pub(super) fn label(self) -> &'static str {
        match self {
            CdSign::Positive => "det(CD) > 0",
            CdSign::Negative => "det(CD) < 0",
        }
    }
}

// The detections for quads, on the tangent plane through the map of scale `scale_arcsec`
// whose determinant has the sign `cd_sign`, brightest first. Pixel x grows to the west in
// both signs, so the sign is that of the y axis.
pub(super) fn stars_to_quad_stars(
    stars: &[DetectedStar],
    (w, h): (usize, usize),
    scale_arcsec: f64,
    params: &MatchParams,
    cd_sign: CdSign,
) -> Vec<QuadStar> {
    let (cx, cy) = image_centre(w, h);
    let rad_per_px = scale_arcsec * ARCSEC_TO_RAD;
    let y_sign = match cd_sign {
        CdSign::Positive => -1.0,
        CdSign::Negative => 1.0,
    };
    quad_detections(stars, (w, h), params)
        .into_iter()
        .map(|i| {
            let s = &stars[i];
            QuadStar {
                source_id: i as i64,
                x: -(s.x - cx) * rad_per_px,
                y: y_sign * (s.y - cy) * rad_per_px,
                mag: -2.5 * libm::log10(s.flux.max(1.0)) as f32,
            }
        })
        .collect()
}

// Indices of the `params.max_stars` detections that are brightest within their cells,
// brightest first. The cells hold about `params.k_neighbors` of them each, as the
// catalog's cells do (`catalog_side::catalog_selection`), so both lists sample the sky at
// the same density when bright stars or false detections crowd one part of the image.
fn quad_detections(
    stars: &[DetectedStar],
    (w, h): (usize, usize),
    params: &MatchParams,
) -> Vec<usize> {
    let mut by_flux: Vec<usize> = (0..stars.len()).collect();
    by_flux.sort_by(|&a, &b| {
        stars[b]
            .flux
            .partial_cmp(&stars[a].flux)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let n = by_flux.len().min(params.max_stars);
    let cell = image_cell_px(w, h, n, params.k_neighbors);
    let points: Vec<(f64, f64)> = by_flux.iter().map(|&i| (stars[i].x, stars[i].y)).collect();
    brightest_in_cells(&points, cell, n, |_| 1.0)
        .into_iter()
        .map(|k| by_flux[k])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn max_stars(max_stars: usize) -> MatchParams {
        MatchParams {
            max_stars,
            ..MatchParams::default()
        }
    }

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

    #[test]
    fn search_radius_returns_override_when_present() {
        let r = search_radius(1000, 1000, 1.5, Some(2.5));
        assert_eq!(r, 2.5);
    }

    #[test]
    fn search_radius_uses_diagonal_with_safety_factor() {
        // diag = sqrt(2) * 1000, scale 3600 "/px → exactly 1 deg/px
        // radius = diag * 1 / 2 * 1.2 = diag * 0.6
        let r = search_radius(1000, 1000, 3600.0, None);
        let expected = libm::sqrt(2.0) * 1000.0 * 0.6;
        assert_eq!(r, expected);
    }

    #[test]
    fn search_radius_scales_linearly_with_plate_scale() {
        let r1 = search_radius(1000, 1000, 1.0, None);
        let r2 = search_radius(1000, 1000, 2.0, None);
        assert_eq!(r2, r1 * 2.0);
    }

    // Four times the width and height give four times the radius exactly: the square
    // root and each later step scale by a power of two without rounding.
    #[test]
    fn search_radius_scales_linearly_with_image_size() {
        let small = search_radius(500, 500, 1.0, None);
        let big = search_radius(2000, 2000, 1.0, None);
        assert_eq!(big, small * 4.0);
    }

    // A 3000 x 4000 image has a 5000 px diagonal, which at 36"/px is 50 degrees.
    #[test]
    fn field_diagonal_deg_is_the_pixel_diagonal_times_the_scale() {
        assert_eq!(field_diagonal_deg(3000, 4000, 36.0), 50.0);
    }

    #[test]
    fn stars_to_quad_stars_sorts_by_descending_flux() {
        let stars = vec![
            det(10.0, 10.0, 100.0),
            det(20.0, 20.0, 500.0),
            det(30.0, 30.0, 300.0),
        ];
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        assert_eq!(out.len(), 3);
        // source_id preserves original index in `stars`
        assert_eq!(out[0].source_id, 1);
        assert_eq!(out[1].source_id, 2);
        assert_eq!(out[2].source_id, 0);
    }

    // 40 bright detections in the top-left quarter of a 1000 px square and 10 fainter
    // ones in each other quarter. 40 stars in cells of 10 make a grid of 2 x 2 cells of
    // 500 px, and each cell gives 10.
    #[test]
    fn stars_to_quad_stars_takes_the_brightest_of_each_cell() {
        let mut stars: Vec<DetectedStar> = (0..40)
            .map(|i| det(10.0 + 10.0 * i as f64, 20.0, 10_000.0 - i as f64))
            .collect();
        for (qx, qy) in [(500.0, 0.0), (0.0, 500.0), (500.0, 500.0)] {
            stars.extend(
                (0..10).map(|i| det(qx + 10.0 + 40.0 * i as f64, qy + 50.0, 100.0 - i as f64)),
            );
        }
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(40), CdSign::Positive);
        let mut per_quarter = [0; 4];
        for q in &out {
            let s = &stars[q.source_id as usize];
            per_quarter[usize::from(s.x >= 500.0) + 2 * usize::from(s.y >= 500.0)] += 1;
        }
        assert_eq!(per_quarter, [10, 10, 10, 10]);
    }

    #[test]
    fn stars_to_quad_stars_truncates_to_max_stars() {
        let stars: Vec<_> = (0..20)
            .map(|i| det(i as f64, i as f64, (20 - i) as f64 * 10.0))
            .collect();
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(5), CdSign::Positive);
        assert_eq!(out.len(), 5);
    }

    #[test]
    fn stars_to_quad_stars_centers_coordinates() {
        let stars = vec![det(500.0, 500.0, 100.0)];
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        // At image center: (500-500)=0 → both coords zero regardless of the CD sign
        assert_eq!(out[0].x, 0.0);
        assert_eq!(out[0].y, 0.0);
    }

    #[test]
    fn stars_to_quad_stars_flips_y_with_the_cd_sign() {
        let stars = vec![det(500.0, 600.0, 100.0)];
        let positive =
            stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        let negative =
            stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Negative);
        // X unchanged; Y sign inverts
        assert_eq!(positive[0].x, negative[0].x);
        assert_eq!(positive[0].y, -negative[0].y);
    }

    #[test]
    fn stars_to_quad_stars_inverts_x_axis() {
        // Positive x offset from center must yield negative output x (RA increases east).
        let stars = vec![det(600.0, 500.0, 100.0)];
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        assert!(out[0].x < 0.0);
    }

    #[test]
    fn stars_to_quad_stars_mag_from_flux_log10() {
        let stars = vec![det(500.0, 500.0, 100.0)];
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        // -2.5 * log10(100) = -5.0
        assert_eq!(out[0].mag, -5.0);
    }

    #[test]
    fn stars_to_quad_stars_mag_clamps_flux_floor() {
        // Flux below 1 must not explode log10; clamped to 1 → mag = 0.
        let stars = vec![det(500.0, 500.0, 0.1)];
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        assert_eq!(out[0].mag, 0.0);
    }

    #[test]
    fn stars_to_quad_stars_plate_scale_scales_offset() {
        let stars = vec![det(600.0, 500.0, 100.0)];
        let a = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        let b = stars_to_quad_stars(&stars, (1000, 1000), 2.0, &max_stars(10), CdSign::Positive);
        // Offset in radians is linear in arcsec/px.
        assert_eq!(b[0].x, a[0].x * 2.0);
    }

    #[test]
    fn stars_to_quad_stars_empty_input_yields_empty() {
        let stars: Vec<DetectedStar> = Vec::new();
        let out = stars_to_quad_stars(&stars, (1000, 1000), 1.0, &max_stars(10), CdSign::Positive);
        assert!(out.is_empty());
    }

    use celestial_catalog::query::catalog::Catalog;
    use celestial_catalog::query::quad::tan_deproject_star;
    use celestial_core::constants::RAD_TO_DEG;
    use celestial_time::julian::JulianDate;

    use crate::sky_scene::{build, SynthStar};

    // The four cells about the cone centre, one per quadrant, each hold a quarter of the
    // cone, so the selection is the `count` brightest stars.
    fn brightest(count: usize) -> QuadSelection {
        QuadSelection {
            count,
            cell_rad: 2.0 * PI,
        }
    }

    // Cells 0.01 rad wide in a cone whose tangent-plane radius is 0.025 rad. Cell (0, 0)
    // lies inside the cone; 0.4316 of the area of cell (2, 0) does. Eight stars in each,
    // alternating in brightness, give the full cell 6 of the 8 quad stars: its sweeps
    // rank 1 to 8, and those of the edge cell rank 2.32, 4.63, 6.95 and so on.
    #[test]
    fn catalog_cone_search_gives_a_cell_the_share_of_its_area_inside_the_cone() {
        let centre = (100.0, 20.0);
        let stars: Vec<SynthStar> = (0..16)
            .map(|i| {
                let k = (i / 2) as f64;
                let (x, y) = if i % 2 == 0 {
                    (0.001 + 0.001 * k, 0.001 + 0.001 * k)
                } else {
                    (0.0201 + 0.0002 * k, 0.001 + 0.0005 * k)
                };
                let (ra, dec) = tan_deproject_star(x, y, centre.0, centre.1);
                SynthStar {
                    source_id: i as i64,
                    ra,
                    dec,
                    mag: 8.0 + 0.1 * i as f32,
                }
            })
            .collect();
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let selection = QuadSelection {
            count: 8,
            cell_rad: 0.01,
        };
        let radius_deg = libm::atan(0.025) * RAD_TO_DEG;
        let epoch = JulianDate::new(2451545.0, 0.0);
        let cone = catalog_cone_search(&catalog, centre, radius_deg, &selection, 0, epoch).unwrap();
        let in_full_cell = cone.quad_stars.iter().filter(|q| q.x < 0.01).count();
        assert_eq!((in_full_cell, cone.quad_stars.len() - in_full_cell), (6, 2));
    }

    #[test]
    fn catalog_cone_search_returns_stars_in_radius_sorted_by_mag() {
        // Cluster of stars within 0.1° of the center; build a catalog and verify
        // the cone search returns them, magnitude-sorted, truncated to max_stars.
        let center_ra = 270.0;
        let center_dec = -30.0;
        let mut stars = Vec::new();
        for i in 0..20 {
            let i = i as f64;
            stars.push(SynthStar {
                source_id: 5000 + i as i64,
                ra: center_ra + 0.02 * (i - 10.0),
                dec: center_dec + 0.01 * (i - 5.0),
                mag: 10.0 - i as f32 * 0.1, // brightest first when sorted
            });
        }
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let ConeStars {
            quad_stars,
            results,
            ..
        } = catalog_cone_search(
            &catalog,
            (center_ra, center_dec),
            1.0,
            &brightest(8),
            0,
            epoch,
        )
        .unwrap();

        // truncate(max_stars) → at most 8 entries.
        assert!(quad_stars.len() <= 8);
        assert_eq!(quad_stars.len(), results.len());
        // Mag-sorted ascending → check monotonicity of results.
        for w in results.windows(2) {
            assert!(w[0].star.mag <= w[1].star.mag);
        }
        // QuadStar source_id mirrors the catalog source_id.
        for (qs, r) in quad_stars.iter().zip(results.iter()) {
            assert_eq!(qs.source_id, r.star.source_id);
        }
    }

    #[test]
    fn catalog_cone_search_empty_far_field_returns_empty() {
        // Stars on one side of the sky, search the other side → no hits.
        let stars = vec![SynthStar {
            source_id: 1,
            ra: 10.0,
            dec: 0.0,
            mag: 9.0,
        }];
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let ConeStars {
            quad_stars: qs,
            results: res,
            ..
        } = catalog_cone_search(&catalog, (200.0, 0.0), 0.5, &brightest(50), 0, epoch).unwrap();
        assert!(qs.is_empty());
        assert!(res.is_empty());
    }

    #[test]
    fn catalog_cone_search_skips_stars_that_fail_tangent_projection() {
        // tan_project_star returns None for stars on the far side of the sphere
        // (denom <= 0). Searching with a huge radius forces such stars into the
        // result list, where they must be filtered out of quad_stars.
        let center_ra = 0.0;
        let center_dec = 0.0;
        let stars = vec![
            SynthStar {
                source_id: 1,
                ra: 0.0,
                dec: 0.0,
                mag: 8.0,
            },
            SynthStar {
                source_id: 2,
                ra: 180.0,
                dec: 0.0,
                mag: 9.0,
            },
        ];
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let ConeStars {
            quad_stars: qs,
            results: res,
            ..
        } = catalog_cone_search(
            &catalog,
            (center_ra, center_dec),
            180.0,
            &brightest(50),
            0,
            epoch,
        )
        .unwrap();
        // Both stars are within the cone, but only the near one projects.
        assert!(!res.is_empty());
        assert!(qs.len() < res.len() || qs.iter().all(|s| s.source_id != 2));
    }
}
