//! Quad-star construction for image detections and catalog cone searches.

use celestial_catalog::query::{
    cone_search, tan_project_star, Catalog, ConeSearchParams, ConeSearchResult, QuadStar,
};
use celestial_core::constants::PI;
use celestial_time::JulianDate;

use crate::detect::DetectedStar;

pub(super) fn search_radius(
    w: usize,
    h: usize,
    scale_arcsec: f64,
    override_deg: Option<f64>,
) -> f64 {
    if let Some(r) = override_deg {
        return r;
    }
    let diag = libm::sqrt((w * w + h * h) as f64);
    diag * scale_arcsec / 3600.0 / 2.0 * 1.2
}

pub(super) fn catalog_cone_search(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    max_stars: usize,
    epoch: JulianDate,
) -> (Vec<QuadStar>, Vec<ConeSearchResult>) {
    let params = ConeSearchParams {
        ra_deg,
        dec_deg,
        radius_deg,
        max_mag: None,
        max_results: None,
        epoch: Some(epoch),
    };
    let mut results = cone_search(catalog, &params);
    results.sort_by(|a, b| {
        a.star.mag.partial_cmp(&b.star.mag).unwrap_or(std::cmp::Ordering::Equal)
    });
    results.truncate(max_stars);

    let stars: Vec<QuadStar> = results
        .iter()
        .filter_map(|r| {
            tan_project_star(r.ra_deg, r.dec_deg, ra_deg, dec_deg)
                .map(|(x, y)| QuadStar {
                    source_id: r.star.source_id,
                    x,
                    y,
                    mag: r.star.mag,
                })
        })
        .collect();

    (stars, results)
}

pub(super) fn stars_to_quad_stars(
    stars: &[DetectedStar],
    w: usize,
    h: usize,
    scale_arcsec: f64,
    max_stars: usize,
    parity_flip: bool,
) -> Vec<QuadStar> {
    let cx = w as f64 / 2.0;
    let cy = h as f64 / 2.0;
    let rad_per_px = scale_arcsec * PI / (180.0 * 3600.0);
    let y_sign = if parity_flip { 1.0 } else { -1.0 };

    let mut sorted: Vec<_> = stars.iter().enumerate().collect();
    sorted.sort_by(|a, b| {
        b.1.flux.partial_cmp(&a.1.flux).unwrap_or(std::cmp::Ordering::Equal)
    });
    sorted.truncate(max_stars);

    sorted
        .into_iter()
        .map(|(i, s)| QuadStar {
            source_id: i as i64,
            x: -(s.x - cx) * rad_per_px,
            y: y_sign * (s.y - cy) * rad_per_px,
            mag: -2.5 * libm::log10(s.flux.max(1.0)) as f32,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn det(x: f64, y: f64, flux: f64) -> DetectedStar {
        DetectedStar {
            x, y, flux,
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

    #[test]
    fn stars_to_quad_stars_sorts_by_descending_flux() {
        let stars = vec![
            det(10.0, 10.0, 100.0),
            det(20.0, 20.0, 500.0),
            det(30.0, 30.0, 300.0),
        ];
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        assert_eq!(out.len(), 3);
        // source_id preserves original index in `stars`
        assert_eq!(out[0].source_id, 1);
        assert_eq!(out[1].source_id, 2);
        assert_eq!(out[2].source_id, 0);
    }

    #[test]
    fn stars_to_quad_stars_truncates_to_max_stars() {
        let stars: Vec<_> = (0..20)
            .map(|i| det(i as f64, i as f64, (20 - i) as f64 * 10.0))
            .collect();
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 5, false);
        assert_eq!(out.len(), 5);
    }

    #[test]
    fn stars_to_quad_stars_centers_coordinates() {
        let stars = vec![det(500.0, 500.0, 100.0)];
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        // At image center: (500-500)=0 → both coords zero regardless of parity
        assert_eq!(out[0].x, 0.0);
        assert_eq!(out[0].y, 0.0);
    }

    #[test]
    fn stars_to_quad_stars_flips_y_with_parity() {
        let stars = vec![det(500.0, 600.0, 100.0)];
        let normal = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        let flipped = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, true);
        // X unchanged; Y sign inverts
        assert_eq!(normal[0].x, flipped[0].x);
        assert_eq!(normal[0].y, -flipped[0].y);
    }

    #[test]
    fn stars_to_quad_stars_inverts_x_axis() {
        // Positive x offset from center must yield negative output x (RA increases east).
        let stars = vec![det(600.0, 500.0, 100.0)];
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        assert!(out[0].x < 0.0);
    }

    #[test]
    fn stars_to_quad_stars_mag_from_flux_log10() {
        let stars = vec![det(500.0, 500.0, 100.0)];
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        // -2.5 * log10(100) = -5.0
        assert_eq!(out[0].mag, -5.0);
    }

    #[test]
    fn stars_to_quad_stars_mag_clamps_flux_floor() {
        // Flux below 1 must not explode log10; clamped to 1 → mag = 0.
        let stars = vec![det(500.0, 500.0, 0.1)];
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        assert_eq!(out[0].mag, 0.0);
    }

    #[test]
    fn stars_to_quad_stars_plate_scale_scales_offset() {
        let stars = vec![det(600.0, 500.0, 100.0)];
        let a = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        let b = stars_to_quad_stars(&stars, 1000, 1000, 2.0, 10, false);
        // Offset in radians is linear in arcsec/px.
        assert_eq!(b[0].x, a[0].x * 2.0);
    }

    #[test]
    fn stars_to_quad_stars_empty_input_yields_empty() {
        let stars: Vec<DetectedStar> = Vec::new();
        let out = stars_to_quad_stars(&stars, 1000, 1000, 1.0, 10, false);
        assert!(out.is_empty());
    }

    use celestial_catalog::query::Catalog;
    use celestial_time::JulianDate;

    use super::super::test_catalog::{build, SynthStar};

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

        let (quad_stars, results) = catalog_cone_search(
            &catalog, center_ra, center_dec, 1.0, 8, epoch,
        );

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
        let stars = vec![SynthStar { source_id: 1, ra: 10.0, dec: 0.0, mag: 9.0 }];
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let (qs, res) = catalog_cone_search(&catalog, 200.0, 0.0, 0.5, 50, epoch);
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
            SynthStar { source_id: 1, ra: 0.0, dec: 0.0, mag: 8.0 },
            SynthStar { source_id: 2, ra: 180.0, dec: 0.0, mag: 9.0 },
        ];
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let (qs, res) = catalog_cone_search(&catalog, center_ra, center_dec, 180.0, 50, epoch);
        // Both stars are within the cone, but only the near one projects.
        assert!(res.len() >= 1);
        assert!(qs.len() < res.len() || qs.iter().all(|s| s.source_id != 2));
    }
}
