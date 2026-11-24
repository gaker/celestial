// The catalog stars the refinement re-matches against.

use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::cone::{cone_search, ConeSearchParams};
use celestial_time::julian::JulianDate;

use crate::error::SolveError;
use crate::match_field::quad_stars::search_radius;

use super::super::solution::WcsSolution;

pub(super) fn prep_catalog_for_refine(
    initial: &WcsSolution,
    catalog: &Catalog,
    epoch: JulianDate,
    target_stars: usize,
) -> Result<(Vec<CatalogEntry>, f64), SolveError> {
    let (w, h) = (initial.width, initial.height);
    let radius_deg = search_radius(w, h, initial.mean_scale_arcsec(), None);
    let (catalog_stars, mag_limit) = search_optimal_magnitude(
        catalog,
        initial.crval1,
        initial.crval2,
        radius_deg,
        epoch,
        target_stars,
    )?;
    log::debug!(
        "refine_wcs: {} catalog stars to mag {:.2} within {:.3} deg (target {})",
        catalog_stars.len(),
        mag_limit,
        radius_deg,
        target_stars,
    );
    Ok((catalog_stars, mag_limit))
}

pub(super) struct CatalogEntry {
    pub(super) ra_deg: f64,
    pub(super) dec_deg: f64,
}

fn cone_search_at_mag(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    max_mag: f64,
    epoch: JulianDate,
) -> Result<Vec<CatalogEntry>, SolveError> {
    let params = ConeSearchParams {
        ra_deg,
        dec_deg,
        radius_deg,
        max_mag: Some(max_mag),
        max_results: None,
        epoch: Some(epoch),
    };
    Ok(cone_search(catalog, &params)?
        .into_iter()
        .map(|r| CatalogEntry {
            ra_deg: r.ra_deg,
            dec_deg: r.dec_deg,
        })
        .collect())
}

fn search_optimal_magnitude(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    epoch: JulianDate,
    target: usize,
) -> Result<(Vec<CatalogEntry>, f64), SolveError> {
    let (mut lo, mut hi) = (8.0_f64, 20.0_f64);
    for _ in 0..8 {
        let mid = (lo + hi) / 2.0;
        let n = cone_search_at_mag(catalog, ra_deg, dec_deg, radius_deg, mid, epoch)?.len();
        log::debug!("  mag search: m={:.2}, {} stars", mid, n);
        if n < target {
            lo = mid;
        } else {
            hi = mid;
        }
        if libm::fabs(hi - lo) < 0.1 {
            break;
        }
    }
    let mag = hi;
    let stars = cone_search_at_mag(catalog, ra_deg, dec_deg, radius_deg, mag, epoch)?;
    Ok((stars, mag))
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::sky_scene::SynthStar;

    use super::super::test_support::synth_catalog;

    #[test]
    fn cone_search_at_mag_filters_by_magnitude() {
        // 3 stars: mag 8, 12, 18. Search with max_mag=10 should return only the
        // bright one; max_mag=20 should return all three.
        let stars = vec![
            SynthStar {
                source_id: 1,
                ra: 0.0,
                dec: 0.0,
                mag: 8.0,
            },
            SynthStar {
                source_id: 2,
                ra: 0.01,
                dec: 0.0,
                mag: 12.0,
            },
            SynthStar {
                source_id: 3,
                ra: 0.0,
                dec: 0.01,
                mag: 18.0,
            },
        ];
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let bright_only = cone_search_at_mag(&c, 0.0, 0.0, 1.0, 10.0, epoch).unwrap();
        assert_eq!(bright_only.len(), 1);
        assert_eq!((bright_only[0].ra_deg, bright_only[0].dec_deg), (0.0, 0.0));

        let all = cone_search_at_mag(&c, 0.0, 0.0, 1.0, 20.0, epoch).unwrap();
        assert_eq!(all.len(), 3);
    }

    #[test]
    fn cone_search_at_mag_returns_empty_when_no_stars_in_cone() {
        let stars = vec![SynthStar {
            source_id: 1,
            ra: 10.0,
            dec: 0.0,
            mag: 8.0,
        }];
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let out = cone_search_at_mag(&c, 200.0, 0.0, 0.5, 20.0, epoch).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn search_optimal_magnitude_converges_when_target_reachable() {
        // 50 stars spread across magnitudes 8-18. Target 10 → search should
        // settle on a magnitude that yields ≥10 stars.
        let stars: Vec<_> = (0..50)
            .map(|i| SynthStar {
                source_id: i,
                ra: (i as f64) * 0.005,
                dec: 0.0,
                mag: 8.0 + (i as f32) * 0.2,
            })
            .collect();
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let (results, mag) = search_optimal_magnitude(&c, 0.0, 0.0, 1.0, epoch, 10).unwrap();
        assert!(
            results.len() >= 10,
            "got {} stars at mag {mag}",
            results.len()
        );
        assert!(
            (8.0..=20.0).contains(&mag),
            "mag {mag} outside expected range"
        );
    }

    #[test]
    fn search_optimal_magnitude_returns_high_cap_when_target_unreachable() {
        // 2 stars in a target-of-100 search → never reaches target; mag should
        // hit the hi cap (20.0).
        let stars = vec![
            SynthStar {
                source_id: 1,
                ra: 0.0,
                dec: 0.0,
                mag: 8.0,
            },
            SynthStar {
                source_id: 2,
                ra: 0.01,
                dec: 0.0,
                mag: 10.0,
            },
        ];
        let (c, _f) = synth_catalog(stars);
        let epoch = JulianDate::new(2451545.0, 0.0);

        let (results, mag) = search_optimal_magnitude(&c, 0.0, 0.0, 1.0, epoch, 100).unwrap();
        assert_eq!(results.len(), 2);
        // The magnitude search stops at its cap.
        assert_eq!(mag, 20.0);
    }
}
