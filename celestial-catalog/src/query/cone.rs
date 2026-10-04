//! Cone search over a HEALPix-indexed star catalog.
//!
//! Given a sky position and radius, [`cone_search`] determines which HEALPix
//! pixels overlap the cone, scans only those pixels, filters by distance and
//! optional magnitude limit, and returns results sorted by angular distance.
//!
//! Proper motion can be propagated from the catalog epoch (J2016.0) to an
//! arbitrary observation epoch before matching.

use celestial_coords::errors::CoordResult;
use celestial_coords::proper_motion::CatalogStar;
use celestial_core::angle::Angle;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tt::TT;

use super::catalog::{Catalog, StarRecord};
use super::healpix::{angular_separation_deg, query_disc_nest};

/// Julian date of epoch J2016.0 (Gaia DR3 reference epoch).
const J2016_JD: f64 = 2457389.0;

/// Parameters for a cone search query.
#[derive(Debug, Clone)]
pub struct ConeSearchParams {
    /// Cone center right ascension, in degrees.
    pub ra_deg: f64,
    /// Cone center declination, in degrees.
    pub dec_deg: f64,
    /// Search radius, in degrees.
    pub radius_deg: f64,
    /// If set, exclude stars fainter than this magnitude.
    pub max_mag: Option<f64>,
    /// If set, return at most this many results (closest first).
    pub max_results: Option<usize>,
    /// If set, propagate proper motion from J2016.0 to this epoch before matching.
    pub epoch: Option<JulianDate>,
}

/// A single star returned from a cone search.
#[derive(Debug, Clone)]
pub struct ConeSearchResult {
    /// The original star record from the catalog.
    pub star: StarRecord,
    /// Right ascension used for matching (propagated if an epoch was given).
    pub ra_deg: f64,
    /// Declination used for matching (propagated if an epoch was given).
    pub dec_deg: f64,
    /// Angular distance from the search center, in degrees.
    pub distance_deg: f64,
}

/// Convenience wrapper that runs a cone search with proper-motion propagation.
///
/// Equivalent to calling [`cone_search`] with `epoch` set and
/// no magnitude or result-count limits.
pub fn cone_search_at_epoch(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    epoch: JulianDate,
) -> CoordResult<Vec<ConeSearchResult>> {
    let params = ConeSearchParams {
        ra_deg,
        dec_deg,
        radius_deg,
        max_mag: None,
        max_results: None,
        epoch: Some(epoch),
    };
    cone_search(catalog, &params)
}

/// Search for stars within a cone on the sky.
///
/// Identifies overlapping HEALPix pixels, scans their star lists, applies
/// optional proper-motion propagation and magnitude filtering, then returns
/// results sorted by angular distance from the cone center.
pub fn cone_search(
    catalog: &Catalog,
    params: &ConeSearchParams,
) -> CoordResult<Vec<ConeSearchResult>> {
    let nside = 1 << catalog.header().order;
    let mut results = Vec::new();
    for pixel in query_disc_nest(nside, params.ra_deg, params.dec_deg, params.radius_deg) {
        for star in catalog.stars_in_pixel(pixel) {
            if let Some(result) = match_star(star, params)? {
                results.push(result);
            }
        }
    }

    results.sort_by(|a, b| {
        a.distance_deg
            .partial_cmp(&b.distance_deg)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    if let Some(max_results) = params.max_results {
        results.truncate(max_results);
    }

    Ok(results)
}

fn match_star(
    star: &StarRecord,
    params: &ConeSearchParams,
) -> CoordResult<Option<ConeSearchResult>> {
    if params
        .max_mag
        .is_some_and(|max_mag| star.mag as f64 > max_mag)
    {
        return Ok(None);
    }
    let (ra_deg, dec_deg) = match params.epoch {
        Some(epoch_jd) => apply_proper_motion(star, epoch_jd)?,
        None => (star.ra, star.dec),
    };
    let distance_deg = angular_separation_deg(params.ra_deg, params.dec_deg, ra_deg, dec_deg);
    if distance_deg > params.radius_deg {
        return Ok(None);
    }
    Ok(Some(ConeSearchResult {
        star: *star,
        ra_deg,
        dec_deg,
        distance_deg,
    }))
}

/// Propagate a star's position from J2016.0 to `epoch_jd`.
fn apply_proper_motion(star: &StarRecord, epoch_jd: JulianDate) -> CoordResult<(f64, f64)> {
    // The catalog carries no radial velocities.
    let catalog_star = CatalogStar {
        ra: Angle::from_degrees(star.ra),
        dec: Angle::from_degrees(star.dec),
        pm_ra_cos_dec_mas_per_year: star.pmra,
        pm_dec_mas_per_year: star.pmdec,
        parallax_mas: star.parallax,
        radial_velocity_km_per_s: 0.0,
    };
    let from = TT::from_julian_date(JulianDate::new(J2016_JD, 0.0));
    let moved = catalog_star.propagate(&from, &TT::from_julian_date(epoch_jd))?;
    Ok((moved.ra.degrees(), moved.dec.degrees()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_angular_distance_same_point() {
        let dist = angular_separation_deg(0.0, 0.0, 0.0, 0.0);
        assert!((dist - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_angular_distance_90_degrees() {
        let dist = angular_separation_deg(0.0, 0.0, 90.0, 0.0);
        assert!((dist - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_angular_distance_pole_to_equator() {
        let dist = angular_separation_deg(0.0, 90.0, 0.0, 0.0);
        assert!((dist - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_angular_distance_antipodes() {
        let dist = angular_separation_deg(0.0, 0.0, 180.0, 0.0);
        assert!((dist - 180.0).abs() < 1e-10);
    }

    fn star_at_100_45(pm: f64) -> StarRecord {
        StarRecord {
            source_id: 1,
            ra: 100.0,
            dec: 45.0,
            pmra: pm,
            pmdec: pm,
            parallax: 0.0,
            mag: 5.0,
            flags: 0,
            _padding: 0,
        }
    }

    fn one_year_on() -> JulianDate {
        JulianDate::new(J2016_JD, 0.0).add_days(365.25)
    }

    // Expected places in these tests are eraPmsafe's. A star at rest comes back where it was,
    // apart from the last bit lost to the round trip through Cartesian coordinates.
    #[test]
    fn test_apply_proper_motion_zero_pm() {
        let place = apply_proper_motion(&star_at_100_45(0.0), one_year_on()).unwrap();
        assert_eq!(place, (100.0, 44.99999999999999));
    }

    #[test]
    fn test_apply_proper_motion_one_year() {
        let place = apply_proper_motion(&star_at_100_45(3600.0), one_year_on()).unwrap();
        assert_eq!(place, (100.0014142382452, 45.00099999127294));
    }

    #[test]
    fn test_apply_proper_motion_non_finite_is_err() {
        let star = star_at_100_45(f64::NAN);
        assert!(apply_proper_motion(&star, one_year_on()).is_err());
    }
}
