//! Cone search over a HEALPix-indexed star catalog.
//!
//! Given a sky position and radius, [`cone_search`] determines which HEALPix
//! pixels overlap the cone, scans only those pixels, filters by distance and
//! optional magnitude limit, and returns results sorted by angular distance.
//!
//! Proper motion can be propagated from the catalog epoch (J2016.0) to an
//! arbitrary observation epoch before matching.

use celestial_time::JulianDate;

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
) -> Vec<ConeSearchResult> {
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
pub fn cone_search(catalog: &Catalog, params: &ConeSearchParams) -> Vec<ConeSearchResult> {
    let order = catalog.header().order;
    let nside = 1 << order;

    let overlapping_pixels =
        query_disc_nest(nside, params.ra_deg, params.dec_deg, params.radius_deg);

    let mut results = Vec::new();

    for pixel in overlapping_pixels {
        let stars = catalog.stars_in_pixel(pixel);

        for star in stars {
            let (ra_obs, dec_obs) = if let Some(epoch_jd) = params.epoch {
                apply_proper_motion(star, epoch_jd)
            } else {
                (star.ra, star.dec)
            };

            let distance_deg =
                angular_separation_deg(params.ra_deg, params.dec_deg, ra_obs, dec_obs);

            if distance_deg > params.radius_deg {
                continue;
            }

            if let Some(max_mag) = params.max_mag {
                if star.mag as f64 > max_mag {
                    continue;
                }
            }

            results.push(ConeSearchResult {
                star: *star,
                ra_deg: ra_obs,
                dec_deg: dec_obs,
                distance_deg,
            });
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

    results
}

/// Linearly propagate a star's position from J2016.0 to `epoch_jd`.
fn apply_proper_motion(star: &StarRecord, epoch_jd: JulianDate) -> (f64, f64) {
    const MAS_PER_DEGREE: f64 = 3_600_000.0;

    let dt_years = (epoch_jd - JulianDate::new(J2016_JD, 0.0)).to_f64() / 365.25;

    let dec_obs = star.dec + star.pmdec * dt_years / MAS_PER_DEGREE;
    let cos_dec = libm::cos(star.dec * celestial_core::constants::PI / 180.0);
    let ra_obs = star.ra + star.pmra * dt_years / MAS_PER_DEGREE / cos_dec;

    (ra_obs, dec_obs)
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

    #[test]
    fn test_apply_proper_motion_zero_pm() {
        let star = StarRecord {
            source_id: 1,
            ra: 100.0,
            dec: 45.0,
            pmra: 0.0,
            pmdec: 0.0,
            parallax: 0.0,
            mag: 5.0,
            flags: 0,
            _padding: 0,
        };

        let (ra, dec) = apply_proper_motion(&star, JulianDate::new(J2016_JD, 0.0).add_days(365.25));
        assert!((ra - 100.0).abs() < 1e-10);
        assert!((dec - 45.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_proper_motion_one_year() {
        let star = StarRecord {
            source_id: 1,
            ra: 100.0,
            dec: 45.0,
            pmra: 3600.0,
            pmdec: 3600.0,
            parallax: 0.0,
            mag: 5.0,
            flags: 0,
            _padding: 0,
        };

        let (ra, dec) = apply_proper_motion(&star, JulianDate::new(J2016_JD, 0.0).add_days(365.25));

        // pmdec is sky rate, converts directly: 3600 mas/yr = 0.001 deg/yr
        let expected_dec = 45.0 + (3600.0 / 3_600_000.0);
        assert!((dec - expected_dec).abs() < 1e-10);

        // pmra is μα* = μα·cos(δ), so ΔRA = μα*/cos(δ) · Δt
        let cos_dec = libm::cos(45.0_f64 * celestial_core::constants::PI / 180.0);
        let expected_ra = 100.0 + (3600.0 / 3_600_000.0) / cos_dec;
        assert!((ra - expected_ra).abs() < 1e-10);
    }
}
