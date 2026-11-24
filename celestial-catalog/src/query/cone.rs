//! Cone search over a HEALPix-indexed star catalog.
//!
//! Given a sky position and radius, [`cone_search`] determines which HEALPix
//! pixels overlap the cone, scans only those pixels, filters by distance and
//! optional magnitude limit, and returns results sorted by angular distance.
//!
//! Proper motion can be propagated from the catalog epoch in the header to an
//! arbitrary observation epoch before matching.

use celestial_coords::proper_motion::CatalogStar;
use celestial_core::angle::Angle;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tt::TT;

use super::catalog::{Catalog, StarRecord};
use super::disc::query_disc_nest;
use super::healpix::{angular_separation_deg, reduce_ra_deg};
use crate::error::{CatalogError, CatalogResult};

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
    /// If set, propagate proper motion from the catalog epoch to this epoch before matching.
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

/// Search for stars within a cone on the sky.
///
/// Identifies overlapping HEALPix pixels, scans their star lists, applies
/// optional proper-motion propagation and magnitude filtering, then returns
/// results sorted by angular distance from the cone center.
pub fn cone_search(
    catalog: &Catalog,
    params: &ConeSearchParams,
) -> CatalogResult<Vec<ConeSearchResult>> {
    check_epoch(params.epoch)?;
    let params = &checked_params(params)?;
    let mut results = stars_in_cone(catalog, params)?;
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

fn stars_in_cone(
    catalog: &Catalog,
    params: &ConeSearchParams,
) -> CatalogResult<Vec<ConeSearchResult>> {
    let header = catalog.header();
    // The header epoch is a Julian year, read as TT. Catalog::open checked that it is finite.
    let propagation = params.epoch.map(|to| {
        let from = JulianDate::from_julian_year(header.epoch);
        (TT::from_julian_date(from), TT::from_julian_date(to))
    });
    let (ra, dec, radius) = (params.ra_deg, params.dec_deg, params.radius_deg);
    let mut results = Vec::new();
    for pixel in query_disc_nest(header.order, ra, dec, radius)? {
        for star in catalog.stars_in_pixel(pixel) {
            results.extend(match_star(star, params, propagation));
        }
    }
    Ok(results)
}

// Returns a copy of `params` with RA reduced into [0, 360) and the radius clamped to 180
// degrees: a larger cone is the whole sky.
fn checked_params(params: &ConeSearchParams) -> CatalogResult<ConeSearchParams> {
    check_finite(params)?;
    if !(-90.0..=90.0).contains(&params.dec_deg) {
        return Err(invalid_cone(format!(
            "Dec {} is outside [-90, 90] degrees",
            params.dec_deg
        )));
    }
    if params.radius_deg < 0.0 {
        return Err(invalid_cone(format!(
            "radius {} is negative",
            params.radius_deg
        )));
    }
    Ok(ConeSearchParams {
        ra_deg: reduce_ra_deg(params.ra_deg),
        radius_deg: params.radius_deg.min(180.0),
        ..params.clone()
    })
}

fn check_finite(params: &ConeSearchParams) -> CatalogResult<()> {
    let values = [params.ra_deg, params.dec_deg, params.radius_deg];
    if values
        .iter()
        .chain(params.max_mag.as_ref())
        .any(|v| !v.is_finite())
    {
        return Err(invalid_cone(format!(
            "RA {}, Dec {}, radius {} and max_mag {:?} must be finite",
            params.ra_deg, params.dec_deg, params.radius_deg, params.max_mag
        )));
    }
    Ok(())
}

// match_star skips a star whose propagation fails, so an epoch that would fail every star must
// be an error here. Propagation fails for any star when the epoch's two parts sum to a value
// that is not finite.
fn check_epoch(epoch: Option<JulianDate>) -> CatalogResult<()> {
    match epoch {
        Some(e) if !(e.jd1() + e.jd2()).is_finite() => Err(invalid_cone(format!(
            "epoch JD {} + {} is not finite",
            e.jd1(),
            e.jd2()
        ))),
        _ => Ok(()),
    }
}

fn invalid_cone(message: String) -> CatalogError {
    CatalogError::invalid_input(format!("cone search: {message}"))
}

// Returns None for a star outside the cone or fainter than max_mag, and for a corrupt record:
// a place that is not finite or has |Dec| above 90, or a motion that cannot be propagated.
// Checking each candidate is cheaper than checking every record at open.
fn match_star(
    star: &StarRecord,
    params: &ConeSearchParams,
    propagation: Option<(TT, TT)>,
) -> Option<ConeSearchResult> {
    let too_faint = params
        .max_mag
        .is_some_and(|max_mag| star.mag as f64 > max_mag);
    if too_faint || !is_valid_place(star.ra, star.dec) {
        return None;
    }
    // check_epoch rejected an epoch that fails every star, so an error here comes from the
    // record: a motion that is not finite, or a space velocity above half the speed of light.
    let (ra_deg, dec_deg) = match propagation {
        Some((from, to)) => apply_proper_motion(star, &from, &to).ok()?,
        None => (star.ra, star.dec),
    };
    let distance_deg = angular_separation_deg(params.ra_deg, params.dec_deg, ra_deg, dec_deg);
    let outside = distance_deg.is_nan() || distance_deg > params.radius_deg;
    if outside || !is_valid_place(ra_deg, dec_deg) {
        return None;
    }
    Some(ConeSearchResult {
        star: *star,
        ra_deg,
        dec_deg,
        distance_deg,
    })
}

fn is_valid_place(ra_deg: f64, dec_deg: f64) -> bool {
    ra_deg.is_finite() && (-90.0..=90.0).contains(&dec_deg)
}

/// Propagate a star's position from the catalog epoch `from` to `to`.
fn apply_proper_motion(star: &StarRecord, from: &TT, to: &TT) -> CatalogResult<(f64, f64)> {
    // The catalog carries no radial velocities.
    let catalog_star = CatalogStar {
        ra: Angle::from_degrees(star.ra),
        dec: Angle::from_degrees(star.dec),
        pm_ra_cos_dec_mas_per_year: star.pmra,
        pm_dec_mas_per_year: star.pmdec,
        parallax_mas: star.parallax,
        radial_velocity_km_per_s: 0.0,
    };
    let moved = catalog_star.propagate(from, to)?;
    Ok((moved.ra.degrees(), moved.dec.degrees()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::catalog::tests::{build_test_catalog, make_star};
    use crate::query::healpix::ang2pix_nest;
    use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};
    use std::collections::BTreeMap;

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

    fn j2016() -> TT {
        TT::from_julian_date(JulianDate::from_julian_year(2016.0))
    }

    fn one_year_on() -> TT {
        TT::from_julian_date(JulianDate::from_julian_year(2016.0).add_days(365.25))
    }

    // Expected places in these tests are eraPmsafe's. A star at rest comes back where it was,
    // apart from the last bit lost to the round trip through Cartesian coordinates.
    #[test]
    fn test_apply_proper_motion_zero_pm() {
        let place = apply_proper_motion(&star_at_100_45(0.0), &j2016(), &one_year_on()).unwrap();
        assert_eq!(place, (100.0, 44.99999999999999));
    }

    #[test]
    fn test_apply_proper_motion_one_year() {
        let place = apply_proper_motion(&star_at_100_45(3600.0), &j2016(), &one_year_on()).unwrap();
        assert_eq!(place, (100.0014142382452, 45.00099999127294));
    }

    #[test]
    fn test_apply_proper_motion_non_finite_is_err() {
        let star = star_at_100_45(f64::NAN);
        assert!(apply_proper_motion(&star, &j2016(), &one_year_on()).is_err());
    }

    // Order 1 NEST pixels 19 and 40 contain (10, 20) and (200, -70).
    fn two_star_catalog() -> tempfile::NamedTempFile {
        build_test_catalog(
            1,
            2016.0,
            21.0,
            &[
                (19, vec![make_star(1, 10.0, 20.0, 5.0, 0)]),
                (40, vec![make_star(2, 200.0, -70.0, 6.0, 0)]),
            ],
        )
    }

    fn params(ra_deg: f64, dec_deg: f64, radius_deg: f64) -> ConeSearchParams {
        ConeSearchParams {
            ra_deg,
            dec_deg,
            radius_deg,
            max_mag: None,
            max_results: None,
            epoch: None,
        }
    }

    #[test]
    fn test_cone_search_rejects_invalid_input() {
        let file = two_star_catalog();
        let catalog = Catalog::open(file.path()).unwrap();
        let nan_mag = ConeSearchParams {
            max_mag: Some(f64::NAN),
            ..params(10.0, 20.0, 1.0)
        };
        let cases = [
            ("NaN RA", params(f64::NAN, 20.0, 1.0)),
            ("infinite RA", params(f64::INFINITY, 20.0, 1.0)),
            ("NaN Dec", params(10.0, f64::NAN, 1.0)),
            ("Dec 95", params(10.0, 95.0, 1.0)),
            ("Dec -90.5", params(10.0, -90.5, 1.0)),
            ("NaN radius", params(10.0, 20.0, f64::NAN)),
            ("infinite radius", params(10.0, 20.0, f64::INFINITY)),
            ("negative radius", params(10.0, 20.0, -1.0)),
            ("NaN max_mag", nan_mag),
        ];
        for (name, p) in cases {
            assert!(cone_search(&catalog, &p).is_err(), "{name} was accepted");
        }
    }

    // Order 0: every record sits in pixel 0, and a 180-degree cone reads every pixel.
    fn catalog_with_corrupt_records() -> tempfile::NamedTempFile {
        let moving = |id, pmra, pmdec, parallax| StarRecord {
            pmra,
            pmdec,
            parallax,
            ..make_star(id, 10.0, 20.0, 5.0, 0)
        };
        let stars = vec![
            make_star(1, 10.0, 20.0, 5.0, 0),
            make_star(2, 200.0, -70.0, 6.0, 0),
            make_star(3, f64::NAN, 20.0, 5.0, 0),
            make_star(4, 10.0, f64::INFINITY, 5.0, 0),
            make_star(5, 10.0, 95.0, 5.0, 0),
            make_star(6, f64::NEG_INFINITY, 20.0, 5.0, 0),
            make_star(7, 10.0, f64::NAN, 5.0, 0),
            moving(8, f64::NAN, 0.0, 1.0),
            moving(9, 0.0, f64::INFINITY, 1.0),
            moving(10, 0.0, 0.0, f64::NAN),
            moving(11, 1e308, 1e308, 1.0),
        ];
        build_test_catalog(0, 2016.0, 21.0, &[(0, stars)])
    }

    fn found_ids(catalog: &Catalog, epoch: Option<JulianDate>) -> Vec<i64> {
        let p = ConeSearchParams {
            epoch,
            ..params(10.0, 20.0, 180.0)
        };
        let mut ids: Vec<i64> = cone_search(catalog, &p)
            .unwrap()
            .iter()
            .map(|r| r.star.source_id)
            .collect();
        ids.sort();
        ids
    }

    #[test]
    fn test_cone_search_skips_records_with_invalid_places() {
        let file = catalog_with_corrupt_records();
        let catalog = Catalog::open(file.path()).unwrap();
        assert_eq!(found_ids(&catalog, None), [1, 2, 8, 9, 10, 11]);
    }

    #[test]
    fn test_cone_search_skips_records_whose_propagated_place_is_invalid() {
        let file = catalog_with_corrupt_records();
        let catalog = Catalog::open(file.path()).unwrap();
        let epoch = JulianDate::from_julian_year(2016.0).add_days(365.25);
        assert_eq!(found_ids(&catalog, Some(epoch)), [1, 2]);
    }

    #[test]
    fn test_cone_search_rejects_non_finite_epoch() {
        let file = two_star_catalog();
        let catalog = Catalog::open(file.path()).unwrap();
        let epochs = [
            JulianDate::new(f64::NAN, 0.0),
            JulianDate::new(2_460_000.5, f64::INFINITY),
            JulianDate::new(f64::MAX, f64::MAX),
        ];
        // The first cone holds a star and the second holds none.
        for (ra, dec) in [(10.0, 20.0), (100.0, 0.0)] {
            for epoch in epochs {
                let p = ConeSearchParams {
                    epoch: Some(epoch),
                    ..params(ra, dec, 0.1)
                };
                let err = cone_search(&catalog, &p)
                    .unwrap_err()
                    .to_string()
                    .to_lowercase();
                assert!(err.contains("epoch"), "({ra}, {dec}) {epoch:?}: {err}");
            }
        }
    }

    #[test]
    fn test_cone_search_reduces_large_ra() {
        let file = two_star_catalog();
        let catalog = Catalog::open(file.path()).unwrap();
        // fmod(3600000000000010, 360) is 10 and fmod(1e300, 360) is 0, both exactly.
        for (ra, reduced, radius) in [(3600000000000010.0, 10.0, 0.01), (1e300, 0.0, 10.0)] {
            let found = cone_search(&catalog, &params(ra, 20.0, radius)).unwrap();
            let expected = cone_search(&catalog, &params(reduced, 20.0, radius)).unwrap();
            let summary = |r: &[ConeSearchResult]| -> Vec<(i64, f64)> {
                r.iter()
                    .map(|r| (r.star.source_id, r.distance_deg))
                    .collect()
            };
            assert_eq!(summary(&found), summary(&expected), "ra {ra}");
            assert_eq!(found.len(), 1, "ra {ra}");
        }
    }

    // Knuth's MMIX linear congruential generator, so the star fields are the same on every run.
    struct Lcg(u64);

    impl Lcg {
        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        fn sphere(&mut self) -> (f64, f64) {
            let ra = self.uniform() * 360.0;
            let z = 2.0 * self.uniform() - 1.0;
            (ra, libm::asin(z) * RAD_TO_DEG)
        }
    }

    // The point `distance_deg` from (ra, dec) at position angle `pa_deg`.
    fn offset(ra: f64, dec: f64, distance_deg: f64, pa_deg: f64) -> (f64, f64) {
        let (d, pa, dec0) = (
            distance_deg * DEG_TO_RAD,
            pa_deg * DEG_TO_RAD,
            dec * DEG_TO_RAD,
        );
        let sin_dec =
            libm::sin(dec0) * libm::cos(d) + libm::cos(dec0) * libm::sin(d) * libm::cos(pa);
        let sin_dec = sin_dec.clamp(-1.0, 1.0);
        let dra = libm::atan2(
            libm::sin(pa) * libm::sin(d) * libm::cos(dec0),
            libm::cos(d) - libm::sin(dec0) * sin_dec,
        );
        (
            reduce_ra_deg(ra + dra * RAD_TO_DEG),
            libm::asin(sin_dec) * RAD_TO_DEG,
        )
    }

    // The centre, rings just inside the cone edge at every degree of position angle, and
    // 20,000 stars spread over the sky.
    fn stars_for_cone(ra: f64, dec: f64, radius: f64) -> Vec<(f64, f64)> {
        let mut points = vec![(reduce_ra_deg(ra), dec)];
        for fraction in [0.5, 0.9, 0.99, 0.999999] {
            for pa in 0..360 {
                points.push(offset(ra, dec, radius * fraction, f64::from(pa)));
            }
        }
        let mut rng = Lcg(12345);
        points.extend((0..20_000).map(|_| rng.sphere()));
        points
    }

    fn catalog_of(order: u32, epoch: f64, stars: &[StarRecord]) -> tempfile::NamedTempFile {
        let mut by_pixel: BTreeMap<u64, Vec<StarRecord>> = BTreeMap::new();
        for star in stars {
            let pixel = ang2pix_nest(order, star.ra, star.dec).unwrap();
            by_pixel.entry(pixel).or_default().push(*star);
        }
        let pixel_stars: Vec<(u64, Vec<StarRecord>)> = by_pixel.into_iter().collect();
        build_test_catalog(order, epoch, 21.0, &pixel_stars)
    }

    // Compares cone_search with a scan of every star, using the distance test match_star uses.
    fn assert_cone_finds_every_star(order: u32, ra: f64, dec: f64, radius: f64) {
        let points = stars_for_cone(ra, dec, radius);
        let stars: Vec<StarRecord> = points
            .iter()
            .enumerate()
            .map(|(i, &(ra, dec))| make_star(i as i64, ra, dec, 10.0, 0))
            .collect();
        let file = catalog_of(order, 2016.0, &stars);
        let catalog = Catalog::open(file.path()).unwrap();
        let found = cone_search(&catalog, &params(ra, dec, radius)).unwrap();
        let mut found: Vec<i64> = found.iter().map(|r| r.star.source_id).collect();
        found.sort_unstable();
        let centre_ra = reduce_ra_deg(ra);
        let inside: Vec<i64> = stars
            .iter()
            .filter(|s| angular_separation_deg(centre_ra, dec, s.ra, s.dec) <= radius)
            .map(|s| s.source_id)
            .collect();
        if found != inside {
            let missed = inside
                .iter()
                .filter(|id| found.binary_search(id).is_err())
                .count();
            panic!(
                "order {order}, cone ({ra}, {dec}) radius {radius}: {missed} of {} stars missed",
                inside.len()
            );
        }
    }

    fn assert_cones_find_every_star(cones: &[(f64, f64, f64)]) {
        for order in [0, 3, 6] {
            for &(ra, dec, radius) in cones {
                assert_cone_finds_every_star(order, ra, dec, radius);
            }
        }
    }

    #[test]
    fn test_cone_search_finds_every_star_near_the_poles() {
        assert_cones_find_every_star(&[
            (0.0, -90.0, 2.0),
            (200.0, 89.5, 1.0),
            (297.83, 89.51, 1.0),
            (0.0, 90.0, 0.5),
            (10.0, -89.9, 0.3),
            (45.0, 88.0, 2.5),
        ]);
    }

    #[test]
    fn test_cone_search_finds_every_star_across_ra_zero() {
        assert_cones_find_every_star(&[
            (0.0, 0.0, 1.0),
            (359.95, 0.0, 0.1),
            (360.0, 10.0, 2.0),
            (0.01, -30.0, 5.0),
            (315.0, 0.0, 5.0),
        ]);
    }

    // Dec 41.810314895778596 is z = 2/3, where the polar caps meet the equatorial belt.
    #[test]
    fn test_cone_search_finds_every_star_at_cap_boundaries() {
        let cap = 41.810314895778596;
        assert_cones_find_every_star(&[
            (45.0, cap, 3.0),
            (0.0, -cap, 1.0),
            (315.0, cap, 0.2),
            (100.0, -cap, 10.0),
        ]);
    }

    #[test]
    fn test_cone_search_finds_every_star_at_radius_extremes() {
        assert_cones_find_every_star(&[
            (10.0, 10.0, 0.0),
            (0.0, 90.0, 0.0),
            (123.4, -45.0, 1e-7),
            (300.0, 89.99999, 1e-5),
            (80.0, 20.0, 90.0),
            (10.0, -60.0, 180.0),
        ]);
    }

    fn fast_star(ra: f64, dec: f64) -> StarRecord {
        StarRecord {
            pmra: 5000.0,
            pmdec: -3000.0,
            parallax: 100.0,
            ..make_star(7, ra, dec, 5.0, 0)
        }
    }

    #[test]
    fn test_cone_search_propagates_from_header_epoch() {
        let star = fast_star(150.0, 30.0);
        let file = catalog_of(1, 1991.25, &[star]);
        let catalog = Catalog::open(file.path()).unwrap();
        let to = JulianDate::from_julian_year(2024.5);
        let search = ConeSearchParams {
            epoch: Some(to),
            ..params(150.0, 30.0, 1.0)
        };
        let found = cone_search(&catalog, &search).unwrap();
        let expected = CatalogStar {
            ra: Angle::from_degrees(star.ra),
            dec: Angle::from_degrees(star.dec),
            pm_ra_cos_dec_mas_per_year: star.pmra,
            pm_dec_mas_per_year: star.pmdec,
            parallax_mas: star.parallax,
            radial_velocity_km_per_s: 0.0,
        }
        .propagate(
            &TT::from_julian_date(JulianDate::from_julian_year(1991.25)),
            &TT::from_julian_date(to),
        )
        .unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(
            (found[0].ra_deg, found[0].dec_deg),
            (expected.ra.degrees(), expected.dec.degrees())
        );
    }

    #[test]
    fn test_cone_search_radius_above_180_is_whole_sky() {
        let file = two_star_catalog();
        let catalog = Catalog::open(file.path()).unwrap();
        for radius in [180.0, 1000.0] {
            let found = cone_search(&catalog, &params(10.0, 20.0, radius)).unwrap();
            let ids: Vec<i64> = found.iter().map(|r| r.star.source_id).collect();
            assert_eq!(ids, vec![1, 2], "radius {radius}");
        }
    }
}
