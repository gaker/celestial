//! Geodetic to geocentric coordinate conversions for Earth-based observers.
//!
//! # Geodetic vs Geocentric Coordinates
//!
//! **Geodetic coordinates** (what GPS gives you) define position relative to the WGS84
//! reference ellipsoid:
//! - Latitude: angle between the equatorial plane and the ellipsoid surface normal
//! - Longitude: angle from the prime meridian
//! - Height: distance above the ellipsoid surface
//!
//! **Geocentric coordinates** define position relative to Earth's center of mass:
//! - The Earth is modeled as an oblate spheroid (equatorial bulge)
//! - At mid-latitudes, geodetic and geocentric latitude differ by up to ~11 arcminutes
//!
//! Topocentric corrections (parallax, aberration, refraction) require knowing the
//! observer's true position in space, not their position on the reference ellipsoid.
//! The geocentric coordinates returned here are the cylindrical components needed for:
//!
//! - **Diurnal parallax**: Moon position shifts by up to 1° depending on observer
//! - **Stellar parallax**: precise baseline for nearby star distances
//! - **Satellite tracking**: ground station positions in Earth-centered frame
//!
//! # WGS84 Ellipsoid Parameters
//!
//! This module uses the WGS84 reference ellipsoid:
//! - Semi-major axis (equatorial radius): 6,378,137.0 m (or 6378.137 km)
//! - Flattening: 1/298.257223563
//! - First eccentricity squared: ~0.00669438
//!
//! # Output Format
//!
//! Both conversion methods return `(u, v)` where:
//! - `u`: distance from Earth's rotation axis (equatorial component)
//! - `v`: distance from equatorial plane (polar component)
//!
//! These are cylindrical coordinates centered on Earth's center of mass.
//! To get Cartesian XYZ, you'd combine with longitude: `x = u*cos(lon)`, `y = u*sin(lon)`, `z = v`.

use crate::constants::{
    WGS84_ECCENTRICITY_SQUARED, WGS84_SEMI_MAJOR_AXIS, WGS84_SEMI_MAJOR_AXIS_KM,
};
use crate::errors::{AstroError, AstroResult, MathErrorKind};

use super::Location;

impl Location {
    /// Converts geodetic coordinates to geocentric cylindrical coordinates in kilometers.
    ///
    /// Uses the WGS84 ellipsoid to compute the observer's position relative to
    /// Earth's center of mass. The result accounts for Earth's equatorial bulge.
    ///
    /// # Returns
    ///
    /// `(u, v)` in kilometers where:
    /// - `u`: perpendicular distance from Earth's rotation axis
    /// - `v`: distance from the equatorial plane (positive north)
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::Location;
    ///
    /// let obs = Location::from_degrees(45.0, 0.0, 0.0)?;
    /// let (u, v) = obs.to_geocentric_km()?;
    ///
    /// // At 45 degrees, u and v are similar but u > v due to Earth's shape
    /// assert!(u > 4500.0 && u < 4600.0);
    /// assert!(v > 4400.0 && v < 4500.0);
    /// # Ok::<(), celestial_core::AstroError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error for degenerate latitude values that would cause division by zero.
    /// In practice, this is unlikely with validated `Location` values.
    pub fn to_geocentric_km(&self) -> AstroResult<(f64, f64)> {
        let lat = self.latitude;
        let height_km = self.height / 1000.0;

        let (sin_lat, cos_lat) = lat.sin_cos();

        let denominator = 1.0 - WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat;
        if denominator <= f64::EPSILON {
            return Err(AstroError::math_error(
                "geocentric_conversion",
                MathErrorKind::DivisionByZero,
                "Latitude too close to critical value causing division by zero",
            ));
        }

        let n = WGS84_SEMI_MAJOR_AXIS_KM / denominator.sqrt();

        let u = (n + height_km) * cos_lat;

        let v = (n * (1.0 - WGS84_ECCENTRICITY_SQUARED) + height_km) * sin_lat;

        Ok((u, v))
    }

    /// Converts geodetic coordinates to geocentric cylindrical coordinates in meters.
    ///
    /// Same algorithm as [`to_geocentric_km`](Self::to_geocentric_km) but returns
    /// results in meters for applications requiring higher precision or SI units.
    ///
    /// This method uses a slightly different formulation internally (computing the
    /// flattening ratio explicitly) but produces equivalent results to the km version.
    ///
    /// # Returns
    ///
    /// `(u, v)` in meters where:
    /// - `u`: perpendicular distance from Earth's rotation axis
    /// - `v`: distance from the equatorial plane (positive north)
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::Location;
    ///
    /// // Equator at sea level
    /// let equator = Location::from_degrees(0.0, 0.0, 0.0)?;
    /// let (u, v) = equator.to_geocentric_meters()?;
    ///
    /// // At equator: u equals semi-major axis, v is zero
    /// assert!((u - 6_378_137.0).abs() < 1.0);
    /// assert!(v.abs() < 1e-9);
    /// # Ok::<(), celestial_core::AstroError>(())
    /// ```
    pub fn to_geocentric_meters(&self) -> AstroResult<(f64, f64)> {
        let height_m = self.height;

        let (phi_sin, phi_cos) = self.latitude.sin_cos();

        let wgs_flattened = 1.0 / 298.257223563;
        let axis_ratio = 1.0 - wgs_flattened;
        let axis_ratio_sq = axis_ratio * axis_ratio;

        let norm_sq = phi_cos * phi_cos + axis_ratio_sq * phi_sin * phi_sin;
        if norm_sq <= f64::EPSILON {
            return Err(AstroError::math_error(
                "geocentric_conversion",
                MathErrorKind::DivisionByZero,
                "Latitude too close to critical value causing division by zero",
            ));
        }

        let prime_vertical_radius = WGS84_SEMI_MAJOR_AXIS / norm_sq.sqrt();
        let as_val = axis_ratio_sq * prime_vertical_radius;

        let equatorial_radius = (prime_vertical_radius + height_m) * phi_cos;
        let z_coordinate = (as_val + height_m) * phi_sin;

        Ok((equatorial_radius, z_coordinate))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geocentric_at_equator() {
        let loc = Location::from_degrees(0.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km().unwrap();

        assert_eq!(
            u,
            crate::constants::WGS84_SEMI_MAJOR_AXIS_KM,
            "u = {} km, expected ~6378.137 km",
            u
        );
        assert_eq!(v, 0.0, "v = {} km, expected ~0 km", v);
    }

    #[test]
    fn test_geocentric_at_north_pole() {
        let loc = Location::from_degrees(90.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km().unwrap();

        assert!(u.abs() < 1e-10, "u = {} km, expected very close to 0 km", u);
        let expected_polar_radius =
            crate::constants::WGS84_SEMI_MAJOR_AXIS_KM * (1.0 - 1.0 / 298.257223563);
        assert_eq!(
            v, expected_polar_radius,
            "v = {} km, expected ~{} km",
            v, expected_polar_radius
        );
    }

    #[test]
    fn test_geocentric_km_rejects_degenerate_latitude() {
        use crate::constants::{HALF_PI, PI};
        let critical_lat = HALF_PI - 1e-10;
        let critical_lat_deg = critical_lat * 180.0 / PI;

        let mut test_lat = critical_lat_deg;
        while test_lat < 90.0 {
            let loc = Location::from_degrees(test_lat, 0.0, 0.0).unwrap();
            if loc.to_geocentric_km().is_err() {
                return;
            }
            test_lat += 1e-12;
        }
    }

    #[test]
    fn test_geocentric_meters_rejects_degenerate_latitude() {
        use crate::constants::{HALF_PI, PI};
        let critical_lat = HALF_PI - 1e-10;
        let critical_lat_deg = critical_lat * 180.0 / PI;

        let mut test_lat = critical_lat_deg;
        while test_lat < 90.0 {
            let loc = Location::from_degrees(test_lat, 0.0, 0.0).unwrap();
            if loc.to_geocentric_meters().is_err() {
                return;
            }
            test_lat += 1e-12;
        }
    }

    #[test]
    fn test_geocentric_meters_handles_equator() {
        let loc = Location::from_degrees(0.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_meters().unwrap();
        crate::test_helpers::assert_float_eq(u, WGS84_SEMI_MAJOR_AXIS, 1);
        crate::test_helpers::assert_float_eq(v, 0.0, 1);
    }

    #[test]
    fn test_geocentric_meters_handles_north_pole() {
        let loc = Location::from_degrees(90.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_meters().unwrap();
        assert!(u.abs() < 1e-9);
        crate::test_helpers::assert_ulp_le(v, 6356752.314245179, 1, "v at pole");
    }

    #[test]
    fn test_geocentric_at_45_degrees() {
        let loc = Location::from_degrees(45.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km().unwrap();

        assert!(u > 4000.0 && u < 5000.0, "u = {} km, expected ~4500 km", u);
        assert!(v > 4000.0 && v < 5000.0, "v = {} km, expected ~4500 km", v);

        assert!(
            u > v,
            "u should be larger than v due to Earth's oblate shape: u={}, v={}",
            u,
            v
        );
        assert!(
            (u - v).abs() < 100.0,
            "At 45°, u and v should be similar: u={}, v={}",
            u,
            v
        );
    }

    #[test]
    fn test_geocentric_with_height() {
        let loc_sea_level = Location::from_degrees(0.0, 0.0, 0.0).unwrap();
        let loc_elevated = Location::from_degrees(0.0, 0.0, 1000.0).unwrap();

        let (u1, v1) = loc_sea_level.to_geocentric_km().unwrap();
        let (u2, v2) = loc_elevated.to_geocentric_km().unwrap();

        assert!(
            (u2 - u1 - 1.0).abs() < 0.001,
            "1km elevation should increase u by ~1km"
        );
        assert!(
            (v2 - v1).abs() < 0.001,
            "At equator, elevation shouldn't affect v much"
        );
    }

    #[test]
    fn test_negative_latitude() {
        let loc = Location::from_degrees(-45.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km().unwrap();

        assert!(u > 0.0, "u should be positive: {}", u);
        assert!(
            v < 0.0,
            "v should be negative in southern hemisphere: {}",
            v
        );
    }

    #[test]
    fn test_geocentric_division_by_zero() {
        let north_pole = Location {
            latitude: crate::constants::HALF_PI,
            longitude: 0.0,
            height: 0.0,
        };

        let result = north_pole.to_geocentric_km();
        assert!(result.is_ok());
    }
}
