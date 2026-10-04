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

use crate::constants::{WGS84_FLATTENING, WGS84_SEMI_MAJOR_AXIS};

use super::Location;

impl Location {
    /// Converts geodetic coordinates to geocentric cylindrical coordinates in kilometers.
    ///
    /// Uses the WGS84 ellipsoid to compute the observer's position relative to
    /// Earth's center of mass. The result accounts for Earth's equatorial bulge.
    /// This is [`to_geocentric_meters`](Self::to_geocentric_meters) divided by 1000.
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
    /// use celestial_core::location::Location;
    ///
    /// let obs = Location::from_degrees(45.0, 0.0, 0.0)?;
    /// let (u, v) = obs.to_geocentric_km();
    ///
    /// // At 45 degrees, u and v are similar but u > v due to Earth's shape
    /// assert!(u > 4500.0 && u < 4600.0);
    /// assert!(v > 4400.0 && v < 4500.0);
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    pub fn to_geocentric_km(&self) -> (f64, f64) {
        let (u, v) = self.to_geocentric_meters();
        (u / 1000.0, v / 1000.0)
    }

    /// Converts geodetic coordinates to geocentric cylindrical coordinates in meters.
    ///
    /// Matches ERFA `gd2gce` on the WGS84 ellipsoid bit for bit.
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
    /// use celestial_core::location::Location;
    ///
    /// // Equator at sea level
    /// let equator = Location::from_degrees(0.0, 0.0, 0.0)?;
    /// let (u, v) = equator.to_geocentric_meters();
    ///
    /// // At equator: u equals semi-major axis, v is zero
    /// assert_eq!(u, 6_378_137.0);
    /// assert_eq!(v, 0.0);
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    pub fn to_geocentric_meters(&self) -> (f64, f64) {
        let height_m = self.height();

        let (phi_sin, phi_cos) = libm::sincos(self.latitude());

        let axis_ratio = 1.0 - WGS84_FLATTENING;
        let axis_ratio_sq = axis_ratio * axis_ratio;

        let norm_sq = phi_cos * phi_cos + axis_ratio_sq * phi_sin * phi_sin;
        let prime_vertical_radius = WGS84_SEMI_MAJOR_AXIS / libm::sqrt(norm_sq);
        let as_val = axis_ratio_sq * prime_vertical_radius;

        let equatorial_radius = (prime_vertical_radius + height_m) * phi_cos;
        let z_coordinate = (as_val + height_m) * phi_sin;

        (equatorial_radius, z_coordinate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{HALF_PI, WGS84_SEMI_MAJOR_AXIS_KM};

    #[test]
    fn test_geocentric_at_equator() {
        let loc = Location::from_degrees(0.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km();

        assert_eq!(
            u, WGS84_SEMI_MAJOR_AXIS_KM,
            "u = {} km, expected ~6378.137 km",
            u
        );
        assert_eq!(v, 0.0, "v = {} km, expected ~0 km", v);
    }

    #[test]
    fn test_geocentric_at_north_pole() {
        let loc = Location::from_degrees(90.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km();

        // f64 π/2 sits just below the true value, so its cosine and u are tiny but not zero.
        // Both are the ERFA gd2gce pole values below, divided by 1000.
        assert_eq!(u, 3.9186209248144715e-13);
        assert_eq!(v, 6356.752314245179);
    }

    // (latitude, u, v) from ERFA gd2gce on WGS84 at zero longitude and height, run with
    // Rust libm for sin/cos.
    const ERFA_GD2GCE: [(f64, f64, f64); 3] = [
        (HALF_PI, 3.9186209248144716e-10, 6356752.314245179),
        (-HALF_PI, 3.9186209248144716e-10, -6356752.314245179),
        (0.7, 4885058.998983235, 4087083.5464733117),
    ];

    #[test]
    fn test_geocentric_meters_matches_erfa_gd2gce() {
        for (latitude, u, v) in ERFA_GD2GCE {
            let loc = Location::new(latitude, 0.0, 0.0).unwrap();
            assert_eq!(loc.to_geocentric_meters(), (u, v), "{latitude}");
        }
    }

    // (latitude, height m, u km, v km): ERFA gd2gce on WGS84 at zero longitude, run with
    // Rust libm for sin/cos, divided by 1000.
    const ERFA_GD2GCE_KM: [(f64, f64, f64, f64); 2] = [
        (
            -0.007087162286871784,
            4407.929206062694,
            6382.385711322742,
            -44.93115779447372,
        ),
        (
            -0.9192720194817078,
            1622.7674418241413,
            3876.8926040765296,
            -5049.676267536883,
        ),
    ];

    #[test]
    fn test_geocentric_km_matches_erfa_gd2gce() {
        for (latitude, height, u, v) in ERFA_GD2GCE_KM {
            let loc = Location::new(latitude, 0.0, height).unwrap();
            assert_eq!(loc.to_geocentric_km(), (u, v), "{latitude}");
        }
    }

    #[test]
    fn test_geocentric_meters_handles_equator() {
        let loc = Location::from_degrees(0.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_meters();
        assert_eq!(u, WGS84_SEMI_MAJOR_AXIS);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_geocentric_meters_handles_north_pole() {
        let loc = Location::from_degrees(90.0, 0.0, 0.0).unwrap();
        let (_, u, v) = ERFA_GD2GCE[0];
        assert_eq!(loc.to_geocentric_meters(), (u, v));
    }

    #[test]
    fn test_geocentric_at_45_degrees() {
        let loc = Location::from_degrees(45.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km();

        assert!(u > 4000.0 && u < 5000.0, "u = {} km, expected ~4500 km", u);
        assert!(v > 4000.0 && v < 5000.0, "v = {} km, expected ~4500 km", v);

        assert!(
            u > v,
            "u should be larger than v due to Earth's oblate shape: u={}, v={}",
            u,
            v
        );
        assert!(
            libm::fabs(u - v) < 100.0,
            "At 45°, u and v should be similar: u={}, v={}",
            u,
            v
        );
    }

    #[test]
    fn test_geocentric_with_height() {
        let loc_elevated = Location::from_degrees(0.0, 0.0, 1000.0).unwrap();
        assert_eq!(
            loc_elevated.to_geocentric_km(),
            (WGS84_SEMI_MAJOR_AXIS_KM + 1.0, 0.0)
        );
    }

    #[test]
    fn test_negative_latitude() {
        let loc = Location::from_degrees(-45.0, 0.0, 0.0).unwrap();
        let (u, v) = loc.to_geocentric_km();

        assert!(u > 0.0, "u should be positive: {}", u);
        assert!(
            v < 0.0,
            "v should be negative in southern hemisphere: {}",
            v
        );
    }
}
