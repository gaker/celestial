//! Observer location on Earth using WGS84 geodetic coordinates.
//!
//! This module provides the [`Location`] type for representing geographic positions.
//! Coordinates are geodetic (latitude/longitude relative to the WGS84 ellipsoid),
//! not geocentric (relative to Earth's center of mass).
//!
//! The distinction matters for precision astronomy: geodetic latitude differs from
//! geocentric latitude by up to ~11 arcminutes at mid-latitudes due to Earth's
//! equatorial bulge.
//!
//! # Coordinate conventions
//!
//! - **Latitude**: North positive, stored in radians, range [-pi/2, pi/2]
//! - **Longitude**: East positive, stored in radians, range [-pi, pi]
//! - **Height**: Meters above the WGS84 ellipsoid (not sea level)
//!
//! # Example
//!
//! ```
//! use astro_core::Location;
//!
//! // Mauna Kea summit
//! let obs = Location::from_degrees(19.8207, -155.4681, 4205.0)?;
//!
//! // Access coordinates
//! assert!((obs.latitude_degrees() - 19.8207).abs() < 1e-10);
//! # Ok::<(), astro_core::AstroError>(())
//! ```

use crate::errors::{AstroError, AstroResult, MathErrorKind};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A geographic location on Earth in WGS84 geodetic coordinates.
///
/// All angular values are stored internally in radians. Use [`Location::from_degrees`]
/// for convenience when working with degree-based coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Location {
    /// Geodetic latitude in radians. North is positive.
    pub latitude: f64,
    /// Geodetic longitude in radians. East is positive.
    pub longitude: f64,
    /// Height above WGS84 ellipsoid in meters.
    pub height: f64,
}

impl Location {
    /// Creates a new location from coordinates in radians.
    ///
    /// # Arguments
    ///
    /// * `latitude` - Geodetic latitude in radians, must be in [-pi/2, pi/2]
    /// * `longitude` - Geodetic longitude in radians, must be in [-pi, pi]
    /// * `height` - Height above WGS84 ellipsoid in meters, must be in [-12000, 100000]
    ///
    /// # Errors
    ///
    /// Returns an error if any coordinate is non-finite or outside its valid range.
    /// The height range covers the Mariana Trench floor to well above aircraft altitude.
    pub fn new(latitude: f64, longitude: f64, height: f64) -> AstroResult<Self> {
        if !latitude.is_finite() {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Latitude must be finite",
            ));
        }
        if !longitude.is_finite() {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Longitude must be finite",
            ));
        }
        if !height.is_finite() {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Height must be finite",
            ));
        }

        if latitude.abs() > crate::constants::HALF_PI {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Latitude outside valid range [-π/2, π/2]",
            ));
        }
        if longitude.abs() > crate::constants::PI {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Longitude outside valid range [-π, π]",
            ));
        }
        if !(-12000.0..=100000.0).contains(&height) {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Height outside reasonable range [-12000, 100000] meters",
            ));
        }

        Ok(Self {
            latitude,
            longitude,
            height,
        })
    }

    /// Creates a new location from coordinates in degrees.
    ///
    /// This is the typical way to create a Location, since most sources
    /// provide coordinates in degrees.
    ///
    /// # Arguments
    ///
    /// * `lat_deg` - Geodetic latitude in degrees, must be in [-90, 90]
    /// * `lon_deg` - Geodetic longitude in degrees, must be in [-180, 180]
    /// * `height_m` - Height above WGS84 ellipsoid in meters
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Location;
    ///
    /// // La Silla Observatory, Chile
    /// let la_silla = Location::from_degrees(-29.2563, -70.7380, 2400.0)?;
    /// # Ok::<(), astro_core::AstroError>(())
    /// ```
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, height_m: f64) -> AstroResult<Self> {
        if !lat_deg.is_finite() {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Latitude degrees must be finite",
            ));
        }
        if !lon_deg.is_finite() {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Longitude degrees must be finite",
            ));
        }
        if lat_deg.abs() > 90.0 {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Latitude outside valid range [-90, 90] degrees",
            ));
        }
        if lon_deg.abs() > 180.0 {
            return Err(AstroError::math_error(
                "location_validation",
                MathErrorKind::InvalidInput,
                "Longitude outside valid range [-180, 180] degrees",
            ));
        }

        Self::new(lat_deg.to_radians(), lon_deg.to_radians(), height_m)
    }

    /// Returns the latitude in degrees.
    pub fn latitude_degrees(&self) -> f64 {
        self.latitude.to_degrees()
    }

    /// Returns the longitude in degrees.
    pub fn longitude_degrees(&self) -> f64 {
        self.longitude.to_degrees()
    }

    /// Returns the latitude as an [`Angle`](crate::Angle).
    pub fn latitude_angle(&self) -> crate::Angle {
        crate::Angle::from_radians(self.latitude)
    }

    /// Returns the longitude as an [`Angle`](crate::Angle).
    pub fn longitude_angle(&self) -> crate::Angle {
        crate::Angle::from_radians(self.longitude)
    }

    /// Returns the Royal Observatory, Greenwich (0, 0, 0).
    ///
    /// Useful as a default or reference location.
    pub fn greenwich() -> Self {
        Self::from_degrees(0.0, 0.0, 0.0).expect("Greenwich coordinates should always be valid")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_location_creation() {
        let loc = Location::new(0.5, 1.0, 100.0).unwrap();
        assert_eq!(loc.latitude, 0.5);
        assert_eq!(loc.longitude, 1.0);
        assert_eq!(loc.height, 100.0);
    }

    #[test]
    fn test_from_degrees() {
        let loc = Location::from_degrees(45.0, 90.0, 1000.0).unwrap();
        assert!((loc.latitude - 45.0_f64.to_radians()).abs() < 1e-15);
        assert!((loc.longitude - 90.0_f64.to_radians()).abs() < 1e-15);
        assert_eq!(loc.height, 1000.0);
    }

    #[test]
    fn test_longitude_degrees_conversion_returns_degrees() {
        let loc = Location::from_degrees(0.0, 180.0, 0.0).unwrap();
        assert_eq!(loc.longitude_degrees(), 180.0);
    }

    #[test]
    fn test_longitude_degrees_conversion_handles_negative() {
        let loc = Location::from_degrees(0.0, -90.0, 0.0).unwrap();
        assert_eq!(loc.longitude_degrees(), -90.0);
    }

    #[test]
    fn test_longitude_angle_returns_angle_object() {
        let loc = Location::from_degrees(0.0, 45.0, 0.0).unwrap();
        let angle = loc.longitude_angle();
        crate::test_helpers::assert_float_eq(angle.degrees(), 45.0, 1);
    }

    #[test]
    fn test_longitude_angle_handles_wraparound() {
        let loc = Location::from_degrees(0.0, -180.0, 0.0).unwrap();
        let angle = loc.longitude_angle();
        crate::test_helpers::assert_float_eq(angle.degrees(), -180.0, 1);
    }

    #[test]
    fn test_location_validation_errors() {
        let result = Location::new(f64::NAN, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Latitude must be finite"));

        let result = Location::new(0.0, f64::NAN, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Longitude must be finite"));

        let result = Location::new(0.0, 0.0, f64::NAN);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Height must be finite"));

        let result = Location::new(f64::INFINITY, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Latitude must be finite"));

        let result = Location::new(0.0, f64::INFINITY, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Longitude must be finite"));

        let result = Location::new(0.0, 0.0, f64::INFINITY);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Height must be finite"));

        let result = Location::new(crate::constants::PI, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range"));

        let result = Location::new(-crate::constants::PI, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range"));

        let result = Location::new(0.0, crate::constants::TWOPI, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range"));

        let result = Location::new(0.0, -crate::constants::TWOPI, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range"));

        let result = Location::new(0.0, 0.0, 200000.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside reasonable range"));

        let result = Location::new(0.0, 0.0, -20000.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside reasonable range"));
    }

    #[test]
    fn test_from_degrees_validation_errors() {
        let result = Location::from_degrees(f64::NAN, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Latitude degrees must be finite"));

        let result = Location::from_degrees(0.0, f64::NAN, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Longitude degrees must be finite"));

        let result = Location::from_degrees(95.0, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range [-90, 90]"));

        let result = Location::from_degrees(-95.0, 0.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range [-90, 90]"));

        let result = Location::from_degrees(0.0, 185.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range [-180, 180]"));

        let result = Location::from_degrees(0.0, -185.0, 0.0);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("outside valid range [-180, 180]"));
    }
}
