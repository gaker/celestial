use super::validation::validate;
use crate::angle::Angle;
use crate::errors::AstroResult;

#[cfg(feature = "serde")]
use crate::errors::AstroError;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A geographic location on Earth in WGS84 geodetic coordinates.
///
/// All angular values are stored internally in radians. Use [`Location::from_degrees`]
/// for convenience when working with degree-based coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(try_from = "UncheckedLocation"))]
pub struct Location {
    /// Geodetic latitude in radians. North is positive.
    latitude: f64,
    /// Geodetic longitude in radians. East is positive.
    longitude: f64,
    /// Height above WGS84 ellipsoid in meters.
    height: f64,
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
        validate(latitude, longitude, height)?;
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
    /// * `lon_deg` - Geodetic longitude in degrees, must be in [-180, 360]. Values above 180
    ///   are East longitudes in 0-360 form and are stored as their [-180, 180] equivalent.
    /// * `height_m` - Height above WGS84 ellipsoid in meters
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::location::Location;
    ///
    /// // La Silla Observatory, Chile
    /// let la_silla = Location::from_degrees(-29.2563, -70.7380, 2400.0)?;
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, height_m: f64) -> AstroResult<Self> {
        let latitude = Angle::from_degrees(lat_deg).radians();
        let longitude = Angle::from_degrees(signed_longitude_degrees(lon_deg)).radians();
        Self::new(latitude, longitude, height_m)
    }

    pub fn latitude(&self) -> f64 {
        self.latitude
    }

    pub fn longitude(&self) -> f64 {
        self.longitude
    }

    pub fn height(&self) -> f64 {
        self.height
    }

    /// Returns the latitude in degrees.
    pub fn latitude_degrees(&self) -> f64 {
        self.latitude_angle().degrees()
    }

    /// Returns the longitude in degrees.
    pub fn longitude_degrees(&self) -> f64 {
        self.longitude_angle().degrees()
    }

    /// Returns the latitude as an [`Angle`].
    pub fn latitude_angle(&self) -> Angle {
        Angle::from_radians(self.latitude)
    }

    /// Returns the longitude as an [`Angle`].
    pub fn longitude_angle(&self) -> Angle {
        Angle::from_radians(self.longitude)
    }

    /// Returns the Royal Observatory, Greenwich (0, 0, 0).
    ///
    /// Useful as a default or reference location.
    pub fn greenwich() -> Self {
        Self {
            latitude: 0.0,
            longitude: 0.0,
            height: 0.0,
        }
    }
}

// East longitudes in (180°, 360°] map to (-180°, 0°]. Subtracting 360 is exact on that
// interval (Sterbenz), so 270°E converts to the same radians as -90°.
fn signed_longitude_degrees(lon_deg: f64) -> f64 {
    if lon_deg > 180.0 && lon_deg <= 360.0 {
        return lon_deg - 360.0;
    }
    lon_deg
}

#[cfg(feature = "serde")]
#[derive(Deserialize)]
struct UncheckedLocation {
    latitude: f64,
    longitude: f64,
    height: f64,
}

#[cfg(feature = "serde")]
impl TryFrom<UncheckedLocation> for Location {
    type Error = AstroError;

    fn try_from(raw: UncheckedLocation) -> AstroResult<Self> {
        Self::new(raw.latitude, raw.longitude, raw.height)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_location_creation() {
        let loc = Location::new(0.5, 1.0, 100.0).unwrap();
        assert_eq!(loc.latitude(), 0.5);
        assert_eq!(loc.longitude(), 1.0);
        assert_eq!(loc.height(), 100.0);
    }

    #[test]
    fn test_from_degrees() {
        let loc = Location::from_degrees(45.0, 90.0, 1000.0).unwrap();
        assert_eq!(loc.latitude, crate::constants::QUARTER_PI);
        assert_eq!(loc.longitude, crate::constants::HALF_PI);
        assert_eq!(loc.height, 1000.0);
    }

    // Inputs where multiplying by the rounded DEG_TO_RAD or RAD_TO_DEG lands 1 ULP off the
    // correctly rounded conversion.
    #[test]
    fn test_from_degrees_is_correctly_rounded() {
        let loc = Location::from_degrees(0.0017, -100.0002, 0.0).unwrap();
        assert_eq!(loc.latitude(), 2.9670597283903603e-5);
        assert_eq!(loc.longitude(), -1.7453327426528338);
    }

    #[test]
    fn test_degree_getters_are_correctly_rounded() {
        let loc = Location::new(0.003, 0.006, 0.0).unwrap();
        assert_eq!(loc.latitude_degrees(), 0.17188733853924695);
        assert_eq!(loc.longitude_degrees(), 0.3437746770784939);
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
        assert_eq!(angle.degrees(), 45.0);
    }

    #[test]
    fn test_longitude_angle_handles_wraparound() {
        let loc = Location::from_degrees(0.0, -180.0, 0.0).unwrap();
        let angle = loc.longitude_angle();
        assert_eq!(angle.degrees(), -180.0);
    }

    // Longitudes above 180° are East longitudes (the 0-360° form used by observatory code
    // lists); they wrap in degrees, which is exact, so 270°E is bit-identical to -90°.
    #[test]
    fn test_from_degrees_accepts_east_longitudes() {
        for (east, signed) in [(270.0, -90.0), (359.0, -1.0), (180.5, -179.5), (360.0, 0.0)] {
            let wrapped = Location::from_degrees(10.0, east, 0.0).unwrap();
            assert_eq!(wrapped, Location::from_degrees(10.0, signed, 0.0).unwrap());
        }
        let antimeridian = Location::from_degrees(0.0, 180.0, 0.0).unwrap();
        assert_eq!(antimeridian.longitude(), crate::constants::PI);
    }

    #[test]
    fn test_greenwich_is_origin() {
        assert_eq!(Location::greenwich(), Location::new(0.0, 0.0, 0.0).unwrap());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_deserialize_validates() {
        let valid: Location =
            serde_json::from_str(r#"{"latitude": 0.5, "longitude": -1.0, "height": 10.0}"#)
                .unwrap();
        assert_eq!(valid, Location::new(0.5, -1.0, 10.0).unwrap());
        let invalid = r#"{"latitude": 100.0, "longitude": 0.0, "height": 0.0}"#;
        assert!(serde_json::from_str::<Location>(invalid).is_err());
    }
}
