//! Observer location on Earth.
//!
//! - [`Location`]: WGS84 geodetic coordinates (latitude, longitude, height)
//! - [`Location::to_geocentric_km`] and [`Location::to_geocentric_meters`]: geodetic-to-geocentric
//!   conversions for parallax corrections
//!
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
//! use celestial_core::location::Location;
//!
//! // Mauna Kea summit
//! let obs = Location::from_degrees(19.8207, -155.4681, 4205.0)?;
//!
//! // Access coordinates
//! assert_eq!(obs.latitude_degrees(), 19.8207);
//! # Ok::<(), celestial_core::errors::AstroError>(())
//! ```

mod geodesy;
mod types;
mod validation;

pub use types::Location;
