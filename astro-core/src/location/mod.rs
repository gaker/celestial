//! Observer location on Earth.
//!
//! - [`Location`]: WGS84 geodetic coordinates (latitude, longitude, height)
//! - [`geodesy`]: geodetic-to-geocentric conversions for parallax corrections

pub mod core;
pub mod geodesy;

pub use core::Location;
