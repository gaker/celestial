//! Angular measurements for astronomical calculations.
//!
//! This module provides [`Angle`], the fundamental angular measurement type used throughout
//! the astronomy library. Angles are stored internally as radians (f64) but can be constructed
//! from and converted to degrees, hours, arcminutes, and arcseconds.
//!
//! # Design Rationale
//!
//! **Why radians internally?** All trigonometric functions in Rust (and most languages) operate
//! on radians. Storing radians avoids repeated conversions during calculations. The degree-based
//! constructors and accessors provide ergonomic APIs for human-readable values.
//!
//! **Why associated constants?** [`Angle::PI`], [`Angle::HALF_PI`], and [`Angle::ZERO`] exist
//! because angles are not just numbers. While `celestial_core::constants::PI` gives you a raw
//! float, `Angle::PI` gives you a typed angle. This prevents accidentally mixing raw radians
//! with Angles and catches unit errors at compile time.
//!
//! # Quick Start
//!
//! ```
//! use celestial_core::angle::Angle;
//! use celestial_core::constants::QUARTER_PI;
//!
//! // Construction - pick the unit that matches your data
//! let from_deg = Angle::from_degrees(45.0);
//! let from_rad = Angle::from_radians(0.785398);
//! let from_hrs = Angle::from_hours(3.0);  // 3h = 45 degrees
//! let from_arcsec = Angle::from_arcseconds(162000.0);  // 45 degrees
//!
//! // Conversion - get any unit you need
//! assert_eq!(from_deg.radians(), QUARTER_PI);
//! assert_eq!(from_deg.hours(), 3.0);
//!
//! // Trigonometry - no conversion needed
//! let (sin, cos) = from_deg.sin_cos();
//! ```
//!
//! # Hour Angles
//!
//! Astronomy uses hours (0-24h) for right ascension. One hour equals 15 degrees:
//!
//! ```
//! use celestial_core::angle::Angle;
//!
//! let ra = Angle::from_hours(6.0);  // 6h RA
//! assert_eq!(ra.degrees(), 90.0);
//! ```
//!
//! # Validation
//!
//! Angles can be validated for specific astronomical contexts:
//!
//! ```
//! use celestial_core::angle::Angle;
//!
//! let dec = Angle::from_degrees(45.0);
//! assert!(dec.validate_declination(false).is_ok());  // -90 to +90
//!
//! let bad_dec = Angle::from_degrees(100.0);
//! assert!(bad_dec.validate_declination(false).is_err());  // Out of range
//!
//! // Right ascension auto-normalizes to [0, 360). 2π has no exact f64 value, so the
//! // wrapped angle lands one ulp above 40°.
//! let ra = Angle::from_degrees(400.0);
//! let normalized = ra.validate_right_ascension().unwrap();
//! assert_eq!(normalized.degrees(), 40.00000000000001);
//! ```
//!
//! # Convenience Functions
//!
//! For terser code, use the free functions [`deg`], [`rad`], [`hours`], [`arcsec`], [`arcmin`]:
//!
//! ```
//! use celestial_core::angle::{deg, hours, arcsec};
//!
//! let a = deg(45.0);
//! let b = hours(3.0);
//! let c = arcsec(162000.0);
//!
//! assert_eq!(a.degrees(), b.degrees());
//! ```
//!
//! # Arithmetic
//!
//! Angles support addition, subtraction, negation, and scalar multiplication/division:
//!
//! ```
//! use celestial_core::angle::Angle;
//!
//! let a = Angle::from_degrees(30.0);
//! let b = Angle::from_degrees(15.0);
//!
//! let sum = a + b;  // 45 degrees
//! let diff = a - b;  // 15 degrees
//! let scaled = a * 2.0;  // 60 degrees
//! let neg = -a;  // -30 degrees
//! ```

mod format;
mod normalize;
mod ops;
mod parse;
#[cfg(feature = "serde")]
mod serde_;
mod shorthand;
mod types;
mod units;
mod validate;

pub use format::{parse_angle, DmsFmt, HmsFmt, ParsedAngle};
pub use normalize::{wrap_0_2pi, wrap_pm_pi};
pub use parse::{parse_dms, parse_hms, AngleUnits, ParseAngle};
pub use shorthand::{arcmin, arcsec, deg, hours, rad};
pub use types::Angle;
pub use validate::{validate_declination, validate_latitude, validate_right_ascension};
