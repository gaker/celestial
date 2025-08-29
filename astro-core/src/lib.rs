//! Core astronomical types and utilities
//!
//! Provides fundamental types used across astronomical calculations, including
//! Earth locations and coordinate system utilities. This crate serves as the 
//! foundation for higher-level astronomical computations requiring precise
//! geographic positioning.
//!
//! # Key Types
//!
//! - [`Location`]: Earth surface positions with WGS84 ellipsoid support
//! - Physical and mathematical constants for astronomical calculations
//!
//! # Precision
//!
//! Geographic coordinates maintain full f64 precision. Geocentric conversions
//! are accurate to micrometers when validated against ERFA/SOFA standards.
//!
//! # Example
//!
//! ```rust
//! use astro_core::Location;
//! 
//! // Create a location for Mauna Kea Observatory
//! let mauna_kea = Location::from_degrees(19.8283, -155.4783, 4207.0)?;
//! 
//! // Get geocentric coordinates for Earth rotation calculations
//! let (u, v) = mauna_kea.to_geocentric_km()?;
//! println!("Distance from Earth's axis: {:.3} km", u);
//! println!("Distance north of equator: {:.3} km", v);
//! # Ok::<(), astro_core::AstroError>(())
//! ```

pub mod location;
pub mod constants;
pub mod errors;

pub use location::Location;
pub use errors::{AstroError, AstroResult, MathErrorKind};
