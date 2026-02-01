//! Low-level astronomical calculations for coordinate transformations.
//!
//! `astro-core` provides the mathematical building blocks for celestial mechanics:
//! rotation matrices, nutation/precession models, angle handling, and geodetic conversions.
//! It implements IAU 2000/2006 standards in pure Rust with no runtime FFI.
//!
//! # Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`angle`] | Angle types, parsing (HMS/DMS), normalization, validation |
//! | [`matrix`] | 3×3 rotation matrices and 3D vectors |
//! | [`nutation`] | IAU 2000A/2000B/2006A nutation models |
//! | [`precession`] | IAU 2000/2006 precession (Fukushima-Williams angles) |
//! | [`cio`] | CIO-based GCRS↔CIRS transformations |
//! | [`obliquity`] | Mean obliquity of the ecliptic (IAU 1980, 2006) |
//! | [`location`] | Observer geodetic coordinates, geocentric conversion |
//! | [`constants`] | Astronomical constants (J2000, WGS84, unit conversions) |
//! | [`errors`] | [`AstroError`] and [`AstroResult`] |
//!
//! # Coordinate Transformation Pipeline
//!
//! GCRS → CIRS transformation (CIO-based):
//!
//! ```ignore
//! // 1. Compute precession-nutation-bias matrix
//! let fw = FukushimaWilliamsAngles::iau2006a(tt_centuries);
//! let nutation = NutationIAU2006A::new().compute(jd1, jd2)?;
//! let npb = fw.build_npb_matrix(nutation.delta_psi, nutation.delta_eps);
//!
//! // 2. Extract CIO quantities
//! let cio = CioSolution::calculate(&npb, tt_centuries)?;
//!
//! // 3. Build GCRS→CIRS matrix
//! let matrix = gcrs_to_cirs_matrix(cio.cip.x, cio.cip.y, cio.s);
//! ```
//!
//! # Re-exports
//!
//! Common types are re-exported at the crate root for convenience:
//!
//! ```
//! use astro_core::{Angle, Vector3, RotationMatrix3, Location};
//! use astro_core::{AstroError, AstroResult, MathErrorKind};
//! ```
//!
//! # Design Notes
//!
//! - **Two-part Julian Dates**: Functions accepting `(jd1, jd2)` preserve precision by
//!   splitting the date. Typically `jd1 = 2451545.0` (J2000.0) and `jd2` is days from epoch.
//!
//! - **Radians internally**: All angular computations use radians. The [`Angle`] type
//!   provides conversion methods for degrees/HMS/DMS display.
//!
//! - **No implicit state**: Models like [`NutationIAU2006A`](nutation::NutationIAU2006A)
//!   are stateless calculators. Call `compute(jd1, jd2)` with any epoch.

pub mod angle;
pub mod cio;
pub mod constants;
pub mod errors;
pub mod location;
pub mod math;
pub mod matrix;
pub mod nutation;
pub mod obliquity;
pub mod precession;
pub mod utils;

pub use angle::Angle;
pub use cio::{gcrs_to_cirs_matrix, CioLocator, CioSolution, CipCoordinates, EquationOfOrigins};
pub use errors::{AstroError, AstroResult, MathErrorKind};
pub use location::Location;
pub use matrix::{RotationMatrix3, Vector3};

pub mod test_helpers;
