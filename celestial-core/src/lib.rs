//! Low-level astronomical calculations for coordinate transformations.
//!
//! `celestial-core` provides the mathematical building blocks for celestial mechanics:
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
//! | [`errors`] | [`AstroError`](errors::AstroError) and [`AstroResult`](errors::AstroResult) |
//!
//! # Coordinate Transformation Pipeline
//!
//! GCRS → CIRS transformation (CIO-based):
//!
//! ```
//! use celestial_core::constants::J2000_JD;
//! use celestial_core::nutation::NutationIAU2006A;
//! use celestial_core::precession::PrecessionIAU2006;
//! use celestial_core::utils::jd_to_centuries;
//! use celestial_core::cio::{gcrs_to_cirs_matrix, CioSolution};
//!
//! // TT as a two-part Julian Date
//! let (jd1, jd2) = (J2000_JD, 9000.0);
//! let tt_centuries = jd_to_centuries(jd1, jd2);
//!
//! // 1. Compute precession-nutation-bias matrix
//! let nutation = NutationIAU2006A::new().compute(jd1, jd2)?;
//! let npb = PrecessionIAU2006::new().npb_matrix_iau2006a(
//!     tt_centuries,
//!     nutation.delta_psi,
//!     nutation.delta_eps,
//! );
//!
//! // 2. Extract CIO quantities
//! let cio = CioSolution::calculate(&npb, tt_centuries)?;
//!
//! // 3. Build GCRS→CIRS matrix
//! let matrix = gcrs_to_cirs_matrix(cio.cip.x, cio.cip.y, cio.s)?;
//! # Ok::<(), celestial_core::errors::AstroError>(())
//! ```
//!
//! # Import paths
//!
//! Items are imported from their module; the crate root re-exports nothing.
//!
//! ```
//! use celestial_core::angle::Angle;
//! use celestial_core::errors::{AstroError, AstroResult, MathErrorKind};
//! use celestial_core::location::Location;
//! use celestial_core::matrix::{RotationMatrix3, Vector3};
//! ```
//!
//! # Design Notes
//!
//! - **Two-part Julian Dates**: Functions accepting `(jd1, jd2)` preserve precision by
//!   splitting the date. Typically `jd1 = 2451545.0` (J2000.0) and `jd2` is days from epoch.
//!
//! - **Radians internally**: All angular computations use radians. The [`Angle`](angle::Angle) type
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

#[cfg(any(test, feature = "test-support"))]
pub mod test_helpers;
