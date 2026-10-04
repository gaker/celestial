//! Nutation models for computing oscillations in Earth's rotational axis.
//!
//! Nutation is the short-period oscillation of Earth's rotational axis about its mean
//! position, superimposed on the longer-term precession. It arises from gravitational
//! torques exerted by the Moon, Sun, and planets on Earth's equatorial bulge. The
//! principal component has a period of 18.6 years (the lunar nodal period) with an
//! amplitude of about 9 arcseconds in obliquity.
//!
//! This module provides implementations of the IAU standard nutation models:
//!
//! | Model | Lunisolar Terms | Planetary Terms | Precision | Use Case |
//! |-------|-----------------|-----------------|-----------|----------|
//! | [`NutationIAU2000A`] | 678 | 687 | ~0.1 µas | High-precision astrometry, VLBI |
//! | [`NutationIAU2000B`] | 77 | 0 (bias only) | ~1 mas | General ephemeris, telescope pointing |
//! | [`NutationIAU2006A`] | 678 | 687 | ~0.1 µas | Use with IAU 2006 precession |
//!
//! # Output
//!
//! All models return [`NutationResult`] containing:
//! - `delta_psi`: Nutation in longitude (radians)
//! - `delta_eps`: Nutation in obliquity (radians)
//!
//! # Time Argument
//!
//! All `compute(jd1, jd2)` methods accept a two-part Julian Date in TT (Terrestrial
//! Time). The split preserves precision: typically `jd1 = 2451545.0` (J2000.0)
//! and `jd2` = days from that epoch.
//!
//! # Example
//!
//! ```
//! use celestial_core::nutation::NutationIAU2006A;
//!
//! let nutation = NutationIAU2006A::new();
//! let result = nutation.compute(2451545.0, 0.0).unwrap();
//!
//! // At J2000.0: Δψ ≈ -13.932 arcsec, Δε ≈ -5.769 arcsec
//! println!("Δψ = {:.6} rad", result.delta_psi);
//! println!("Δε = {:.6} rad", result.delta_eps);
//! ```
//!
//! # Contents
//!
//! - [`NutationIAU2000A`]: Full IAU 2000A model (678 lunisolar + 687 planetary terms)
//! - [`NutationIAU2000B`]: Truncated IAU 2000B model (77 terms + planetary bias)
//! - [`NutationIAU2006A`]: IAU 2000A with J2 corrections for IAU 2006 precession compatibility
//! - `fundamental_args`: Delaunay arguments and planetary mean longitudes
//! - `lunisolar_terms`: Coefficient table for lunisolar nutation series
//! - `planetary_terms`: Coefficient table for planetary nutation series
//! - [`NutationResult`]: nutation in longitude and obliquity

#[cfg(feature = "test-support")]
pub mod fundamental_args;
#[cfg(not(feature = "test-support"))]
pub(crate) mod fundamental_args;

#[cfg(feature = "test-support")]
pub mod lunisolar_terms;
#[cfg(not(feature = "test-support"))]
mod lunisolar_terms;

#[cfg(feature = "test-support")]
pub mod planetary_terms;
#[cfg(not(feature = "test-support"))]
mod planetary_terms;

mod iau2000a;
mod iau2000b;
mod iau2006a;
mod types;

pub use iau2000a::NutationIAU2000A;
pub use iau2000b::NutationIAU2000B;
pub use iau2006a::NutationIAU2006A;
pub use types::NutationResult;
