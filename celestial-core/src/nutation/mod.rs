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
//! All `compute(jd1, jd2)` methods accept a two-part Julian Date in TDB (Barycentric
//! Dynamical Time). The split preserves precision: typically `jd1 = 2451545.0` (J2000.0)
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
//! // At J2000.0: Δψ ≈ -0.04 arcsec, Δε ≈ -0.007 arcsec
//! println!("Δψ = {:.6} rad", result.delta_psi);
//! println!("Δε = {:.6} rad", result.delta_eps);
//! ```
//!
//! # Sub-modules
//!
//! - [`iau2000a`]: Full IAU 2000A model (678 lunisolar + 687 planetary terms)
//! - [`iau2000b`]: Truncated IAU 2000B model (77 terms + planetary bias)
//! - [`iau2006a`]: IAU 2000A with J2 corrections for IAU 2006 precession compatibility
//! - [`fundamental_args`]: Delaunay arguments and planetary mean longitudes
//! - [`lunisolar_terms`]: Coefficient table for lunisolar nutation series
//! - [`planetary_terms`]: Coefficient table for planetary nutation series
//! - [`types`]: [`NutationResult`] and [`NutationModel`] wrapper

#[cfg(feature = "erfa-tests")]
pub mod fundamental_args;

#[cfg(not(feature = "erfa-tests"))]
mod fundamental_args;

pub mod iau2000a;
pub mod iau2000b;
pub mod iau2006a;
pub mod lunisolar_terms;
pub mod planetary_terms;
pub mod types;

pub use fundamental_args::{IERS2010FundamentalArgs, MHB2000FundamentalArgs};
pub use iau2000a::NutationIAU2000A;
pub use iau2000b::NutationIAU2000B;
pub use iau2006a::NutationIAU2006A;
pub use types::{NutationModel, NutationResult};
