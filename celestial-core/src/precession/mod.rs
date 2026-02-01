//! Precession models for transforming coordinates between epochs.
//!
//! Precession is the slow, continuous change in the orientation of Earth's rotational
//! axis caused by gravitational torques from the Sun and Moon acting on Earth's
//! equatorial bulge. The axis traces a cone in space with a period of approximately
//! 26,000 years (the "Platonic year"), causing the celestial pole to drift among
//! the stars and the vernal equinox to move westward along the ecliptic.
//!
//! This module provides rotation matrices that transform celestial coordinates from
//! one epoch to another, accounting for the accumulated precession between epochs.
//!
//! # Available Models
//!
//! ## IAU 2000
//!
//! The IAU 2000 precession model ([`PrecessionIAU2000`]) uses the Lieske (1977)
//! precession angles with corrections from the IAU 2000A nutation model. It computes
//! three rotation angles (psi_A, omega_A, chi_A) to construct the precession matrix.
//!
//! The frame bias matrix accounts for the offset between the J2000.0 dynamical frame
//! and the ICRS (International Celestial Reference System), which amounts to a few
//! milliarcseconds.
//!
//! ## IAU 2006
//!
//! The IAU 2006 precession model ([`PrecessionIAU2006`]) uses the Fukushima-Williams
//! four-angle formulation (gamma_bar, phi_bar, psi_bar, epsilon_A). This parameterization
//! provides improved numerical stability and separates the frame bias from the
//! precession proper.
//!
//! The Fukushima-Williams angles represent:
//! - **gamma_bar**: Frame bias in right ascension
//! - **phi_bar**: Obliquity of the ecliptic at J2000.0
//! - **psi_bar**: Precession in longitude
//! - **epsilon_A**: Mean obliquity of date
//!
//! IAU 2006 is the current standard and should be preferred for new applications.
//!
//! # Output Matrices
//!
//! Both models produce a [`PrecessionResult`] containing:
//!
//! - **bias_matrix**: Transforms from GCRS (Geocentric Celestial Reference System)
//!   to the mean equator and equinox of J2000.0. This is constant for a given model.
//!
//! - **precession_matrix**: Transforms from the mean equator and equinox of J2000.0
//!   to the mean equator and equinox of date. At J2000.0 (t=0), this is the identity.
//!
//! - **bias_precession_matrix**: The combined transformation from GCRS to the mean
//!   equator and equinox of date (bias_precession = precession * bias for IAU 2000,
//!   computed directly via Fukushima-Williams for IAU 2006).
//!
//! # Time Argument
//!
//! IAU 2000 takes time as Julian centuries of TT (Terrestrial Time) since J2000.0.
//! IAU 2006 takes a two-part Julian Date in TT for improved numerical precision.

pub mod iau2000;
pub mod iau2006;
pub mod types;

pub use iau2000::PrecessionIAU2000;
pub use iau2006::PrecessionIAU2006;
pub use types::{BiasMatrix, PrecessionMatrix, PrecessionResult};

/// Trait for types that can compute precession matrices.
///
/// Implementors provide access to both IAU 2000 and IAU 2006 precession models,
/// allowing consistent precession calculations across different time representations.
pub trait PrecessionCalculator {
    /// Computes precession using the IAU 2000 model.
    ///
    /// # Arguments
    ///
    /// * `tt_centuries` - Julian centuries of TT since J2000.0
    fn precession_iau2000(&self, tt_centuries: f64) -> crate::AstroResult<PrecessionResult>;

    /// Computes precession using the IAU 2006 model.
    ///
    /// # Arguments
    ///
    /// * `tt_centuries` - Julian centuries of TT since J2000.0
    fn precession_iau2006(&self, tt_centuries: f64) -> crate::AstroResult<PrecessionResult>;
}
