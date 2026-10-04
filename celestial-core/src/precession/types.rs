//! Types for representing precession computation results.
//!
//! The three rotation matrices in a precession result (bias, precession, and
//! combined) transform between different celestial reference frames.
//!
//! # Frame Relationships
//!
//! - **Bias matrix**: Transforms from GCRS (Geocentric Celestial Reference System)
//!   to mean equator and equinox of J2000.0, accounting for the small offset
//!   between the ICRS pole and the mean celestial pole at J2000.0.
//!
//! - **Precession matrix**: Transforms from mean equator and equinox of J2000.0
//!   to the mean equator and equinox of the target date.
//!
//! - **Bias-precession matrix**: The combined transformation from GCRS directly
//!   to the mean equator and equinox of the target date.

use crate::matrix::RotationMatrix3;

/// Complete result of a precession computation.
///
/// Contains all three rotation matrices needed for transformations between
/// GCRS and the mean equator/equinox of a target date. The matrices are
/// computed together because they share intermediate calculations.
///
/// # Usage
///
/// For most transformations, use `bias_precession_matrix` directly. The
/// individual `bias_matrix` and `precession_matrix` are provided for cases
/// where only one component is needed, or for debugging and validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecessionResult {
    /// The frame bias matrix (GCRS to mean J2000.0).
    pub bias_matrix: RotationMatrix3,

    /// The precession matrix (mean J2000.0 to mean of date).
    pub precession_matrix: RotationMatrix3,

    /// The combined bias-precession matrix (GCRS to mean of date).
    pub bias_precession_matrix: RotationMatrix3,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FukushimaWilliamsAngles {
    pub gamma_bar: f64,
    pub phi_bar: f64,
    pub psi_bar: f64,
    pub epsilon_a: f64,
}
