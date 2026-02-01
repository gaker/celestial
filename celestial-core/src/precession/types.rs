//! Types for representing precession computation results.
//!
//! This module provides type aliases and result structures for precession
//! calculations. The three rotation matrices (bias, precession, and combined)
//! are used to transform between different celestial reference frames.
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

/// Rotation matrix accounting for frame bias between GCRS and mean J2000.0.
///
/// The bias arises because the ICRS axes are defined kinematically rather than
/// dynamically, resulting in a small angular offset from the mean equator and
/// equinox of J2000.0. This offset is approximately 23 milliarcseconds in the
/// equator (dx) and 7 milliarcseconds in the ecliptic (dy).
pub type BiasMatrix = RotationMatrix3;

/// Rotation matrix for precession from J2000.0 to a target epoch.
///
/// Precession is the slow, gravity-induced wobble of Earth's rotational axis,
/// with a period of approximately 26,000 years. This matrix transforms
/// coordinates from the mean equator and equinox of J2000.0 to the mean
/// equator and equinox of the specified date.
pub type PrecessionMatrix = RotationMatrix3;

/// Combined bias and precession rotation matrix.
///
/// This matrix is the product of the bias and precession matrices, providing
/// a single transformation from GCRS to the mean equator and equinox of the
/// target date. Using the combined matrix is more efficient than applying
/// bias and precession separately when both corrections are needed.
pub type BiasPrecessionMatrix = RotationMatrix3;

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
#[derive(Debug, Clone)]
pub struct PrecessionResult {
    /// The frame bias matrix (GCRS to mean J2000.0).
    pub bias_matrix: BiasMatrix,

    /// The precession matrix (mean J2000.0 to mean of date).
    pub precession_matrix: PrecessionMatrix,

    /// The combined bias-precession matrix (GCRS to mean of date).
    pub bias_precession_matrix: BiasPrecessionMatrix,
}
