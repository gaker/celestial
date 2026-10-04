//! Equation of Origins computation.
//!
//! The Equation of Origins (EO) is the arc on the CIP equator between the Celestial
//! Intermediate Origin (CIO) and the equinox. It connects the modern CIO-based system
//! to the classical equinox-based system.
//!
//! In practice, EO is the difference between ERA (Earth Rotation Angle, measured from
//! the CIO) and GAST (Greenwich Apparent Sidereal Time, measured from the equinox):
//!
//! ```text
//! GAST = ERA - EO
//! ```
//!
//! Typical values are a few arcseconds, slowly drifting due to precession.
//!
//! # When to use this
//!
//! - Converting between CIO-based and equinox-based right ascension
//! - Relating ERA to sidereal time
//! - Legacy interoperability with equinox-based catalogs and software
//!
//! For purely CIO-based work (GCRS to CIRS), you don't need EO directly — use
//! [`CioSolution`](crate::cio::CioSolution) instead.

use crate::errors::{AstroError, AstroResult, MathErrorKind};
use crate::matrix::RotationMatrix3;

/// Computes EO from the NPB (nutation-precession-bias) matrix and CIO locator.
///
/// This is the rigorous method. It extracts the equinox position from
/// the NPB matrix and computes the arc to the CIO.
///
/// # Arguments
///
/// * `npb_matrix` - Combined frame bias, precession, and nutation rotation matrix
/// * `s` - CIO locator in radians (from [`CioLocator::calculate`](crate::cio::CioLocator::calculate))
///
/// # Returns
///
/// Equation of Origins in radians. Positive when the equinox is west of the CIO.
pub fn equation_of_origins(npb_matrix: &RotationMatrix3, s: f64) -> AstroResult<f64> {
    let eo = equinox_to_cio_arc(npb_matrix.elements(), s);
    if eo.is_finite() {
        return Ok(eo);
    }
    Err(AstroError::math_error(
        "equation_of_origins",
        MathErrorKind::NotFinite,
        "Equation of origins is undefined for this matrix and locator",
    ))
}

fn equinox_to_cio_arc(matrix: &[[f64; 3]; 3], s: f64) -> f64 {
    let x = matrix[2][0];
    let ax = x / (1.0 + matrix[2][2]);
    let xs = 1.0 - ax * x;
    let ys = -ax * matrix[2][1];
    let zs = -x;

    let p = matrix[0][0] * xs + matrix[0][1] * ys + matrix[0][2] * zs;
    let q = matrix[1][0] * xs + matrix[1][1] * ys + matrix[1][2] * zs;

    if p != 0.0 || q != 0.0 {
        s - libm::atan2(q, p)
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::RotationMatrix3;

    #[test]
    fn test_equation_of_origins_identity_matrix() {
        let identity = RotationMatrix3::identity();
        let s = 0.0;

        let eo = equation_of_origins(&identity, s).unwrap();

        assert_eq!(eo, 0.0);
    }

    // Expected values are ERFA eors outputs with Rust libm for atan2.
    #[test]
    fn test_equation_of_origins_small_rotation_matches_erfa_eors() {
        let mut small_rotation = RotationMatrix3::identity();
        small_rotation.rotate_z(1e-6);

        let eo = equation_of_origins(&small_rotation, 1e-7).unwrap();

        assert_eq!(eo, 1.1e-6);
    }

    #[test]
    fn test_equation_of_origins_matches_erfa_eors() {
        let npb = RotationMatrix3::from_array([
            [
                0.9999989440476104,
                -0.0013328817612400115,
                -0.0005790767434730085,
            ],
            [
                0.0013328582543089545,
                0.9999991109044506,
                -4.097782710401556e-5,
            ],
            [
                0.0005791308472168153,
                4.0205956615939944e-5,
                0.9999998314954572,
            ],
        ])
        .unwrap();

        let eo = equation_of_origins(&npb, -1.220040848472272e-8).unwrap();

        assert_eq!(eo, -0.0013328827151307448);
    }

    #[test]
    fn test_rejects_undefined_equation_of_origins() {
        let mut antipodal_pole = RotationMatrix3::identity();
        antipodal_pole.rotate_x(crate::constants::PI);
        assert!(equation_of_origins(&antipodal_pole, 0.0).is_err());
        let identity = RotationMatrix3::identity();
        assert!(equation_of_origins(&identity, f64::NAN).is_err());
    }
}
