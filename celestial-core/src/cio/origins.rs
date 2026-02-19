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

use crate::errors::AstroResult;
use crate::matrix::RotationMatrix3;

/// Computes the Equation of Origins from precession-nutation quantities.
///
/// This is a stateless utility type — all methods are associated functions.
pub struct EquationOfOrigins;

impl EquationOfOrigins {
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
    pub fn from_npb_and_locator(npb_matrix: &RotationMatrix3, s: f64) -> AstroResult<f64> {
        let matrix = npb_matrix.elements();

        let x = matrix[2][0];
        let ax = x / (1.0 + matrix[2][2]);
        let xs = 1.0 - ax * x;
        let ys = -ax * matrix[2][1];
        let zs = -x;

        let p = matrix[0][0] * xs + matrix[0][1] * ys + matrix[0][2] * zs;
        let q = matrix[1][0] * xs + matrix[1][1] * ys + matrix[1][2] * zs;

        let eo = if p != 0.0 || q != 0.0 {
            s - libm::atan2(q, p)
        } else {
            s
        };

        Ok(eo)
    }

    /// Approximates EO from the equation of equinoxes and CIO locator.
    ///
    /// A simpler but less accurate method: `EO ≈ EE - s`, where EE is the equation
    /// of equinoxes (nutation in right ascension). Useful when you already have EE
    /// but not the full NPB matrix.
    ///
    /// For sub-milliarcsecond accuracy, use [`from_npb_and_locator`](Self::from_npb_and_locator).
    ///
    /// # Arguments
    ///
    /// * `equation_of_equinoxes` - Nutation in right ascension (radians)
    /// * `s` - CIO locator (radians)
    pub fn from_equation_of_equinoxes_approximation(
        equation_of_equinoxes: f64,
        s: f64,
    ) -> AstroResult<f64> {
        let eo_approx = equation_of_equinoxes - s;

        Ok(eo_approx)
    }

    /// Converts EO from radians to arcseconds.
    pub fn to_arcseconds(eo_radians: f64) -> f64 {
        eo_radians * crate::constants::RAD_TO_DEG * 3600.0
    }

    /// Converts EO from radians to milliarcseconds.
    pub fn to_milliarcseconds(eo_radians: f64) -> f64 {
        eo_radians * crate::constants::RAD_TO_DEG * 3600.0 * 1000.0
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

        let eo = EquationOfOrigins::from_npb_and_locator(&identity, s).unwrap();

        assert!(
            eo.abs() < 1e-15,
            "EO should be zero for identity matrix: {}",
            eo
        );
    }

    #[test]
    fn test_equation_of_origins_small_rotation() {
        let mut small_rotation = RotationMatrix3::identity();
        small_rotation.rotate_z(1e-6);
        let s = 1e-7;

        let eo = EquationOfOrigins::from_npb_and_locator(&small_rotation, s).unwrap();

        assert!(eo.abs() < 1e-3, "EO should be small: {}", eo);
        assert!(
            eo.abs() > 1e-12,
            "EO should be non-zero for rotation: {}",
            eo
        );
    }

    #[test]
    fn test_equation_of_origins_approximation() {
        let equation_of_equinoxes = 1e-6;
        let s = 5e-7;

        let eo =
            EquationOfOrigins::from_equation_of_equinoxes_approximation(equation_of_equinoxes, s)
                .unwrap();

        let expected = equation_of_equinoxes - s;
        assert!(
            (eo - expected).abs() < 1e-15,
            "Approximation should match expected value"
        );
    }

    #[test]
    fn test_unit_conversions() {
        let eo_rad = 1e-6;

        let eo_arcsec = EquationOfOrigins::to_arcseconds(eo_rad);
        let eo_mas = EquationOfOrigins::to_milliarcseconds(eo_rad);

        assert!((eo_arcsec * 1000.0 - eo_mas).abs() < 1e-10);
        assert!(eo_arcsec > 0.0 && eo_arcsec < 1.0);
    }

    #[test]
    fn test_large_eo_validation() {
        let mut normal_matrix = RotationMatrix3::identity();
        normal_matrix.rotate_z(1e-6);
        let normal_s = 1e-7;

        let eo = EquationOfOrigins::from_npb_and_locator(&normal_matrix, normal_s).unwrap();

        assert!(eo.abs() < 1e-3);
    }
}
