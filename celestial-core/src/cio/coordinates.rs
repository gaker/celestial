//! CIP (Celestial Intermediate Pole) coordinates for the IAU 2000/2006 precession-nutation model.
//!
//! The Celestial Intermediate Pole defines the axis around which Earth rotates in the
//! Celestial Intermediate Reference System (CIRS). Its position relative to the GCRS
//! is described by two small angles X and Y, which encode the combined effects of
//! precession and nutation.
//!
//! # When to use this
//!
//! CIP coordinates are needed when:
//! - Converting between GCRS (geocentric celestial) and CIRS (intermediate) frames
//! - Computing Earth Rotation Angle for sidereal time
//! - Implementing the full IAU 2000/2006 transformation chain
//!
//! For most high-precision applications, you'll extract CIP coordinates from a
//! precession-nutation matrix rather than computing them directly.
//!
//! # Coordinate ranges
//!
//! X and Y are stored in radians. Current values are on the order of 10^-7 radians
//! (~0.02 arcseconds). The validation threshold of 0.2 radians (~11 degrees) catches
//! obviously invalid matrices while allowing for long-term secular drift.

use crate::errors::{AstroError, AstroResult};
use crate::matrix::RotationMatrix3;

/// Position of the Celestial Intermediate Pole in the GCRS.
///
/// X and Y are the direction cosines of the CIP unit vector projected onto
/// the GCRS equatorial plane. They're extracted from elements [2][0] and [2][1]
/// of the NPB (nutation-precession-bias) matrix.
///
/// Units are radians internally. Use [`to_arcseconds`](Self::to_arcseconds) for display.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CipCoordinates {
    pub x: f64,
    pub y: f64,
}

impl CipCoordinates {
    /// Creates CIP coordinates from X and Y values in radians.
    ///
    /// No validation is performed. For coordinates extracted from a real
    /// precession-nutation matrix, use [`from_npb_matrix`](Self::from_npb_matrix).
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Extracts CIP coordinates from a nutation-precession-bias matrix.
    ///
    /// The NPB matrix transforms GCRS to mean-of-date coordinates. X and Y
    /// are taken from the third row (elements [2][0] and [2][1]), which
    /// represents the CIP direction in GCRS.
    ///
    /// # Errors
    ///
    /// Returns an error if X or Y exceeds 0.2 radians, indicating the matrix
    /// doesn't represent a valid Earth orientation.
    pub fn from_npb_matrix(npb_matrix: &RotationMatrix3) -> AstroResult<Self> {
        let matrix = npb_matrix.elements();

        let x = matrix[2][0];
        let y = matrix[2][1];

        if x.abs() > 0.2 || y.abs() > 0.2 {
            return Err(AstroError::math_error(
                "CIP coordinate extraction",
                crate::errors::MathErrorKind::InvalidInput,
                &format!(
                    "CIP coordinates out of reasonable range: X={:.6}, Y={:.6}",
                    x, y
                ),
            ));
        }

        Ok(Self { x, y })
    }

    /// Distance of the CIP from the GCRS pole, in radians.
    ///
    /// This is sqrt(X^2 + Y^2), useful for gauging the total pole offset.
    pub fn magnitude(&self) -> f64 {
        libm::sqrt(self.x * self.x + self.y * self.y)
    }

    /// Returns (X, Y) converted to degrees.
    pub fn to_degrees(&self) -> (f64, f64) {
        (
            self.x * crate::constants::RAD_TO_DEG,
            self.y * crate::constants::RAD_TO_DEG,
        )
    }

    /// Returns (X, Y) converted to arcseconds.
    ///
    /// Arcseconds are the conventional unit for reporting CIP coordinates
    /// in publications and IERS bulletins.
    pub fn to_arcseconds(&self) -> (f64, f64) {
        (
            self.x * crate::constants::RAD_TO_DEG * 3600.0,
            self.y * crate::constants::RAD_TO_DEG * 3600.0,
        )
    }
}

impl std::fmt::Display for CipCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (x_as, y_as) = self.to_arcseconds();
        write!(f, "CIP(X={:.3}\", Y={:.3}\")", x_as, y_as)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::RotationMatrix3;

    #[test]
    fn test_cip_coordinates_creation() {
        let cip = CipCoordinates::new(1e-6, -2e-6);
        assert_eq!(cip.x, 1e-6);
        assert_eq!(cip.y, -2e-6);
    }

    #[test]
    fn test_cip_magnitude() {
        let cip = CipCoordinates::new(3e-6, 4e-6);
        assert!((cip.magnitude() - 5e-6).abs() < 1e-12);
    }

    #[test]
    fn test_cip_display() {
        let cip = CipCoordinates::new(1e-6, -1e-6);
        let display = format!("{}", cip);
        assert!(display.contains("CIP"));
        assert!(display.contains("0.206")); // ~1e-6 radians in arcseconds
    }

    #[test]
    fn test_cip_from_identity_matrix() {
        // Identity matrix should give CIP coordinates of exactly zero
        let identity = RotationMatrix3::identity();
        let cip = CipCoordinates::from_npb_matrix(&identity).unwrap();

        // For identity matrix, third column is [0, 0, 1]
        // So X = 0, Y = 0 exactly
        assert_eq!(cip.x, 0.0);
        assert_eq!(cip.y, 0.0);
    }

    #[test]
    fn test_cip_validation_error() {
        // Create a matrix with unreasonably large CIP coordinates in third row
        // CIP X,Y are extracted from matrix[2][0] and matrix[2][1]
        // This should trigger the validation error (threshold is 0.2 radians ~= 11.5 degrees)
        let invalid_matrix = RotationMatrix3::from_array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.25, 0.05, 1.0], // X = 0.25 (too large - exceeds 0.2 rad threshold)
        ]);

        let result = CipCoordinates::from_npb_matrix(&invalid_matrix);

        // Should return an error due to coordinates being out of reasonable range
        assert!(result.is_err());

        let error_message = format!("{}", result.unwrap_err());
        assert!(error_message.contains("CIP coordinates out of reasonable range"));
    }

    #[test]
    fn test_cip_reasonable_values() {
        // Test that reasonable CIP coordinate values work correctly
        let reasonable_cip = CipCoordinates::new(1e-6, 1e-6);

        // These are the exact computed values for our test case
        assert_eq!(reasonable_cip.x, 1e-6);
        assert_eq!(reasonable_cip.y, 1e-6);
        assert_eq!(reasonable_cip.magnitude(), libm::sqrt(2e-12_f64)); // sqrt(x² + y²)
    }
}
