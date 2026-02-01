//! Celestial Intermediate Origin (CIO) based celestial-to-terrestrial transformations.
//!
//! This module implements the IAU 2000/2006 CIO-based transformation from GCRS (Geocentric
//! Celestial Reference System) to CIRS (Celestial Intermediate Reference System). The CIO
//! approach is the modern replacement for the classical equinox-based method.
//!
//! # Components
//!
//! - [`CipCoordinates`]: X/Y coordinates of the Celestial Intermediate Pole
//! - [`CioLocator`]: The CIO locator `s`, positioning the origin on the CIP equator
//! - [`EquationOfOrigins`]: Relates CIO-based and equinox-based right ascension
//! - [`CioSolution`]: Bundles all CIO quantities for a given epoch
//!
//! # Usage
//!
//! For most use cases, compute a [`CioSolution`] from the NPB (nutation-precession-bias) matrix:
//!
//! ```ignore
//! let solution = CioSolution::calculate(&npb_matrix, tt_centuries)?;
//! let cirs_matrix = gcrs_to_cirs_matrix(solution.cip.x, solution.cip.y, solution.s);
//! ```

pub mod coordinates;
pub mod locator;
pub mod origins;

pub use coordinates::CipCoordinates;
pub use locator::CioLocator;
pub use origins::EquationOfOrigins;

use crate::errors::AstroResult;
use crate::matrix::RotationMatrix3;

/// Builds the GCRS-to-CIRS rotation matrix from CIP coordinates and CIO locator.
///
/// This implements the IAU 2006 CIO-based transformation using three rotations:
/// R₃(E) · R₂(d) · R₃(-(E+s)) where E = atan2(Y, X) and d = atan(sqrt(X²+Y²/(1-X²-Y²))).
pub fn gcrs_to_cirs_matrix(x: f64, y: f64, s: f64) -> RotationMatrix3 {
    let r2 = x * x + y * y;
    let e = if r2 > 0.0 { y.atan2(x) } else { 0.0 };
    let d = (r2 / (1.0 - r2)).sqrt().atan();

    let mut matrix = RotationMatrix3::identity();
    matrix.rotate_z(e);
    matrix.rotate_y(d);
    matrix.rotate_z(-(e + s));

    matrix
}

/// All CIO-based quantities for a given epoch.
///
/// Bundles CIP coordinates, CIO locator, and equation of origins — everything needed
/// for the GCRS↔CIRS transformation.
#[derive(Debug, Clone, PartialEq)]
pub struct CioSolution {
    /// CIP X/Y coordinates (radians)
    pub cip: CipCoordinates,
    /// CIO locator s (radians)
    pub s: f64,
    /// Equation of origins (radians) — difference between CIO-based and equinox-based RA
    pub equation_of_origins: f64,
}

impl CioSolution {
    /// Computes all CIO quantities from an NPB matrix and TT centuries since J2000.
    pub fn calculate(
        npb_matrix: &crate::matrix::RotationMatrix3,
        tt_centuries: f64,
    ) -> AstroResult<Self> {
        let cip = CipCoordinates::from_npb_matrix(npb_matrix)?;

        let locator = CioLocator::iau2006a(tt_centuries);
        let s = locator.calculate(cip.x, cip.y)?;

        let equation_of_origins = EquationOfOrigins::from_npb_and_locator(npb_matrix, s)?;

        Ok(Self {
            cip,
            s,
            equation_of_origins,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cio_identity_matrix_returns_zero_components() {
        let identity = crate::matrix::RotationMatrix3::identity();
        let solution = CioSolution::calculate(&identity, 0.0).unwrap();

        assert!(solution.cip.x.abs() < 1e-15);
        assert!(solution.cip.y.abs() < 1e-15);
        assert!(solution.s.abs() < 1e-8);
        assert!(solution.equation_of_origins.abs() < 1e-8);
    }

    #[test]
    fn gcrs_to_cirs_matrix_with_zero_inputs_returns_identity() {
        let matrix = gcrs_to_cirs_matrix(0.0, 0.0, 0.0);
        let identity = RotationMatrix3::identity();
        assert!(matrix.max_difference(&identity) < 1e-15);
    }

    #[test]
    fn gcrs_to_cirs_matrix_is_rotation_matrix() {
        let matrix = gcrs_to_cirs_matrix(1e-6, 2e-6, 5e-9);
        assert!(matrix.is_rotation_matrix(1e-14));
    }

    #[test]
    fn gcrs_to_cirs_matrix_small_cip_produces_near_identity() {
        let x = 1e-6;
        let y = 1e-6;
        let s = 1e-9;
        let matrix = gcrs_to_cirs_matrix(x, y, s);
        let identity = RotationMatrix3::identity();
        assert!(matrix.max_difference(&identity) < 1e-5);
    }
}
