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
//! - [`equation_of_origins`]: Relates CIO-based and equinox-based right ascension
//! - [`CioSolution`]: Bundles all CIO quantities for a given epoch
//!
//! # Usage
//!
//! For most use cases, compute a [`CioSolution`] from the NPB (nutation-precession-bias) matrix:
//!
//! ```
//! use celestial_core::cio::{gcrs_to_cirs_matrix, CioSolution};
//! use celestial_core::constants::J2000_JD;
//! use celestial_core::nutation::NutationIAU2006A;
//! use celestial_core::precession::PrecessionIAU2006;
//! use celestial_core::utils::jd_to_centuries;
//!
//! let (jd1, jd2) = (J2000_JD, 9000.0); // TT
//! let tt_centuries = jd_to_centuries(jd1, jd2);
//! let nut = NutationIAU2006A::new().compute(jd1, jd2)?;
//! let npb_matrix =
//!     PrecessionIAU2006::new().npb_matrix_iau2006a(tt_centuries, nut.delta_psi, nut.delta_eps);
//!
//! let solution = CioSolution::calculate(&npb_matrix, tt_centuries)?;
//! let cirs_matrix = gcrs_to_cirs_matrix(solution.cip.x, solution.cip.y, solution.s)?;
//! # Ok::<(), celestial_core::errors::AstroError>(())
//! ```

mod coordinates;
mod locator;
mod locator_terms;
mod origins;

pub use coordinates::CipCoordinates;
pub use locator::CioLocator;
pub use origins::equation_of_origins;

use crate::errors::{AstroError, AstroResult, MathErrorKind};
use crate::matrix::RotationMatrix3;

/// Builds the GCRS-to-CIRS rotation matrix from CIP coordinates and CIO locator.
///
/// This implements the IAU 2006 CIO-based transformation using three rotations:
/// R₃(-(E+s)) · R₂(d) · R₃(E) where E = atan2(Y, X) and d = atan(sqrt((X²+Y²)/(1-X²-Y²))).
pub fn gcrs_to_cirs_matrix(x: f64, y: f64, s: f64) -> AstroResult<RotationMatrix3> {
    let r2 = x * x + y * y;
    let defined = r2 <= 1.0 && s.is_finite();
    if !defined {
        return Err(AstroError::math_error(
            "gcrs_to_cirs_matrix",
            MathErrorKind::InvalidInput,
            "CIP must be finite and within the unit circle, and s finite",
        ));
    }
    let e = if r2 > 0.0 { libm::atan2(y, x) } else { 0.0 };
    let d = libm::atan(libm::sqrt(r2 / (1.0 - r2)));

    let mut matrix = RotationMatrix3::identity();
    matrix.rotate_z(e);
    matrix.rotate_y(d);
    matrix.rotate_z(-(e + s));

    Ok(matrix)
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

        let equation_of_origins = equation_of_origins(npb_matrix, s)?;

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

    // s is ERFA s06 at J2000.0 with X = Y = 0. ERFA eors returns s unchanged for the
    // identity matrix.
    #[test]
    fn cio_identity_matrix_gives_zero_cip_and_eo_equal_to_s() {
        let identity = crate::matrix::RotationMatrix3::identity();
        let solution = CioSolution::calculate(&identity, 0.0).unwrap();
        let s = -9.756652246326891e-9;
        let expected = CioSolution {
            cip: CipCoordinates::new(0.0, 0.0),
            s,
            equation_of_origins: s,
        };
        assert_eq!(solution, expected);
    }

    #[test]
    fn gcrs_to_cirs_matrix_with_zero_inputs_returns_identity() {
        let matrix = gcrs_to_cirs_matrix(0.0, 0.0, 0.0).unwrap();
        assert_eq!(matrix, RotationMatrix3::identity());
    }

    #[test]
    fn gcrs_to_cirs_matrix_is_rotation_matrix() {
        let matrix = gcrs_to_cirs_matrix(1e-6, 2e-6, 5e-9).unwrap();
        assert!(matrix.is_rotation_matrix(1e-14));
    }

    // Expected values are ERFA c2ixys outputs with Rust libm for atan2/atan/sin/cos.
    #[test]
    fn gcrs_to_cirs_matrix_matches_erfa_c2ixys() {
        let matrix = gcrs_to_cirs_matrix(1e-6, 1e-6, 1e-9).unwrap();
        let expected = RotationMatrix3::from_array([
            [0.9999999999995, -1.000500071679511e-9, -9.99999999e-7],
            [9.995000382900798e-10, 0.9999999999995, -1.000000001e-6],
            [1.0000000000000002e-6, 1e-6, 0.999999999999],
        ])
        .unwrap();
        assert_eq!(matrix, expected);
    }

    #[test]
    fn cio_solution_rejects_non_finite_inputs() {
        let mut nan = crate::matrix::RotationMatrix3::identity();
        nan.rotate_y(f64::NAN);
        assert!(CioSolution::calculate(&nan, 0.0).is_err());
        let identity = crate::matrix::RotationMatrix3::identity();
        assert!(CioSolution::calculate(&identity, f64::NAN).is_err());
    }

    #[test]
    fn gcrs_to_cirs_matrix_rejects_undefined_cip() {
        assert!(gcrs_to_cirs_matrix(1.0, 0.5, 0.0).is_err());
        assert!(gcrs_to_cirs_matrix(f64::NAN, 0.0, 0.0).is_err());
        assert!(gcrs_to_cirs_matrix(0.0, 0.0, f64::INFINITY).is_err());
        assert!(gcrs_to_cirs_matrix(1.0, 0.0, 0.0).is_ok());
    }
}
