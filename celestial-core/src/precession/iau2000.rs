//! IAU 2000 precession model.
//!
//! This module implements the IAU 2000A precession model, which describes
//! the gradual shift of Earth's rotational axis and equatorial plane due to
//! gravitational torques from the Sun and Moon acting on Earth's equatorial
//! bulge.
//!
//! # Background
//!
//! Precession causes the celestial pole to trace a circle around the ecliptic
//! pole over approximately 26,000 years. The IAU 2000 model computes precession
//! as a correction to the earlier IAU 1976 (Lieske) precession, applying
//! small adjustments derived from VLBI observations.
//!
//! The model separates two components:
//!
//! - **Frame bias**: A small fixed rotation accounting for the offset between
//!   the dynamical mean equator and equinox of J2000.0 and the ICRS origin.
//!
//! - **Precession**: Time-dependent rotation from J2000.0 to the mean equator
//!   and equinox of date.
//!
//! # Reference Frame
//!
//! The frame bias rotates from the Geocentric Celestial Reference System (GCRS)
//! to the mean equator and equinox of J2000.0. The bias parameters are:
//!
//! - Right ascension of the pole: -14.6 mas
//! - Longitude of the pole: -41.775 mas
//! - Obliquity of the pole: -6.8192 mas
//!
//! These values are from the IERS Conventions (2010), Table 5.1.
//!
//! # Algorithm
//!
//! The precession matrix is constructed using the Lieske (1979) angles with
//! IAU 2000 corrections:
//!
//! 1. Compute the Lieske precession angles (psi_A, omega_A, chi_A) using
//!    polynomial expressions in Julian centuries from J2000.0
//!
//! 2. Apply the IAU 2000 precession-rate corrections:
//!    - dpsi_pr = -0.29965"/century
//!    - deps_pr = -0.02524"/century
//!
//! 3. Construct the rotation matrix R_x(eps_0) R_z(-psi_A) R_x(-omega_A) R_z(chi_A)
//!
//! # Accuracy
//!
//! The IAU 2000 precession model is accurate to approximately 0.3 mas/century
//! over several centuries around J2000.0. For higher accuracy over longer time
//! spans, see the IAU 2006 precession model which uses improved polynomial
//! expressions.
//!
//! # References
//!
//! - IERS Conventions (2010), Chapter 5
//! - Lieske et al. (1977), A&A 58, 1-16
//! - Mathews, Herring & Buffett (2002), J. Geophys. Res. 107

use super::types::PrecessionResult;
use crate::constants::ARCSEC_TO_RAD;
use crate::matrix::RotationMatrix3;

/// IAU 2000 precession computation.
///
/// Computes the precession and frame bias matrices for transforming
/// coordinates between J2000.0 and the mean equator and equinox of date.
///
/// # Example
///
/// ```
/// use celestial_core::precession::PrecessionIAU2000;
///
/// let precession = PrecessionIAU2000::new();
///
/// // Compute for 0.5 Julian centuries (50 years) after J2000.0
/// let result = precession.compute(0.5).unwrap();
///
/// // Access individual matrices
/// let _bias = &result.bias_matrix;           // GCRS to mean J2000.0
/// let _prec = &result.precession_matrix;     // Mean J2000.0 to mean of date
/// let _combined = &result.bias_precession_matrix; // GCRS to mean of date
/// ```
pub struct PrecessionIAU2000;

impl PrecessionIAU2000 {
    /// Creates a new IAU 2000 precession calculator.
    pub fn new() -> Self {
        Self
    }

    /// Computes precession matrices for the given time.
    ///
    /// # Arguments
    ///
    /// * `tt_centuries` - Julian centuries of TT from J2000.0 (JD 2451545.0).
    ///   Positive values are after J2000.0, negative values before.
    ///
    /// # Returns
    ///
    /// A [`PrecessionResult`] containing:
    /// - `bias_matrix`: Rotation from GCRS to mean equator/equinox of J2000.0
    /// - `precession_matrix`: Rotation from mean J2000.0 to mean of date
    /// - `bias_precession_matrix`: Combined rotation from GCRS to mean of date
    ///
    /// The combined matrix is computed as `precession_matrix * bias_matrix`.
    ///
    /// # Notes
    ///
    /// The bias matrix is constant (independent of time) and accounts for
    /// the small misalignment between the ICRS and the dynamical frame.
    /// The precession matrix is identity at t=0 and diverges with time.
    pub fn compute(&self, tt_centuries: f64) -> crate::AstroResult<PrecessionResult> {
        let bias_matrix = self.frame_bias_matrix_iau2000();

        let precession_matrix = self.precession_matrix_iau2000(tt_centuries)?;

        let bias_precession_matrix = precession_matrix.multiply(&bias_matrix);

        Ok(PrecessionResult {
            bias_matrix,
            precession_matrix,
            bias_precession_matrix,
        })
    }

    /// Computes the frame bias matrix for IAU 2000.
    ///
    /// The frame bias accounts for the offset between the GCRS (defined by
    /// extragalactic radio sources) and the mean dynamical frame of J2000.0.
    /// This is a small, constant rotation of order tens of milliarcseconds.
    ///
    /// The rotation sequence is R_z(dRA) R_y(dLon * sin(eps0)) R_x(-dObl),
    /// where the bias parameters are defined at the J2000.0 epoch.
    fn frame_bias_matrix_iau2000(&self) -> RotationMatrix3 {
        // Mean obliquity at J2000.0 (arcseconds converted to radians)
        const EPS0: f64 = 84381.448 * ARCSEC_TO_RAD;

        // Frame bias parameters from IERS Conventions (2010), Table 5.1
        // All values in milliarcseconds, converted to radians
        const FRAME_BIAS_LONGITUDE: f64 = -0.041775 / 1000.0 * ARCSEC_TO_RAD;
        const FRAME_BIAS_OBLIQUITY: f64 = -0.0068192 / 1000.0 * ARCSEC_TO_RAD;
        const FRAME_BIAS_RA_OFFSET: f64 = -0.0146 / 1000.0 * ARCSEC_TO_RAD;

        let mut rb = RotationMatrix3::identity();
        rb.rotate_z(FRAME_BIAS_RA_OFFSET);
        rb.rotate_y(FRAME_BIAS_LONGITUDE * EPS0.sin());
        rb.rotate_x(-FRAME_BIAS_OBLIQUITY);

        rb
    }

    /// Computes the precession matrix from mean J2000.0 to mean of date.
    ///
    /// Uses the Lieske et al. (1977) precession angles with IAU 2000 corrections.
    /// The angles psi_A (precession in longitude), omega_A (obliquity of the
    /// ecliptic), and chi_A (planetary precession) are computed as polynomials
    /// in time, then small corrections are applied for the IAU 2000 model.
    ///
    /// The rotation sequence R_x(eps0) R_z(-psi_A) R_x(-omega_A) R_z(chi_A)
    /// transforms vectors from the mean equator/equinox of J2000.0 to the
    /// mean equator/equinox of date.
    fn precession_matrix_iau2000(&self, tt_centuries: f64) -> crate::AstroResult<RotationMatrix3> {
        // Mean obliquity at J2000.0
        const EPS0: f64 = 84381.448 * ARCSEC_TO_RAD;

        let t = tt_centuries;

        // Lieske (1979) precession angles (arcseconds, converted to radians)
        // psi_A: precession in longitude
        // omega_A: obliquity of the ecliptic of date
        // chi_A: planetary precession
        let psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * ARCSEC_TO_RAD;
        let oma77 = EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * ARCSEC_TO_RAD;
        let chia = (10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * ARCSEC_TO_RAD;

        // IAU 2000 precession-rate corrections (milliarcseconds/century)
        // These adjust the Lieske angles based on VLBI observations
        let dpsipr = -0.29965 / 1000.0 * ARCSEC_TO_RAD * t;
        let depspr = -0.02524 / 1000.0 * ARCSEC_TO_RAD * t;

        // Apply corrections to Lieske angles
        let psia = psia77 + dpsipr;
        let oma = oma77 + depspr;

        // Build precession matrix using four rotations
        let mut rp = RotationMatrix3::identity();

        rp.rotate_x(EPS0); // Rotate to ecliptic of J2000.0
        rp.rotate_z(-psia); // Precess along ecliptic
        rp.rotate_x(-oma); // Rotate to equator of date
        rp.rotate_z(chia); // Account for planetary precession

        Ok(rp)
    }
}

impl Default for PrecessionIAU2000 {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::assert_ulp_le;

    #[test]
    fn test_new_and_default() {
        let p1 = PrecessionIAU2000::new();
        let p2 = PrecessionIAU2000::default();
        let r1 = p1.compute(0.0).unwrap();
        let r2 = p2.compute(0.0).unwrap();
        assert_eq!(r1.bias_matrix, r2.bias_matrix);
    }

    #[test]
    fn test_compute_returns_rotation_matrices() {
        let p = PrecessionIAU2000::new();
        let result = p.compute(0.5).unwrap();
        assert!(result.bias_matrix.is_rotation_matrix(1e-14));
        assert!(result.precession_matrix.is_rotation_matrix(1e-14));
        assert!(result.bias_precession_matrix.is_rotation_matrix(1e-14));
    }

    #[test]
    fn test_bias_matrix_is_constant() {
        let p = PrecessionIAU2000::new();
        let r1 = p.compute(0.0).unwrap();
        let r2 = p.compute(1.0).unwrap();
        assert_eq!(r1.bias_matrix, r2.bias_matrix);
    }

    #[test]
    fn test_precession_at_j2000_is_identity() {
        let p = PrecessionIAU2000::new();
        let result = p.compute(0.0).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_ulp_le(
                    result.precession_matrix.get(i, j),
                    expected,
                    1,
                    &format!("precession[{},{}] at t=0", i, j),
                );
            }
        }
    }

    #[test]
    fn test_precession_changes_with_time() {
        let p = PrecessionIAU2000::new();
        let r0 = p.compute(0.0).unwrap();
        let r1 = p.compute(1.0).unwrap();
        // Verify precession matrix is NOT identity after 1 century
        // At least one element must differ significantly from identity
        let mut differs = false;
        for i in 0..3 {
            for j in 0..3 {
                let identity_val = if i == j { 1.0 } else { 0.0 };
                if (r1.precession_matrix.get(i, j) - identity_val).abs() > 1e-4 {
                    differs = true;
                }
            }
        }
        assert!(
            differs,
            "precession matrix should differ from identity after 1 century"
        );
        // Also verify the two matrices are not equal
        assert_ne!(r0.precession_matrix, r1.precession_matrix);
    }

    #[test]
    fn test_combined_matrix_consistency() {
        let p = PrecessionIAU2000::new();
        let result = p.compute(0.5).unwrap();
        let expected = result.precession_matrix.multiply(&result.bias_matrix);
        for i in 0..3 {
            for j in 0..3 {
                assert_ulp_le(
                    result.bias_precession_matrix.get(i, j),
                    expected.get(i, j),
                    1,
                    &format!("bias_precession[{},{}]", i, j),
                );
            }
        }
    }
}
