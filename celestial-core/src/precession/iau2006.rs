//! IAU 2006 precession model.
//!
//! This module implements the IAU 2006 precession model, which supersedes the
//! IAU 2000 precession. The key improvement is the adoption of a new precession
//! rate in longitude (the "P03" solution) derived by Capitaine et al. (2003),
//! which corrected the IAU 2000 value by approximately -0.3 milliarcseconds per
//! year. This correction addressed a long-standing discrepancy with VLBI
//! observations.
//!
//! The model uses the Fukushima-Williams parameterization (gamb, phib, psib,
//! epsa), which provides a more stable numerical formulation than the classical
//! Euler angles. These four angles define the orientation of the mean equator
//! and equinox of date relative to J2000.0.
//!
//! # Background
//!
//! The Earth's rotation axis precesses around the ecliptic pole with a period
//! of approximately 26,000 years. The IAU 2006 precession model describes this
//! motion through polynomial expressions in Julian centuries from J2000.0,
//! derived from the most accurate celestial reference frame observations
//! available at the time.
//!
//! The Fukushima-Williams angles are:
//! - **gamb**: Frame bias in the longitude direction
//! - **phib**: Frame bias in the obliquity direction
//! - **psib**: Precession in longitude (luni-solar + planetary)
//! - **epsa**: Mean obliquity of the ecliptic
//!
//! # References
//!
//! - IERS Conventions (2010), Chapter 5
//! - Capitaine, N., Wallace, P.T., & Chapront, J. (2003), A&A 412, 567-586
//! - Hilton, J.L., et al. (2006), Celest. Mech. Dyn. Astron. 94, 351-367

use super::types::PrecessionResult;
use crate::constants::ARCSEC_TO_RAD;
use crate::constants::J2000_JD;
use crate::matrix::RotationMatrix3;

/// IAU 2006 precession model using the Fukushima-Williams parameterization.
///
/// This struct provides methods to compute precession matrices that transform
/// coordinates between the J2000.0 reference frame and the mean equator and
/// equinox of a given date.
///
/// # Example
///
/// ```
/// use celestial_core::precession::PrecessionIAU2006;
/// use celestial_core::constants::J2000_JD;
///
/// let precession = PrecessionIAU2006::new();
/// let result = precession.compute(J2000_JD, 3652.5).unwrap(); // ~10 years after J2000
///
/// // The result contains three matrices:
/// // - bias_matrix: transforms from GCRS to mean J2000.0 frame
/// // - precession_matrix: transforms from mean J2000.0 to mean of date
/// // - bias_precession_matrix: combined transform from GCRS to mean of date
/// ```
pub struct PrecessionIAU2006;

impl Default for PrecessionIAU2006 {
    fn default() -> Self {
        Self::new()
    }
}

impl PrecessionIAU2006 {
    /// Creates a new IAU 2006 precession model instance.
    pub fn new() -> Self {
        Self
    }

    /// Computes precession matrices for the given Julian Date.
    ///
    /// The date is specified as a two-part Julian Date for maximum precision:
    /// `jd = date1 + date2`. A common convention is `date1 = J2000_JD` and
    /// `date2 = days since J2000.0`.
    ///
    /// # Returns
    ///
    /// A [`PrecessionResult`] containing:
    /// - `bias_matrix`: The frame bias matrix at J2000.0
    /// - `precession_matrix`: Pure precession from mean J2000.0 to mean of date
    /// - `bias_precession_matrix`: Combined bias + precession matrix
    ///
    /// # Algorithm
    ///
    /// 1. Compute Fukushima-Williams angles at the target date
    /// 2. Build the combined bias-precession matrix from these angles
    /// 3. Compute the J2000.0 bias matrix (FW angles at t=0)
    /// 4. Extract pure precession by removing the bias: P = BP * B^T
    pub fn compute(&self, date1: f64, date2: f64) -> crate::AstroResult<PrecessionResult> {
        let t = ((date1 - J2000_JD) + date2) / crate::constants::DAYS_PER_JULIAN_CENTURY;

        let (gamb, phib, psib, _epsa_unused) = self.fukushima_williams_angles(t);
        let epsa = self.obliquity_from_t(t);

        let bias_precession_matrix = self.fw_angles_to_matrix(gamb, phib, psib, epsa);

        let (gamb_j2000, phib_j2000, psib_j2000, _epsa_j2000_unused) =
            self.fukushima_williams_angles(0.0);
        let epsa_j2000 = self.obliquity_from_t(0.0);
        let bias_matrix = self.fw_angles_to_matrix(gamb_j2000, phib_j2000, psib_j2000, epsa_j2000);

        let bias_inverse = bias_matrix.transpose();
        let precession_matrix = bias_precession_matrix.multiply(&bias_inverse);

        Ok(PrecessionResult {
            bias_matrix,
            precession_matrix,
            bias_precession_matrix,
        })
    }

    /// Computes the four Fukushima-Williams angles for a given time.
    ///
    /// These angles parameterize the orientation of the mean equator and equinox
    /// of date relative to the GCRS. The polynomial coefficients are from
    /// Hilton et al. (2006) and implement the IAU 2006 precession.
    ///
    /// # Arguments
    ///
    /// * `t` - Julian centuries of TDB (or TT, the difference is negligible)
    ///   since J2000.0
    ///
    /// # Returns
    ///
    /// A tuple `(gamb, phib, psib, epsa)` in radians:
    /// - `gamb`: F-W angle gamma_bar (related to frame bias in longitude)
    /// - `phib`: F-W angle phi_bar (related to frame bias in obliquity)
    /// - `psib`: F-W angle psi_bar (precession in longitude)
    /// - `epsa`: Mean obliquity of the ecliptic
    ///
    /// # Note
    ///
    /// The polynomials use Horner's method for numerical stability.
    pub fn fukushima_williams_angles(&self, t: f64) -> (f64, f64, f64, f64) {
        let gamb = (-0.052928
            + (10.556378
                + (0.4932044 + (-0.00031238 + (-0.000002788 + (0.0000000260) * t) * t) * t) * t)
                * t)
            * ARCSEC_TO_RAD;

        let phib = (84381.412819
            + (-46.811016
                + (0.0511268 + (0.00053289 + (-0.000000440 + (-0.0000000176) * t) * t) * t) * t)
                * t)
            * ARCSEC_TO_RAD;

        let psib = (-0.041775
            + (5038.481484
                + (1.5584175 + (-0.00018522 + (-0.000026452 + (-0.0000000148) * t) * t) * t) * t)
                * t)
            * ARCSEC_TO_RAD;

        let epsa = (84381.406
            + (-46.836769
                + (-0.0001831 + (0.00200340 + (-0.000000576 + (-0.0000000434) * t) * t) * t) * t)
                * t)
            * ARCSEC_TO_RAD;

        (gamb, phib, psib, epsa)
    }

    /// Computes the IAU 2006 mean obliquity of the ecliptic.
    ///
    /// This is the angle between the mean equator and the ecliptic plane,
    /// accounting for the secular decrease due to planetary perturbations.
    /// The polynomial is from Capitaine et al. (2003), Eq. 37.
    ///
    /// At J2000.0, the mean obliquity is approximately 84381.406 arcseconds
    /// (about 23.4 degrees), decreasing by roughly 47 arcseconds per century.
    fn obliquity_from_t(&self, t: f64) -> f64 {
        (84381.406
            + (-46.836769
                + (-0.0001831 + (0.00200340 + (-0.000000576 + (-0.0000000434) * t) * t) * t) * t)
                * t)
            * ARCSEC_TO_RAD
    }

    /// Constructs a rotation matrix from the four Fukushima-Williams angles.
    ///
    /// The matrix is built as a sequence of four rotations:
    /// 1. Rz(gamb) - rotation about the z-axis by gamma_bar
    /// 2. Rx(phib) - rotation about the new x-axis by phi_bar
    /// 3. Rz(-psib) - rotation about the new z-axis by -psi_bar
    /// 4. Rx(-epsa) - rotation about the new x-axis by -epsilon_A
    ///
    /// This sequence transforms vectors from the GCRS (or the mean equator
    /// and equinox of J2000.0 when using bias-free angles) to the mean
    /// equator and equinox of date.
    ///
    /// # Arguments
    ///
    /// * `gamb` - F-W angle gamma_bar in radians
    /// * `phib` - F-W angle phi_bar in radians
    /// * `psib` - F-W angle psi_bar in radians
    /// * `epsa` - Mean obliquity epsilon_A in radians
    pub fn fw_angles_to_matrix(
        &self,
        gamb: f64,
        phib: f64,
        psib: f64,
        epsa: f64,
    ) -> RotationMatrix3 {
        let mut matrix = RotationMatrix3::identity();

        matrix.rotate_z(gamb);

        matrix.rotate_x(phib);

        matrix.rotate_z(-psib);

        matrix.rotate_x(-epsa);

        matrix
    }

    /// Computes the combined nutation-precession-bias (NPB) matrix for IAU 2006/2000A.
    ///
    /// This method combines the IAU 2006 precession with nutation corrections
    /// (typically from the IAU 2000A nutation model) to produce a single rotation
    /// matrix that transforms GCRS coordinates to the true equator and equinox
    /// of date.
    ///
    /// The nutation corrections `dpsi` (nutation in longitude) and `deps`
    /// (nutation in obliquity) are added to the precession angles psi_bar and
    /// epsilon_A respectively before constructing the F-W matrix.
    ///
    /// # Arguments
    ///
    /// * `tt_centuries` - Julian centuries of TT since J2000.0
    /// * `dpsi` - Nutation in longitude in radians
    /// * `deps` - Nutation in obliquity in radians
    ///
    /// # Returns
    ///
    /// A rotation matrix that transforms from GCRS to the true equator and
    /// equinox of date, incorporating frame bias, precession, and nutation.
    pub fn npb_matrix_iau2006a(&self, tt_centuries: f64, dpsi: f64, deps: f64) -> RotationMatrix3 {
        let (gamb, phib, psib, epsa) = self.fukushima_williams_angles(tt_centuries);

        self.fw_angles_to_matrix(gamb, phib, psib + dpsi, epsa + deps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::assert_ulp_le;

    #[test]
    fn test_new_and_default() {
        let p1 = PrecessionIAU2006::new();
        let p2 = PrecessionIAU2006::default();
        let r1 = p1.compute(J2000_JD, 0.0).unwrap();
        let r2 = p2.compute(J2000_JD, 0.0).unwrap();
        assert_eq!(r1.bias_matrix, r2.bias_matrix);
    }

    #[test]
    fn test_compute_returns_rotation_matrices() {
        let p = PrecessionIAU2006::new();
        let result = p
            .compute(J2000_JD, 0.5 * crate::constants::DAYS_PER_JULIAN_CENTURY)
            .unwrap();
        assert!(result.bias_matrix.is_rotation_matrix(1e-14));
        assert!(result.precession_matrix.is_rotation_matrix(1e-14));
        assert!(result.bias_precession_matrix.is_rotation_matrix(1e-14));
    }

    #[test]
    fn test_bias_matrix_is_constant() {
        let p = PrecessionIAU2006::new();
        let r1 = p.compute(J2000_JD, 0.0).unwrap();
        let r2 = p
            .compute(J2000_JD, crate::constants::DAYS_PER_JULIAN_CENTURY)
            .unwrap();
        assert_eq!(r1.bias_matrix, r2.bias_matrix);
    }

    #[test]
    fn test_precession_at_j2000_is_identity() {
        let p = PrecessionIAU2006::new();
        let result = p.compute(J2000_JD, 0.0).unwrap();
        // Precession = bias_precession * bias_inverse, so numerical noise accumulates
        // Diagonal elements: check ULP against 1.0
        // Off-diagonal elements: check they're tiny (numerical noise ~1e-22)
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_ulp_le(
                        result.precession_matrix.get(i, j),
                        1.0,
                        4,
                        &format!("precession[{},{}] at t=0", i, j),
                    );
                } else {
                    assert!(
                        result.precession_matrix.get(i, j).abs() < 1e-15,
                        "precession[{},{}] at t=0 should be ~0, got {}",
                        i,
                        j,
                        result.precession_matrix.get(i, j)
                    );
                }
            }
        }
    }

    #[test]
    fn test_precession_changes_with_time() {
        let p = PrecessionIAU2006::new();
        let r0 = p.compute(J2000_JD, 0.0).unwrap();
        let r1 = p
            .compute(J2000_JD, crate::constants::DAYS_PER_JULIAN_CENTURY)
            .unwrap();
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
        assert_ne!(r0.precession_matrix, r1.precession_matrix);
    }

    #[test]
    fn test_combined_matrix_consistency() {
        let p = PrecessionIAU2006::new();
        let result = p
            .compute(J2000_JD, 0.5 * crate::constants::DAYS_PER_JULIAN_CENTURY)
            .unwrap();
        let bias_inverse = result.bias_matrix.transpose();
        let expected = result.bias_precession_matrix.multiply(&bias_inverse);
        for i in 0..3 {
            for j in 0..3 {
                assert_ulp_le(
                    result.precession_matrix.get(i, j),
                    expected.get(i, j),
                    2,
                    &format!("precession[{},{}]", i, j),
                );
            }
        }
    }

    #[test]
    fn test_fukushima_williams_angles_at_j2000() {
        let p = PrecessionIAU2006::new();
        let (gamb, phib, psib, epsa) = p.fukushima_williams_angles(0.0);
        assert_ulp_le(gamb, -0.052928 * ARCSEC_TO_RAD, 1, "gamb at t=0");
        assert_ulp_le(phib, 84381.412819 * ARCSEC_TO_RAD, 1, "phib at t=0");
        assert_ulp_le(psib, -0.041775 * ARCSEC_TO_RAD, 1, "psib at t=0");
        assert_ulp_le(epsa, 84381.406 * ARCSEC_TO_RAD, 1, "epsa at t=0");
    }

    #[test]
    fn test_fukushima_williams_angles_change_with_time() {
        let p = PrecessionIAU2006::new();
        let (gamb0, phib0, psib0, epsa0) = p.fukushima_williams_angles(0.0);
        let (gamb1, phib1, psib1, epsa1) = p.fukushima_williams_angles(1.0);
        assert_ne!(gamb0, gamb1);
        assert_ne!(phib0, phib1);
        assert_ne!(psib0, psib1);
        assert_ne!(epsa0, epsa1);
    }

    #[test]
    fn test_fw_angles_to_matrix_returns_rotation() {
        let p = PrecessionIAU2006::new();
        let (gamb, phib, psib, epsa) = p.fukushima_williams_angles(0.5);
        let matrix = p.fw_angles_to_matrix(gamb, phib, psib, epsa);
        assert!(matrix.is_rotation_matrix(1e-14));
    }

    #[test]
    fn test_npb_matrix_returns_rotation() {
        let p = PrecessionIAU2006::new();
        let matrix = p.npb_matrix_iau2006a(0.5, 0.001 * ARCSEC_TO_RAD, 0.0005 * ARCSEC_TO_RAD);
        assert!(matrix.is_rotation_matrix(1e-14));
    }

    #[test]
    fn test_npb_matrix_with_zero_nutation() {
        let p = PrecessionIAU2006::new();
        let (gamb, phib, psib, epsa) = p.fukushima_williams_angles(0.5);
        let fw_matrix = p.fw_angles_to_matrix(gamb, phib, psib, epsa);
        let npb_matrix = p.npb_matrix_iau2006a(0.5, 0.0, 0.0);
        for i in 0..3 {
            for j in 0..3 {
                assert_ulp_le(
                    npb_matrix.get(i, j),
                    fw_matrix.get(i, j),
                    1,
                    &format!("npb[{},{}] with zero nutation", i, j),
                );
            }
        }
    }

    #[test]
    fn test_two_part_date_equivalence() {
        let p = PrecessionIAU2006::new();
        let r1 = p.compute(J2000_JD, 1000.0).unwrap();
        let r2 = p.compute(J2000_JD + 500.0, 500.0).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert_ulp_le(
                    r1.bias_precession_matrix.get(i, j),
                    r2.bias_precession_matrix.get(i, j),
                    1,
                    &format!("two-part date[{},{}]", i, j),
                );
            }
        }
    }
}
