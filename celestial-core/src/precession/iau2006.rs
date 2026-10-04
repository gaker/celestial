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
//! - **gamb**: GCRS right ascension of the intersection of the ecliptic of date
//!   with the GCRS equator
//! - **phib**: Obliquity of the ecliptic of date on the GCRS equator
//! - **psib**: Precession angle along the ecliptic of date
//! - **epsa**: Mean obliquity of the ecliptic of date
//!
//! # References
//!
//! - IERS Conventions (2010), Chapter 5
//! - Capitaine, N., Wallace, P.T., & Chapront, J. (2003), A&A 412, 567-586
//! - Hilton, J.L., et al. (2006), Celest. Mech. Dyn. Astron. 94, 351-367

use super::types::{FukushimaWilliamsAngles, PrecessionResult};
use crate::constants::ARCSEC_TO_RAD;
use crate::math::polynomial;
use crate::matrix::RotationMatrix3;
use crate::obliquity::iau_2006_obliquity_at;

// Fukushima-Williams angle polynomials in arcseconds, constant term first.
const GAMMA_BAR_ARCSEC: [f64; 6] = [
    -0.052928,
    10.556378,
    0.4932044,
    -0.00031238,
    -0.000002788,
    0.0000000260,
];
const PHI_BAR_ARCSEC: [f64; 6] = [
    84381.412819,
    -46.811016,
    0.0511268,
    0.00053289,
    -0.000000440,
    -0.0000000176,
];
const PSI_BAR_ARCSEC: [f64; 6] = [
    -0.041775,
    5038.481484,
    1.5584175,
    -0.00018522,
    -0.000026452,
    -0.0000000148,
];

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
#[derive(Debug, Clone, Copy, Default)]
pub struct PrecessionIAU2006;

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
    pub fn compute(&self, date1: f64, date2: f64) -> crate::errors::AstroResult<PrecessionResult> {
        let t = crate::utils::checked_jd_to_centuries(date1, date2)?;
        let bias_precession_matrix = self.bias_precession_at(t);
        let bias_matrix = self.bias_precession_at(0.0);
        let precession_matrix = bias_precession_matrix.multiply(&bias_matrix.transpose());

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
    /// * `t` - Julian centuries of TT since J2000.0
    ///
    /// # Returns
    ///
    /// [`FukushimaWilliamsAngles`] in radians:
    /// - `gamma_bar`: GCRS right ascension of the intersection of the ecliptic of
    ///   date with the GCRS equator
    /// - `phi_bar`: obliquity of the ecliptic of date on the GCRS equator
    /// - `psi_bar`: precession in longitude
    /// - `epsilon_a`: mean obliquity of the ecliptic
    ///
    /// # Note
    ///
    /// The polynomials use Horner's method for numerical stability.
    pub fn fukushima_williams_angles(&self, t: f64) -> FukushimaWilliamsAngles {
        FukushimaWilliamsAngles {
            gamma_bar: polynomial(&GAMMA_BAR_ARCSEC, t) * ARCSEC_TO_RAD,
            phi_bar: polynomial(&PHI_BAR_ARCSEC, t) * ARCSEC_TO_RAD,
            psi_bar: polynomial(&PSI_BAR_ARCSEC, t) * ARCSEC_TO_RAD,
            epsilon_a: iau_2006_obliquity_at(t),
        }
    }

    fn bias_precession_at(&self, t: f64) -> RotationMatrix3 {
        fw_matrix(&self.fukushima_williams_angles(t))
    }

    /// Computes the combined nutation-precession-bias (NPB) matrix for IAU 2006/2000A.
    ///
    /// This method combines the IAU 2006 precession with nutation corrections
    /// from the IAU 2006A nutation model
    /// ([`NutationIAU2006A`](crate::nutation::NutationIAU2006A)) to produce a single
    /// rotation matrix that transforms GCRS coordinates to the true equator and
    /// equinox of date.
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
        let fw = self.fukushima_williams_angles(tt_centuries);
        fw_matrix(&FukushimaWilliamsAngles {
            psi_bar: fw.psi_bar + dpsi,
            epsilon_a: fw.epsilon_a + deps,
            ..fw
        })
    }
}

// Maps GCRS (or mean J2000.0, for bias-free angles) to the mean equator and equinox of date.
pub(super) fn fw_matrix(fw: &FukushimaWilliamsAngles) -> RotationMatrix3 {
    let mut matrix = RotationMatrix3::identity();
    matrix.rotate_z(fw.gamma_bar);
    matrix.rotate_x(fw.phi_bar);
    matrix.rotate_z(-fw.psi_bar);
    matrix.rotate_x(-fw.epsilon_a);
    matrix
}
