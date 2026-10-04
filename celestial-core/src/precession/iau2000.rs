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

// Mean obliquity at J2000.0, shared by the frame bias and the precession rotations.
const EPS0: f64 = 84381.448 * ARCSEC_TO_RAD;

/// IAU 2000 precession computation.
///
/// Computes the precession and frame bias matrices for transforming
/// coordinates between J2000.0 and the mean equator and equinox of date.
///
/// # Example
///
/// ```
/// use celestial_core::constants::J2000_JD;
/// use celestial_core::precession::PrecessionIAU2000;
///
/// let precession = PrecessionIAU2000::new();
///
/// // Compute for 18262.5 days (50 years) after J2000.0
/// let result = precession.compute(J2000_JD, 18262.5).unwrap();
///
/// // Access individual matrices
/// let _bias = &result.bias_matrix;           // GCRS to mean J2000.0
/// let _prec = &result.precession_matrix;     // Mean J2000.0 to mean of date
/// let _combined = &result.bias_precession_matrix; // GCRS to mean of date
/// ```
#[derive(Debug, Clone, Copy, Default)]
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
    /// * `date1`, `date2` - TT as a two-part Julian Date, `jd = date1 + date2`.
    ///   A common convention is `date1 = J2000_JD` and `date2 = days since J2000.0`.
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
    pub fn compute(&self, date1: f64, date2: f64) -> crate::errors::AstroResult<PrecessionResult> {
        let t = crate::utils::checked_jd_to_centuries(date1, date2)?;
        let bias_matrix = self.frame_bias_matrix_iau2000();

        let precession_matrix = self.precession_matrix_iau2000(t);

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
    fn frame_bias_matrix_iau2000(&self) -> RotationMatrix3 {
        // Frame bias parameters from IERS Conventions (2010), Table 5.1
        const FRAME_BIAS_LONGITUDE: f64 = -0.041775 * ARCSEC_TO_RAD;
        const FRAME_BIAS_OBLIQUITY: f64 = -0.0068192 * ARCSEC_TO_RAD;
        const FRAME_BIAS_RA_OFFSET: f64 = -0.0146 * ARCSEC_TO_RAD;

        let mut rb = RotationMatrix3::identity();
        rb.rotate_z(FRAME_BIAS_RA_OFFSET);
        rb.rotate_y(FRAME_BIAS_LONGITUDE * libm::sin(EPS0));
        rb.rotate_x(-FRAME_BIAS_OBLIQUITY);

        rb
    }

    /// Computes the precession matrix from mean J2000.0 to mean of date.
    ///
    /// Uses the Lieske et al. (1977) precession angles with IAU 2000 corrections:
    /// psi_A (precession in longitude), omega_A (obliquity of the ecliptic), and
    /// chi_A (planetary precession).
    fn precession_matrix_iau2000(&self, tt_centuries: f64) -> RotationMatrix3 {
        let (psia, oma, chia) = precession_angles_iau2000(tt_centuries);
        let mut rp = RotationMatrix3::identity();
        rp.rotate_x(EPS0);
        rp.rotate_z(-psia);
        rp.rotate_x(-oma);
        rp.rotate_z(chia);
        rp
    }
}

// Lieske et al. (1977) psi_A, omega_A and chi_A in radians, with the IAU 2000
// precession-rate corrections applied to psi_A and omega_A.
fn precession_angles_iau2000(t: f64) -> (f64, f64, f64) {
    let psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * ARCSEC_TO_RAD;
    let oma77 = EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * ARCSEC_TO_RAD;
    let chia = (10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * ARCSEC_TO_RAD;

    let dpsipr = -0.29965 * ARCSEC_TO_RAD * t;
    let depspr = -0.02524 * ARCSEC_TO_RAD * t;

    (psia77 + dpsipr, oma77 + depspr, chia)
}
