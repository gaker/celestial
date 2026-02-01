//! IAU 2006A nutation model.
//!
//! This module implements the IAU 2006A nutation model, which combines the
//! IAU 2000A nutation series with corrections for compatibility with the
//! IAU 2006 precession model.
//!
//! The IAU 2000A nutation was originally developed alongside the IAU 2000
//! precession model. When the IAU adopted the improved IAU 2006 precession
//! in 2006 (Capitaine et al. 2003), small adjustments to the nutation
//! angles became necessary to maintain consistency. The IAU 2006A model
//! applies these adjustments via the frame bias correction factor J2.
//!
//! The correction is applied as:
//!
//! ```text
//! Δψ_2006A = Δψ_2000A × (1 + 0.4697×10⁻⁶ + fJ2)
//! Δε_2006A = Δε_2000A × (1 + fJ2)
//!
//! where fJ2 = -2.7774×10⁻⁶ × t
//! and t is Julian centuries from J2000.0 TT
//! ```
//!
//! The 0.4697 µas factor corrects for the change in the dynamical ellipticity
//! of the Earth between the IAU 2000 and IAU 2006 precession models.
//!
//! Reference: IERS Conventions (2010), Chapter 5, Section 5.5.4

use super::iau2000a::NutationIAU2000A;
use super::types::NutationResult;
use crate::errors::AstroResult;

/// IAU 2006A nutation calculator.
///
/// Wraps [`NutationIAU2000A`] and applies the J2 frame bias corrections
/// required for use with IAU 2006 precession. This is the recommended
/// nutation model for high-accuracy applications using IAU 2006 precession.
///
/// # Example
///
/// ```
/// use astro_core::nutation::NutationIAU2006A;
///
/// let nutation = NutationIAU2006A::new();
///
/// // Compute nutation at J2000.0
/// let result = nutation.compute(2451545.0, 0.0).unwrap();
///
/// // delta_psi and delta_eps are in radians
/// println!("Δψ = {} rad", result.delta_psi);
/// println!("Δε = {} rad", result.delta_eps);
/// ```
pub struct NutationIAU2006A {
    iau2000a: NutationIAU2000A,
}

impl Default for NutationIAU2006A {
    fn default() -> Self {
        Self::new()
    }
}

impl NutationIAU2006A {
    /// Creates a new IAU 2006A nutation calculator.
    pub fn new() -> Self {
        Self {
            iau2000a: NutationIAU2000A::new(),
        }
    }

    /// Computes nutation angles Δψ (nutation in longitude) and Δε (nutation in obliquity).
    ///
    /// The computation follows IERS Conventions (2010):
    /// 1. Compute IAU 2000A nutation angles (lunisolar + planetary terms)
    /// 2. Apply the J2 frame bias correction for IAU 2006 precession compatibility
    ///
    /// # Arguments
    ///
    /// * `jd1` - First part of two-part Julian Date (TT scale)
    /// * `jd2` - Second part of two-part Julian Date (TT scale)
    ///
    /// The two-part JD should satisfy `jd1 + jd2 = JD`. For best precision,
    /// set `jd1` to 2451545.0 (J2000.0) and `jd2` to the offset from J2000.
    ///
    /// # Returns
    ///
    /// [`NutationResult`] containing:
    /// - `delta_psi`: Nutation in longitude (radians)
    /// - `delta_eps`: Nutation in obliquity (radians)
    ///
    /// # Accuracy
    ///
    /// The IAU 2000A series includes 1365 terms (678 lunisolar + 687 planetary)
    /// providing sub-milliarcsecond accuracy. The J2 correction is at the
    /// microarcsecond level.
    pub fn compute(&self, jd1: f64, jd2: f64) -> AstroResult<NutationResult> {
        let t =
            ((jd1 - crate::constants::J2000_JD) + jd2) / crate::constants::DAYS_PER_JULIAN_CENTURY;

        // J2 frame bias correction factor
        let fj2 = -2.7774e-6 * t;

        let res = self.iau2000a.compute(jd1, jd2)?;
        let dp = res.delta_psi;
        let de = res.delta_eps;

        // Apply corrections: 0.4697e-6 is the fixed J2 rate correction,
        // fj2 is the time-dependent part
        Ok(NutationResult {
            delta_psi: dp + dp * (0.4697e-6 + fj2),
            delta_eps: de + de * fj2,
        })
    }
}
