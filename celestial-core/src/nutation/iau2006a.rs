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
//! applies these adjustments as scale factors on the IAU 2000A angles.
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
//! The constant 0.4697×10⁻⁶ factor on Δψ is the P03 precession adjustment: it
//! accounts for the change in the Earth's dynamical ellipticity implied by the
//! IAU 2006 precession rate. fJ2 accounts for the secular change in the Earth's
//! dynamical form factor J2.
//!
//! Reference: IERS Conventions (2010), Chapter 5, Section 5.5.4

use super::iau2000a::NutationIAU2000A;
use super::types::NutationResult;
use crate::errors::AstroResult;

/// IAU 2006A nutation calculator.
///
/// Wraps [`NutationIAU2000A`] and applies the P03 and secular-J2 adjustments
/// required for use with IAU 2006 precession. This is the recommended
/// nutation model for high-accuracy applications using IAU 2006 precession.
///
/// # Example
///
/// ```
/// use celestial_core::nutation::NutationIAU2006A;
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
#[derive(Debug, Clone, Copy, Default)]
pub struct NutationIAU2006A {
    iau2000a: NutationIAU2000A,
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
    /// 2. Apply the P03 and secular-J2 adjustments for IAU 2006 precession compatibility
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
        let t = crate::utils::checked_jd_to_centuries(jd1, jd2)?;
        let fj2 = -2.7774e-6 * t;

        let res = self.iau2000a.nutation_at(t);
        let dp = res.delta_psi;
        let de = res.delta_eps;

        // 0.4697e-6 is the P03 dynamical-ellipticity factor (Δψ only);
        // fj2 is the secular J2 change (both angles)
        Ok(NutationResult {
            delta_psi: dp + dp * (0.4697e-6 + fj2),
            delta_eps: de + de * fj2,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // (jd1, jd2, dpsi, deps) from ERFA nut06a, run with Rust libm for sin/cos/fmod.
    const ERFA_NUT06A: [(f64, f64, f64, f64); 6] = [
        (
            2400000.5,
            53736.0,
            -9.630912025821214e-6,
            4.063238496887236e-5,
        ),
        (
            2451545.0,
            0.0,
            -6.754425598969512e-5,
            -2.7970831192374137e-5,
        ),
        (
            2400000.5,
            60000.0,
            -4.4963372912472335e-5,
            3.753544209495469e-5,
        ),
        (
            2451545.0,
            -219150.0,
            4.3899113687041e-5,
            -4.124811886127149e-5,
        ),
        (
            2451545.0,
            219150.0,
            -5.088630519406993e-5,
            3.108378096964397e-5,
        ),
        // One of the rare dates where the order of operations in the general
        // precession argument changes the result.
        (
            2451545.0,
            -630023.9141714298,
            -5.351930654489544e-6,
            4.295298093314962e-5,
        ),
    ];

    #[test]
    fn test_matches_erfa_nut06a() {
        for (jd1, jd2, dpsi, deps) in ERFA_NUT06A {
            let result = NutationIAU2006A::new().compute(jd1, jd2).unwrap();
            let got = (result.delta_psi, result.delta_eps);
            assert_eq!(got, (dpsi, deps), "{jd1} + {jd2}");
        }
    }

    #[test]
    fn test_rejects_epoch_outside_model_range() {
        let model = NutationIAU2006A::new();
        assert!(model.compute(f64::NAN, 0.0).is_err());
        assert!(model.compute(2451545.0, f64::INFINITY).is_err());
        assert!(model.compute(f64::MAX, f64::MAX).is_err());
        assert!(model.compute(2451545.0, 1e300).is_err());
        assert!(model.compute(2451545.0, -730501.0).is_err());
    }
}
