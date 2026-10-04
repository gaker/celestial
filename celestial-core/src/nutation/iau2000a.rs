//! IAU 2000A nutation model.
//!
//! Implements the IAU 2000A nutation model as defined by the International
//! Astronomical Union. This model computes the nutation in longitude (delta_psi)
//! and nutation in obliquity (delta_eps) for a given epoch.
//!
//! ## Model Specification
//!
//! IAU 2000A is a nutation model based on the MHB2000 (Mathews,
//! Herring, Buffett 2002) rigid-Earth series with:
//!
//! - **678 lunisolar terms**: Trigonometric series based on 5 fundamental arguments
//!   (Moon's mean anomaly, Sun's mean anomaly, Moon's argument of latitude,
//!   mean elongation of Moon from Sun, longitude of Moon's ascending node)
//! - **687 planetary terms**: Additional terms involving planetary mean longitudes
//!   (Mercury through Neptune) and general precession in longitude
//!
//! ## Precision
//!
//! - Formal precision: ~0.1 microarcsecond (μas) for epochs near J2000.0
//! - Accuracy degrades for epochs far from J2000.0 due to polynomial approximations
//! - Suitable for applications requiring sub-milliarcsecond precision
//!
//! ## Reference
//!
//! - IERS Conventions (2010), Chapter 5
//! - Mathews, Herring & Buffett (2002), J. Geophys. Res. 107, B4

use super::fundamental_args::{IERS2010FundamentalArgs, MHB2000FundamentalArgs};
use super::lunisolar_terms::LUNISOLAR_TERMS_F64;
use super::planetary_terms::PLANETARY_TERMS_F64;
use super::types::NutationResult;
use crate::constants::{TENTH_MICROARCSEC_TO_RAD, TWOPI};
use crate::errors::AstroResult;
use crate::math::fmod;

/// IAU 2000A nutation calculator.
///
/// Computes nutation angles using the full IAU 2000A model with 678 lunisolar
/// terms and 687 planetary terms. The computation follows the MHB2000 formulation
/// with coefficients expressed in microarcseconds.
///
/// # Example
///
/// ```
/// use celestial_core::nutation::NutationIAU2000A;
///
/// let nut = NutationIAU2000A::new();
///
/// // J2000.0 epoch (two-part JD for precision)
/// let jd1 = 2451545.0;
/// let jd2 = 0.0;
///
/// let result = nut.compute(jd1, jd2).unwrap();
/// // result.delta_psi: nutation in longitude (radians)
/// // result.delta_eps: nutation in obliquity (radians)
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub struct NutationIAU2000A;

impl NutationIAU2000A {
    /// Creates a new IAU 2000A nutation calculator.
    pub fn new() -> Self {
        Self
    }

    /// Computes nutation for the given epoch.
    ///
    /// Evaluates both lunisolar and planetary nutation series at the specified
    /// Julian Date, returning nutation in longitude (delta_psi) and obliquity
    /// (delta_eps) in radians.
    ///
    /// # Arguments
    ///
    /// * `jd1` - First part of two-part Julian Date (typically the integer day)
    /// * `jd2` - Second part of two-part Julian Date (typically the fractional day)
    ///
    /// The epoch is computed as `jd1 + jd2`. The two-part representation preserves
    /// precision when the epoch is far from J2000.0.
    ///
    /// # Returns
    ///
    /// [`NutationResult`] containing:
    /// - `delta_psi`: Nutation in longitude (radians)
    /// - `delta_eps`: Nutation in obliquity (radians)
    pub fn compute(&self, jd1: f64, jd2: f64) -> AstroResult<NutationResult> {
        let t = crate::utils::checked_jd_to_centuries(jd1, jd2)?;
        Ok(self.nutation_at(t))
    }

    pub(super) fn nutation_at(&self, t: f64) -> NutationResult {
        let lunisolar_args = [
            t.moon_mean_anomaly(),
            t.sun_mean_anomaly_mhb(),
            t.mean_argument_of_latitude(),
            t.mean_elongation_mhb(),
            t.moon_ascending_node_longitude(),
        ];
        let (delta_psi_ls, delta_eps_ls) = self.compute_lunisolar(&lunisolar_args, t);

        let (delta_psi_planetary, delta_eps_planetary) = self.compute_planetary(t);

        NutationResult {
            delta_psi: delta_psi_planetary + delta_psi_ls,
            delta_eps: delta_eps_planetary + delta_eps_ls,
        }
    }

    /// Computes the lunisolar nutation contribution.
    ///
    /// Evaluates 678 terms of the lunisolar nutation series. Each term is a
    /// trigonometric function of a linear combination of the five fundamental
    /// arguments of lunisolar motion.
    ///
    /// The series has the form:
    /// ```text
    /// delta_psi = sum_i (A_i + A'_i * t) * sin(arg_i) + A''_i * cos(arg_i)
    /// delta_eps = sum_i (B_i + B'_i * t) * cos(arg_i) + B''_i * sin(arg_i)
    /// ```
    ///
    /// where `arg_i = n_l * l + n_lp * l' + n_F * F + n_D * D + n_Om * Om` and
    /// coefficients are in microarcseconds.
    ///
    /// # Arguments
    ///
    /// * `args` - Five fundamental arguments in radians: \[l, l', F, D, Om\]
    ///   - l: Moon's mean anomaly
    ///   - l': Sun's mean anomaly
    ///   - F: Moon's argument of latitude
    ///   - D: Mean elongation of Moon from Sun
    ///   - Om: Longitude of Moon's ascending node
    /// * `t` - Julian centuries from J2000.0 (TT)
    ///
    /// # Returns
    ///
    /// Tuple of (delta_psi, delta_eps) in radians.
    pub fn compute_lunisolar(&self, args: &[f64; 5], t: f64) -> (f64, f64) {
        lunisolar_series(&LUNISOLAR_TERMS_F64, args, t)
    }

    /// Computes the planetary nutation contribution.
    ///
    /// Evaluates 687 terms of the planetary nutation series. Each term depends on
    /// the mean longitudes of the planets (Mercury through Neptune) plus the
    /// general precession in longitude.
    ///
    /// The planetary series is smaller in amplitude than the lunisolar series
    /// but essential for sub-milliarcsecond accuracy. The largest planetary terms
    /// arise from resonances between planetary and lunar orbital periods.
    ///
    /// # Arguments
    ///
    /// * `t` - Julian centuries from J2000.0 (TT)
    ///
    /// # Returns
    ///
    /// Tuple of (delta_psi, delta_eps) in radians.
    pub fn compute_planetary(&self, t: f64) -> (f64, f64) {
        let args = planetary_args(t);
        let mut dpsi = 0.0;
        let mut deps = 0.0;
        for [multipliers @ .., sp, cp, se, ce] in PLANETARY_TERMS_F64.iter().rev() {
            let (sarg, carg) = libm::sincos(series_argument(multipliers, &args));
            dpsi += sp * sarg + cp * carg;
            deps += se * sarg + ce * carg;
        }
        (
            dpsi * TENTH_MICROARCSEC_TO_RAD,
            deps * TENTH_MICROARCSEC_TO_RAD,
        )
    }
}

// Shared with IAU 2000B, which sums only the leading terms of the same table.
pub(super) fn lunisolar_series(terms: &[[f64; 11]], args: &[f64; 5], t: f64) -> (f64, f64) {
    let mut dpsi = 0.0;
    let mut deps = 0.0;
    for [multipliers @ .., sp, spt, cp, ce, cet, se] in terms.iter().rev() {
        let (sarg, carg) = libm::sincos(series_argument(multipliers, args));
        dpsi += (sp + spt * t) * sarg + cp * carg;
        deps += (ce + cet * t) * carg + se * sarg;
    }
    (
        dpsi * TENTH_MICROARCSEC_TO_RAD,
        deps * TENTH_MICROARCSEC_TO_RAD,
    )
}

// Summed left to right in table order; the ERFA pins depend on that order.
fn series_argument<const N: usize>(multipliers: &[f64; N], args: &[f64; N]) -> f64 {
    let sum = multipliers
        .iter()
        .zip(args)
        .map(|(n, a)| n * a)
        .reduce(|sum, product| sum + product);
    fmod(sum.unwrap_or(0.0), TWOPI)
}

// The lunar arguments (l, F, D, Ω) are MHB2000's linear forms, not the IERS polynomials
// the lunisolar series uses. Order matches the planetary table's multiplier columns.
fn planetary_args(t: f64) -> [f64; 13] {
    [
        fmod(2.35555598 + 8328.6914269554 * t, TWOPI),
        fmod(1.627905234 + 8433.466158131 * t, TWOPI),
        fmod(5.198466741 + 7771.3771468121 * t, TWOPI),
        fmod(2.18243920 - 33.757045 * t, TWOPI),
        t.mercury_lng(),
        t.venus_lng(),
        t.earth_lng(),
        t.mars_lng(),
        t.jupiter_lng(),
        t.saturn_lng(),
        t.uranus_lng(),
        t.neptune_longitude_mhb(),
        t.precession(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    // (jd1, jd2, dpsi, deps) from ERFA nut00a, run with Rust libm for sin/cos/fmod.
    const ERFA_NUT00A: [(f64, f64, f64, f64); 6] = [
        (
            2400000.5,
            53736.0,
            -9.630909107116424e-6,
            4.0632391740016646e-5,
        ),
        (
            2451545.0,
            0.0,
            -6.754422426417298e-5,
            -2.7970831192374137e-5,
        ),
        (
            2400000.5,
            60000.0,
            -4.496338070306366e-5,
            3.753546622895506e-5,
        ),
        (
            2451545.0,
            -219150.0,
            4.389836152812474e-5,
            -4.124743149757404e-5,
        ),
        (
            2451545.0,
            219150.0,
            -5.0887129295862736e-5,
            3.108429897083574e-5,
        ),
        // One of the rare dates where the order of operations in the general
        // precession argument changes the result.
        (
            2451545.0,
            -630023.9141714298,
            -5.351671754571772e-6,
            4.295092325358382e-5,
        ),
    ];

    #[test]
    fn test_matches_erfa_nut00a() {
        for (jd1, jd2, dpsi, deps) in ERFA_NUT00A {
            let result = NutationIAU2000A::new().compute(jd1, jd2).unwrap();
            let got = (result.delta_psi, result.delta_eps);
            assert_eq!(got, (dpsi, deps), "{jd1} + {jd2}");
        }
    }

    #[test]
    fn test_rejects_epoch_outside_model_range() {
        let model = NutationIAU2000A::new();
        assert!(model.compute(f64::NAN, 0.0).is_err());
        assert!(model.compute(2451545.0, f64::INFINITY).is_err());
        assert!(model.compute(f64::MAX, f64::MAX).is_err());
        assert!(model.compute(2451545.0, 1e300).is_err());
        assert!(model.compute(2451545.0, -730501.0).is_err());
    }
}
