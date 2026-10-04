//! IAU 2000B Nutation Model
//!
//! This module implements the IAU 2000B nutation model, a truncated version of the full
//! IAU 2000A model designed for applications where sub-milliarcsecond precision is not required.
//!
//! # Model Description
//!
//! IAU 2000B reduces computational cost by:
//! - Using only the first 77 lunisolar terms (out of 678 in IAU 2000A)
//! - Omitting the 687 planetary terms entirely
//! - Applying fixed bias corrections to approximate the omitted planetary effects
//!
//! The planetary bias corrections are:
//! - Longitude (Δψ): -0.135 milliarcseconds
//! - Obliquity (Δε): +0.388 milliarcseconds
//!
//! # Accuracy
//!
//! The IAU 2000B model achieves accuracy of approximately 1 milliarcsecond over the
//! period 1995-2050. This is sufficient for many practical applications including
//! amateur telescope pointing and general ephemeris work, but insufficient for
//! high-precision astrometry, VLBI, or pulsar timing.
//!
//! # Fundamental Arguments
//!
//! The model uses five Delaunay arguments computed from polynomial expressions
//! in Julian centuries from J2000.0 (TT):
//!
//! - `l` (mean anomaly of the Moon)
//! - `l'` (mean anomaly of the Sun)
//! - `F` (mean argument of latitude of the Moon)
//! - `D` (mean elongation of the Moon from the Sun)
//! - `Ω` (mean longitude of the Moon's ascending node)
//!
//! # References
//!
//! - McCarthy, D. D. & Luzum, B. J., "An Abridged Model of the Precession-Nutation
//!   of the Celestial Pole", Celestial Mechanics and Dynamical Astronomy, 2003
//! - IERS Conventions (2003), Chapter 5
//! - SOFA Library: `iauNut00b`

use super::iau2000a::lunisolar_series;
use super::lunisolar_terms::LUNISOLAR_TERMS_F64;
use super::types::NutationResult;
use crate::constants::{ARCSEC_TO_RAD, CIRCULAR_ARCSECONDS, MILLIARCSEC_TO_RAD};
use crate::errors::AstroResult;
use crate::math::fmod;

// The 2000B lunisolar series is the leading 77 terms of the 2000A table.
const LUNISOLAR_TERM_COUNT: usize = 77;
const PLANETARY_BIAS_LONGITUDE: f64 = -0.135 * MILLIARCSEC_TO_RAD;
const PLANETARY_BIAS_OBLIQUITY: f64 = 0.388 * MILLIARCSEC_TO_RAD;

/// IAU 2000B nutation calculator.
///
/// A simplified nutation model using 77 lunisolar terms plus fixed planetary bias
/// corrections. Provides ~1 mas accuracy, suitable for applications not requiring
/// the full precision of [`NutationIAU2000A`](super::NutationIAU2000A).
///
/// # Example
///
/// ```
/// use celestial_core::nutation::NutationIAU2000B;
///
/// let nut = NutationIAU2000B::new();
/// // Compute nutation for J2000.0 (two-part JD: 2451545.0 + 0.0)
/// let result = nut.compute(2451545.0, 0.0).unwrap();
///
/// // delta_psi and delta_eps are in radians
/// println!("Δψ = {} rad", result.delta_psi);
/// println!("Δε = {} rad", result.delta_eps);
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub struct NutationIAU2000B;

impl NutationIAU2000B {
    /// Creates a new IAU 2000B nutation calculator.
    pub fn new() -> Self {
        Self
    }

    /// Computes nutation angles for a given Julian Date.
    ///
    /// # Arguments
    ///
    /// * `jd1` - First part of two-part Julian Date (TT). Typically the integer part
    ///   or J2000 epoch (2451545.0).
    /// * `jd2` - Second part of two-part Julian Date (TT). Typically the fractional
    ///   part or offset from `jd1`.
    ///
    /// The two-part representation preserves precision. The split is arbitrary;
    /// `jd1 + jd2` must equal the desired Julian Date.
    ///
    /// # Returns
    ///
    /// Returns a [`NutationResult`] containing:
    /// - `delta_psi`: Nutation in longitude (radians)
    /// - `delta_eps`: Nutation in obliquity (radians)
    ///
    /// Both values are IAU 2000B approximations with ~1 mas accuracy.
    pub fn compute(&self, jd1: f64, jd2: f64) -> AstroResult<NutationResult> {
        let t = crate::utils::checked_jd_to_centuries(jd1, jd2)?;

        let terms = &LUNISOLAR_TERMS_F64[..LUNISOLAR_TERM_COUNT];
        let (delta_psi_ls, delta_eps_ls) = lunisolar_series(terms, &delaunay_args(t), t);
        let delta_psi = delta_psi_ls + PLANETARY_BIAS_LONGITUDE;
        let delta_eps = delta_eps_ls + PLANETARY_BIAS_OBLIQUITY;
        Ok(NutationResult {
            delta_psi,
            delta_eps,
        })
    }
}

// Linear forms of the Delaunay arguments (l, l', F, D, Ω), in arcseconds before the
// conversion to radians; 2000B drops the higher-order terms the 2000A arguments carry.
fn delaunay_args(t: f64) -> [f64; 5] {
    [
        fmod(485868.249036 + 1717915923.2178 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD,
        fmod(1287104.79305 + 129596581.0481 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD,
        fmod(335779.526232 + 1739527262.8478 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD,
        fmod(1072260.70369 + 1602961601.2090 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD,
        fmod(450160.398036 + -6962890.5431 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    // (jd1, jd2, dpsi, deps) from ERFA nut00b, run with Rust libm for sin/cos/fmod.
    const ERFA_NUT00B: [(f64, f64, f64, f64); 5] = [
        (
            2400000.5,
            53736.0,
            -9.632552291148318e-6,
            4.063197106621162e-5,
        ),
        (
            2451545.0,
            0.0,
            -6.754261253992235e-5,
            -2.7970923310985653e-5,
        ),
        (
            2400000.5,
            60000.0,
            -4.496465903077654e-5,
            3.753571642681161e-5,
        ),
        (
            2451545.0,
            -219150.0,
            4.379669022155589e-5,
            -4.1276856183100264e-5,
        ),
        (
            2451545.0,
            219150.0,
            -5.08020843173988e-5,
            3.112354967405106e-5,
        ),
    ];

    #[test]
    fn test_matches_erfa_nut00b() {
        for (jd1, jd2, dpsi, deps) in ERFA_NUT00B {
            let result = NutationIAU2000B::new().compute(jd1, jd2).unwrap();
            let got = (result.delta_psi, result.delta_eps);
            assert_eq!(got, (dpsi, deps), "{jd1} + {jd2}");
        }
    }

    #[test]
    fn test_rejects_epoch_outside_model_range() {
        let model = NutationIAU2000B::new();
        assert!(model.compute(f64::NAN, 0.0).is_err());
        assert!(model.compute(2451545.0, f64::INFINITY).is_err());
        assert!(model.compute(f64::MAX, f64::MAX).is_err());
        assert!(model.compute(2451545.0, 1e300).is_err());
        assert!(model.compute(2451545.0, -730501.0).is_err());
    }
}
