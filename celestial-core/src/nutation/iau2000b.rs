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
//! in Julian centuries from J2000.0 (TDB):
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

use super::lunisolar_terms::LUNISOLAR_TERMS;
use super::types::NutationResult;
use crate::constants::{
    ARCSEC_TO_RAD, CIRCULAR_ARCSECONDS, MICROARCSEC_TO_RAD, MILLIARCSEC_TO_RAD, TWOPI,
};
use crate::errors::AstroResult;
use crate::math::fmod;

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
pub struct NutationIAU2000B;

impl Default for NutationIAU2000B {
    fn default() -> Self {
        Self::new()
    }
}

impl NutationIAU2000B {
    /// Creates a new IAU 2000B nutation calculator.
    pub fn new() -> Self {
        Self
    }

    /// Computes nutation angles for a given Julian Date.
    ///
    /// # Arguments
    ///
    /// * `jd1` - First part of two-part Julian Date (TDB). Typically the integer part
    ///   or J2000 epoch (2451545.0).
    /// * `jd2` - Second part of two-part Julian Date (TDB). Typically the fractional
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
        let t = crate::utils::jd_to_centuries(jd1, jd2);

        let (delta_psi_ls, delta_eps_ls) = self.compute_lunisolar(t);

        const PLANETARY_BIAS_LONGITUDE: f64 = -0.135 * MILLIARCSEC_TO_RAD;
        const PLANETARY_BIAS_OBLIQUITY: f64 = 0.388 * MILLIARCSEC_TO_RAD;

        let delta_psi = delta_psi_ls + PLANETARY_BIAS_LONGITUDE;
        let delta_eps = delta_eps_ls + PLANETARY_BIAS_OBLIQUITY;
        Ok(NutationResult {
            delta_psi,
            delta_eps,
        })
    }

    /// Computes the lunisolar nutation contribution.
    ///
    /// Evaluates the first 77 terms of the lunisolar nutation series using
    /// the five Delaunay fundamental arguments. Each term contributes
    /// sine and cosine components to both longitude and obliquity.
    ///
    /// # Arguments
    ///
    /// * `t` - Julian centuries from J2000.0 (TDB)
    ///
    /// # Returns
    ///
    /// Tuple of (Δψ, Δε) in radians representing the lunisolar contribution
    /// to nutation in longitude and obliquity.
    fn compute_lunisolar(&self, t: f64) -> (f64, f64) {
        // Delaunay arguments (arcseconds, then converted to radians)
        // l: Mean anomaly of the Moon
        let el = fmod(485868.249036 + 1717915923.2178 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD;
        // l': Mean anomaly of the Sun
        let elp = fmod(1287104.79305 + 129596581.0481 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD;
        // F: Mean argument of latitude of the Moon
        let f = fmod(335779.526232 + 1739527262.8478 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD;
        // D: Mean elongation of the Moon from the Sun
        let d = fmod(1072260.70369 + 1602961601.2090 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD;
        // Ω: Mean longitude of the Moon's ascending node
        let om = fmod(450160.398036 + -6962890.5431 * t, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD;

        let mut dpsi = 0.0;
        let mut deps = 0.0;

        for &(nl, nlp, nf, nd, nom, sp, spt, cp, ce, cet, se) in
            LUNISOLAR_TERMS.iter().take(77).rev()
        {
            let arg = fmod(
                (nl as f64) * el
                    + (nlp as f64) * elp
                    + (nf as f64) * f
                    + (nd as f64) * d
                    + (nom as f64) * om,
                TWOPI,
            );

            let (sarg, carg) = libm::sincos(arg);

            dpsi += (sp + spt * t) * sarg + cp * carg;
            deps += (ce + cet * t) * carg + se * sarg;
        }

        (dpsi * MICROARCSEC_TO_RAD, deps * MICROARCSEC_TO_RAD)
    }
}
