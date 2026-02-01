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
use super::lunisolar_terms::LUNISOLAR_TERMS;
use super::planetary_terms::PLANETARY_TERMS;
use super::types::NutationResult;
use crate::constants::{MICROARCSEC_TO_RAD, TWOPI};
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
/// use celestial_core::nutation::iau2000a::NutationIAU2000A;
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
#[derive(Debug, Clone, Copy)]
pub struct NutationIAU2000A;

impl Default for NutationIAU2000A {
    fn default() -> Self {
        Self::new()
    }
}

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
        let t = crate::utils::jd_to_centuries(jd1, jd2);

        let lunisolar_args = [
            t.moon_mean_anomaly(),
            t.sun_mean_anomaly_mhb(),
            t.mean_argument_of_latitude(),
            t.mean_elongation_mhb(),
            t.moon_ascending_node_longitude(),
        ];
        let (delta_psi_ls, delta_eps_ls) = self.compute_lunisolar(&lunisolar_args, t);

        let (delta_psi_planetary, delta_eps_planetary) = self.compute_planetary(t);

        Ok(NutationResult {
            delta_psi: delta_psi_planetary + delta_psi_ls,
            delta_eps: delta_eps_planetary + delta_eps_ls,
        })
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
        let mut dpsi = 0.0;
        let mut deps = 0.0;

        for term in LUNISOLAR_TERMS.iter().rev() {
            let arg = fmod(
                (term.0 as f64) * args[0]
                    + (term.1 as f64) * args[1]
                    + (term.2 as f64) * args[2]
                    + (term.3 as f64) * args[3]
                    + (term.4 as f64) * args[4],
                TWOPI,
            );

            let (sarg, carg) = arg.sin_cos();

            dpsi += (term.5 + term.6 * t) * sarg + term.7 * carg;
            deps += (term.8 + term.9 * t) * carg + term.10 * sarg;
        }

        (dpsi * MICROARCSEC_TO_RAD, deps * MICROARCSEC_TO_RAD)
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
        let al = fmod(2.35555598 + 8328.6914269554 * t, TWOPI);
        let af = fmod(1.627905234 + 8433.466158131 * t, TWOPI);
        let ad = fmod(5.198466741 + 7771.3771468121 * t, TWOPI);
        let aom = fmod(2.18243920 - 33.757045 * t, TWOPI);
        let apa = t.precession();

        let alme = t.mercury_lng();
        let alve = t.venus_lng();
        let alea = t.earth_lng();
        let alma = t.mars_lng();
        let alju = t.jupiter_lng();
        let alsa = t.saturn_lng();
        let alur = t.uranus_lng();
        let alne = t.neptune_longitude_mhb();

        let mut dpsi = 0.0;
        let mut deps = 0.0;

        for &(nl, nf, nd, nom, nme, nve, nea, nma, nju, nsa, nur, nne, npa, sp, cp, se, ce) in
            PLANETARY_TERMS.iter().rev()
        {
            let arg = fmod(
                (nl as f64) * al
                    + (nf as f64) * af
                    + (nd as f64) * ad
                    + (nom as f64) * aom
                    + (nme as f64) * alme
                    + (nve as f64) * alve
                    + (nea as f64) * alea
                    + (nma as f64) * alma
                    + (nju as f64) * alju
                    + (nsa as f64) * alsa
                    + (nur as f64) * alur
                    + (nne as f64) * alne
                    + (npa as f64) * apa,
                TWOPI,
            );

            let (sarg, carg) = arg.sin_cos();

            dpsi += (sp as f64) * sarg + (cp as f64) * carg;
            deps += (se as f64) * sarg + (ce as f64) * carg;
        }

        (dpsi * MICROARCSEC_TO_RAD, deps * MICROARCSEC_TO_RAD)
    }
}
