//! CIO locator (s) for the IAU 2006/2000A precession-nutation model.
//!
//! The CIO locator `s` positions the Celestial Intermediate Origin on the CIP equator.
//! It's the arc length from the GCRS x-axis intersection to the CIO, measured along
//! the CIP equator. This small angle (microarcseconds) completes the transformation
//! from GCRS to CIRS coordinates.
//!
//! # The transformation chain
//!
//! GCRS -> CIRS requires three pieces:
//! 1. CIP coordinates (X, Y) - where the pole is
//! 2. CIO locator (s) - where the origin is
//! 3. Earth Rotation Angle - how much Earth has rotated
//!
//! This module provides piece #2.
//!
//! # When to use this
//!
//! You need the CIO locator when:
//! - Building the GCRS-to-CIRS rotation matrix via `gcrs_to_cirs_matrix(x, y, s)`
//! - Computing equation of the origins (difference between CEO-based and equinox-based sidereal time)
//! - Implementing IAU 2000/2006 compliant coordinate transformations
//!
//! For most uses, [`CioSolution::calculate`](super::CioSolution::calculate) handles this automatically.
//!
//! # Algorithm
//!
//! Uses the IAU 2006/2000A series expansion with 66 periodic terms across 5 polynomial
//! orders, plus a polynomial part. The full expression is:
//!
//! ```text
//! s = series(t) - X*Y/2
//! ```
//!
//! where `t` is TT centuries from J2000.0. The `-X*Y/2` term accounts for the
//! frame rotation induced by the CIP motion.
//!
//! # References
//!
//! - Capitaine et al. (2003), A&A 400, 1145-1154
//! - IERS Conventions (2010), Chapter 5
//! - SOFA library: `iauS06` function

use crate::constants::ARCSEC_TO_RAD;
use crate::errors::{AstroError, AstroResult};
use crate::math::polynomial;
use crate::nutation::fundamental_args::IERS2010FundamentalArgs;

use super::locator_terms::{SeriesTerm, S0, S1, S2, S3, S4, SP};

/// Computes the CIO locator angle `s` for a given epoch.
///
/// The locator is model-dependent. Currently only IAU 2006A is implemented,
/// which uses the IAU 2006 precession with IAU 2000A nutation.
///
/// # Example
///
/// ```
/// use celestial_core::cio::CioLocator;
///
/// // Compute s for J2000.0 + 0.5 centuries (year ~2050)
/// let locator = CioLocator::iau2006a(0.5);
///
/// // X, Y from CIP coordinates (typically from precession-nutation matrix)
/// let x = 1.0e-7;  // radians
/// let y = 2.0e-7;  // radians
///
/// let s = locator.calculate(x, y).unwrap();
/// // s is in radians, typically on the order of 10^-8
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CioLocator {
    tt_centuries: f64,
}

#[inline]
fn series_sum(constant: f64, terms: &[SeriesTerm], fa: &[f64; 8]) -> f64 {
    let mut w = constant;
    for term in terms.iter().rev() {
        let mut arg = 0.0;
        for (i, item) in fa.iter().enumerate() {
            arg += f64::from(term.coeffs[i]) * item;
        }
        w += term.sine * libm::sin(arg) + term.cosine * libm::cos(arg);
    }
    w
}

fn cio_series_s_plus_xy_half(t: f64) -> f64 {
    let fa = fundamental_arguments(t);
    let coefficients = [
        series_sum(SP[0], &S0, &fa),
        series_sum(SP[1], &S1, &fa),
        series_sum(SP[2], &S2, &fa),
        series_sum(SP[3], &S3, &fa),
        series_sum(SP[4], &S4, &fa),
        SP[5],
    ];
    polynomial(&coefficients, t) * ARCSEC_TO_RAD
}

fn fundamental_arguments(t: f64) -> [f64; 8] {
    [
        t.moon_mean_anomaly(),
        t.sun_mean_anomaly(),
        t.mean_argument_of_latitude(),
        t.mean_elongation(),
        t.moon_ascending_node_longitude(),
        t.venus_lng(),
        t.earth_lng(),
        t.precession(),
    ]
}

impl CioLocator {
    /// Creates a CIO locator using the IAU 2006/2000A model.
    ///
    /// # Parameters
    ///
    /// * `tt_centuries` - TT (Terrestrial Time) as Julian centuries from J2000.0.
    ///   Computed as `(JD_TT - 2451545.0) / celestial_core::constants::DAYS_PER_JULIAN_CENTURY`.
    pub fn iau2006a(tt_centuries: f64) -> Self {
        Self { tt_centuries }
    }

    /// Computes the CIO locator `s` given the CIP coordinates.
    ///
    /// # Parameters
    ///
    /// * `x` - CIP X coordinate in radians (from NPB matrix element `[2][0]`)
    /// * `y` - CIP Y coordinate in radians (from NPB matrix element `[2][1]`)
    ///
    /// # Returns
    ///
    /// The CIO locator `s` in radians.
    ///
    /// # Errors
    ///
    /// Returns an error if the epoch is more than 20 centuries from J2000.0,
    /// where the series expansion becomes unreliable.
    pub fn calculate(&self, x: f64, y: f64) -> AstroResult<f64> {
        let t = self.tt_centuries;
        validate_locator_inputs(t, x, y)?;

        let s_series = cio_series_s_plus_xy_half(t);
        let s = s_series - 0.5 * x * y;

        Ok(s)
    }
}

fn validate_locator_inputs(t: f64, x: f64, y: f64) -> AstroResult<()> {
    crate::utils::check_model_epoch(t)?;
    if x.is_finite() && y.is_finite() {
        return Ok(());
    }
    Err(AstroError::math_error(
        "CIO locator calculation",
        crate::errors::MathErrorKind::NotFinite,
        "CIP coordinates must be finite",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    // (jd1, jd2, x, y, s) from ERFA s06, run with Rust libm for sin/cos.
    const ERFA_S06: [(f64, f64, f64, f64, f64); 2] = [
        (
            2400000.5,
            53736.0,
            0.0005791308486706011,
            4.020579816732961e-5,
            -1.2200322130764642e-8,
        ),
        (2451545.0, 0.0, 0.0, 0.0, -9.756652246326891e-9),
    ];

    #[test]
    fn test_cio_locator_matches_erfa_s06() {
        for (jd1, jd2, x, y, expected) in ERFA_S06 {
            let t = ((jd1 - crate::constants::J2000_JD) + jd2)
                / crate::constants::DAYS_PER_JULIAN_CENTURY;
            let s = CioLocator::iau2006a(t).calculate(x, y).unwrap();
            assert_eq!(s, expected, "{jd1} + {jd2}");
        }
    }

    #[test]
    fn test_cio_locator_time_dependency() {
        let locator_past = CioLocator::iau2006a(-1.0);
        let locator_future = CioLocator::iau2006a(1.0);

        let s_past = locator_past.calculate(0.0, 0.0).unwrap();
        let s_future = locator_future.calculate(0.0, 0.0).unwrap();

        assert_ne!(s_past, s_future);
        assert!(
            libm::fabs(s_future - s_past) > 1e-8,
            "CIO locator should show time dependence"
        );
    }

    #[test]
    fn test_cio_locator_cip_dependency() {
        let locator = CioLocator::iau2006a(0.0);

        let s_zero = locator.calculate(0.0, 0.0).unwrap();
        let s_offset = locator.calculate(1e-6, 1e-6).unwrap();

        assert_ne!(s_zero, s_offset);
    }

    #[test]
    fn test_extreme_time_validation() {
        let locator = CioLocator::iau2006a(25.0);
        let result = locator.calculate(0.0, 0.0);

        assert!(result.is_err());
    }

    #[test]
    fn test_rejects_non_finite_inputs() {
        assert!(CioLocator::iau2006a(f64::NAN).calculate(0.0, 0.0).is_err());
        let locator = CioLocator::iau2006a(0.0);
        assert!(locator.calculate(f64::NAN, 0.0).is_err());
        assert!(locator.calculate(0.0, f64::INFINITY).is_err());
    }
}
