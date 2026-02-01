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

use crate::constants::{ARCSEC_TO_RAD, TWOPI};
use crate::errors::{AstroError, AstroResult};

/// Computes the CIO locator angle `s` for a given epoch.
///
/// The locator is model-dependent. Currently only IAU 2006A is implemented,
/// which uses the IAU 2006 precession with IAU 2000A nutation.
///
/// # Example
///
/// ```
/// use astro_core::cio::CioLocator;
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
#[derive(Debug, Clone)]
pub struct CioLocator {
    tt_centuries: f64,
    model: CioModel,
}

#[derive(Clone, Copy)]
struct SeriesTerm {
    coeffs: [i8; 8],
    sine: f64,
    cosine: f64,
}

#[allow(clippy::excessive_precision)]
const SP: [f64; 6] = [
    94.00e-6,
    3808.65e-6,
    -122.68e-6,
    -72574.11e-6,
    27.98e-6,
    15.62e-6,
];

#[allow(clippy::excessive_precision)]
const S0: [SeriesTerm; 33] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: -2640.73e-6,
        cosine: 0.39e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: -63.53e-6,
        cosine: 0.02e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 3, 0, 0, 0],
        sine: -11.75e-6,
        cosine: -0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 1, 0, 0, 0],
        sine: -11.21e-6,
        cosine: -0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 2, 0, 0, 0],
        sine: 4.57e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 3, 0, 0, 0],
        sine: -2.02e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 1, 0, 0, 0],
        sine: -1.98e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 3, 0, 0, 0],
        sine: 1.72e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 0, 0, 1, 0, 0, 0],
        sine: 1.41e-6,
        cosine: 0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 0, 0, -1, 0, 0, 0],
        sine: 1.26e-6,
        cosine: 0.01e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, -1, 0, 0, 0],
        sine: 0.63e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, 1, 0, 0, 0],
        sine: 0.63e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 3, 0, 0, 0],
        sine: -0.46e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 1, 0, 0, 0],
        sine: -0.45e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 4, -4, 4, 0, 0, 0],
        sine: -0.36e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 1, -1, 1, -8, 12, 0],
        sine: 0.24e-6,
        cosine: 0.12e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 0, 0, 0, 0],
        sine: -0.32e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 2, 0, 0, 0],
        sine: -0.28e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 3, 0, 0, 0],
        sine: -0.27e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 1, 0, 0, 0],
        sine: -0.26e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 0, 0, 0, 0],
        sine: 0.21e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, -2, 2, -3, 0, 0, 0],
        sine: -0.19e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, -2, 2, -1, 0, 0, 0],
        sine: -0.18e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 0, 8, -13, -1],
        sine: 0.10e-6,
        cosine: -0.05e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 2, 0, 0, 0, 0],
        sine: -0.15e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, -2, 0, -1, 0, 0, 0],
        sine: 0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 2, 0, 0, 0],
        sine: 0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, -2, 1, 0, 0, 0],
        sine: -0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, -2, -1, 0, 0, 0],
        sine: -0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 4, -2, 4, 0, 0, 0],
        sine: -0.13e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 4, 0, 0, 0],
        sine: 0.11e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, 0, -3, 0, 0, 0],
        sine: -0.11e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, 0, -1, 0, 0, 0],
        sine: -0.11e-6,
        cosine: 0.00e-6,
    },
];

#[allow(clippy::excessive_precision)]
const S1: [SeriesTerm; 3] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: -0.07e-6,
        cosine: 3.57e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: 1.73e-6,
        cosine: -0.03e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 3, 0, 0, 0],
        sine: 0.00e-6,
        cosine: 0.48e-6,
    },
];

#[allow(clippy::excessive_precision)]
const S2: [SeriesTerm; 25] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: 743.52e-6,
        cosine: -0.17e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 2, 0, 0, 0],
        sine: 56.91e-6,
        cosine: 0.06e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 2, 0, 0, 0],
        sine: 9.84e-6,
        cosine: -0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: -8.85e-6,
        cosine: 0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 0, 0, 0, 0, 0, 0],
        sine: -6.38e-6,
        cosine: -0.05e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, 0, 0, 0, 0],
        sine: -3.07e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 2, 0, 0, 0],
        sine: 2.23e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 1, 0, 0, 0],
        sine: 1.67e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 2, 0, 0, 0],
        sine: 1.30e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, -2, 2, -2, 0, 0, 0],
        sine: 0.93e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, -2, 0, 0, 0, 0],
        sine: 0.68e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 1, 0, 0, 0],
        sine: -0.55e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, 0, -2, 0, 0, 0],
        sine: 0.53e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 2, 0, 0, 0, 0],
        sine: -0.27e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, 1, 0, 0, 0],
        sine: -0.27e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, -2, -2, 0, 0, 0],
        sine: -0.26e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, -1, 0, 0, 0],
        sine: -0.25e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 1, 0, 0, 0],
        sine: 0.22e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, 0, -2, 0, 0, 0, 0],
        sine: -0.21e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, -2, 0, -1, 0, 0, 0],
        sine: 0.20e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 2, 2, 0, 0, 0],
        sine: 0.17e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, 2, 0, 2, 0, 0, 0],
        sine: 0.13e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, 0, 0, 0, 0, 0, 0],
        sine: -0.13e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, -2, 2, 0, 0, 0],
        sine: -0.12e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 0, 0, 0, 0],
        sine: -0.11e-6,
        cosine: 0.00e-6,
    },
];

#[allow(clippy::excessive_precision)]
const S3: [SeriesTerm; 4] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: 0.30e-6,
        cosine: -23.42e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 2, 0, 0, 0],
        sine: -0.03e-6,
        cosine: -1.46e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 2, 0, 0, 0],
        sine: -0.01e-6,
        cosine: -0.25e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: 0.00e-6,
        cosine: 0.23e-6,
    },
];

#[allow(clippy::excessive_precision)]
const S4: [SeriesTerm; 1] = [SeriesTerm {
    coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
    sine: -0.26e-6,
    cosine: -0.01e-6,
}];

#[inline]
fn sum_terms(w: &mut f64, terms: &[SeriesTerm], fa: &[f64; 8]) {
    for term in terms.iter().rev() {
        let mut arg = 0.0;
        for (i, item) in fa.iter().enumerate() {
            arg += f64::from(term.coeffs[i]) * item;
        }
        *w += term.sine * arg.sin() + term.cosine * arg.cos();
    }
}

fn cio_series_s_plus_xy_half(t: f64) -> f64 {
    let fa = fundamental_arguments(t);

    let mut w0 = SP[0];
    sum_terms(&mut w0, &S0, &fa);

    let mut w1 = SP[1];
    sum_terms(&mut w1, &S1, &fa);

    let mut w2 = SP[2];
    sum_terms(&mut w2, &S2, &fa);

    let mut w3 = SP[3];
    sum_terms(&mut w3, &S3, &fa);

    let mut w4 = SP[4];
    sum_terms(&mut w4, &S4, &fa);

    let w5 = SP[5];

    let series_arcsec = w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t;
    series_arcsec * ARCSEC_TO_RAD
}

fn fundamental_arguments(t: f64) -> [f64; 8] {
    [
        mean_anomaly_moon(t),
        mean_anomaly_sun(t),
        mean_longitude_moon_minus_node(t),
        mean_elongation_moon_sun(t),
        mean_longitude_ascending_node_moon(t),
        mean_longitude_venus(t),
        mean_longitude_earth(t),
        general_precession_longitude(t),
    ]
}

#[allow(clippy::excessive_precision)]
fn mean_anomaly_moon(t: f64) -> f64 {
    (485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470)))))
        % crate::constants::CIRCULAR_ARCSECONDS
        * ARCSEC_TO_RAD
}

#[allow(clippy::excessive_precision)]
fn mean_anomaly_sun(t: f64) -> f64 {
    (1287104.793048 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))))
        % crate::constants::CIRCULAR_ARCSECONDS
        * ARCSEC_TO_RAD
}

#[allow(clippy::excessive_precision)]
fn mean_longitude_moon_minus_node(t: f64) -> f64 {
    (335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))))
        % crate::constants::CIRCULAR_ARCSECONDS
        * ARCSEC_TO_RAD
}

#[allow(clippy::excessive_precision)]
fn mean_elongation_moon_sun(t: f64) -> f64 {
    (1072260.703692 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))))
        % crate::constants::CIRCULAR_ARCSECONDS
        * ARCSEC_TO_RAD
}

#[allow(clippy::excessive_precision)]
fn mean_longitude_ascending_node_moon(t: f64) -> f64 {
    (450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939)))))
        % crate::constants::CIRCULAR_ARCSECONDS
        * ARCSEC_TO_RAD
}

#[allow(clippy::excessive_precision)]
fn mean_longitude_venus(t: f64) -> f64 {
    (3.176146697 + 1021.3285546211 * t) % TWOPI
}

#[allow(clippy::excessive_precision)]
fn mean_longitude_earth(t: f64) -> f64 {
    (1.753470314 + 628.3075849991 * t) % TWOPI
}

#[allow(clippy::excessive_precision)]
fn general_precession_longitude(t: f64) -> f64 {
    (0.024381750 + 0.00000538691 * t) * t
}

/// The precession-nutation model used for CIO locator computation.
#[derive(Debug, Clone, Copy, PartialEq)]
enum CioModel {
    /// IAU 2006 precession with IAU 2000A nutation (66 terms).
    IAU2006A,
}

impl CioLocator {
    /// Creates a CIO locator using the IAU 2006/2000A model.
    ///
    /// # Parameters
    ///
    /// * `tt_centuries` - TT (Terrestrial Time) as Julian centuries from J2000.0.
    ///   Computed as `(JD_TT - 2451545.0) / astro_core::constants::DAYS_PER_JULIAN_CENTURY`.
    pub fn iau2006a(tt_centuries: f64) -> Self {
        Self {
            tt_centuries,
            model: CioModel::IAU2006A,
        }
    }

    /// Computes the CIO locator `s` given the CIP coordinates.
    ///
    /// # Parameters
    ///
    /// * `x` - CIP X coordinate in radians (from NPB matrix element [2][0])
    /// * `y` - CIP Y coordinate in radians (from NPB matrix element [2][1])
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
        match self.model {
            CioModel::IAU2006A => self.calculate_iau2006a(x, y),
        }
    }

    fn calculate_iau2006a(&self, x: f64, y: f64) -> AstroResult<f64> {
        let t = self.tt_centuries;

        if t.abs() > 20.0 {
            return Err(AstroError::math_error(
                "CIO locator calculation",
                crate::errors::MathErrorKind::InvalidInput,
                &format!(
                    "Time too far from J2000.0 for CIO locator: {:.1} centuries",
                    t
                ),
            ));
        }

        let s_series = cio_series_s_plus_xy_half(t);
        let s = s_series - 0.5 * x * y;

        Ok(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cio_locator_at_j2000() {
        let locator = CioLocator::iau2006a(0.0);
        let s = locator.calculate(0.0, 0.0).unwrap();

        assert!(
            s.abs() < 1e-6,
            "CIO locator at J2000.0 should be small: {}",
            s
        );
    }

    #[test]
    fn test_cio_locator_time_dependency() {
        let locator_past = CioLocator::iau2006a(-1.0);
        let locator_future = CioLocator::iau2006a(1.0);

        let s_past = locator_past.calculate(0.0, 0.0).unwrap();
        let s_future = locator_future.calculate(0.0, 0.0).unwrap();

        assert_ne!(s_past, s_future);
        assert!(
            (s_future - s_past).abs() > 1e-8,
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
}
