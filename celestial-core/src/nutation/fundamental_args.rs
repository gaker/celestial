//! Fundamental arguments for nutation and precession models.
//!
//! Fundamental arguments are periodic angular quantities derived from the orbital
//! mechanics of the Earth-Moon-Sun system and planetary motions. They form the
//! basis for computing nutation, precession, and other Earth orientation effects.
//!
//! This module provides two trait-based interfaces:
//!
//! - [`IERS2010FundamentalArgs`]: Arguments from IERS Conventions (2010), including
//!   planetary longitudes, lunar arguments, and the general precession.
//!
//! - [`MHB2000FundamentalArgs`]: Additional arguments from the Mathews-Herring-Buffett
//!   2000 nutation model (MHB2000), used in IAU 2000A nutation.
//!
//! All methods are implemented on `f64` representing time in Julian centuries (TDB)
//! from J2000.0. Results are in radians, normalized to [0, 2π) where applicable.
//!
//! # References
//!
//! - IERS Conventions (2010), Chapter 5
//! - Mathews, Herring, & Buffett 2002, J. Geophys. Res. 107(B4)

use crate::constants::{ARCSEC_TO_RAD, CIRCULAR_ARCSECONDS, TWOPI};
use crate::math::fmod;

/// Fundamental arguments from IERS Conventions (2010).
///
/// These arguments are functions of TDB Julian centuries from J2000.0.
/// They include mean planetary longitudes and lunar orbital elements
/// required for computing nutation and precession.
///
/// # Usage
///
/// ```
/// use celestial_core::nutation::IERS2010FundamentalArgs;
///
/// let t: f64 = 0.1; // Julian centuries from J2000.0
/// let l = t.moon_mean_anomaly();
/// let f = t.mean_argument_of_latitude();
/// ```
pub trait IERS2010FundamentalArgs {
    /// Mean longitude of Mercury (radians).
    fn mercury_lng(&self) -> f64;

    /// Mean longitude of Venus (radians).
    fn venus_lng(&self) -> f64;

    /// Mean longitude of Earth (radians).
    fn earth_lng(&self) -> f64;

    /// Mean longitude of Mars (radians).
    fn mars_lng(&self) -> f64;

    /// Mean longitude of Jupiter (radians).
    fn jupiter_lng(&self) -> f64;

    /// Mean longitude of Saturn (radians).
    fn saturn_lng(&self) -> f64;

    /// Mean longitude of Uranus (radians).
    fn uranus_lng(&self) -> f64;

    /// General accumulated precession in longitude (radians).
    ///
    /// This is the precession of the ecliptic along the equator, not normalized
    /// to [0, 2π).
    fn precession(&self) -> f64;

    /// Mean anomaly of the Moon (radians), denoted l.
    ///
    /// Computed using a 4th-order polynomial in arcseconds, then converted
    /// to radians and normalized.
    fn moon_mean_anomaly(&self) -> f64;

    /// Mean argument of latitude of the Moon (radians), denoted F.
    ///
    /// The angular distance from the ascending node to the Moon, measured
    /// along the lunar orbit.
    fn mean_argument_of_latitude(&self) -> f64;

    /// Mean longitude of the Moon's ascending node (radians), denoted Ω.
    ///
    /// The point where the Moon's orbit crosses the ecliptic from south to
    /// north.
    fn moon_ascending_node_longitude(&self) -> f64;
}

impl IERS2010FundamentalArgs for f64 {
    #[inline]
    fn mercury_lng(&self) -> f64 {
        fmod(4.402608842 + 2608.7903141574 * self, TWOPI)
    }

    #[inline]
    fn venus_lng(&self) -> f64 {
        fmod(3.176146697 + 1021.3285546211 * self, TWOPI)
    }

    #[inline]
    fn earth_lng(&self) -> f64 {
        fmod(1.753470314 + 628.3075849991 * self, TWOPI)
    }

    #[inline]
    fn mars_lng(&self) -> f64 {
        fmod(6.203480913 + 334.0612426700 * self, TWOPI)
    }

    #[inline]
    fn jupiter_lng(&self) -> f64 {
        fmod(0.599546497 + 52.9690962641 * self, TWOPI)
    }

    #[inline]
    fn saturn_lng(&self) -> f64 {
        fmod(0.874016757 + 21.3299104960 * self, TWOPI)
    }

    #[inline]
    fn uranus_lng(&self) -> f64 {
        fmod(5.481293872 + 7.4781598567 * self, TWOPI)
    }

    #[inline]
    fn precession(&self) -> f64 {
        0.024381750 * self + 0.00000538691 * self * self
    }

    #[inline]
    fn moon_mean_anomaly(&self) -> f64 {
        let l = 485868.249036
            + self * (1717915923.2178 + self * (31.8792 + self * (0.051635 - self * 0.00024470)));
        fmod(l, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD
    }

    #[inline]
    fn mean_argument_of_latitude(&self) -> f64 {
        let f = 335779.526232
            + self * (1739527262.8478 + self * (-12.7512 + self * (-0.001037 + self * 0.00000417)));
        fmod(f, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD
    }

    #[inline]
    fn moon_ascending_node_longitude(&self) -> f64 {
        let om = 450160.398036
            + self * (-6962890.5431 + self * (7.4722 + self * (0.007702 - self * 0.00005939)));
        fmod(om, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD
    }
}

/// Additional fundamental arguments from the MHB2000 nutation model.
///
/// These arguments supplement [`IERS2010FundamentalArgs`] with expressions
/// specific to the Mathews-Herring-Buffett 2000 nutation series. The `_mhb`
/// suffix distinguishes these from IERS 2010 versions where expressions differ.
///
/// # Usage
///
/// ```
/// use celestial_core::nutation::MHB2000FundamentalArgs;
///
/// let t: f64 = 0.1; // Julian centuries from J2000.0
/// let lp = t.sun_mean_anomaly_mhb();
/// let d = t.mean_elongation_mhb();
/// ```
pub trait MHB2000FundamentalArgs {
    /// Mean anomaly of the Sun (radians), denoted l'.
    ///
    /// Uses the MHB2000 polynomial coefficients.
    fn sun_mean_anomaly_mhb(&self) -> f64;

    /// Mean elongation of the Moon from the Sun (radians), denoted D.
    ///
    /// The angular separation between the Moon and Sun as seen from Earth.
    fn mean_elongation_mhb(&self) -> f64;

    /// Mean longitude of Neptune (radians).
    fn neptune_longitude_mhb(&self) -> f64;
}

impl MHB2000FundamentalArgs for f64 {
    #[inline]
    fn sun_mean_anomaly_mhb(&self) -> f64 {
        let lp = 1287104.79305
            + self * (129596581.0481 + self * (-0.5532 + self * (0.000136 - self * 0.00001149)));
        fmod(lp, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD
    }

    #[inline]
    fn mean_elongation_mhb(&self) -> f64 {
        let d = 1072260.70369
            + self * (1602961601.2090 + self * (-6.3706 + self * (0.006593 - self * 0.00003169)));
        fmod(d, CIRCULAR_ARCSECONDS) * ARCSEC_TO_RAD
    }

    #[inline]
    fn neptune_longitude_mhb(&self) -> f64 {
        fmod(5.321159000 + 3.8127774000 * self, TWOPI)
    }
}
