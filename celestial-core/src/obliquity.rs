//! Mean obliquity of the ecliptic.
//!
//! The obliquity is the angle between Earth's equatorial plane and the ecliptic
//! (the plane of Earth's orbit around the Sun). It's approximately 23.4° and
//! decreases slowly due to gravitational perturbations from other planets.
//!
//! This module provides two IAU models:
//!
//! | Function | Model | J2000.0 Value | Polynomial Order |
//! |----------|-------|---------------|------------------|
//! | [`iau_2006_mean_obliquity`] | IAU 2006 | 84381.406″ | 5th order |
//! | [`iau_1980_mean_obliquity`] | IAU 1980 | 84381.448″ | 3rd order |
//!
//! Both return the *mean* obliquity — the smoothly varying component without
//! short-period nutation oscillations. For the *true* obliquity (mean + nutation
//! in obliquity), add [`NutationResult::delta_eps`](crate::nutation::NutationResult::delta_eps).
//!
//! # Time Argument
//!
//! Both functions accept a two-part Julian Date in TDB. Split as `(jd1, jd2)`
//! where typically `jd1 = 2451545.0` (J2000.0) and `jd2` is days from that epoch.
//!
//! # Example
//!
//! ```
//! use celestial_core::obliquity::iau_2006_mean_obliquity;
//! use celestial_core::constants::J2000_JD;
//!
//! // At J2000.0
//! let eps = iau_2006_mean_obliquity(J2000_JD, 0.0);
//! let eps_deg = eps.to_degrees();
//! assert!((eps_deg - 23.4392794).abs() < 1e-6);
//! ```

use crate::constants::J2000_JD;

/// Mean obliquity of the ecliptic using the IAU 2006 precession model.
///
/// Returns the mean obliquity in radians. This is a 5th-order polynomial
/// valid for several centuries around J2000.0.
///
/// At J2000.0: ε₀ = 84381.406″ ≈ 23°26′21.406″
pub fn iau_2006_mean_obliquity(date1: f64, date2: f64) -> f64 {
    let t = ((date1 - J2000_JD) + date2) / crate::constants::DAYS_PER_JULIAN_CENTURY;

    let obliquity_arcsec = 84381.406
        + (-46.836769
            + (-0.0001831 + (0.00200340 + (-0.000000576 + (-0.0000000434) * t) * t) * t) * t)
            * t;

    obliquity_arcsec * (crate::constants::PI / (180.0 * 3600.0))
}

/// Mean obliquity of the ecliptic using the IAU 1980 model.
///
/// Returns the mean obliquity in radians. This is a 3rd-order polynomial,
/// less accurate than the IAU 2006 model but still used with IAU 1980 nutation.
///
/// At J2000.0: ε₀ = 84381.448″ ≈ 23°26′21.448″
pub fn iau_1980_mean_obliquity(date1: f64, date2: f64) -> f64 {
    let t = ((date1 - J2000_JD) + date2) / crate::constants::DAYS_PER_JULIAN_CENTURY;

    let obliquity_arcsec = 84381.448 + (-46.8150 + (-0.00059 + (0.001813) * t) * t) * t;

    obliquity_arcsec * (crate::constants::PI / (180.0 * 3600.0))
}
