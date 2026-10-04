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
//! Both functions accept a two-part Julian Date in TT. Split as `(jd1, jd2)`
//! where typically `jd1 = 2451545.0` (J2000.0) and `jd2` is days from that epoch.
//! Dates that are not finite, or more than
//! [`MAX_CENTURIES_FROM_J2000`](crate::constants::MAX_CENTURIES_FROM_J2000) from J2000.0,
//! return an error.
//!
//! # Example
//!
//! ```
//! use celestial_core::obliquity::iau_2006_mean_obliquity;
//! use celestial_core::constants::J2000_JD;
//!
//! // At J2000.0: 84381.406″, about 23.4392794°
//! let eps = iau_2006_mean_obliquity(J2000_JD, 0.0).unwrap();
//! assert_eq!(eps, 0.4090926006005829);
//! ```

use crate::constants::ARCSEC_TO_RAD;
use crate::errors::AstroResult;
use crate::utils::checked_jd_to_centuries;

/// Mean obliquity of the ecliptic using the IAU 2006 precession model.
///
/// Returns the mean obliquity in radians. This is a 5th-order polynomial
/// valid for several centuries around J2000.0.
///
/// At J2000.0: ε₀ = 84381.406″ ≈ 23°26′21.406″
pub fn iau_2006_mean_obliquity(date1: f64, date2: f64) -> AstroResult<f64> {
    Ok(iau_2006_obliquity_at(checked_jd_to_centuries(
        date1, date2,
    )?))
}

pub(crate) fn iau_2006_obliquity_at(t: f64) -> f64 {
    (84381.406
        + (-46.836769
            + (-0.0001831 + (0.00200340 + (-0.000000576 + (-0.0000000434) * t) * t) * t) * t)
            * t)
        * ARCSEC_TO_RAD
}

/// Mean obliquity of the ecliptic using the IAU 1980 model.
///
/// Returns the mean obliquity in radians. This is a 3rd-order polynomial,
/// less accurate than the IAU 2006 model but still used with IAU 1980 nutation.
///
/// At J2000.0: ε₀ = 84381.448″ ≈ 23°26′21.448″
pub fn iau_1980_mean_obliquity(date1: f64, date2: f64) -> AstroResult<f64> {
    let t = checked_jd_to_centuries(date1, date2)?;

    let obliquity_arcsec = 84381.448 + (-46.8150 + (-0.00059 + (0.001813) * t) * t) * t;

    Ok(obliquity_arcsec * ARCSEC_TO_RAD)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD};

    // (jd1, jd2, obliquity) from ERFA obl06 and obl80.
    const ERFA_OBL06: [(f64, f64, f64); 5] = [
        (2400000.5, 54388.0, 0.4090749229387258),
        (2451545.0, 0.0, 0.4090926006005829),
        (2400000.5, 60000.0, 0.4090400339553416),
        (2451545.0, -219150.0, 0.4104528950884669),
        (2451545.0, 219150.0, 0.40773223496051214),
    ];

    const ERFA_OBL80: [(f64, f64, f64); 5] = [
        (2400000.5, 54388.0, 0.4090751347643816),
        (2451545.0, 0.0, 0.40909280422232897),
        (2400000.5, 60000.0, 0.40904026189211384),
        (2451545.0, -219150.0, 0.41045259582761134),
        (2451545.0, 219150.0, 0.40773280666819484),
    ];

    #[test]
    fn test_iau_2006_mean_obliquity_matches_erfa_obl06() {
        for (jd1, jd2, expected) in ERFA_OBL06 {
            assert_eq!(
                iau_2006_mean_obliquity(jd1, jd2).unwrap(),
                expected,
                "{jd1} + {jd2}"
            );
        }
    }

    #[test]
    fn test_iau_1980_mean_obliquity_matches_erfa_obl80() {
        for (jd1, jd2, expected) in ERFA_OBL80 {
            assert_eq!(
                iau_1980_mean_obliquity(jd1, jd2).unwrap(),
                expected,
                "{jd1} + {jd2}"
            );
        }
    }

    #[test]
    fn test_mean_obliquity_rejects_epochs_outside_the_model_range() {
        let beyond = 20.5 * DAYS_PER_JULIAN_CENTURY;
        for model in [iau_2006_mean_obliquity, iau_1980_mean_obliquity] {
            assert!(model(J2000_JD, beyond).is_err());
            assert!(model(J2000_JD, -beyond).is_err());
            assert!(model(f64::NAN, 0.0).is_err());
        }
    }
}
