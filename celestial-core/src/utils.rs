//! Utility functions for time conversions.
//!
//! [`jd_to_centuries`] converts a two-part Julian Date to Julian centuries from J2000.0,
//! the time unit used by most IAU precession/nutation models.
//!
//! Angle normalization lives in [`crate::angle`] ([`wrap_pm_pi`](crate::angle::wrap_pm_pi),
//! [`wrap_0_2pi`](crate::angle::wrap_0_2pi)).

use crate::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD, MAX_CENTURIES_FROM_J2000};
use crate::errors::{AstroError, AstroResult, MathErrorKind};

/// Converts a two-part Julian Date to Julian centuries from J2000.0.
///
/// The two-part split preserves precision. Typically:
/// - `jd1 = 2451545.0` (J2000.0 epoch)
/// - `jd2` = days from that epoch
///
/// One Julian century = 36525 days.
///
/// # Example
///
/// ```
/// use celestial_core::utils::jd_to_centuries;
/// use celestial_core::constants::J2000_JD;
///
/// // At J2000.0 → t = 0
/// assert_eq!(jd_to_centuries(J2000_JD, 0.0), 0.0);
///
/// // One century later → t = 1
/// assert_eq!(jd_to_centuries(J2000_JD, celestial_core::constants::DAYS_PER_JULIAN_CENTURY), 1.0);
/// ```
#[inline]
pub fn jd_to_centuries(jd1: f64, jd2: f64) -> f64 {
    ((jd1 - J2000_JD) + jd2) / DAYS_PER_JULIAN_CENTURY
}

pub(crate) fn checked_jd_to_centuries(jd1: f64, jd2: f64) -> AstroResult<f64> {
    check_model_epoch(jd_to_centuries(jd1, jd2))
}

pub(crate) fn check_model_epoch(t: f64) -> AstroResult<f64> {
    if !t.is_finite() {
        return Err(AstroError::math_error(
            "model epoch",
            MathErrorKind::NotFinite,
            "Julian date must be finite",
        ));
    }
    if libm::fabs(t) <= MAX_CENTURIES_FROM_J2000 {
        return Ok(t);
    }
    Err(AstroError::math_error(
        "model epoch",
        MathErrorKind::OutOfRange,
        &format!(
            "Epoch {:.1} centuries from J2000.0 is outside the model range",
            t
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jd_to_centuries_j2000() {
        let t = jd_to_centuries(J2000_JD, 0.0);
        assert_eq!(t, 0.0);
    }

    #[test]
    fn test_jd_to_centuries_one_century() {
        let t = jd_to_centuries(J2000_JD, crate::constants::DAYS_PER_JULIAN_CENTURY);
        assert_eq!(t, 1.0);
    }

    #[test]
    fn test_jd_to_centuries_negative() {
        let t = jd_to_centuries(J2000_JD, -crate::constants::DAYS_PER_JULIAN_CENTURY);
        assert_eq!(t, -1.0);
    }

    #[test]
    fn test_jd_to_centuries_two_part() {
        let t = jd_to_centuries(crate::constants::MJD_ZERO_POINT, 51544.5);
        assert_eq!(t, 0.0);
    }

    #[test]
    fn test_jd_to_centuries_precision() {
        let jd2 = 0.123456789;
        let t = jd_to_centuries(J2000_JD, jd2);
        assert_eq!(t, 0.123456789 / crate::constants::DAYS_PER_JULIAN_CENTURY);
    }

    #[test]
    fn test_checked_jd_to_centuries_enforces_model_range() {
        let limit_days = MAX_CENTURIES_FROM_J2000 * DAYS_PER_JULIAN_CENTURY;
        assert_eq!(checked_jd_to_centuries(J2000_JD, limit_days).unwrap(), 20.0);
        assert_eq!(
            checked_jd_to_centuries(J2000_JD, -limit_days).unwrap(),
            -20.0
        );
        assert!(checked_jd_to_centuries(J2000_JD, limit_days + 1.0).is_err());
        assert!(checked_jd_to_centuries(J2000_JD, 1e300).is_err());
        assert!(checked_jd_to_centuries(f64::NAN, 0.0).is_err());
    }
}
