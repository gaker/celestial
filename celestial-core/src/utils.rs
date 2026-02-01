//! Utility functions for time and angle conversions.
//!
//! Helper functions for common operations: Julian Date to centuries conversion,
//! angle normalization, and angular differences. These are building blocks used
//! throughout the library.
//!
//! # Time Conversion
//!
//! [`jd_to_centuries`] converts a two-part Julian Date to Julian centuries from J2000.0,
//! the time unit used by most IAU precession/nutation models.
//!
//! # Angle Normalization
//!
//! | Function | Input | Output Range |
//! |----------|-------|--------------|
//! | [`normalize_longitude`] | degrees | (-180°, 180°] |
//! | [`normalize_latitude`] | degrees | [-90°, 90°] (clamped) |
//! | [`normalize_angle_rad`] | radians | (-π, π] |
//!
//! # Angular Difference
//!
//! [`angular_difference`] computes the shortest signed difference between two angles
//! in degrees, handling the wraparound at ±180°.

use crate::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD, PI, TWOPI};

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

/// Normalizes longitude to the range (-180°, 180°].
///
/// Wraps values outside the range by adding/subtracting 360°.
#[inline]
pub fn normalize_longitude(lon: f64) -> f64 {
    let mut normalized = lon % 360.0;
    if normalized > 180.0 {
        normalized -= 360.0;
    } else if normalized < -180.0 {
        normalized += 360.0;
    }
    normalized
}

/// Clamps latitude to the valid range [-90°, 90°].
///
/// Values outside the range are clamped to the nearest pole.
#[inline]
pub fn normalize_latitude(lat: f64) -> f64 {
    lat.clamp(-90.0, 90.0)
}

/// Normalizes an angle in radians to the range (-π, π].
#[inline]
pub fn normalize_angle_rad(angle: f64) -> f64 {
    let mut normalized = angle % TWOPI;
    if normalized > PI {
        normalized -= TWOPI;
    } else if normalized < -PI {
        normalized += TWOPI;
    }
    normalized
}

/// Normalizes an angle in radians to the range [0, 2π).
#[inline]
pub fn normalize_angle_to_positive(angle: f64) -> f64 {
    let mut a = angle % TWOPI;
    if a < 0.0 {
        a += TWOPI;
    }
    a
}

/// Computes the shortest signed angular difference `a - b` in degrees.
///
/// Handles wraparound at ±180°. The result is in the range (-180°, 180°].
///
/// # Example
///
/// ```
/// use celestial_core::utils::angular_difference;
///
/// // Simple case
/// assert_eq!(angular_difference(90.0, 45.0), 45.0);
///
/// // Across the 0°/360° boundary: 10° is 20° ahead of 350°
/// assert!((angular_difference(10.0, 350.0) - 20.0).abs() < 1e-12);
/// ```
#[inline]
pub fn angular_difference(a: f64, b: f64) -> f64 {
    let mut diff = a - b;
    if diff > 180.0 {
        diff -= 360.0;
    } else if diff < -180.0 {
        diff += 360.0;
    }
    diff
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
        let expected = 0.123456789 / crate::constants::DAYS_PER_JULIAN_CENTURY;
        assert!((t - expected).abs() < 1e-15);
    }

    #[test]
    fn test_normalize_longitude() {
        assert_eq!(normalize_longitude(0.0), 0.0);
        assert_eq!(normalize_longitude(180.0), 180.0);
        assert_eq!(normalize_longitude(-180.0), -180.0);
        assert_eq!(normalize_longitude(181.0), -179.0);
        assert_eq!(normalize_longitude(-181.0), 179.0);
        assert_eq!(normalize_longitude(360.0), 0.0);
        assert_eq!(normalize_longitude(720.0), 0.0);
        assert_eq!(normalize_longitude(450.0), 90.0);
    }

    #[test]
    fn test_normalize_latitude() {
        assert_eq!(normalize_latitude(0.0), 0.0);
        assert_eq!(normalize_latitude(45.0), 45.0);
        assert_eq!(normalize_latitude(-45.0), -45.0);
        assert_eq!(normalize_latitude(90.0), 90.0);
        assert_eq!(normalize_latitude(-90.0), -90.0);
        assert_eq!(normalize_latitude(100.0), 90.0);
        assert_eq!(normalize_latitude(-100.0), -90.0);
    }

    #[test]
    fn test_normalize_angle_rad() {
        assert_eq!(normalize_angle_rad(0.0), 0.0);
        assert!((normalize_angle_rad(PI) - PI).abs() < 1e-15);
        assert!((normalize_angle_rad(-PI) - (-PI)).abs() < 1e-15);
        assert!((normalize_angle_rad(TWOPI)).abs() < 1e-15);
        assert!((normalize_angle_rad(3.0 * PI) - PI).abs() < 1e-15);
    }

    #[test]
    fn test_angular_difference() {
        assert_eq!(angular_difference(0.0, 0.0), 0.0);
        assert_eq!(angular_difference(90.0, 45.0), 45.0);
        assert_eq!(angular_difference(45.0, 90.0), -45.0);
        assert!((angular_difference(10.0, 350.0) - 20.0).abs() < 1e-12);
        assert!((angular_difference(-170.0, 170.0) - 20.0).abs() < 1e-12);
        assert!((angular_difference(350.0, 10.0) + 20.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_angle_to_positive() {
        assert_eq!(normalize_angle_to_positive(0.0), 0.0);
        assert!((normalize_angle_to_positive(TWOPI)).abs() < 1e-15);
        assert!((normalize_angle_to_positive(-PI) - PI).abs() < 1e-15);
        assert!((normalize_angle_to_positive(3.0 * PI) - PI).abs() < 1e-15);
        assert!(normalize_angle_to_positive(1.0) >= 0.0);
        assert!(normalize_angle_to_positive(-1.0) >= 0.0);
        assert!(normalize_angle_to_positive(-1.0) < TWOPI);
    }
}
