//! Angle normalization for astronomical coordinate systems.
//!
//! Different astronomical quantities require different angular ranges:
//!
//! | Quantity | Range | Function |
//! |----------|-------|----------|
//! | Right Ascension | [0, 2pi) | [`wrap_0_2pi`] |
//! | Hour Angle | [-pi, +pi) | [`wrap_pm_pi`] |
//! | Longitude (celestial) | [-pi, +pi) | [`wrap_pm_pi`] |
//!
//! Wrapping preserves the direction on the sphere. An angle of 370 degrees
//! represents the same direction as 10 degrees, so `wrap_0_2pi` returns 10 degrees.
//!
//! # Why Two Wrapping Functions?
//!
//! Right ascension and hour angle are both cyclic, but their conventions differ:
//!
//! - **Right Ascension** uses [0, 24h) or [0, 360 deg) because negative RA makes no sense.
//!   Stars at RA = 23h 59m are close to RA = 0h 01m on the sky.
//!
//! - **Hour Angle** uses [-12h, +12h) because it represents "hours from meridian."
//!   Negative means east of meridian (not yet crossed), positive means west (already crossed).
//!   The discontinuity at +/-180 degrees is at the anti-meridian, far from the observing position.
//!
//! - **Celestial Longitude** (e.g., in galactic or ecliptic coordinates) typically uses
//!   [-180, +180) to center the "interesting" region (galactic center, vernal equinox)
//!   at zero, with the discontinuity 180 degrees away.
//!
//! # Example
//!
//! ```
//! use celestial_core::angle::{wrap_0_2pi, wrap_pm_pi};
//! use celestial_core::constants::PI;
//!
//! // Right ascension: always positive
//! let ra = wrap_0_2pi(-0.5)?;  // -0.5 rad -> ~5.78 rad
//! assert!(ra > 0.0 && ra < 2.0 * PI);
//!
//! // Hour angle: centered on zero
//! let ha = wrap_pm_pi(3.5)?;  // 3.5 rad -> ~-2.78 rad (wrapped)
//! assert!(ha >= -PI && ha < PI);
//! # Ok::<(), celestial_core::errors::AstroError>(())
//! ```
//!
//! # Algorithm Notes
//!
//! The wrapping functions use `libm::fmod` (via [`crate::math::fmod`]) because the
//! crate uses `libm` for all floating-point math, for deterministic results across
//! platforms. Like Rust's `%` on `f64`, `fmod` is a remainder that keeps the sign of
//! the dividend, so it alone does not produce the target range:
//!
//! - `-1.0 % 360.0` = `-1.0`
//! - `fmod(-1.0, 360.0)` = `-1.0`
//!
//! After `fmod`, we adjust for the desired range.

use crate::constants::{PI, TWOPI};
use crate::errors::{AstroError, AstroResult, MathErrorKind};
use crate::math::fmod;

/// Wraps an angle to [-pi, +pi) radians.
///
/// Use for quantities where the discontinuity should be at +/-180 degrees
/// (the "back" of the circle), not at 0/360 degrees.
///
/// # Arguments
///
/// * `x` - Angle in radians (any value, including negative or > 2pi)
///
/// # Returns
///
/// The equivalent angle in [-pi, +pi).
///
/// # When to Use
///
/// - **Hour angle**: hours from meridian, negative = east, positive = west
/// - **Longitude differences**: shortest arc between two longitudes
/// - **Galactic/ecliptic longitude**: if you want galactic center at l=0
/// - **Position angle differences**: relative rotation between two frames
///
/// # Examples
///
/// ```
/// use celestial_core::angle::wrap_pm_pi;
/// use celestial_core::constants::PI;
///
/// // 270 degrees -> -90 degrees
/// let x = wrap_pm_pi(3.0 * PI / 2.0)?;
/// assert_eq!(x, -PI / 2.0);
///
/// // -270 degrees -> +90 degrees
/// let y = wrap_pm_pi(-3.0 * PI / 2.0)?;
/// assert_eq!(y, PI / 2.0);
///
/// // Already in range: unchanged
/// let z = wrap_pm_pi(1.0)?;
/// assert_eq!(z, 1.0);
/// # Ok::<(), celestial_core::errors::AstroError>(())
/// ```
///
/// # Algorithm
///
/// 1. Reduce to [-2pi, +2pi) via `fmod(x, 2pi)`
/// 2. If result is >= pi or <= -pi, subtract/add 2pi to bring into range
#[inline]
pub fn wrap_pm_pi(x: f64) -> AstroResult<f64> {
    let w = fmod(finite_angle(x, "wrap_pm_pi")?, TWOPI);
    if libm::fabs(w) >= PI {
        return Ok(w - libm::copysign(TWOPI, x));
    }

    Ok(w)
}

/// Wraps an angle to [0, 2pi) radians.
///
/// Use for quantities that are conventionally non-negative, with the
/// discontinuity at 0/360 degrees (midnight/noon for time-like quantities).
///
/// # Arguments
///
/// * `x` - Angle in radians (any value, including negative or > 2pi)
///
/// # Returns
///
/// The equivalent angle in [0, 2pi).
///
/// # When to Use
///
/// - **Right ascension**: 0h to 24h, never negative
/// - **Azimuth**: 0 to 360 degrees, measured from north through east
/// - **Sidereal time**: 0h to 24h
/// - **Mean anomaly, true anomaly**: orbital angles
///
/// # Examples
///
/// ```
/// use celestial_core::angle::wrap_0_2pi;
/// use celestial_core::constants::PI;
///
/// // Negative angle -> positive equivalent
/// let x = wrap_0_2pi(-PI / 2.0)?;  // -90 deg -> 270 deg
/// assert_eq!(x, 3.0 * PI / 2.0);
///
/// // Angle > 2pi -> reduced
/// let y = wrap_0_2pi(5.0 * PI)?;  // 900 deg -> 180 deg
/// assert_eq!(y, PI);
///
/// // Already in range: unchanged
/// let z = wrap_0_2pi(1.0)?;
/// assert_eq!(z, 1.0);
/// # Ok::<(), celestial_core::errors::AstroError>(())
/// ```
///
/// # Algorithm
///
/// 1. Reduce to (-2pi, +2pi) via `fmod(x, 2pi)`
/// 2. If result is negative, add 2pi to make it positive
#[inline]
pub fn wrap_0_2pi(x: f64) -> AstroResult<f64> {
    let w = fmod(finite_angle(x, "wrap_0_2pi")?, TWOPI);
    if w < 0.0 {
        Ok(w + TWOPI)
    } else {
        Ok(w)
    }
}

pub(super) fn finite_angle(x: f64, operation: &str) -> AstroResult<f64> {
    if x.is_finite() {
        return Ok(x);
    }
    Err(AstroError::math_error(
        operation,
        MathErrorKind::NotFinite,
        "Angle must be finite",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wrap_pm_pi() {
        assert_eq!(wrap_pm_pi(1.0).unwrap(), 1.0);
        assert_eq!(wrap_pm_pi(3.0 * PI / 2.0).unwrap(), -PI / 2.0);
        assert_eq!(wrap_pm_pi(-3.0 * PI / 2.0).unwrap(), PI / 2.0);
        // ERFA anpm convention: +π maps to -π.
        assert_eq!(wrap_pm_pi(PI).unwrap(), -PI);
    }

    #[test]
    fn test_wrap_0_2pi() {
        assert_eq!(wrap_0_2pi(1.0).unwrap(), 1.0);
        assert_eq!(wrap_0_2pi(-PI / 2.0).unwrap(), 3.0 * PI / 2.0);
        assert_eq!(wrap_0_2pi(3.0 * PI).unwrap(), PI);
        assert_eq!(wrap_0_2pi(TWOPI).unwrap(), 0.0);
    }

    #[test]
    fn test_wrap_rejects_non_finite() {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(wrap_pm_pi(x).is_err());
            assert!(wrap_0_2pi(x).is_err());
        }
    }
}
