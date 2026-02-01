//! Angle normalization for astronomical coordinate systems.
//!
//! Different astronomical quantities require different angular ranges:
//!
//! | Quantity | Range | Function |
//! |----------|-------|----------|
//! | Right Ascension | [0, 2pi) | [`wrap_0_2pi`] |
//! | Hour Angle | [-pi, +pi) | [`wrap_pm_pi`] |
//! | Longitude (celestial) | [-pi, +pi) | [`wrap_pm_pi`] |
//! | Declination | [-pi/2, +pi/2] | [`clamp_dec`] |
//! | Latitude | [-pi/2, +pi/2] | [`clamp_dec`] |
//!
//! # Wrapping vs Clamping
//!
//! **Wrapping** preserves the direction on the sphere. An angle of 370 degrees
//! represents the same direction as 10 degrees, so `wrap_0_2pi` returns 10 degrees.
//!
//! **Clamping** enforces physical limits. Declination cannot exceed +/-90 degrees
//! because you cannot go "past" the pole. [`clamp_dec`] enforces this by saturating
//! at the limits rather than wrapping.
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
//! use celestial_core::angle::{wrap_0_2pi, wrap_pm_pi, clamp_dec};
//! use std::f64::consts::PI;
//!
//! // Right ascension: always positive
//! let ra = wrap_0_2pi(-0.5);  // -0.5 rad -> ~5.78 rad
//! assert!(ra > 0.0 && ra < 2.0 * PI);
//!
//! // Hour angle: centered on zero
//! let ha = wrap_pm_pi(3.5);  // 3.5 rad -> ~-2.78 rad (wrapped)
//! assert!(ha >= -PI && ha < PI);
//!
//! // Declination: cannot exceed poles
//! let dec = clamp_dec(2.0);  // 2.0 rad -> pi/2 (clamped to +90 deg)
//! assert!((dec - PI / 2.0).abs() < 1e-10);
//! ```
//!
//! # Algorithm Notes
//!
//! The wrapping functions use `libm::fmod` (via [`crate::math::fmod`]) rather than
//! the `%` operator because Rust's `%` is a remainder, not a modulo. For negative
//! numbers, remainder and modulo differ:
//!
//! - `-1.0 % 360.0` = `-1.0` (remainder, keeps sign of dividend)
//! - `fmod(-1.0, 360.0)` = `-1.0` (same as %, but well-defined for floats)
//!
//! After `fmod`, we adjust for the desired range.

use crate::constants::{HALF_PI, PI, TWOPI};
use crate::math::fmod;

/// Specifies which normalization convention to apply.
///
/// Used when a function needs to normalize angles but the appropriate range
/// depends on what the angle represents.
#[derive(Copy, Clone, Debug)]
pub enum NormalizeMode {
    /// Right ascension: wrap to [0, 2pi).
    Ra0To2Pi,
    /// Longitude or hour angle: wrap to [-pi, +pi).
    LonMinusPiToPi,
    /// Declination or latitude: clamp to [-pi/2, +pi/2].
    DecClamped,
}

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
/// use std::f64::consts::PI;
///
/// // 270 degrees -> -90 degrees
/// let x = wrap_pm_pi(3.0 * PI / 2.0);
/// assert!((x - (-PI / 2.0)).abs() < 1e-10);
///
/// // -270 degrees -> +90 degrees
/// let y = wrap_pm_pi(-3.0 * PI / 2.0);
/// assert!((y - (PI / 2.0)).abs() < 1e-10);
///
/// // Already in range: unchanged
/// let z = wrap_pm_pi(1.0);
/// assert!((z - 1.0).abs() < 1e-10);
/// ```
///
/// # Algorithm
///
/// 1. Reduce to [-2pi, +2pi) via `fmod(x, 2pi)`
/// 2. If result is >= pi or <= -pi, subtract/add 2pi to bring into range
#[inline]
pub fn wrap_pm_pi(x: f64) -> f64 {
    let w = fmod(x, TWOPI);
    if w.abs() >= PI {
        return w - TWOPI.copysign(x);
    }

    w
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
/// use std::f64::consts::PI;
///
/// // Negative angle -> positive equivalent
/// let x = wrap_0_2pi(-PI / 2.0);  // -90 deg -> 270 deg
/// assert!((x - 3.0 * PI / 2.0).abs() < 1e-10);
///
/// // Angle > 2pi -> reduced
/// let y = wrap_0_2pi(5.0 * PI);  // 900 deg -> 180 deg
/// assert!((y - PI).abs() < 1e-10);
///
/// // Already in range: unchanged
/// let z = wrap_0_2pi(1.0);
/// assert!((z - 1.0).abs() < 1e-10);
/// ```
///
/// # Algorithm
///
/// 1. Reduce to (-2pi, +2pi) via `fmod(x, 2pi)`
/// 2. If result is negative, add 2pi to make it positive
#[inline]
pub fn wrap_0_2pi(x: f64) -> f64 {
    let w = fmod(x, TWOPI);
    if w < 0.0 {
        w + TWOPI
    } else {
        w
    }
}

/// Clamps an angle to [-pi/2, +pi/2] radians (i.e., [-90, +90] degrees).
///
/// Unlike wrapping, clamping saturates at the limits. This is appropriate for
/// quantities that have physical bounds you cannot exceed.
///
/// # Arguments
///
/// * `x` - Angle in radians
///
/// # Returns
///
/// The angle clamped to [-pi/2, +pi/2].
///
/// # When to Use
///
/// - **Declination**: celestial latitude, poles are at +/-90 degrees
/// - **Geographic latitude**: -90 (south pole) to +90 (north pole)
/// - **Altitude**: horizon at 0, zenith at +90 (though altitude can be negative)
///
/// # Why Clamp Instead of Wrap?
///
/// You cannot go "past" the north pole by walking north. If you try, you end up
/// walking south on the other side. This is fundamentally different from longitude,
/// where walking east forever eventually brings you back to where you started.
///
/// Clamping is a safety mechanism. If your calculation produces declination = 100 degrees,
/// something is wrong upstream. Clamping to 90 degrees prevents downstream code from
/// breaking, but you should investigate why the input was out of range.
///
/// # Examples
///
/// ```
/// use celestial_core::angle::clamp_dec;
/// use std::f64::consts::FRAC_PI_2;
///
/// // Within range: unchanged
/// let x = clamp_dec(0.5);
/// assert!((x - 0.5).abs() < 1e-10);
///
/// // Above +90 degrees: clamped to +90
/// let y = clamp_dec(2.0);
/// assert!((y - FRAC_PI_2).abs() < 1e-10);
///
/// // Below -90 degrees: clamped to -90
/// let z = clamp_dec(-2.0);
/// assert!((z - (-FRAC_PI_2)).abs() < 1e-10);
/// ```
///
/// # See Also
///
/// For declinations that might legitimately exceed +/-90 degrees (e.g., pier-flipped
/// telescope positions), see the validation functions in the [`validate`](super::validate)
/// module which offer extended range options.
#[inline]
pub fn clamp_dec(x: f64) -> f64 {
    x.clamp(-HALF_PI, HALF_PI)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wrap_pm_pi() {
        // In range: unchanged
        assert_eq!(wrap_pm_pi(1.0), 1.0);
        // Positive overflow: 270° -> -90°
        assert!((wrap_pm_pi(3.0 * PI / 2.0) - (-PI / 2.0)).abs() < 1e-15);
        // Negative overflow: -270° -> +90°
        assert!((wrap_pm_pi(-3.0 * PI / 2.0) - (PI / 2.0)).abs() < 1e-15);
        // At boundary: ±π both wrap (abs >= PI triggers adjustment)
        assert!((wrap_pm_pi(PI) - (-PI)).abs() < 1e-15);
    }

    #[test]
    fn test_wrap_0_2pi() {
        // In range: unchanged
        assert_eq!(wrap_0_2pi(1.0), 1.0);
        // Negative becomes positive: -90° -> 270°
        assert!((wrap_0_2pi(-PI / 2.0) - (3.0 * PI / 2.0)).abs() < 1e-15);
        // Overflow: 3π -> π
        assert!((wrap_0_2pi(3.0 * PI) - PI).abs() < 1e-15);
        // At 2π: wraps to 0
        assert!(wrap_0_2pi(TWOPI).abs() < 1e-15);
    }

    #[test]
    fn test_clamp_dec() {
        // In range: unchanged
        assert_eq!(clamp_dec(0.5), 0.5);
        // At boundary: unchanged
        assert_eq!(clamp_dec(HALF_PI), HALF_PI);
        // Overflow: clamped to ±π/2
        assert_eq!(clamp_dec(2.0), HALF_PI);
        assert_eq!(clamp_dec(-2.0), -HALF_PI);
    }
}
