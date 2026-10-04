use crate::julian::finite_arg;
use crate::TimeResult;
use celestial_core::constants::TWOPI;
use celestial_core::math::fmod;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SiderealAngle {
    angle_hours: f64,
    exact_radians: Option<f64>,
}

impl SiderealAngle {
    pub fn from_hours(hours: f64) -> TimeResult<Self> {
        Ok(Self {
            angle_hours: Self::normalize_hours(finite_arg("hours", hours)?),
            exact_radians: None,
        })
    }

    pub fn from_degrees(degrees: f64) -> TimeResult<Self> {
        Self::from_hours(finite_arg("degrees", degrees)? / 15.0)
    }

    pub fn from_radians(radians: f64) -> TimeResult<Self> {
        let radians = finite_arg("radians", radians)?;
        // An angle already in [0, 2π] is kept as given: wrap_0_2pi, like ERFA's
        // anp, can return 2π itself, and wrapping again would turn that into 0.
        let radians = if (0.0..=TWOPI).contains(&radians) {
            radians
        } else {
            wrap_to_period(radians, TWOPI)
        };
        Ok(Self {
            angle_hours: Self::normalize_hours(radians * 12.0 / celestial_core::constants::PI),
            exact_radians: Some(radians),
        })
    }

    pub fn hours(&self) -> f64 {
        self.angle_hours
    }

    pub fn degrees(&self) -> f64 {
        self.angle_hours * 15.0
    }

    pub fn radians(&self) -> f64 {
        if let Some(exact) = self.exact_radians {
            exact
        } else {
            self.angle_hours * celestial_core::constants::PI / 12.0
        }
    }

    // A tiny negative input wraps to 24 - ε, which rounds to 24 itself.
    fn normalize_hours(hours: f64) -> f64 {
        let wrapped = wrap_to_period(hours, 24.0);
        if wrapped == 24.0 {
            0.0
        } else {
            wrapped
        }
    }

    pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> TimeResult<f64> {
        let target_ra_hours = finite_arg("target_ra_hours", target_ra_hours)?;
        let hour_angle = fmod(self.hours() - target_ra_hours, 24.0);
        Ok(if hour_angle >= 12.0 {
            hour_angle - 24.0
        } else if hour_angle < -12.0 {
            hour_angle + 24.0
        } else {
            hour_angle
        })
    }
}

fn wrap_to_period(x: f64, period: f64) -> f64 {
    let wrapped = fmod(x, period);
    if wrapped < 0.0 {
        wrapped + period
    } else {
        wrapped
    }
}

// Whether the radians were given or derived from hours is not part of the
// value: two angles are equal when every accessor agrees.
impl PartialEq for SiderealAngle {
    fn eq(&self, other: &Self) -> bool {
        self.hours() == other.hours() && self.radians() == other.radians()
    }
}

impl fmt::Display for SiderealAngle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.6}h", self.angle_hours)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TimeError;

    fn not_finite(name: &str, shown: &str) -> TimeError {
        TimeError::ConversionError(format!("{} must be finite, got {}", name, shown))
    }

    #[test]
    fn test_non_finite_input_is_rejected() {
        let lst = SiderealAngle::from_hours(6.0).unwrap();
        for (bad, shown) in [
            (f64::NAN, "NaN"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
        ] {
            let rejected = |name| Err(not_finite(name, shown));
            assert_eq!(SiderealAngle::from_hours(bad), rejected("hours"));
            assert_eq!(SiderealAngle::from_degrees(bad), rejected("degrees"));
            assert_eq!(SiderealAngle::from_radians(bad), rejected("radians"));
            assert_eq!(
                lst.hour_angle_to_target(bad),
                Err(not_finite("target_ra_hours", shown))
            );
        }
    }

    #[test]
    fn test_angle_conversions() {
        let angle = SiderealAngle::from_hours(6.0).unwrap();

        assert_eq!(angle.hours(), 6.0);
        assert_eq!(angle.degrees(), 90.0);
        assert_eq!(angle.radians(), celestial_core::constants::HALF_PI);
    }

    #[test]
    fn test_normalization() {
        let angle1 = SiderealAngle::from_hours(25.5).unwrap();
        assert_eq!(angle1.hours(), 1.5);

        let angle2 = SiderealAngle::from_hours(-1.5).unwrap();
        assert_eq!(angle2.hours(), 22.5);
    }

    #[test]
    fn test_tiny_negative_angles_wrap_to_zero() {
        assert_eq!(SiderealAngle::from_hours(-1e-17).unwrap().hours(), 0.0);
        assert_eq!(SiderealAngle::from_degrees(-1e-16).unwrap().degrees(), 0.0);
    }

    #[test]
    fn test_equality_ignores_how_the_angle_was_given() {
        let given_hours = SiderealAngle::from_hours(6.0).unwrap();
        let given_radians =
            SiderealAngle::from_radians(celestial_core::constants::HALF_PI).unwrap();
        assert_eq!(given_hours, given_radians);
        assert_ne!(
            given_hours,
            SiderealAngle::from_hours(6.000000000000001).unwrap()
        );
    }

    #[test]
    fn test_hour_angle_calculation() {
        let lst = SiderealAngle::from_hours(12.0).unwrap();
        let target_ra = 6.0;
        let hour_angle = lst.hour_angle_to_target(target_ra).unwrap();
        assert_eq!(hour_angle, 6.0);
    }

    #[test]
    fn test_hour_angle_wraps_to_plus_minus_12() {
        let cases = [
            (1.0, 23.0, 2.0),
            (23.0, 1.0, -2.0),
            (6.0, 18.0, -12.0),
            (18.0, 6.0, -12.0),
            (0.25, 12.0, -11.75),
            (12.0, 0.25, 11.75),
            (5.0, 29.0, 0.0),
            (5.0, -30.5, 11.5),
            (0.0, 0.0, 0.0),
        ];
        for (lst_hours, ra_hours, expected) in cases {
            let lst = SiderealAngle::from_hours(lst_hours).unwrap();
            assert_eq!(
                lst.hour_angle_to_target(ra_hours).unwrap(),
                expected,
                "LST {lst_hours} h, RA {ra_hours} h"
            );
        }
    }
}
