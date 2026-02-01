use astro_core::math::fmod;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SiderealAngle {
    angle_hours: f64,
    exact_radians: Option<f64>,
}

impl SiderealAngle {
    pub fn from_hours(hours: f64) -> Self {
        Self {
            angle_hours: Self::normalize_hours(hours),
            exact_radians: None,
        }
    }

    pub fn from_degrees(degrees: f64) -> Self {
        Self::from_hours(degrees / 15.0)
    }

    pub fn from_radians(radians: f64) -> Self {
        Self::from_hours(radians * 12.0 / astro_core::constants::PI)
    }

    pub(crate) fn from_radians_exact(radians: f64) -> Self {
        let hours = radians * 12.0 / astro_core::constants::PI;
        Self {
            angle_hours: Self::normalize_hours(hours),
            exact_radians: Some(radians),
        }
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
            self.angle_hours * astro_core::constants::PI / 12.0
        }
    }

    fn normalize_hours(hours: f64) -> f64 {
        let mut normalized = fmod(hours, 24.0);
        if normalized < 0.0 {
            normalized += 24.0;
        }
        normalized
    }

    pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> f64 {
        self.hours() - target_ra_hours
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

    #[test]
    fn test_angle_conversions() {
        let angle = SiderealAngle::from_hours(6.0);

        assert_eq!(angle.hours(), 6.0);
        assert_eq!(angle.degrees(), 90.0);
        assert!((angle.radians() - astro_core::constants::HALF_PI).abs() < 1e-15);
    }

    #[test]
    fn test_normalization() {
        let angle1 = SiderealAngle::from_hours(25.5);
        assert_eq!(angle1.hours(), 1.5);

        let angle2 = SiderealAngle::from_hours(-1.5);
        assert_eq!(angle2.hours(), 22.5);
    }

    #[test]
    fn test_hour_angle_calculation() {
        let lst = SiderealAngle::from_hours(12.0);
        let target_ra = 6.0;
        let hour_angle = lst.hour_angle_to_target(target_ra);
        assert_eq!(hour_angle, 6.0);
    }
}
