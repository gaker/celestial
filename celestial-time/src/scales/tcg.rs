//! Geocentric Coordinate Time (TCG) time scale.
//!
//! TCG is the proper time of a clock at rest at the geocenter, free from Earth's
//! gravitational potential. It runs faster than TT by approximately 22 microseconds
//! per year due to gravitational time dilation.
//!
//! # Background
//!
//! TCG was introduced by the IAU in 1991 as the coordinate time for the Geocentric
//! Celestial Reference System (GCRS). While TT is adjusted to match the rate of
//! proper time on Earth's geoid, TCG ticks at the rate of a clock experiencing
//! no gravitational potential.
//!
//! The relationship between TCG and TT is defined by IAU Resolution B1.9 (2000):
//!
//! ```text
//! TCG - TT = L_G * (JD_TT - T_0) * 86400
//!
//! where:
//!   L_G = 6.969290134e-10 (defining constant)
//!   T_0 = 2443144.5003725 (TCG/TT coincidence epoch, 1977-01-01 00:00:32.184)
//! ```
//!
//! # Usage
//!
//! ```
//! use celestial_time::{JulianDate, TCG};
//! use celestial_time::scales::tcg_from_calendar;
//!
//! let tcg = TCG::j2000();
//! let jd = tcg.to_julian_date();
//!
//! let tcg_cal = tcg_from_calendar(2000, 1, 1, 12, 0, 0.0);
//! ```
//!
//! # Precision
//!
//! TCG values use split Julian Date storage internally. The struct methods preserve
//! full f64 precision through all arithmetic operations. Conversions to/from TT
//! maintain nanosecond accuracy.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::SECONDS_PER_DAY_F64;
use std::fmt;
use std::str::FromStr;

/// Geocentric Coordinate Time.
///
/// Wraps a `JulianDate` representing an instant in the TCG time scale.
/// TCG is the coordinate time for the Geocentric Celestial Reference System,
/// running ~6.97e-10 faster than TT (about 22 microseconds per year).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TCG(JulianDate);

impl TCG {
    /// Creates a TCG instant from Unix timestamp components.
    ///
    /// Converts seconds and nanoseconds since 1970-01-01 00:00:00 to TCG.
    /// Note: This assumes the Unix timestamp is already in the TCG scale.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let total_seconds =
            seconds as f64 + nanos as f64 / celestial_core::constants::NANOSECONDS_PER_SECOND_F64;
        let jd = JulianDate::from_f64(UNIX_EPOCH_JD + total_seconds / SECONDS_PER_DAY_F64);
        Self(jd)
    }

    /// Creates a TCG instant from a Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Returns the J2000.0 epoch (2000-01-01 12:00:00) in TCG.
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Adds seconds to this TCG instant, returning a new TCG.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Adds days to this TCG instant, returning a new TCG.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }
}

impl fmt::Display for TCG {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TCG {}", self.0)
    }
}

/// Conversion from JulianDate to TCG.
impl From<JulianDate> for TCG {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

/// Parses ISO 8601 formatted strings into TCG.
///
/// Accepts standard date-time formats like "2000-01-01T12:00:00".
/// Fractional seconds are supported.
impl FromStr for TCG {
    type Err = TimeError;

    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

/// Creates a TCG instant from calendar date components.
///
/// Uses proleptic Gregorian calendar. No leap second or time zone handling;
/// the values are interpreted directly as TCG coordinates.
pub fn tcg_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> TCG {
    let jd = JulianDate::from_calendar(year, month, day, hour, minute, second);
    TCG::from_julian_date(jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_tcg_constructors() {
        assert_eq!(TCG::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(TCG::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            tcg_from_calendar(2000, 1, 1, 12, 0, 0.0)
                .to_julian_date()
                .to_f64(),
            J2000_JD
        );

        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let tcg_direct = TCG::from_julian_date(jd);
        let tcg_from_trait: TCG = jd.into();
        assert_eq!(
            tcg_direct.to_julian_date().jd1(),
            tcg_from_trait.to_julian_date().jd1()
        );
        assert_eq!(
            tcg_direct.to_julian_date().jd2(),
            tcg_from_trait.to_julian_date().jd2()
        );
    }

    #[test]
    fn test_tcg_arithmetic() {
        let tcg = TCG::j2000();
        assert_eq!(tcg.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(
            tcg.add_seconds(3600.0).to_julian_date().to_f64(),
            J2000_JD + 1.0 / 24.0
        );
    }

    #[test]
    fn test_tcg_string_parsing() {
        assert_eq!(
            TCG::from_str("2000-01-01T12:00:00")
                .unwrap()
                .to_julian_date()
                .to_f64(),
            TCG::j2000().to_julian_date().to_f64()
        );

        let result = TCG::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);

        assert!(TCG::from_str("invalid-date").is_err());
    }

    #[test]
    fn test_tcg_display() {
        let display_str = format!("{}", TCG::from_julian_date(JulianDate::new(J2000_JD, 0.5)));
        assert!(display_str.starts_with("TCG"));
        assert!(display_str.contains("2451545"));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_tcg_serde_round_trip() {
        let test_cases = [
            ("J2000", TCG::j2000()),
            ("Unix epoch", TCG::new(0, 0)),
            (
                "Modern date",
                tcg_from_calendar(2024, 6, 15, 14, 30, 45.123),
            ),
            (
                "High precision",
                tcg_from_calendar(1990, 12, 31, 23, 59, 59.999999999),
            ),
        ];

        for (name, original) in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: TCG = serde_json::from_str(&json).unwrap();

            let jd1_diff =
                (original.to_julian_date().jd1() - deserialized.to_julian_date().jd1()).abs();
            let jd2_diff =
                (original.to_julian_date().jd2() - deserialized.to_julian_date().jd2()).abs();
            let total_diff =
                (original.to_julian_date().to_f64() - deserialized.to_julian_date().to_f64()).abs();

            assert!(jd1_diff < 1e-14, "{}: jd1 diff {:.2e}", name, jd1_diff);
            assert!(jd2_diff < 1e-14, "{}: jd2 diff {:.2e}", name, jd2_diff);
            assert!(
                total_diff < 1e-14,
                "{}: total diff {:.2e}",
                name,
                total_diff
            );
        }
    }
}
