//! Barycentric Coordinate Time (TCB) representation.
//!
//! TCB is the coordinate time for the barycentric reference frame, as defined by the IAU.
//! It ticks faster than TDB by approximately 1.55e-8 due to gravitational time dilation
//! (Earth sits in the Sun's gravitational well).
//!
//! # Relationship to TDB
//!
//! TCB and TDB are related by a linear transformation plus periodic terms:
//!
//! ```text
//! TCB - TDB = L_B * (JD_TCB - T_0) * 86400
//! ```
//!
//! Where:
//! - L_B = 1.550519768e-8 (IAU 2006 Resolution B3)
//! - T_0 = 2443144.5003725 (TCB-TDB epoch, 1977 Jan 1.0 TAI)
//!
//! # Usage
//!
//! ```
//! use astro_time::{JulianDate, TCB};
//!
//! // Create from Julian Date
//! let tcb = TCB::from_julian_date(JulianDate::j2000());
//!
//! // Create from calendar
//! use astro_time::scales::tcb::tcb_from_calendar;
//! let tcb = tcb_from_calendar(2000, 1, 1, 12, 0, 0.0);
//!
//! // Parse from ISO 8601
//! let tcb: TCB = "2000-01-01T12:00:00".parse().unwrap();
//! ```
//!
//! # When to Use TCB
//!
//! TCB is the natural time coordinate for barycentric calculations (solar system dynamics,
//! pulsar timing, VLBI). For most terrestrial applications, TDB is more practical since
//! it stays close to TT.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use astro_core::constants::SECONDS_PER_DAY_F64;
use std::fmt;
use std::str::FromStr;

/// Barycentric Coordinate Time.
///
/// Wraps a Julian Date interpreted in the TCB time scale. TCB is the proper time
/// for a clock at the solar system barycenter, far from gravitational sources.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TCB(JulianDate);

impl TCB {
    /// Creates TCB from Unix timestamp components.
    ///
    /// Interprets the given seconds and nanoseconds as elapsed since Unix epoch
    /// (1970-01-01T00:00:00) in the TCB scale.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let total_seconds =
            seconds as f64 + nanos as f64 / astro_core::constants::NANOSECONDS_PER_SECOND_F64;
        let jd = JulianDate::from_f64(UNIX_EPOCH_JD + total_seconds / SECONDS_PER_DAY_F64);
        Self(jd)
    }

    /// Creates TCB from a Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Returns TCB at the J2000.0 epoch (2000-01-01T12:00:00).
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Returns a new TCB advanced by the given seconds.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Returns a new TCB advanced by the given days.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }
}

impl fmt::Display for TCB {
    /// Formats as "TCB <julian_date>".
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TCB {}", self.0)
    }
}

impl From<JulianDate> for TCB {
    /// Converts a Julian Date to TCB.
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

impl FromStr for TCB {
    type Err = TimeError;

    /// Parses an ISO 8601 datetime string as TCB.
    ///
    /// Accepts formats like "2000-01-01T12:00:00" or "2000-01-01T12:00:00.123".
    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

/// Creates TCB from calendar components.
///
/// Converts the given Gregorian calendar date and time to TCB. No time scale
/// corrections are applied; the calendar date is interpreted directly as TCB.
pub fn tcb_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> TCB {
    let jd = JulianDate::from_calendar(year, month, day, hour, minute, second);
    TCB::from_julian_date(jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;
    use astro_core::constants::J2000_JD;

    #[test]
    fn test_tcb_construction() {
        let test_cases: [(&str, TCB, f64); 3] = [
            ("new(0, 0) -> Unix epoch", TCB::new(0, 0), UNIX_EPOCH_JD),
            ("j2000() -> J2000_JD", TCB::j2000(), J2000_JD),
            (
                "calendar J2000 -> J2000_JD",
                tcb_from_calendar(2000, 1, 1, 12, 0, 0.0),
                J2000_JD,
            ),
        ];

        for (name, tcb, expected_jd) in test_cases {
            assert_eq!(tcb.to_julian_date().to_f64(), expected_jd, "{}", name);
        }
    }

    #[test]
    fn test_tcb_arithmetic() {
        let tcb = TCB::j2000();

        assert_eq!(tcb.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(
            tcb.add_seconds(3600.0).to_julian_date().to_f64(),
            J2000_JD + 1.0 / 24.0
        );
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_tcb_serde_round_trip() {
        let test_cases = [
            TCB::j2000(),
            TCB::new(0, 0),
            tcb_from_calendar(2024, 6, 15, 14, 30, 45.123),
            tcb_from_calendar(1990, 12, 31, 23, 59, 59.999999999),
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: TCB = serde_json::from_str(&json).unwrap();

            let diff =
                (original.to_julian_date().to_f64() - deserialized.to_julian_date().to_f64()).abs();
            assert!(diff < 1e-14, "serde precision loss: {:.2e}", diff);
        }
    }

    #[test]
    fn test_tcb_display() {
        let tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let s = format!("{}", tcb);

        assert!(s.starts_with("TCB"));
        assert!(s.contains("2451545"));
    }

    #[test]
    fn test_tcb_from_julian_date_trait() {
        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let tcb_direct = TCB::from_julian_date(jd);
        let tcb_trait: TCB = jd.into();

        assert_eq!(tcb_direct.to_julian_date(), tcb_trait.to_julian_date());
    }

    #[test]
    fn test_tcb_string_parsing() {
        let result = TCB::from_str("2000-01-01T12:00:00").unwrap();
        assert_eq!(result.to_julian_date().to_f64(), J2000_JD);

        let result = TCB::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);

        assert!(TCB::from_str("invalid-date").is_err());
    }
}
