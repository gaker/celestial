//! International Atomic Time (TAI) scale.
//!
//! TAI is the reference time scale for astronomical time conversions. It is maintained
//! by the Bureau International des Poids et Mesures (BIPM) as a weighted average of
//! over 400 atomic clocks worldwide.
//!
//! # Background
//!
//! TAI runs continuously without leap seconds. Its epoch is January 1, 1958, when
//! TAI and UT1 were approximately synchronized. TAI now leads UTC by 37 seconds
//! (as of 2017), with the difference increasing each time a leap second is added.
//!
//! Key relationships:
//!
//! ```text
//! TT  = TAI + 32.184 seconds (fixed)
//! GPS = TAI - 19 seconds (fixed)
//! UTC = TAI - leap_seconds (variable, table-based)
//! ```
//!
//! # Usage
//!
//! ```
//! use astro_time::{JulianDate, TAI};
//! use astro_time::scales::tai::tai_from_calendar;
//!
//! // From Julian Date
//! let tai = TAI::j2000();
//!
//! // From calendar date
//! let tai = tai_from_calendar(2024, 6, 15, 12, 30, 0.0);
//!
//! // Arithmetic
//! let later = tai.add_seconds(3600.0);
//! let next_day = tai.add_days(1.0);
//! ```
//!
//! # Precision
//!
//! TAI stores time as a split Julian Date (jd1, jd2) to preserve full f64 precision.
//! The split representation avoids precision loss when adding small time increments
//! to large Julian Date values.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use astro_core::constants::SECONDS_PER_DAY_F64;
use std::fmt;
use std::str::FromStr;

/// International Atomic Time representation.
///
/// Wraps a `JulianDate` to provide TAI-specific semantics. TAI serves as the
/// hub for conversions between other time scales (UTC, TT, GPS, TDB, etc.).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TAI(JulianDate);

impl TAI {
    /// Creates TAI from Unix timestamp components.
    ///
    /// Converts seconds since 1970-01-01 00:00:00 plus nanoseconds to TAI.
    /// Note: This assumes the input is already in TAI, not UTC.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let total_seconds =
            seconds as f64 + nanos as f64 / astro_core::constants::NANOSECONDS_PER_SECOND_F64;
        let jd = JulianDate::from_f64(UNIX_EPOCH_JD + total_seconds / SECONDS_PER_DAY_F64);
        Self(jd)
    }

    /// Creates TAI from a JulianDate.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Creates TAI from raw Julian Date components.
    ///
    /// Useful when you already have the split JD representation and want to
    /// avoid the overhead of creating a JulianDate first.
    pub fn from_julian_date_raw(jd1: f64, jd2: f64) -> Self {
        Self(JulianDate::new(jd1, jd2))
    }

    /// Returns TAI at the J2000.0 epoch (2000-01-01 12:00:00 TT).
    ///
    /// JD = 2451545.0
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Returns a new TAI offset by the given number of seconds.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Returns a new TAI offset by the given number of days.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }
}

impl fmt::Display for TAI {
    /// Formats as "TAI {julian_date}".
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TAI {}", self.0)
    }
}

impl From<JulianDate> for TAI {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

impl FromStr for TAI {
    type Err = TimeError;

    /// Parses an ISO 8601 datetime string as TAI.
    ///
    /// The input is interpreted directly as TAI with no UTC-to-TAI conversion.
    /// For UTC input, parse as UTC first, then convert to TAI.
    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

/// Creates TAI from calendar components.
///
/// Converts a Gregorian calendar date and time directly to TAI. No leap second
/// or timezone corrections are applied. The input is assumed to already be TAI.
///
/// # Arguments
///
/// * `year` - Gregorian year (negative for BCE)
/// * `month` - Month (1-12)
/// * `day` - Day of month (1-31)
/// * `hour` - Hour (0-23)
/// * `minute` - Minute (0-59)
/// * `second` - Second with optional fractional part (0.0-60.0)
pub fn tai_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> TAI {
    let jd = JulianDate::from_calendar(year, month, day, hour, minute, second);
    TAI::from_julian_date(jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;
    use astro_core::constants::J2000_JD;

    #[test]
    fn test_tai_constructors() {
        assert_eq!(TAI::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(TAI::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            tai_from_calendar(2000, 1, 1, 12, 0, 0.0).to_julian_date().to_f64(),
            J2000_JD
        );

        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let tai_direct = TAI::from_julian_date(jd);
        let tai_from_trait: TAI = jd.into();
        assert_eq!(tai_direct.to_julian_date().jd1(), tai_from_trait.to_julian_date().jd1());
        assert_eq!(tai_direct.to_julian_date().jd2(), tai_from_trait.to_julian_date().jd2());
    }

    #[test]
    fn test_tai_arithmetic() {
        let tai = TAI::j2000();
        assert_eq!(tai.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(tai.add_seconds(3600.0).to_julian_date().to_f64(), J2000_JD + 1.0 / 24.0);
    }

    #[test]
    fn test_tai_display() {
        let display_str = format!("{}", TAI::from_julian_date(JulianDate::new(J2000_JD, 0.5)));
        assert!(display_str.starts_with("TAI"));
        assert!(display_str.contains("2451545"));
    }

    #[test]
    fn test_tai_string_parsing() {
        assert_eq!(
            TAI::from_str("2000-01-01T12:00:00").unwrap().to_julian_date().to_f64(),
            TAI::j2000().to_julian_date().to_f64()
        );

        let result = TAI::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);

        assert!(TAI::from_str("invalid-date").is_err());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_tai_serde_round_trip() {
        let test_cases = [
            TAI::j2000(),
            TAI::new(0, 0),
            tai_from_calendar(2024, 6, 15, 14, 30, 45.123),
            tai_from_calendar(1990, 12, 31, 23, 59, 59.999999999),
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: TAI = serde_json::from_str(&json).unwrap();

            let total_diff =
                (original.to_julian_date().to_f64() - deserialized.to_julian_date().to_f64()).abs();
            assert!(total_diff < 1e-14, "serde precision loss: {:.2e}", total_diff);
        }
    }
}
