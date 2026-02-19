//! Coordinated Universal Time (UTC) representation.
//!
//! UTC is the primary civil time standard. It tracks TAI but is adjusted with leap seconds
//! to stay within 0.9 seconds of UT1 (Earth rotation time). This module provides the UTC
//! time scale type and calendar-based construction.
//!
//! # Background
//!
//! UTC was introduced in 1960 and has used its current leap second system since 1972.
//! The offset TAI-UTC grows by 1 second each time a leap second is inserted (typically
//! June 30 or December 31 at 23:59:60 UTC). As of 2024, TAI-UTC = 37 seconds.
//!
//! ```text
//! TAI = UTC + (TAI-UTC offset from leap second table)
//! UTC day length = 86400s (normal) or 86401s (positive leap second)
//! ```
//!
//! # Usage
//!
//! ```
//! use celestial_time::{JulianDate, UTC};
//! use celestial_time::scales::utc::utc_from_calendar;
//!
//! // From Unix timestamp
//! let utc = UTC::new(1704067200, 0); // 2024-01-01 00:00:00 UTC
//!
//! // From calendar components
//! let utc = utc_from_calendar(2024, 1, 1, 12, 30, 45.5);
//!
//! // From Julian Date
//! let utc = UTC::from_julian_date(JulianDate::j2000());
//! ```
//!
//! # Leap Second Handling
//!
//! The `utc_from_calendar` function adjusts day length when a leap second occurs.
//! It queries the TAI-UTC offset at multiple points within the day to detect
//! the discontinuity and scales the time fraction accordingly.
//!
//! # Precision
//!
//! Internally stores time as a split Julian Date for nanosecond-level precision.
//! The `new()` constructor separates days from sub-day time to preserve all
//! significant digits in the fractional portion.

use super::common::{get_tai_utc_offset, next_calendar_day};
use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::SECONDS_PER_DAY_F64;
use std::fmt;
use std::str::FromStr;

/// UTC time scale backed by a split Julian Date.
///
/// Wraps `JulianDate` to represent Coordinated Universal Time. Supports
/// construction from Unix timestamps, calendar components, or raw Julian Dates.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct UTC(JulianDate);

impl UTC {
    /// Creates UTC from Unix timestamp (seconds and nanoseconds since 1970-01-01 00:00:00).
    ///
    /// Days are computed separately from sub-day time to preserve precision.
    /// The resulting Julian Date uses jd1 for whole days and jd2 for the fractional part.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let days = seconds / celestial_core::constants::SECONDS_PER_DAY;
        let remainder_seconds = seconds % celestial_core::constants::SECONDS_PER_DAY;
        let jd1 = UNIX_EPOCH_JD + days as f64;
        let jd2 = (remainder_seconds as f64
            + nanos as f64 / celestial_core::constants::NANOSECONDS_PER_SECOND_F64)
            / SECONDS_PER_DAY_F64;
        Self(JulianDate::new(jd1, jd2))
    }

    /// Creates UTC from a Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Returns UTC at the J2000.0 epoch (2000-01-01 12:00:00 TT, JD 2451545.0).
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Returns a new UTC offset by the given seconds.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Returns a new UTC offset by the given days.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }

    /// Returns the current UTC time from the system clock.
    pub fn now() -> Self {
        use std::time::{SystemTime, UNIX_EPOCH};
        let duration = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default();
        Self::new(duration.as_secs() as i64, duration.subsec_nanos())
    }

    /// Formats as ISO 8601 string (YYYY-MM-DDTHH:MM:SS.sss).
    ///
    /// Falls back to "JD{value}" if calendar conversion fails.
    pub fn to_iso8601(&self) -> String {
        use crate::scales::conversions::utc_tai::julian_to_calendar;
        let jd = self.to_julian_date();
        if let Ok((year, month, day, frac)) = julian_to_calendar(jd.jd1(), jd.jd2()) {
            let total_seconds = frac * SECONDS_PER_DAY_F64;
            let hour = (total_seconds / 3600.0) as u8;
            let minute = ((total_seconds % 3600.0) / 60.0) as u8;
            let second = total_seconds % 60.0;
            format!(
                "{:04}-{:02}-{:02}T{:02}:{:02}:{:06.3}",
                year, month, day, hour, minute, second
            )
        } else {
            format!("JD{:.6}", jd.jd1() + jd.jd2())
        }
    }
}

/// Creates UTC from calendar components, handling leap seconds.
///
/// Computes the TAI-UTC offset at the start, middle, and end of the day
/// to detect leap second insertions. If a leap second occurs, the day
/// is treated as 86401 seconds instead of 86400.
///
/// # Panics
///
/// Panics if the month is invalid (not 1-12).
pub fn utc_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UTC {
    let base_jd = JulianDate::from_calendar(year, month, day, 0, 0, 0.0);

    let mut day_length = SECONDS_PER_DAY_F64;

    let dat0 = get_tai_utc_offset(year, month as i32, day as i32, 0.0);
    let dat12 = get_tai_utc_offset(year, month as i32, day as i32, 0.5);

    let (next_year, next_month, next_day) = next_calendar_day(year, month as i32, day as i32)
        .expect("Invalid month in UTC calendar conversion");
    let dat24 = get_tai_utc_offset(next_year, next_month, next_day, 0.0);

    let dleap = dat24 - (2.0 * dat12 - dat0);
    day_length += dleap;

    let time_fraction = (60.0 * (60 * hour as i32 + minute as i32) as f64 + second) / day_length;

    UTC::from_julian_date(JulianDate::new(
        base_jd.jd1(),
        base_jd.jd2() + time_fraction,
    ))
}

/// Displays as "UTC {julian_date}".
impl fmt::Display for UTC {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "UTC {}", self.0)
    }
}

/// Converts JulianDate to UTC.
impl From<JulianDate> for UTC {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

/// Parses ISO 8601 formatted strings into UTC.
impl FromStr for UTC {
    type Err = TimeError;

    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_utc_constructors() {
        assert_eq!(UTC::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(UTC::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            utc_from_calendar(2000, 1, 1, 12, 0, 0.0)
                .to_julian_date()
                .to_f64(),
            J2000_JD
        );

        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let utc_direct = UTC::from_julian_date(jd);
        let utc_from_trait: UTC = jd.into();
        assert_eq!(utc_direct, utc_from_trait);
    }

    #[test]
    fn test_utc_arithmetic() {
        let utc = UTC::j2000();
        assert_eq!(utc.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(
            utc.add_seconds(3600.0).to_julian_date().to_f64(),
            J2000_JD + 1.0 / 24.0
        );
    }

    #[test]
    fn test_utc_display() {
        let display_str = format!("{}", UTC::from_julian_date(JulianDate::new(J2000_JD, 0.5)));
        assert!(display_str.starts_with("UTC"));
        assert!(display_str.contains("2451545"));
    }

    #[test]
    fn test_utc_string_parsing() {
        assert_eq!(
            UTC::from_str("2000-01-01T12:00:00")
                .unwrap()
                .to_julian_date()
                .to_f64(),
            J2000_JD
        );

        let result = UTC::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);

        assert!(UTC::from_str("invalid-date").is_err());
    }

    #[test]
    fn test_utc_new_precision_preservation() {
        let seconds_50_years = 50 * 365 * celestial_core::constants::SECONDS_PER_DAY as u32;
        let nanos = 123456789u32;

        let utc = UTC::new(seconds_50_years as i64, nanos);
        let jd = utc.to_julian_date();

        let expected_days = seconds_50_years / celestial_core::constants::SECONDS_PER_DAY as u32;
        let remainder_secs = seconds_50_years % celestial_core::constants::SECONDS_PER_DAY as u32;
        let expected_jd1 = UNIX_EPOCH_JD + expected_days as f64;
        let expected_jd2 = (remainder_secs as f64
            + nanos as f64 / celestial_core::constants::NANOSECONDS_PER_SECOND_F64)
            / SECONDS_PER_DAY_F64;

        assert_eq!(jd.jd1(), expected_jd1);
        assert_eq!(jd.jd2(), expected_jd2);
    }

    #[test]
    fn test_tai_utc_offset_edge_cases() {
        assert_eq!(get_tai_utc_offset(2000, 1, 1, -0.5), 0.0);
        assert_eq!(get_tai_utc_offset(2000, 1, 1, 1.5), 0.0);
        assert_eq!(get_tai_utc_offset(1950, 6, 15, 0.5), 0.0);
        assert!(get_tai_utc_offset(1960, 1, 1, 0.0) > 0.0);
    }

    #[test]
    fn test_next_calendar_day() {
        assert!(next_calendar_day(2000, 13, 15).is_err());

        let cases: &[(i32, i32, i32, (i32, i32, i32))] = &[
            (2000, 2, 28, (2000, 2, 29)),
            (1999, 2, 28, (1999, 3, 1)),
            (2000, 4, 30, (2000, 5, 1)),
            (2000, 12, 31, (2001, 1, 1)),
        ];

        for &(y, m, d, expected) in cases {
            assert_eq!(next_calendar_day(y, m, d).unwrap(), expected);
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_utc_serde_round_trip() {
        let test_cases = [
            UTC::j2000(),
            UTC::new(0, 0),
            utc_from_calendar(2024, 6, 15, 14, 30, 45.123),
            utc_from_calendar(1990, 12, 31, 23, 59, 59.0),
            utc_from_calendar(2015, 6, 30, 23, 59, 59.999),
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: UTC = serde_json::from_str(&json).unwrap();
            assert_eq!(original, deserialized);
        }
    }
}
