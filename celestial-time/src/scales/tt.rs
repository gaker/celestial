//! Terrestrial Time (TT) time scale.
//!
//! TT is the modern successor to Ephemeris Time (ET) and Terrestrial Dynamical Time (TDT).
//! It provides a uniform time scale for geocentric ephemerides and is the basis for
//! planetary position calculations referenced to Earth's center.
//!
//! # Relationship to TAI
//!
//! TT differs from TAI by a fixed offset:
//!
//! ```text
//! TT = TAI + 32.184 seconds
//! ```
//!
//! The 32.184s offset was chosen to maintain continuity with ET at the 1977 epoch.
//!
//! # Usage
//!
//! ```
//! use celestial_time::{JulianDate, TT};
//!
//! // Create TT at J2000.0 epoch
//! let tt = TT::j2000();
//!
//! // From calendar date
//! use celestial_time::scales::tt::tt_from_calendar;
//! let tt = tt_from_calendar(2000, 1, 1, 12, 0, 0.0);
//!
//! // Parse from ISO 8601
//! let tt: TT = "2000-01-01T12:00:00".parse().unwrap();
//!
//! // Julian centuries since J2000.0 (for precession/nutation)
//! let centuries = tt.centuries_since_j2000();
//! ```
//!
//! # Precision
//!
//! TT uses split Julian Date storage internally, preserving microsecond accuracy
//! across the full date range. The `centuries_since_j2000()` method provides
//! the T parameter used in IAU precession and nutation models.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::{J2000_JD, SECONDS_PER_DAY_F64};
use std::fmt;
use std::str::FromStr;

/// Terrestrial Time representation.
///
/// Wraps a split Julian Date for high-precision time storage.
/// TT is the primary time scale for geocentric ephemeris calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TT(JulianDate);

impl TT {
    /// Creates TT from Unix timestamp components.
    ///
    /// Converts seconds and nanoseconds since Unix epoch (1970-01-01T00:00:00)
    /// to TT Julian Date representation.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let total_seconds =
            seconds as f64 + nanos as f64 / celestial_core::constants::NANOSECONDS_PER_SECOND_F64;
        let jd = JulianDate::from_f64(UNIX_EPOCH_JD + total_seconds / SECONDS_PER_DAY_F64);
        Self(jd)
    }

    /// Creates TT from a split Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Creates TT from raw Julian Date components.
    ///
    /// Use when you have separate jd1/jd2 values and want to avoid
    /// intermediate JulianDate construction.
    pub fn from_julian_date_raw(jd1: f64, jd2: f64) -> Self {
        Self(JulianDate::new(jd1, jd2))
    }

    /// Returns TT at J2000.0 epoch (2000-01-01T12:00:00 TT).
    ///
    /// This is the fundamental epoch for modern astronomical calculations.
    /// JD = 2451545.0.
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Returns a new TT offset by the given number of seconds.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Returns a new TT offset by the given number of days.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }

    /// Creates TT from a single-value Julian Date.
    ///
    /// For high-precision work, prefer `from_julian_date` with split values.
    pub fn from_jd(jd: f64) -> TimeResult<Self> {
        Ok(Self(JulianDate::from_f64(jd)))
    }

    /// Returns the Julian year corresponding to this TT instant.
    ///
    /// Julian year = 2000.0 + (JD - J2000_JD) / 365.25
    pub fn julian_year(&self) -> f64 {
        2000.0 + (self.0.to_f64() - J2000_JD) / 365.25
    }

    /// Returns Julian centuries since J2000.0 (the T parameter).
    ///
    /// This is the time argument used in IAU precession and nutation series.
    /// One Julian century = 36525 days.
    pub fn centuries_since_j2000(&self) -> f64 {
        (self.0.to_f64() - J2000_JD) / celestial_core::constants::DAYS_PER_JULIAN_CENTURY
    }
}

/// Formats as "TT {julian_date}".
impl fmt::Display for TT {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TT {}", self.0)
    }
}

/// Converts JulianDate to TT directly.
impl From<JulianDate> for TT {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

/// Parses TT from ISO 8601 format (e.g., "2000-01-01T12:00:00").
///
/// Assumes the input string represents a TT instant directly.
impl FromStr for TT {
    type Err = TimeError;

    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

/// Creates TT from calendar components.
///
/// Converts Gregorian calendar date and time to TT. The input is interpreted
/// directly as TT with no UTC or leap second corrections applied.
///
/// # Arguments
///
/// * `year` - Gregorian year (negative for BCE)
/// * `month` - Month (1-12)
/// * `day` - Day of month (1-31)
/// * `hour` - Hour (0-23)
/// * `minute` - Minute (0-59)
/// * `second` - Second with fractional part (0.0 to <61.0)
pub fn tt_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> TT {
    let jd = JulianDate::from_calendar(year, month, day, hour, minute, second);
    TT::from_julian_date(jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;

    #[test]
    fn test_tt_constructors() {
        assert_eq!(TT::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(TT::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            tt_from_calendar(2000, 1, 1, 12, 0, 0.0)
                .to_julian_date()
                .to_f64(),
            J2000_JD
        );
        assert_eq!(
            TT::from_jd(J2000_JD).unwrap().to_julian_date().to_f64(),
            J2000_JD
        );
    }

    #[test]
    fn test_tt_from_julian_date_raw() {
        let tt = TT::from_julian_date_raw(J2000_JD, 0.5);
        assert_eq!(tt.to_julian_date().jd1(), J2000_JD);
        assert_eq!(tt.to_julian_date().jd2(), 0.5);
    }

    #[test]
    fn test_tt_arithmetic() {
        let tt = TT::j2000();
        assert_eq!(tt.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(
            tt.add_seconds(3600.0).to_julian_date().to_f64(),
            J2000_JD + 1.0 / 24.0
        );
    }

    #[test]
    fn test_tt_julian_year_and_centuries() {
        let tt = TT::j2000();
        assert_eq!(tt.julian_year(), 2000.0);
        assert_eq!(tt.centuries_since_j2000(), 0.0);

        let tt_plus_century = tt.add_days(celestial_core::constants::DAYS_PER_JULIAN_CENTURY);
        assert_eq!(tt_plus_century.centuries_since_j2000(), 1.0);
    }

    #[test]
    fn test_tt_from_julian_date_trait() {
        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let tt_direct = TT::from_julian_date(jd);
        let tt_from_trait: TT = jd.into();

        assert_eq!(
            tt_direct.to_julian_date().jd1(),
            tt_from_trait.to_julian_date().jd1()
        );
        assert_eq!(
            tt_direct.to_julian_date().jd2(),
            tt_from_trait.to_julian_date().jd2()
        );
    }

    #[test]
    fn test_tt_display() {
        let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let display_str = format!("{}", tt);
        assert!(display_str.starts_with("TT"));
        assert!(display_str.contains("2451545"));
    }

    #[test]
    fn test_tt_string_parsing() {
        assert_eq!(
            TT::from_str("2000-01-01T12:00:00")
                .unwrap()
                .to_julian_date()
                .to_f64(),
            TT::j2000().to_julian_date().to_f64()
        );
        assert!(TT::from_str("invalid-date").is_err());
    }

    #[test]
    fn test_tt_string_parsing_fractional_seconds() {
        let result = TT::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected = tt_from_calendar(2000, 1, 1, 12, 0, 0.123);
        assert_eq!(
            result.to_julian_date().to_f64(),
            expected.to_julian_date().to_f64()
        );
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_tt_serde_round_trip() {
        let test_cases = [
            TT::j2000(),
            TT::new(0, 0),
            tt_from_calendar(2024, 6, 15, 14, 30, 45.123),
            tt_from_calendar(1990, 12, 31, 23, 59, 59.999999999),
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: TT = serde_json::from_str(&json).unwrap();

            let total_diff =
                (original.to_julian_date().to_f64() - deserialized.to_julian_date().to_f64()).abs();
            assert!(
                total_diff < 1e-14,
                "serde precision loss: {:.2e}",
                total_diff
            );
        }
    }
}
