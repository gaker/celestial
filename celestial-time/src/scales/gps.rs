//! GPS Time scale.
//!
//! GPS Time is the time standard used by GPS satellites. It is synchronized with TAI
//! but offset by exactly 19 seconds: `TAI = GPS + 19s`.
//!
//! # Background
//!
//! GPS Time started on January 6, 1980 at 00:00:00 UTC. At that moment, GPS and UTC
//! were synchronized, and TAI was already 19 seconds ahead of UTC. GPS does not
//! include leap seconds, so the TAI-GPS offset remains constant while UTC-GPS
//! diverges with each new leap second.
//!
//! As of 2024, UTC is 18 seconds behind GPS (37 seconds behind TAI).
//!
//! # Representation
//!
//! Internally stored as a split Julian Date for nanosecond-level precision.
//! See [`JulianDate`] for details on the two-part representation.
//!
//! # Usage
//!
//! ```
//! use celestial_time::{GPS, JulianDate};
//!
//! // From calendar date
//! let gps = celestial_time::scales::gps::gps_from_calendar(2024, 3, 15, 12, 0, 0.0);
//!
//! // From Julian Date
//! let gps = GPS::from_julian_date(JulianDate::j2000());
//!
//! // Arithmetic
//! let later = gps.add_seconds(3600.0);
//! let next_day = gps.add_days(1.0);
//! ```
//!
//! # Conversions
//!
//! GPS converts to/from TAI via a fixed 19-second offset. See
//! [`scales::conversions::gps_tai`](crate::scales::conversions) for the conversion traits.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::SECONDS_PER_DAY_F64;
use std::fmt;
use std::str::FromStr;

/// GPS Time representation.
///
/// Wraps a [`JulianDate`] to provide type safety and GPS-specific operations.
/// The inner Julian Date uses split storage (jd1 + jd2) to preserve precision.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GPS(JulianDate);

impl GPS {
    /// Creates GPS time from Unix timestamp components.
    ///
    /// Converts seconds and nanoseconds since Unix epoch (1970-01-01 00:00:00)
    /// to a Julian Date representation.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let total_seconds =
            seconds as f64 + nanos as f64 / celestial_core::constants::NANOSECONDS_PER_SECOND_F64;
        let jd = JulianDate::from_f64(UNIX_EPOCH_JD + total_seconds / SECONDS_PER_DAY_F64);
        Self(jd)
    }

    /// Creates GPS time from a Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Creates GPS time from raw Julian Date components.
    ///
    /// Prefer this over `from_julian_date(JulianDate::new(...))` when you already
    /// have the split components, as it avoids intermediate allocations.
    pub fn from_julian_date_raw(jd1: f64, jd2: f64) -> Self {
        Self(JulianDate::new(jd1, jd2))
    }

    /// Returns GPS time at J2000.0 epoch (2000-01-01 12:00:00 TT).
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Adds seconds to this GPS time, returning a new instance.
    ///
    /// Precision is preserved by adding to the smaller-magnitude JD component.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Adds days to this GPS time, returning a new instance.
    ///
    /// Precision is preserved by adding to the smaller-magnitude JD component.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }
}

/// Formats as "GPS <julian_date>".
impl fmt::Display for GPS {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GPS {}", self.0)
    }
}

/// Converts a Julian Date to GPS time.
impl From<JulianDate> for GPS {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

/// Parses an ISO 8601 datetime string as GPS time.
///
/// The parsed time is interpreted directly as GPS (no leap second handling).
impl FromStr for GPS {
    type Err = TimeError;

    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

/// Creates GPS time from calendar components.
///
/// Uses direct calendar-to-JD conversion with no leap second corrections.
/// For UTC calendar dates that need leap second handling, parse as UTC first
/// then convert to GPS via TAI.
pub fn gps_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> GPS {
    let jd = JulianDate::from_calendar(year, month, day, hour, minute, second);
    GPS::from_julian_date(jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_gps_constructors() {
        assert_eq!(GPS::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(GPS::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            gps_from_calendar(2000, 1, 1, 12, 0, 0.0)
                .to_julian_date()
                .to_f64(),
            J2000_JD
        );

        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let gps_direct = GPS::from_julian_date(jd);
        let gps_from_trait: GPS = jd.into();
        assert_eq!(
            gps_direct.to_julian_date().jd1(),
            gps_from_trait.to_julian_date().jd1()
        );
        assert_eq!(
            gps_direct.to_julian_date().jd2(),
            gps_from_trait.to_julian_date().jd2()
        );
    }

    #[test]
    fn test_gps_arithmetic() {
        let gps = GPS::j2000();
        assert_eq!(gps.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(
            gps.add_seconds(3600.0).to_julian_date().to_f64(),
            J2000_JD + 1.0 / 24.0
        );
    }

    #[test]
    fn test_gps_display() {
        let display_str = format!("{}", GPS::from_julian_date(JulianDate::new(J2000_JD, 0.5)));
        assert!(display_str.starts_with("GPS"));
        assert!(display_str.contains("2451545"));
    }

    #[test]
    fn test_gps_string_parsing() {
        assert_eq!(
            GPS::from_str("2000-01-01T12:00:00")
                .unwrap()
                .to_julian_date()
                .to_f64(),
            GPS::j2000().to_julian_date().to_f64()
        );

        let result = GPS::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);

        assert!(GPS::from_str("invalid-date").is_err());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_gps_serde_round_trip() {
        let test_cases = [
            GPS::j2000(),
            GPS::new(0, 0),
            gps_from_calendar(2024, 6, 15, 14, 30, 45.123),
            gps_from_calendar(1990, 12, 31, 23, 59, 59.999999999),
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: GPS = serde_json::from_str(&json).unwrap();

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
