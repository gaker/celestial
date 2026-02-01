//! Universal Time UT1 time scale.
//!
//! UT1 is the principal form of Universal Time, defined by Earth's rotation angle.
//! Unlike atomic time scales (TAI, TT), UT1 tracks the actual rotational position
//! of the Earth, making it essential for astronomical observations and coordinate
//! transformations.
//!
//! # Background
//!
//! Earth's rotation is irregular due to tidal friction, core-mantle coupling, and
//! atmospheric effects. UT1 accumulates an unpredictable offset from UTC (typically
//! |DUT1| < 0.9s). The offset UT1-UTC is published by IERS in Bulletin A/B.
//!
//! ```text
//! UT1 = UTC + DUT1    (where DUT1 from IERS observations)
//! ```
//!
//! # Usage
//!
//! ```
//! use celestial_time::{JulianDate, UT1};
//! use celestial_time::scales::ut1::ut1_from_calendar;
//!
//! // From Unix timestamp components
//! let ut1 = UT1::new(0, 0);  // Unix epoch in UT1
//!
//! // From calendar date
//! let ut1 = ut1_from_calendar(2000, 1, 1, 12, 0, 0.0);
//!
//! // From Julian Date
//! let ut1 = UT1::from_julian_date(JulianDate::j2000());
//! ```
//!
//! # Relationship to Other Scales
//!
//! UT1 is required for:
//! - Sidereal time calculations (GMST, GAST, ERA)
//! - Earth orientation parameters
//! - Topocentric coordinate transformations
//!
//! Conversion from UTC requires external Earth Orientation Parameters (EOP) data.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use celestial_core::constants::MJD_ZERO_POINT;

use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::{NANOSECONDS_PER_SECOND_F64, SECONDS_PER_DAY, SECONDS_PER_DAY_F64};
use std::fmt;
use std::str::FromStr;

/// Universal Time UT1, based on Earth's rotation angle.
///
/// Internally stores time as a split Julian Date for full precision.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct UT1(JulianDate);

impl UT1 {
    /// Creates UT1 from Unix timestamp components.
    ///
    /// Converts seconds and nanoseconds since Unix epoch (1970-01-01 00:00:00)
    /// to a split Julian Date representation.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let days = seconds / SECONDS_PER_DAY;
        let remainder_seconds = seconds % SECONDS_PER_DAY;
        let jd1 = UNIX_EPOCH_JD + days as f64;
        let jd2 = (remainder_seconds as f64 + nanos as f64 / NANOSECONDS_PER_SECOND_F64)
            / SECONDS_PER_DAY_F64;
        Self(JulianDate::new(jd1, jd2))
    }

    /// Creates UT1 from a Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Returns UT1 at the J2000.0 epoch (2000-01-01T12:00:00).
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Adds seconds to this UT1 instant. Negative values subtract.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Adds days to this UT1 instant. Negative values subtract.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }
}

/// Creates UT1 from Gregorian calendar components.
///
/// Uses a proleptic Gregorian calendar algorithm. The date is converted to
/// Modified Julian Date, then to split Julian Date with the time fraction
/// stored separately for precision.
///
/// # Arguments
///
/// * `year` - Gregorian year (negative for BCE)
/// * `month` - Month 1-12
/// * `day` - Day of month 1-31
/// * `hour` - Hour 0-23
/// * `minute` - Minute 0-59
/// * `second` - Second with fractional part
pub fn ut1_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UT1 {
    let my = (month as i32 - 14) / 12;
    let iypmy = year + my;

    let mjd_zero = MJD_ZERO_POINT;

    let modified_jd = ((1461 * (iypmy + 4800)) / 4 + (367 * (month as i32 - 2 - 12 * my)) / 12
        - (3 * ((iypmy + 4900) / 100)) / 4
        + day as i32
        - 2432076) as f64;

    let time_fraction =
        (60.0 * (60 * hour as i32 + minute as i32) as f64 + second) / SECONDS_PER_DAY_F64;
    let jd1 = mjd_zero + modified_jd;
    let jd2 = time_fraction;

    UT1::from_julian_date(JulianDate::new(jd1, jd2))
}

/// Formats as "UT1 {julian_date}".
impl fmt::Display for UT1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "UT1 {}", self.0)
    }
}

/// Converts a Julian Date to UT1.
impl From<JulianDate> for UT1 {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

/// Parses UT1 from an ISO 8601 string.
///
/// Accepts formats like "2000-01-01T12:00:00" or "2000-01-01T12:00:00.123".
impl FromStr for UT1 {
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
    fn test_ut1_constructors() {
        assert_eq!(UT1::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(UT1::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            ut1_from_calendar(2000, 1, 1, 12, 0, 0.0).to_julian_date().to_f64(),
            J2000_JD
        );

        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let ut1_direct = UT1::from_julian_date(jd);
        let ut1_from_trait: UT1 = jd.into();
        assert_eq!(ut1_direct, ut1_from_trait);
    }

    #[test]
    fn test_ut1_arithmetic() {
        let ut1 = UT1::j2000();
        assert_eq!(ut1.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(ut1.add_seconds(3600.0).to_julian_date().to_f64(), J2000_JD + 1.0 / 24.0);
    }

    #[test]
    fn test_ut1_display() {
        let display_str = format!("{}", UT1::from_julian_date(JulianDate::new(J2000_JD, 0.5)));
        assert!(display_str.starts_with("UT1"));
        assert!(display_str.contains("2451545"));
    }

    #[test]
    fn test_ut1_string_parsing() {
        assert_eq!(
            UT1::from_str("2000-01-01T12:00:00").unwrap().to_julian_date().to_f64(),
            UT1::j2000().to_julian_date().to_f64()
        );

        let result = UT1::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);

        assert!(UT1::from_str("invalid-date").is_err());
    }
}
