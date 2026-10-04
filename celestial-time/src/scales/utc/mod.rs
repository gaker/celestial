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
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::utc::UTC;
//! use celestial_time::scales::utc::utc_from_calendar;
//!
//! // From Unix timestamp
//! let utc = UTC::new(1704067200, 0).unwrap(); // 2024-01-01 00:00:00 UTC
//!
//! // From calendar components
//! let utc = utc_from_calendar(2024, 1, 1, 12, 30, 45.5).unwrap();
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
use super::conversions::utc_tai::julian_to_calendar;
use crate::constants::UNIX_EPOCH_JD;
use crate::julian::{day_fraction, JulianDate};
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::{
    NANOSECONDS_PER_SECOND, NANOSECONDS_PER_SECOND_F64, SECONDS_PER_DAY, SECONDS_PER_DAY_F64,
};
use std::str::FromStr;
use std::time::{SystemTime, UNIX_EPOCH};

time_scale! {
    /// UTC time scale backed by a split Julian Date.
    ///
    /// Wraps `JulianDate` to represent Coordinated Universal Time. Supports
    /// construction from Unix timestamps, calendar components, or raw Julian Dates.
    UTC,
    j2000_note = "This is 64.184 s after the J2000.0 epoch, which is defined in TT."
}

impl UTC {
    /// Creates UTC from Unix timestamp (seconds and nanoseconds since 1970-01-01 00:00:00).
    ///
    /// The result is the same as `utc_from_calendar` on the Unix time's date and
    /// time of day, so a day that ends in a leap second is 86,401 s long. Unix time
    /// has no value for the leap second itself (23:59:60); use `utc_from_calendar`
    /// or `FromStr` for that instant.
    ///
    /// # Errors
    ///
    /// Returns an error if `nanos` is 1,000,000,000 or more, or the date is outside
    /// the calendar's range.
    pub fn new(seconds: i64, nanos: u32) -> TimeResult<Self> {
        if nanos >= NANOSECONDS_PER_SECOND {
            return Err(TimeError::ConversionError(format!(
                "nanos must be less than {}, got {}",
                NANOSECONDS_PER_SECOND, nanos
            )));
        }
        let days = seconds.div_euclid(SECONDS_PER_DAY);
        let (year, month, day, _) = julian_to_calendar(UNIX_EPOCH_JD, days as f64)?;
        let second_of_day = seconds.rem_euclid(SECONDS_PER_DAY);
        let hour = (second_of_day / 3600) as u8;
        let minute = (second_of_day % 3600 / 60) as u8;
        let second = (second_of_day % 60) as f64 + f64::from(nanos) / NANOSECONDS_PER_SECOND_F64;
        utc_from_calendar(year, month as u8, day as u8, hour, minute, second)
    }

    /// Returns the current UTC time from the system clock.
    pub fn now() -> TimeResult<Self> {
        Self::from_system_time(SystemTime::now())
    }

    // A clock set before 1970 comes back from duration_since as an Err that
    // holds the distance to the epoch.
    fn from_system_time(time: SystemTime) -> TimeResult<Self> {
        let nanos = match time.duration_since(UNIX_EPOCH) {
            Ok(since) => since.as_nanos() as i128,
            Err(before) => -(before.duration().as_nanos() as i128),
        };
        let per_second = i128::from(NANOSECONDS_PER_SECOND);
        // SystemTime keeps its whole seconds in an i64, so the quotient fits.
        let seconds = nanos.div_euclid(per_second) as i64;
        Self::new(seconds, nanos.rem_euclid(per_second) as u32)
    }

    /// Formats as ISO 8601 string (YYYY-MM-DDTHH:MM:SS.sss), rounded to the
    /// millisecond. During a leap second the seconds field reads 60.
    ///
    /// Returns an error if the date is outside the calendar's range.
    pub fn to_iso8601(&self) -> TimeResult<String> {
        let jd = self.to_julian_date();
        let (year, month, day, fraction) = julian_to_calendar(jd.jd1(), jd.jd2())?;
        // Only a whole leap second stretches the day here. The pre-1972 steps
        // are formatted on an 86,400 s day, as eraD2dtf does.
        let leap = leap_at_end_of_day(year, month, day)?;
        let leap = if libm::fabs(leap) > 0.5 { leap } else { 0.0 };
        let seconds = SECONDS_PER_DAY_F64 * (fraction + fraction * leap / SECONDS_PER_DAY_F64);
        let millis = libm::round(1000.0 * seconds) as i64;

        let day_end = if leap > 0.0 {
            MS_PER_DAY + 1000
        } else {
            MS_PER_DAY
        };
        if millis < day_end {
            return Ok(iso_string(year, month, day, millis));
        }
        let (year, month, day) = next_calendar_day(year, month, day)?;
        Ok(iso_string(year, month, day, 0))
    }
}

const MS_PER_MINUTE: i64 = 60_000;
const MS_PER_DAY: i64 = 86_400_000;
const LAST_MINUTE_OF_DAY: i64 = 24 * 60 - 1;

// A leap second belongs to the day's last minute, which then runs to 61 s.
fn iso_string(year: i32, month: i32, day: i32, millis: i64) -> String {
    let minute_of_day = (millis / MS_PER_MINUTE).min(LAST_MINUTE_OF_DAY);
    let millis_of_minute = millis - minute_of_day * MS_PER_MINUTE;
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}.{:03}",
        year,
        month,
        day,
        minute_of_day / 60,
        minute_of_day % 60,
        millis_of_minute / 1000,
        millis_of_minute % 1000
    )
}

// How much longer than 86,400 s the UTC day is: a leap second, or a pre-1972
// step. eraDtf2d and eraD2dtf both take TAI-UTC in this order.
fn leap_at_end_of_day(year: i32, month: i32, day: i32) -> TimeResult<f64> {
    let dat0 = get_tai_utc_offset(year, month, day, 0.0)?;
    let dat12 = get_tai_utc_offset(year, month, day, 0.5)?;
    let (next_year, next_month, next_day) = next_calendar_day(year, month, day)?;
    let dat24 = get_tai_utc_offset(next_year, next_month, next_day, 0.0)?;
    Ok(dat24 - (2.0 * dat12 - dat0))
}

/// Creates UTC from calendar components, handling leap seconds.
///
/// Computes the TAI-UTC offset at the start, middle, and end of the day
/// to detect leap second insertions. If a leap second occurs, the day
/// is treated as 86401 seconds instead of 86400.
///
/// # Errors
///
/// Returns an error if the date doesn't exist, or the time of day is out of
/// range. A seconds value of 60 or more is accepted only in the last minute of
/// a day that ends with a leap second.
pub fn utc_from_calendar(
    year: i32,
    month: u8,
    day: u8,
    hour: u8,
    minute: u8,
    second: f64,
) -> TimeResult<UTC> {
    let base_jd = JulianDate::from_calendar(year, month, day, 0, 0, 0.0)?;
    let dleap = leap_at_end_of_day(year, month.into(), day.into())?;
    let time_fraction = day_fraction(hour, minute, second, dleap)?;
    Ok(UTC::from_julian_date(JulianDate::new(
        base_jd.jd1(),
        base_jd.jd2() + time_fraction,
    )))
}

/// Parses ISO 8601 formatted strings into UTC.
impl FromStr for UTC {
    type Err = TimeError;

    fn from_str(s: &str) -> TimeResult<Self> {
        let s = s.trim();
        let parsed = parse_iso8601(s.strip_suffix('Z').unwrap_or(s))?;
        utc_from_calendar(
            parsed.year,
            parsed.month,
            parsed.day,
            parsed.hour,
            parsed.minute,
            parsed.second,
        )
    }
}

#[cfg(test)]
mod tests;
