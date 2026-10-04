//! Conversions between Coordinated Universal Time (UTC) and International Atomic Time (TAI).
//!
//! UTC and TAI are both atomic time scales, but UTC includes leap seconds to stay within
//! 0.9 seconds of UT1 (Earth rotation time). TAI runs continuously without adjustments.
//!
//! # The UTC-TAI Relationship
//!
//! TAI is always ahead of UTC by an integer number of seconds (since 1972). The offset
//! started at 10 seconds on 1972-01-01 and has grown to 37 seconds as of 2017-01-01:
//!
//! ```text
//! TAI = UTC + (leap seconds accumulated)
//! ```
//!
//! Before 1972, the relationship was more complex, involving both step offsets and
//! continuous drift corrections.
//!
//! # Leap Seconds
//!
//! Leap seconds keep UTC synchronized with Earth's rotation:
//!
//! - The IERS monitors the difference between UT1 (Earth rotation) and UTC
//! - When |UT1 - UTC| approaches 0.9 seconds, a leap second is announced
//! - Leap seconds are inserted at the end of June 30 or December 31
//! - Only positive leap seconds have occurred (Earth rotation is slowing)
//! - Since 1972, leap seconds have been exactly 1 second adjustments
//!
//! The leap second table (`TAI_UTC_OFFSETS` in constants.rs) records all adjustments:
//!
//! | Date       | TAI-UTC (seconds) |
//! |------------|-------------------|
//! | 1972-01-01 | 10.0              |
//! | 1972-07-01 | 11.0              |
//! | ...        | ...               |
//! | 2017-01-01 | 37.0              |
//!
//! # Pre-1972 Handling
//!
//! Before the modern leap second system, UTC used a different adjustment model:
//!
//! - Step offsets at irregular intervals (not exactly 1 second)
//! - Continuous drift corrections between steps
//! - The `UTC_DRIFT_CORRECTIONS` table provides (MJD reference, drift rate) pairs
//!
//! For pre-1972 dates, the offset is computed as:
//!
//! ```text
//! TAI - UTC = base_offset + (MJD - reference_MJD) * drift_rate
//! ```
//!
//! The first 14 entries in `TAI_UTC_OFFSETS` (indices 0-13) use this drift model.
//!
//! # TAI to UTC Conversion Algorithm
//!
//! Converting TAI to UTC requires finding which UTC day corresponds to a given TAI instant.
//! This is non-trivial because leap seconds create discontinuities. The algorithm uses
//! iterative refinement:
//!
//! 1. Start with a UTC guess equal to TAI
//! 2. Convert the UTC guess to TAI
//! 3. Compute the difference from the target TAI
//! 4. Adjust the UTC guess by this difference
//! 5. Repeat for 3 iterations (lands within 20 ps)
//!
//! From 1972 on, TAI-UTC is a whole number of seconds and one step suffices. Before
//! 1972 the drift term makes each step shrink the error by a factor of ~1e-8.
//!
//! # UTC to TAI Conversion Algorithm
//!
//! The forward conversion (UTC to TAI) handles leap seconds and drift corrections:
//!
//! 1. Convert Julian Date to calendar date (year, month, day, fraction)
//! 2. Look up the TAI-UTC offset at the start of the day (0h)
//! 3. Look up the offset at mid-day (12h) to detect drift (pre-1972)
//! 4. Look up the offset at the start of the next day to detect leap seconds
//! 5. Apply drift and leap second corrections to the day fraction
//! 6. Add the base offset to get TAI
//!
//! The three-point lookup (0h, 12h, next day 0h) correctly handles both:
//! - Pre-1972 linear drift (detected by 0h vs 12h difference)
//! - Leap seconds (detected by comparing end of day to start of next day)
//!
//! # Precision
//!
//! Round-trip conversions (UTC -> TAI -> UTC or TAI -> UTC -> TAI) come back within
//! 20 ps. The worst case is just after 0h, where the other scale falls on the previous
//! day and its fraction, just under 1, is rounded twice. The iterative refinement and
//! careful handling of Julian Date components preserve floating-point precision.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::utc::UTC;
//! use celestial_time::scales::conversions::{ToTAI, ToUTC};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! // At 2000-01-01 12:00 UTC (JD 2451545.0 in UTC), TAI-UTC = 32 seconds
//! let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tai = utc.to_tai().unwrap();
//!
//! // 32 s as a fraction of a day, rounded as ERFA's eraUtctai rounds it
//! assert_eq!(tai.to_julian_date(), JulianDate::new(J2000_JD, 0.00037037037037035425));
//! assert_eq!(tai.to_utc().unwrap(), utc);
//! ```
//!
//! # Helper Functions
//!
//! This module also provides calendar conversion utilities used by the UTC-TAI algorithms:
//!
//! - [`julian_to_calendar`]: Convert Julian Date to (year, month, day, day_fraction)
//! - [`calendar_to_julian`]: Convert (year, month, day) to Julian Date
//!
//! These functions handle the full range of historical dates. `julian_to_calendar` uses
//! compensated summation (Kahan algorithm) to preserve precision when combining Julian
//! Date components.
//!
//! # References
//!
//! - IERS Bulletins: Leap second announcements
//! - USNO: History of leap seconds and TAI-UTC differences
//! - ITU-R TF.460-6: Standard-frequency and time-signal emissions
//! - Explanatory Supplement to the Astronomical Almanac, 3rd ed., Chapter 3

mod calendar;

use super::super::common::{get_tai_utc_offset, next_calendar_day, validate_calendar_date};
use super::{ToTAI, ToUTC};
use crate::julian::JulianDate;
use crate::scales::tai::TAI;
use crate::scales::utc::UTC;
use crate::TimeResult;
use calendar::{check_julian_range, day_number_and_fraction, gregorian_date};
use celestial_core::constants::{MJD_ZERO_POINT, SECONDS_PER_DAY_F64};

impl ToTAI for UTC {
    /// Convert UTC to TAI by adding the accumulated leap seconds.
    ///
    /// Looks up the TAI-UTC offset for the given date and applies drift corrections
    /// for pre-1972 dates. Handles leap second boundaries correctly.
    fn to_tai(&self) -> TimeResult<TAI> {
        utc_to_tai(self.to_julian_date())
    }
}

impl ToUTC for TAI {
    /// Convert TAI to UTC using iterative refinement.
    ///
    /// Uses 3 iterations to converge on the correct UTC instant, handling
    /// leap second boundaries where a single TAI instant may map to the
    /// leap second itself.
    fn to_utc(&self) -> TimeResult<UTC> {
        tai_to_utc(self.to_julian_date())
    }
}

/// Convert a UTC Julian Date to TAI.
///
/// This function handles both the modern leap second era (1972+) and the pre-1972
/// drift correction era. The algorithm:
///
/// 1. Orders jd1 and jd2 by magnitude for precision preservation.
///
/// 2. Converts to calendar date to look up the appropriate TAI-UTC offset.
///
/// 3. Samples the offset at three points to detect drift and leap seconds:
///    - Start of day (0h): base offset
///    - Mid-day (12h): detects linear drift (pre-1972)
///    - Start of next day: detects leap seconds
///
/// 4. Computes drift rate from the 0h/12h difference (zero for post-1972 dates).
///
/// 5. Computes leap second amount from the end-of-day discontinuity.
///
/// 6. Scales the day fraction to account for drift and leap seconds, then adds
///    the base offset.
///
/// The correction is applied to the smaller-magnitude JD component to preserve
/// floating-point precision.
fn utc_to_tai(utc_jd: JulianDate) -> TimeResult<TAI> {
    let (utc_int, utc_frac, big1) = if libm::fabs(utc_jd.jd1()) >= libm::fabs(utc_jd.jd2()) {
        (utc_jd.jd1(), utc_jd.jd2(), true)
    } else {
        (utc_jd.jd2(), utc_jd.jd1(), false)
    };

    let (year, month, day, day_fraction) = julian_to_calendar(utc_int, utc_frac)?;
    let (offset_0h, day_fraction) = stretch_day_fraction(year, month, day, day_fraction)?;

    let (z1, z2) = calendar_to_julian(year, month, day)?;

    let mut tai_frac = z1 - utc_int;
    tai_frac += z2;
    tai_frac += day_fraction + offset_0h / SECONDS_PER_DAY_F64;

    let (tai_jd1, tai_jd2) = if big1 {
        (utc_int, tai_frac)
    } else {
        (tai_frac, utc_int)
    };

    Ok(TAI::from_julian_date(JulianDate::new(tai_jd1, tai_jd2)))
}

// Returns TAI-UTC at 0h and the day fraction rescaled to the day's length:
// 86,400 s plus the pre-1972 drift across the day plus any leap second at its end.
fn stretch_day_fraction(
    year: i32,
    month: i32,
    day: i32,
    day_fraction: f64,
) -> TimeResult<(f64, f64)> {
    let offset_0h = get_tai_utc_offset(year, month, day, 0.0)?;
    let offset_12h = get_tai_utc_offset(year, month, day, 0.5)?;
    let (next_year, next_month, next_day) = next_calendar_day(year, month, day)?;
    let offset_24h = get_tai_utc_offset(next_year, next_month, next_day, 0.0)?;

    let drift_rate = 2.0 * (offset_12h - offset_0h);
    let leap_amount = offset_24h - (offset_0h + drift_rate);

    let mut stretched = day_fraction * ((SECONDS_PER_DAY_F64 + leap_amount) / SECONDS_PER_DAY_F64);
    stretched *= (SECONDS_PER_DAY_F64 + drift_rate) / SECONDS_PER_DAY_F64;
    Ok((offset_0h, stretched))
}

/// Convert a TAI Julian Date to UTC using iterative refinement.
///
/// The inverse conversion (TAI to UTC) cannot be done with a simple table lookup
/// because leap seconds create discontinuities: during a leap second, UTC stays
/// at 23:59:60 while TAI advances. The algorithm uses Newton-like iteration:
///
/// 1. Initialize UTC guess = TAI (close enough to converge quickly).
///
/// 2. For each iteration:
///    - Convert the UTC guess to TAI using `utc_to_tai`
///    - Compute the residual: target_TAI - computed_TAI
///    - Add the residual to the UTC guess
///
/// 3. After 3 iterations, the UTC value converts back to within 20 ps of the target.
///
/// The iteration count of 3 is sufficient because:
/// - The initial guess is within ~37 seconds of the answer
/// - From 1972 on the offset is a whole number of seconds, so one step lands
///   within 1 ulp
/// - Before 1972 the drift term shrinks the error by a factor of ~1e-8 per step
///
/// The algorithm correctly handles leap second boundaries where the UTC day
/// "stretches" to include 86401 seconds.
fn tai_to_utc(tai_jd: JulianDate) -> TimeResult<UTC> {
    let (tai_int, tai_frac, big1) = if libm::fabs(tai_jd.jd1()) >= libm::fabs(tai_jd.jd2()) {
        (tai_jd.jd1(), tai_jd.jd2(), true)
    } else {
        (tai_jd.jd2(), tai_jd.jd1(), false)
    };

    let utc_int = tai_int;
    let utc_frac = refine_utc_fraction(tai_int, tai_frac)?;

    let (utc_jd1, utc_jd2) = if big1 {
        (utc_int, utc_frac)
    } else {
        (utc_frac, utc_int)
    };

    Ok(UTC::from_julian_date(JulianDate::new(utc_jd1, utc_jd2)))
}

// The UTC guess starts at TAI and keeps TAI's larger part, so only the smaller
// part moves.
fn refine_utc_fraction(tai_int: f64, tai_frac: f64) -> TimeResult<f64> {
    const TAI_TO_UTC_ITERATIONS: usize = 3;
    let mut utc_frac = tai_frac;
    for _ in 0..TAI_TO_UTC_ITERATIONS {
        let guess_tai = utc_to_tai_jd(tai_int, utc_frac)?;
        utc_frac += tai_int - guess_tai.jd1();
        utc_frac += tai_frac - guess_tai.jd2();
    }
    Ok(utc_frac)
}

/// Helper for iterative TAI->UTC conversion.
///
/// Wraps `utc_to_tai` to work with separate JD components, preserving the
/// split-JD precision during iteration.
fn utc_to_tai_jd(utc_int: f64, utc_frac: f64) -> TimeResult<JulianDate> {
    let utc = UTC::from_julian_date(JulianDate::new(utc_int, utc_frac));
    let tai = utc.to_tai()?;
    Ok(tai.to_julian_date())
}

/// Convert a two-part Julian Date to calendar date with day fraction.
///
/// Returns `(year, month, day, day_fraction)` where:
/// - `year`: Astronomical year (0 is 1 BCE, -1 is 2 BCE)
/// - `month`: 1-12
/// - `day`: 1-31
/// - `day_fraction`: 0.0 to 1.0 (fraction of day from midnight)
///
/// # Algorithm
///
/// The conversion uses integer arithmetic for the calendar calculation (avoiding
/// floating-point error accumulation) and Kahan compensated summation for the
/// fractional day.
///
/// Steps:
/// 1. Round each JD component to the nearest integer, keeping the fractional parts
/// 2. Sum the fractional parts using Kahan summation (adds 0.5 to shift from noon to midnight)
/// 3. Handle edge cases where the fraction overflows [0, 1)
/// 4. Apply the standard algorithm for Julian Day Number to Gregorian calendar
///
/// The compensated summation is critical: without it, adding two fractional parts
/// can lose precision when they have opposite signs or very different magnitudes.
///
/// # Valid Range
///
/// - Minimum: JD -68569.5 (-4900-03-01, which is 4901 BCE)
/// - Maximum: JD 1e9 (far future)
///
/// # Errors
///
/// Returns `TimeError::ConversionError` if the Julian Date is outside the valid range.
///
/// # Example
///
/// ```
/// use celestial_time::scales::conversions::utc_tai::julian_to_calendar;
///
/// let (year, month, day, frac) = julian_to_calendar(2451545.0, 0.0).unwrap();
/// assert_eq!((year, month, day, frac), (2000, 1, 1, 0.5)); // J2000.0 is noon
/// ```
pub fn julian_to_calendar(jd1: f64, jd2: f64) -> TimeResult<(i32, i32, i32, f64)> {
    check_julian_range(jd1 + jd2)?;
    let (jd, fraction) = day_number_and_fraction(jd1, jd2);
    let (year, month, day) = gregorian_date(jd);
    Ok((year, month, day, fraction))
}

/// Convert a calendar date to a two-part Julian Date.
///
/// Returns `(jd1, jd2)` where:
/// - `jd1` = MJD zero point (2400000.5)
/// - `jd2` = Modified Julian Date for the given calendar date at midnight
///
/// This split preserves precision: the large constant is in jd1, and the
/// smaller date-dependent value is in jd2.
///
/// # Algorithm
///
/// Uses the standard formula for Gregorian calendar to Julian Day Number,
/// expressed as a Modified Julian Date for precision.
///
/// # Example
///
/// ```
/// use celestial_time::scales::conversions::utc_tai::calendar_to_julian;
///
/// // Midnight at the start of 2000-01-01, JD 2451544.5
/// assert_eq!(calendar_to_julian(2000, 1, 1).unwrap(), (2400000.5, 51544.0));
/// ```
pub fn calendar_to_julian(year: i32, month: i32, day: i32) -> TimeResult<(f64, f64)> {
    validate_calendar_date(year, month, day)?;
    let (year, month, day) = (i64::from(year), i64::from(month), i64::from(day));
    let my = (month - 14) / 12;
    let iypmy = year + my;

    let modified_jd = (1461 * (iypmy + 4800)) / 4 + (367 * (month - 2 - 12 * my)) / 12
        - (3 * ((iypmy + 4900) / 100)) / 4
        + day
        - 2432076;

    Ok((MJD_ZERO_POINT, modified_jd as f64))
}

#[cfg(test)]
mod tests;
