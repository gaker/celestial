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
//! 5. Repeat for 3 iterations (converges to sub-picosecond accuracy)
//!
//! Three iterations suffice because the leap second table lookup is stable once
//! we're within a few seconds of the correct UTC.
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
//! Round-trip conversions (UTC -> TAI -> UTC or TAI -> UTC -> TAI) achieve ~1 picosecond
//! accuracy. The iterative refinement and careful handling of Julian Date components
//! preserve floating-point precision.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::{UTC, TAI};
//! use celestial_time::scales::conversions::{ToTAI, ToUTC};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! // At J2000.0 (2000-01-01 12:00 TT), TAI-UTC = 32 seconds
//! let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tai = utc.to_tai().unwrap();
//!
//! let offset_days = tai.to_julian_date().to_f64() - utc.to_julian_date().to_f64();
//! let offset_seconds = offset_days * 86400.0;
//! assert!((offset_seconds - 32.0).abs() < 0.001);
//! ```
//!
//! # Helper Functions
//!
//! This module also provides calendar conversion utilities used by the UTC-TAI algorithms:
//!
//! - [`julian_to_calendar`]: Convert Julian Date to (year, month, day, day_fraction)
//! - [`calendar_to_julian`]: Convert (year, month, day) to Julian Date
//!
//! These functions handle the full range of historical dates and use compensated
//! summation (Kahan algorithm) to preserve precision when combining Julian Date components.
//!
//! # References
//!
//! - IERS Bulletins: Leap second announcements
//! - USNO: History of leap seconds and TAI-UTC differences
//! - ITU-R TF.460-6: Standard-frequency and time-signal emissions
//! - Explanatory Supplement to the Astronomical Almanac, 3rd ed., Chapter 3

use super::super::common::{get_tai_utc_offset, next_calendar_day};
use super::{ToTAI, ToUTC};
use crate::julian::JulianDate;
use crate::scales::{TAI, UTC};
use crate::{TimeError, TimeResult};
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

impl ToUTC for UTC {
    /// Identity conversion. Returns self unchanged.
    fn to_utc(&self) -> TimeResult<UTC> {
        Ok(*self)
    }
}

/// Convert a UTC Julian Date to TAI.
///
/// This function handles both the modern leap second era (1972+) and the pre-1972
/// drift correction era. The algorithm:
///
/// 1. Separates the Julian Date into integer and fractional parts, tracking which
///    component has larger magnitude for precision preservation.
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
pub fn utc_to_tai(utc_jd: JulianDate) -> TimeResult<TAI> {
    let (utc_int, utc_frac, big1) = if utc_jd.jd1().abs() >= utc_jd.jd2().abs() {
        (utc_jd.jd1(), utc_jd.jd2(), true)
    } else {
        (utc_jd.jd2(), utc_jd.jd1(), false)
    };

    let (year, month, day, mut day_fraction) = julian_to_calendar(utc_int, utc_frac)?;

    let offset_0h = get_tai_utc_offset(year, month, day, 0.0);

    let offset_12h = get_tai_utc_offset(year, month, day, 0.5);

    let (next_year, next_month, next_day) = next_calendar_day(year, month, day)?;
    let offset_24h = get_tai_utc_offset(next_year, next_month, next_day, 0.0);

    let drift_rate = 2.0 * (offset_12h - offset_0h);
    let leap_amount = offset_24h - (offset_0h + drift_rate);

    day_fraction *= (SECONDS_PER_DAY_F64 + leap_amount) / SECONDS_PER_DAY_F64;
    day_fraction *= (SECONDS_PER_DAY_F64 + drift_rate) / SECONDS_PER_DAY_F64;

    let (z1, z2) = calendar_to_julian(year, month, day);

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
/// 3. After 3 iterations, the UTC value is accurate to ~1 picosecond.
///
/// The iteration count of 3 is sufficient because:
/// - The initial guess is within ~37 seconds of the answer
/// - Each iteration reduces the error by a factor of ~10^14
/// - Floating-point precision limits further improvement
///
/// The algorithm correctly handles leap second boundaries where the UTC day
/// "stretches" to include 86401 seconds.
pub fn tai_to_utc(tai_jd: JulianDate) -> TimeResult<UTC> {
    const TAI_TO_UTC_ITERATIONS: usize = 3;
    let (tai_int, tai_frac, big1) = if tai_jd.jd1().abs() >= tai_jd.jd2().abs() {
        (tai_jd.jd1(), tai_jd.jd2(), true)
    } else {
        (tai_jd.jd2(), tai_jd.jd1(), false)
    };

    let utc_int = tai_int;
    let mut utc_frac = tai_frac;

    for _ in 0..TAI_TO_UTC_ITERATIONS {
        let guess_tai = utc_to_tai_jd(utc_int, utc_frac)?;
        utc_frac += tai_int - guess_tai.jd1();
        utc_frac += tai_frac - guess_tai.jd2();
    }

    let (utc_jd1, utc_jd2) = if big1 {
        (utc_int, utc_frac)
    } else {
        (utc_frac, utc_int)
    };

    Ok(UTC::from_julian_date(JulianDate::new(utc_jd1, utc_jd2)))
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
/// - `year`: Gregorian year (can be negative for BCE dates)
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
/// - Minimum: JD -68569.5 (around 4713 BCE, near the Julian epoch)
/// - Maximum: JD 1e9 (far future)
///
/// # Errors
///
/// Returns `TimeError::ConversionError` if the Julian Date is outside the valid range.
///
/// # Example
///
/// ```ignore
/// let (year, month, day, frac) = julian_to_calendar(2451545.0, 0.0)?;
/// assert_eq!((year, month, day), (2000, 1, 1));  // J2000.0 epoch
/// assert!((frac - 0.5).abs() < 1e-10);  // Noon
/// ```
pub fn julian_to_calendar(jd1: f64, jd2: f64) -> TimeResult<(i32, i32, i32, f64)> {
    let dj = jd1 + jd2;
    const DJMIN: f64 = -68569.5;
    const DJMAX: f64 = 1e9;

    if !(DJMIN..=DJMAX).contains(&dj) {
        return Err(TimeError::ConversionError(format!(
            "Julian Date {} out of valid range [{}, {}]",
            dj, DJMIN, DJMAX
        )));
    }

    fn nearest_int(a: f64) -> f64 {
        if a.abs() < 0.5 {
            0.0
        } else if a < 0.0 {
            (a - 0.5).ceil()
        } else {
            (a + 0.5).floor()
        }
    }

    let day_int_1 = nearest_int(jd1);
    let frac_1 = jd1 - day_int_1;
    let mut jd = day_int_1 as i64;

    let day_int_2 = nearest_int(jd2);
    let frac_2 = jd2 - day_int_2;
    jd += day_int_2 as i64;

    let mut sum = 0.5;
    let mut correction = 0.0;
    let fractions = [frac_1, frac_2];

    for frac in fractions.iter() {
        let temp = sum + frac;
        correction += if sum.abs() >= frac.abs() {
            (sum - temp) + frac
        } else {
            (frac - temp) + sum
        };
        sum = temp;

        if sum >= 1.0 {
            jd += 1;
            sum -= 1.0;
        }
    }
    let mut fraction = sum + correction;
    correction = fraction - sum;

    if fraction < 0.0 {
        fraction = sum + 1.0;
        correction += (1.0 - fraction) + sum;
        sum = fraction;
        fraction = sum + correction;
        correction = fraction - sum;
        jd -= 1
    }

    #[allow(clippy::excessive_precision)]
    const DBL_EPSILON: f64 = 2.220_446_049_250_313_1e-16;
    if (fraction - 1.0) >= -DBL_EPSILON / 4.0 {
        let temp = sum - 1.0;
        correction += (sum - temp) - 1.0;
        sum = temp;
        fraction = sum + correction;

        if (-DBL_EPSILON / 2.0) < fraction {
            jd += 1;
            fraction = if fraction > 0.0 { fraction } else { 0.0 };
        }
    }

    let mut l = jd + 68569_i64;
    let n = (4_i64 * l) / 146097_i64;
    l -= (146097_i64 * n + 3_i64) / 4_i64;
    let i = (4000_i64 * (l + 1_i64)) / 1461001_i64;
    l -= (1461_i64 * i) / 4_i64 - 31_i64;
    let k = (80_i64 * l) / 2447_i64;
    let day = (l - (2447_i64 * k) / 80_i64) as i32;
    let l_final = k / 11_i64;
    let month = (k + 2_i64 - 12_i64 * l_final) as i32;
    let year = (100_i64 * (n - 49_i64) + i + l_final) as i32;

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
/// ```ignore
/// let (jd1, jd2) = calendar_to_julian(2000, 1, 1);
/// // jd1 + jd2 = JD 2451544.5 (midnight starting J2000.0)
/// ```
pub fn calendar_to_julian(year: i32, month: i32, day: i32) -> (f64, f64) {
    let my = (month - 14) / 12;
    let iypmy = year + my;

    let modified_jd = ((1461 * (iypmy + 4800)) / 4 + (367 * (month - 2 - 12 * my)) / 12
        - (3 * ((iypmy + 4900) / 100)) / 4
        + day
        - 2432076) as f64;

    (MJD_ZERO_POINT, modified_jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_identity_conversions() {
        let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let result = utc.to_utc().unwrap();
        assert_eq!(utc.to_julian_date().jd1(), result.to_julian_date().jd1());
        assert_eq!(utc.to_julian_date().jd2(), result.to_julian_date().jd2());
    }

    #[test]
    fn test_utc_tai_leap_second_offset() {
        // Known leap second values at key dates (TAI-UTC in seconds)
        assert_eq!(get_tai_utc_offset(1972, 1, 1, 0.0), 10.0); // First leap second era
        assert_eq!(get_tai_utc_offset(1980, 1, 1, 0.0), 19.0);
        assert_eq!(get_tai_utc_offset(1999, 1, 1, 0.0), 32.0);
        assert_eq!(get_tai_utc_offset(2017, 1, 1, 0.0), 37.0); // Most recent

        // Pre-1972 drift era (before discrete leap seconds)
        let offset_1970 = get_tai_utc_offset(1970, 6, 15, 0.5);
        assert!(offset_1970 > 0.0 && offset_1970 < 15.0);

        // Edge cases returning 0.0
        assert_eq!(get_tai_utc_offset(1959, 12, 31, 0.0), 0.0); // Before 1960
        assert_eq!(get_tai_utc_offset(2000, 1, 1, 1.5), 0.0); // Invalid fraction > 1
        assert_eq!(get_tai_utc_offset(2000, 1, 1, -0.5), 0.0); // Invalid fraction < 0
        assert_eq!(get_tai_utc_offset(1960, 0, 1, 0.0), 0.0); // Loop exhaustion
    }

    #[test]
    fn test_utc_tai_round_trip_precision() {
        // Tolerance for iterative refinement (3 iterations in tai_to_utc)
        // 1e-14 days ~ 1 picosecond, acceptable for iterative algorithm
        const TOLERANCE: f64 = 1e-14;

        let test_cases: &[(f64, f64)] = &[
            (J2000_JD, 0.123456789),
            (J2000_JD, 0.0),
            (J2000_JD, 0.999999),
            (0.5, J2000_JD),  // jd2 > jd1 (alternate split in utc_to_tai)
            (0.1, celestial_core::constants::J2000_JD), // jd2 > jd1 (alternate split in tai_to_utc)
        ];

        for &(jd1, jd2) in test_cases {
            // UTC -> TAI -> UTC round trip
            let utc = UTC::from_julian_date(JulianDate::new(jd1, jd2));
            let tai = utc.to_tai().unwrap();
            let utc_back = tai.to_utc().unwrap();
            let diff_utc =
                (utc.to_julian_date().to_f64() - utc_back.to_julian_date().to_f64()).abs();
            assert!(
                diff_utc < TOLERANCE,
                "UTC->TAI->UTC failed for ({}, {}): diff={:.2e}",
                jd1,
                jd2,
                diff_utc
            );

            // TAI -> UTC -> TAI round trip
            let tai = TAI::from_julian_date(JulianDate::new(jd1, jd2));
            let utc = tai.to_utc().unwrap();
            let tai_back = utc.to_tai().unwrap();
            let diff_tai =
                (tai.to_julian_date().to_f64() - tai_back.to_julian_date().to_f64()).abs();
            assert!(
                diff_tai < TOLERANCE,
                "TAI->UTC->TAI failed for ({}, {}): diff={:.2e}",
                jd1,
                jd2,
                diff_tai
            );
        }
    }

    #[test]
    fn test_calendar_helper_functions() {
        // next_calendar_day: month transitions
        assert_eq!(next_calendar_day(2000, 1, 31).unwrap(), (2000, 2, 1));
        assert_eq!(next_calendar_day(2000, 12, 31).unwrap(), (2001, 1, 1));
        assert_eq!(next_calendar_day(2000, 2, 29).unwrap(), (2000, 3, 1));

        // next_calendar_day: leap year February
        assert_eq!(next_calendar_day(2000, 2, 28).unwrap(), (2000, 2, 29));
        assert_eq!(next_calendar_day(1900, 2, 28).unwrap(), (1900, 3, 1));

        // next_calendar_day: 30-day months
        assert_eq!(next_calendar_day(2000, 4, 30).unwrap(), (2000, 5, 1));
        assert_eq!(next_calendar_day(2000, 6, 30).unwrap(), (2000, 7, 1));
        assert_eq!(next_calendar_day(2000, 9, 30).unwrap(), (2000, 10, 1));
        assert_eq!(next_calendar_day(2000, 11, 30).unwrap(), (2000, 12, 1));

        // next_calendar_day: error cases
        assert!(next_calendar_day(2000, 0, 1).is_err());
        assert!(next_calendar_day(2000, 13, 1).is_err());
        assert!(next_calendar_day(2000, -1, 1).is_err());

        // julian_to_calendar: out of range errors
        assert!(julian_to_calendar(1e10, 0.0).is_err());
        assert!(julian_to_calendar(-1e6, 0.0).is_err());

        // julian_to_calendar: negative fraction correction path
        let (y, m, d, frac) = julian_to_calendar(celestial_core::constants::J2000_JD, -0.6).unwrap();
        assert!(y > 0 && (1..=12).contains(&m) && (1..=31).contains(&d) && frac >= 0.0);

        // julian_to_calendar: Kahan summation else branch (|frac_2| > |sum|)
        let (y, m, d, frac) = julian_to_calendar(2451544.6, 0.2).unwrap();
        assert!(y > 0 && (1..=12).contains(&m) && (1..=31).contains(&d));
        assert!(frac >= 0.0 && frac <= 1.0);

        // julian_to_calendar: near-1.0 fraction correction
        let (y, m, d, frac) = julian_to_calendar(2451544.75, 0.75).unwrap();
        assert!(y > 0 && (1..=12).contains(&m) && (1..=31).contains(&d));
        assert!(frac >= 0.0 && frac <= 1.0);
    }
}
