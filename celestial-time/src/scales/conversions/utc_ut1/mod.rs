//! Conversions between Coordinated Universal Time (UTC) and Universal Time (UT1).
//!
//! UT1 is the principal form of Universal Time, directly tied to Earth's rotation angle.
//! UTC is the civil time standard maintained by atomic clocks. The difference between
//! them, called DUT1 (= UT1 - UTC), is published by the IERS in Bulletin A.
//!
//! # The DUT1 Offset
//!
//! DUT1 measures how much Earth's actual rotation deviates from the uniform UTC clock:
//!
//! ```text
//! UT1 = UTC + DUT1
//! ```
//!
//! The IERS keeps |DUT1| < 0.9 seconds by inserting leap seconds into UTC. DUT1 values
//! are published weekly in IERS Bulletin A at 10 µs resolution. For high-precision
//! applications (astrometry, VLBI, satellite tracking), the correct DUT1 must be
//! obtained from IERS data for the specific date.
//!
//! # Conversion Path
//!
//! UT1 and UTC convert directly through the DUT1 offset. `ToUTCWithDUT1` handles
//! leap second boundaries correctly by adjusting the effective DUT1 value near
//! discontinuities.
//!
//! # Leap Second Handling
//!
//! The UT1→UTC conversion is complicated by leap seconds: when UTC inserts a leap
//! second, there's a discontinuity in the UTC-TAI offset. The `adjust_dut1_for_leap_second`
//! function scans nearby days for offset changes and smoothly interpolates the
//! correction across the leap second boundary.
//!
//! # Precision
//!
//! From 1972 on, round-trip conversions (UTC → UT1 → UTC or UT1 → UTC → UT1) come back
//! within 1 ulp of jd2 (at most ~10 ps) when using consistent DUT1 values. The two-part
//! Julian Date representation preserves full f64 precision throughout.
//!
//! Before 1972, TAI-UTC drifts through each day and some days end in a step of a
//! fraction of a second. As in ERFA, UTC → UT1 goes through TAI with TAI-UTC taken at
//! 0h, while UT1 → UTC subtracts DUT1 directly. A round trip there can be off by up to
//! 2.6 ms, or up to 0.11 s on a day that ends in a step.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::ut1::UT1;
//! use celestial_time::scales::utc::UTC;
//! use celestial_time::scales::conversions::utc_ut1::{ToUT1WithDUT1, ToUTCWithDUT1};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! // DUT1 for 2000-01-01 was approximately +0.3 seconds (from IERS Bulletin A)
//! let dut1 = 0.3;
//!
//! let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let ut1 = utc.to_ut1_with_dut1(dut1).unwrap();
//!
//! // 0.3 s as a fraction of a day, rounded as ERFA's eraUtcut1 rounds it
//! assert_eq!(ut1.to_julian_date(), JulianDate::new(J2000_JD, 3.472222222206101e-6));
//!
//! // The way back lands 1.4 ps short of 0h, as eraUt1utc does
//! let utc_back = ut1.to_utc_with_dut1(dut1).unwrap();
//! assert_eq!(utc_back.to_julian_date(), JulianDate::new(J2000_JD, -1.612115456861747e-17));
//! ```
//!
//! # References
//!
//! - IERS Bulletin A: Weekly publication of UT1-UTC values
//! - IERS Conventions (2010): Chapter 5, Earth Rotation
//! - USNO Earth Orientation Parameters

use super::super::common::get_tai_utc_offset; // Direct import from common
use super::ut1_tai::ToUT1WithOffset;
use super::utc_tai::{calendar_to_julian, julian_to_calendar};
use super::ToTAI;
use crate::julian::{finite_arg, finite_jd, JulianDate};
use crate::scales::ut1::UT1;
use crate::scales::utc::UTC;
use crate::TimeResult;
use celestial_core::constants::SECONDS_PER_DAY_F64;

/// Convert to UT1 using a known DUT1 (UT1-UTC) offset.
///
/// DUT1 values must be obtained from IERS Bulletin A for the specific date.
/// The offset is typically in the range -0.9 to +0.9 seconds.
pub trait ToUT1WithDUT1 {
    /// Convert to UT1 given the DUT1 offset in seconds.
    ///
    /// # Arguments
    ///
    /// * `dut1_seconds` - The UT1-UTC offset in seconds (from IERS Bulletin A)
    ///
    /// # Returns
    ///
    /// The corresponding UT1 instant. The conversion chains through TAI:
    /// UTC → TAI → UT1, computing the UT1-TAI offset from DUT1 and the
    /// TAI-UTC offset for the date.
    fn to_ut1_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UT1>;
}

/// Convert to UTC using a known DUT1 (UT1-UTC) offset.
///
/// This is the inverse of `ToUT1WithDUT1`. Given a UT1 instant and the
/// DUT1 offset, computes the corresponding UTC instant.
pub trait ToUTCWithDUT1 {
    /// Convert to UTC given the DUT1 offset in seconds.
    ///
    /// # Arguments
    ///
    /// * `dut1_seconds` - The UT1-UTC offset in seconds (from IERS Bulletin A)
    ///
    /// # Returns
    ///
    /// The corresponding UTC instant. Handles leap second boundaries by
    /// adjusting the effective DUT1 value near discontinuities.
    fn to_utc_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UTC>;
}

impl ToUT1WithDUT1 for UTC {
    /// Convert UTC to UT1 by computing UT1-TAI from DUT1 and TAI-UTC.
    ///
    /// The conversion uses the relationship:
    ///
    /// ```text
    /// UT1 - TAI = DUT1 - (TAI - UTC) = DUT1 - TAI_UTC_offset
    /// ```
    ///
    /// The TAI-UTC offset is looked up from the leap second table at 0h on
    /// the UTC date. The result chains: UTC → TAI → UT1.
    fn to_ut1_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UT1> {
        finite_arg("dut1_seconds", dut1_seconds)?;
        let tai = self.to_tai()?;

        // TAI-UTC at 0h, as ERFA's utcut1 takes it. Before 1972 it drifts through
        // the day, so this is up to 2.6 ms from the value at the time itself.
        let utc_jd = self.to_julian_date();
        let (year, month, day, _) = julian_to_calendar(utc_jd.jd1(), utc_jd.jd2())?;
        let tai_utc_seconds = get_tai_utc_offset(year, month, day, 0.0)?;
        let ut1_tai_offset = dut1_seconds - tai_utc_seconds;

        tai.to_ut1_with_offset(ut1_tai_offset)
    }
}

/// Adjust DUT1 for leap second discontinuities near the given Julian Date.
///
/// When converting UT1 to UTC near a leap second boundary, the naive subtraction
/// of DUT1 can place the result on the wrong side of the discontinuity. This
/// function detects nearby leap seconds and adjusts the effective DUT1 value
/// to smoothly interpolate across the boundary.
///
/// # Algorithm
///
/// 1. Scan days from (JD - 1) to (JD + 3) looking for TAI-UTC offset changes
/// 2. If a change > 0.5 seconds is found, a leap second occurred
/// 3. If the leap second and DUT1 have the same sign, subtract the leap from DUT1
/// 4. Compute the fraction of the way past the leap second boundary
/// 5. Gradually add back the leap second contribution based on that fraction
///
/// The range [-1, +3] days ensures leap seconds are detected whether the input
/// time is just before, during, or just after the discontinuity.
///
/// # Arguments
///
/// * `jd_big` - Larger magnitude component of the Julian Date
/// * `jd_small` - Smaller magnitude component of the Julian Date
/// * `dut1` - The raw DUT1 offset in seconds
///
/// # Returns
///
/// The adjusted DUT1 value that accounts for any nearby leap second.
fn adjust_dut1_for_leap_second(jd_big: f64, jd_small: f64, dut1: f64) -> TimeResult<f64> {
    let mut duts = dut1;
    let mut prev_offset = 0.0;

    for i in -1..=3 {
        let jd_frac = jd_small + i as f64;
        let (year, month, day, _) = julian_to_calendar(jd_big, jd_frac)?;
        let curr_offset = get_tai_utc_offset(year, month, day, 0.0)?;

        if i == -1 {
            prev_offset = curr_offset;
            continue;
        }

        let delta = curr_offset - prev_offset;
        if libm::fabs(delta) < 0.5 {
            prev_offset = curr_offset;
            continue;
        }

        if delta * duts >= 0.0 {
            duts -= delta;
        }

        let (leap_d1, leap_d2) = calendar_to_julian(year, month, day)?;
        let time_past_leap =
            (jd_big - leap_d1) + (jd_small - (leap_d2 - 1.0 + duts / SECONDS_PER_DAY_F64));

        if time_past_leap > 0.0 {
            let fraction = libm::fmin(
                time_past_leap * SECONDS_PER_DAY_F64 / (SECONDS_PER_DAY_F64 + delta),
                1.0,
            );
            duts += delta * fraction;
        }
        break;
    }

    Ok(duts)
}

impl ToUTCWithDUT1 for UT1 {
    /// Convert UT1 to UTC by subtracting the adjusted DUT1 offset.
    ///
    /// The conversion:
    ///
    /// 1. Determines which JD component has larger magnitude (for precision)
    /// 2. Adjusts DUT1 for any nearby leap second boundaries
    /// 3. Subtracts the adjusted DUT1 from the smaller-magnitude component
    /// 4. Preserves the original JD component ordering
    ///
    /// The leap second adjustment ensures correct behavior at discontinuities
    /// where the UTC scale gains an extra second.
    fn to_utc_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UTC> {
        finite_arg("dut1_seconds", dut1_seconds)?;
        let ut1_jd = finite_jd(self.to_julian_date())?;
        let (big, small, big_first) = if libm::fabs(ut1_jd.jd1()) >= libm::fabs(ut1_jd.jd2()) {
            (ut1_jd.jd1(), ut1_jd.jd2(), true)
        } else {
            (ut1_jd.jd2(), ut1_jd.jd1(), false)
        };

        let adjusted_dut1 = adjust_dut1_for_leap_second(big, small, dut1_seconds)?;
        let small_corrected = small - adjusted_dut1 / SECONDS_PER_DAY_F64;

        let (utc_jd1, utc_jd2) = if big_first {
            (big, small_corrected)
        } else {
            (small_corrected, big)
        };
        Ok(UTC::from_julian_date(JulianDate::new(utc_jd1, utc_jd2)))
    }
}

#[cfg(test)]
mod tests;
