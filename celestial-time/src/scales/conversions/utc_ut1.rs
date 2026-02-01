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
//! are published weekly in IERS Bulletin A with ~1 ms precision. For high-precision
//! applications (astrometry, VLBI, satellite tracking), the correct DUT1 must be
//! obtained from IERS data for the specific date.
//!
//! # Conversion Paths
//!
//! This module provides two paths from UT1 to UTC:
//!
//! ```text
//! Direct:      UT1 ←→ UTC    (via DUT1 offset)
//! Via TAI:     UT1 → TAI → UTC    (for verification)
//! ```
//!
//! The direct path (`ToUTCWithDUT1`) handles leap second boundaries correctly by
//! adjusting the effective DUT1 value near discontinuities. The TAI path
//! (`ToUTCViaTAI`) chains through intermediate time scales and serves as a
//! verification mechanism.
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
//! Round-trip conversions (UTC → UT1 → UTC or UT1 → UTC → UT1) achieve ~1 picosecond
//! accuracy when using consistent DUT1 values. The two-part Julian Date representation
//! preserves full f64 precision throughout.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::{UTC, UT1};
//! use celestial_time::scales::conversions::{ToUT1WithDUT1, ToUTCWithDUT1};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! // DUT1 for 2000-01-01 was approximately +0.3 seconds (from IERS Bulletin A)
//! let dut1 = 0.3;
//!
//! let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let ut1 = utc.to_ut1_with_dut1(dut1).unwrap();
//!
//! // UT1 should be DUT1 seconds ahead of UTC
//! let diff_days = ut1.to_julian_date().to_f64() - utc.to_julian_date().to_f64();
//! let diff_seconds = diff_days * 86400.0;
//! assert!((diff_seconds - dut1).abs() < 0.01);
//!
//! // Round-trip preserves the original time
//! let utc_back = ut1.to_utc_with_dut1(dut1).unwrap();
//! let round_trip_diff = (utc.to_julian_date().to_f64()
//!     - utc_back.to_julian_date().to_f64()).abs();
//! assert!(round_trip_diff < 1e-14);  // ~1 picosecond
//! ```
//!
//! # References
//!
//! - IERS Bulletin A: Weekly publication of UT1-UTC values
//! - IERS Conventions (2010): Chapter 5, Earth Rotation
//! - USNO Earth Orientation Parameters

use super::super::common::get_tai_utc_offset; // Direct import from common
use super::ut1_tai::{ToTAIWithOffset, ToUT1WithOffset};
use super::utc_tai::{calendar_to_julian, julian_to_calendar};
use super::{ToTAI, ToUTC};
use crate::julian::JulianDate;
use crate::scales::{UT1, UTC};
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
    /// The TAI-UTC offset is looked up from the leap second table for the
    /// specific date. The result chains: UTC → TAI → UT1.
    fn to_ut1_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UT1> {
        let tai = self.to_tai()?;

        let utc_jd = self.to_julian_date();
        let (year, month, day, day_fraction) = julian_to_calendar(utc_jd.jd1(), utc_jd.jd2())?;
        let tai_utc_seconds = get_tai_utc_offset(year, month, day, day_fraction);
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
        let curr_offset = get_tai_utc_offset(year, month, day, 0.0);

        if i == -1 {
            prev_offset = curr_offset;
            continue;
        }

        let delta = curr_offset - prev_offset;
        if delta.abs() < 0.5 {
            prev_offset = curr_offset;
            continue;
        }

        // Found leap second boundary
        if delta * duts >= 0.0 {
            duts -= delta;
        }

        let (leap_d1, leap_d2) = calendar_to_julian(year, month, day);
        let time_past_leap =
            (jd_big - leap_d1) + (jd_small - (leap_d2 - 1.0 + duts / SECONDS_PER_DAY_F64));

        if time_past_leap > 0.0 {
            let fraction =
                (time_past_leap * SECONDS_PER_DAY_F64 / (SECONDS_PER_DAY_F64 + delta)).min(1.0);
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
        let ut1_jd = self.to_julian_date();
        let (big, small, big_first) = if ut1_jd.jd1().abs() >= ut1_jd.jd2().abs() {
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

/// Alternative UT1 to UTC conversion path that chains through TAI.
///
/// This trait provides a verification mechanism: the direct path (`ToUTCWithDUT1`)
/// and the TAI path should produce identical results. Any discrepancy indicates
/// a bug in one of the conversion implementations.
pub trait ToUTCViaTAI {
    /// Convert UT1 to UTC by chaining: UT1 → TAI → UTC.
    ///
    /// # Arguments
    ///
    /// * `dut1_seconds` - The UT1-UTC offset in seconds (from IERS Bulletin A)
    ///
    /// # Algorithm
    ///
    /// 1. Compute UT1-TAI offset from DUT1 and TAI-UTC for the date
    /// 2. Convert UT1 to TAI using the computed offset
    /// 3. Convert TAI to UTC using the leap second table
    ///
    /// This path uses the same TAI-UTC lookup as the direct conversion but
    /// exercises different code paths, making it useful for testing.
    fn to_utc_via_tai_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UTC>;
}

impl ToUTCViaTAI for UT1 {
    fn to_utc_via_tai_with_dut1(&self, dut1_seconds: f64) -> TimeResult<UTC> {
        let ut1_jd = self.to_julian_date();
        let (year, month, day, day_fraction) = julian_to_calendar(ut1_jd.jd1(), ut1_jd.jd2())?;
        let tai_utc_seconds = get_tai_utc_offset(year, month, day, day_fraction);
        let ut1_tai_offset = dut1_seconds - tai_utc_seconds;

        let tai = self.to_tai_with_offset(ut1_tai_offset)?;

        tai.to_utc()
    }
}

#[cfg(test)]
mod tests {
    use super::super::{ToUT1, ToUTC};
    use super::*;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_identity_conversions() {
        let ut1 = UT1::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let identity_ut1 = ut1.to_ut1().unwrap();
        assert_eq!(
            ut1.to_julian_date().jd1(),
            identity_ut1.to_julian_date().jd1()
        );
        assert_eq!(
            ut1.to_julian_date().jd2(),
            identity_ut1.to_julian_date().jd2()
        );

        let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let identity_utc = utc.to_utc().unwrap();
        assert_eq!(
            utc.to_julian_date().jd1(),
            identity_utc.to_julian_date().jd1()
        );
        assert_eq!(
            utc.to_julian_date().jd2(),
            identity_utc.to_julian_date().jd2()
        );
    }

    #[test]
    fn test_dut1_offset_relationship() {
        let dut1_values = [-0.9, 0.0, 0.9];

        for dut1 in dut1_values {
            let utc = UTC::from_julian_date(JulianDate::new(J2000_JD, 0.0));
            let ut1 = utc.to_ut1_with_dut1(dut1).unwrap();

            let ut1_jd = ut1.to_julian_date().to_f64();
            let diff_seconds = (ut1_jd - J2000_JD) * SECONDS_PER_DAY_F64;

            assert!(
                diff_seconds > -1.0 && diff_seconds < 1.0,
                "DUT1={}: UT1-UTC difference should be within 1 second: {} seconds",
                dut1,
                diff_seconds
            );

            let ut1_reverse = UT1::from_julian_date(JulianDate::new(J2000_JD, 0.0));
            let utc_reverse = ut1_reverse.to_utc_with_dut1(dut1).unwrap();
            let utc_jd = utc_reverse.to_julian_date().to_f64();
            let reverse_diff = (J2000_JD - utc_jd) * SECONDS_PER_DAY_F64;

            assert!(
                (reverse_diff - dut1).abs() < 0.1,
                "DUT1={}: UTC should be behind UT1 by ~DUT1: {} seconds",
                dut1,
                reverse_diff
            );
        }

        let ut1_normal = UT1::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let utc_normal = ut1_normal.to_utc_with_dut1(0.3).unwrap();
        assert!(
            utc_normal.to_julian_date().jd1().abs() > utc_normal.to_julian_date().jd2().abs(),
            "Should preserve larger JD1 component"
        );

        let ut1_flipped = UT1::from_julian_date(JulianDate::new(0.1, J2000_JD));
        let utc_flipped = ut1_flipped.to_utc_with_dut1(0.3).unwrap();
        assert!(
            utc_flipped.to_julian_date().jd2().abs() > utc_flipped.to_julian_date().jd1().abs(),
            "Should preserve larger JD2 component"
        );
    }

    #[test]
    fn test_utc_ut1_round_trip_precision() {
        let tolerance = 1e-14; // ~1 picosecond

        let jd_splits = [
            (J2000_JD, 0.123456789),
            (J2000_JD, 0.5),
            (0.1, J2000_JD),
            (J2000_JD, 0.999999999),
            (J2000_JD, 0.25),
        ];

        let dut1_values = [-0.9, 0.0, 0.9];

        for (jd1, jd2) in jd_splits {
            for dut1 in dut1_values {
                let original_utc = UTC::from_julian_date(JulianDate::new(jd1, jd2));
                let ut1 = original_utc.to_ut1_with_dut1(dut1).unwrap();
                let round_trip_utc = ut1.to_utc_with_dut1(dut1).unwrap();

                let diff = (original_utc.to_julian_date().to_f64()
                    - round_trip_utc.to_julian_date().to_f64())
                .abs();
                assert!(
                    diff < tolerance,
                    "UTC->UT1->UTC round trip (jd1={}, jd2={}, dut1={}): {:.2e} days exceeds {:.0e}",
                    jd1,
                    jd2,
                    dut1,
                    diff,
                    tolerance
                );

                let original_ut1 = UT1::from_julian_date(JulianDate::new(jd1, jd2));
                let utc = original_ut1.to_utc_with_dut1(dut1).unwrap();
                let round_trip_ut1 = utc.to_ut1_with_dut1(dut1).unwrap();

                let diff_reverse = (original_ut1.to_julian_date().to_f64()
                    - round_trip_ut1.to_julian_date().to_f64())
                .abs();
                assert!(
                    diff_reverse < tolerance,
                    "UT1->UTC->UT1 round trip (jd1={}, jd2={}, dut1={}): {:.2e} days exceeds {:.0e}",
                    jd1,
                    jd2,
                    dut1,
                    diff_reverse,
                    tolerance
                );
            }
        }
    }

    #[test]
    fn test_leap_second_boundary_handling() {
        let leap_dates = [
            (2441499.5, "1972-07-01"),
            (2441683.5, "1973-01-01"),
            (2442048.5, "1974-01-01"),
        ];

        for (jd, description) in leap_dates {
            let ut1_at_leap = UT1::from_julian_date(JulianDate::new(jd, 0.0));
            let utc = ut1_at_leap.to_utc_with_dut1(0.0).unwrap();
            assert!(
                utc.to_julian_date().to_f64() > 0.0,
                "{}: Should produce valid UTC at leap second",
                description
            );

            let ut1_after_leap = UT1::from_julian_date(JulianDate::new(jd, 0.001));
            let utc_after = ut1_after_leap.to_utc_with_dut1(0.0).unwrap();
            assert!(
                utc_after.to_julian_date().to_f64() > jd,
                "{}: UTC should be after leap second start",
                description
            );

            let utc_direct = ut1_at_leap.to_utc_with_dut1(0.0).unwrap();
            let utc_via_tai = ut1_at_leap.to_utc_via_tai_with_dut1(0.0).unwrap();
            let diff = (utc_direct.to_julian_date().to_f64()
                - utc_via_tai.to_julian_date().to_f64())
            .abs();
            assert!(
                diff < 1e-14,
                "{}: Direct vs TAI-intermediate should match: {:.2e} days",
                description,
                diff
            );
        }
    }
}
