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
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::tai::tai_from_calendar;
//!
//! // From Julian Date
//! let tai = TAI::j2000();
//!
//! // From calendar date
//! let tai = tai_from_calendar(2024, 6, 15, 12, 30, 0.0).unwrap();
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

time_scale! {
    /// International Atomic Time representation.
    ///
    /// Wraps a `JulianDate` to provide TAI-specific semantics. TAI serves as the
    /// hub for conversions between other time scales (UTC, TT, GPS, UT1, etc.).
    TAI,
    j2000_note = "This is 32.184 s after the J2000.0 epoch, which is defined in TT."
}

uniform_day_scale!(TAI, tai_from_calendar);
