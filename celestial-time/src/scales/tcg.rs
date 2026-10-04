//! Geocentric Coordinate Time (TCG) time scale.
//!
//! TCG is the proper time of a clock at rest at the geocenter, free from Earth's
//! gravitational potential. It runs faster than TT by approximately 22 ms per year
//! due to gravitational time dilation.
//!
//! # Background
//!
//! TCG was introduced by the IAU in 1991 as the coordinate time for the Geocentric
//! Celestial Reference System (GCRS). While TT is adjusted to match the rate of
//! proper time on Earth's geoid, TCG ticks at the rate of a clock experiencing
//! no gravitational potential.
//!
//! The relationship between TCG and TT is defined by IAU Resolution B1.9 (2000):
//!
//! ```text
//! TCG - TT = L_G * (JD_TCG - T_0) * 86400
//!
//! where:
//!   L_G = 6.969290134e-10 (defining constant)
//!   T_0 = 2443144.5003725 (TCG/TT coincidence epoch, 1977-01-01 00:00:32.184)
//! ```
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tcg::TCG;
//! use celestial_time::scales::tcg::tcg_from_calendar;
//!
//! let tcg = TCG::j2000();
//! let jd = tcg.to_julian_date();
//!
//! let tcg_cal = tcg_from_calendar(2000, 1, 1, 12, 0, 0.0).unwrap();
//! ```
//!
//! # Precision
//!
//! TCG values use split Julian Date storage internally. The struct methods preserve
//! full f64 precision through all arithmetic operations. Conversions to/from TT
//! maintain nanosecond accuracy.

time_scale! {
    /// Geocentric Coordinate Time.
    ///
    /// Wraps a `JulianDate` representing an instant in the TCG time scale.
    /// TCG is the coordinate time for the Geocentric Celestial Reference System,
    /// running ~6.97e-10 faster than TT (about 22 ms per year).
    TCG
}

uniform_day_scale!(TCG, tcg_from_calendar);
