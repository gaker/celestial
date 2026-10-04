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
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::ut1::UT1;
//! use celestial_time::scales::ut1::ut1_from_calendar;
//!
//! // From Unix timestamp components
//! let ut1 = UT1::new(0, 0);  // Unix epoch in UT1
//!
//! // From calendar date
//! let ut1 = ut1_from_calendar(2000, 1, 1, 12, 0, 0.0).unwrap();
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

time_scale! {
    /// Universal Time UT1, based on Earth's rotation angle.
    ///
    /// Internally stores time as a split Julian Date for full precision.
    UT1
}

uniform_day_scale!(UT1, ut1_from_calendar);
