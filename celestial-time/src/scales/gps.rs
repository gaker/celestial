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
//! See [`JulianDate`](crate::julian::JulianDate) for details on the two-part representation.
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::gps::GPS;
//!
//! // From calendar date
//! let gps = celestial_time::scales::gps::gps_from_calendar(2024, 3, 15, 12, 0, 0.0).unwrap();
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
//! [`scales::conversions::gps_tai`](crate::scales::conversions::gps_tai) for the conversion traits.

time_scale! {
    /// GPS Time representation.
    ///
    /// Wraps a [`JulianDate`](crate::julian::JulianDate) to provide type safety and
    /// GPS-specific operations. The inner Julian Date uses split storage (jd1 + jd2)
    /// to preserve precision.
    GPS,
    j2000_note = "This is 51.184 s after the J2000.0 epoch, which is defined in TT."
}

uniform_day_scale!(GPS, gps_from_calendar);
