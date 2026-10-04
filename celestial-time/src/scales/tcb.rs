//! Barycentric Coordinate Time (TCB) representation.
//!
//! TCB is the coordinate time for the barycentric reference frame, as defined by the IAU.
//! It ticks faster than TDB by approximately 1.55e-8. That rate is L_C + L_G: Earth's
//! orbital velocity and the Sun's and planets' potential (L_C), plus Earth's own
//! potential at the geoid (L_G).
//!
//! # Relationship to TDB
//!
//! TCB and TDB are related by a linear transformation:
//!
//! ```text
//! TCB - TDB = L_B * (JD_TCB - T_0) * 86400 - TDB_0
//! ```
//!
//! Where:
//! - L_B = 1.550519768e-8 (IAU 2006 Resolution B3)
//! - T_0 = 2443144.5003725 (TCB-TDB epoch, 1977 Jan 1.0 TAI)
//! - TDB_0 = -6.55e-5 s
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tcb::TCB;
//!
//! // Create from Julian Date
//! let tcb = TCB::from_julian_date(JulianDate::j2000());
//!
//! // Create from calendar
//! use celestial_time::scales::tcb::tcb_from_calendar;
//! let tcb = tcb_from_calendar(2000, 1, 1, 12, 0, 0.0).unwrap();
//!
//! // Parse from ISO 8601
//! let tcb: TCB = "2000-01-01T12:00:00".parse().unwrap();
//! ```
//!
//! # When to Use TCB
//!
//! TCB is the natural time coordinate for barycentric calculations (solar system dynamics,
//! pulsar timing, VLBI). For most terrestrial applications, TDB is more practical since
//! it stays close to TT.

time_scale! {
    /// Barycentric Coordinate Time.
    ///
    /// Wraps a Julian Date interpreted in the TCB time scale. TCB is the coordinate
    /// time of the BCRS, equivalent to a clock at rest far from the solar system.
    TCB
}

uniform_day_scale!(TCB, tcb_from_calendar);
