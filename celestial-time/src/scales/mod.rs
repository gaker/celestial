//! Astronomical time scales.
//!
//! Provides implementations of the eight primary time scales used in astronomical
//! calculations: UTC, TAI, TT, UT1, GPS, TDB, TCB, and TCG.
//!
//! # Time Scale Overview
//!
//! | Scale | Description | TAI Relationship |
//! |-------|-------------|------------------|
//! | TAI | International Atomic Time | Reference |
//! | UTC | Coordinated Universal Time | TAI - leap seconds |
//! | TT | Terrestrial Time | TAI + 32.184s |
//! | UT1 | Earth rotation time | Requires EOP data |
//! | GPS | GPS satellite time | TAI - 19s |
//! | TCG | Geocentric Coordinate Time | Linear scale from TT |
//! | TDB | Barycentric Dynamical Time | TT + periodic terms |
//! | TCB | Barycentric Coordinate Time | Linear scale from TDB |
//!
//! # Usage
//!
//! Each time scale is a newtype wrapping a Julian Date. Create instances via
//! `from_julian_date` or `*_from_calendar` helper functions:
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::utc::UTC;
//! use celestial_time::scales::tai::tai_from_calendar;
//! use celestial_time::scales::tt::tt_from_calendar;
//!
//! // From Julian Date
//! let tai = TAI::from_julian_date(JulianDate::new(2451545.0, 0.0));
//!
//! // From calendar components
//! let tt = tt_from_calendar(2000, 1, 1, 12, 0, 0.0).unwrap();
//! ```
//!
//! # Conversions
//!
//! Convert between scales using traits from the [`conversions`] submodule:
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::gps::GPS;
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::conversions::{ToTAI, ToTT, ToGPS};
//!
//! let tai = TAI::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let tt = tai.to_tt().unwrap();
//! let gps = tai.to_gps().unwrap();
//! ```
//!
//! Scales without a direct conversion go through an intermediate one. For
//! example, GPS to TT is `gps.to_tai()?.to_tt()`.
//!
//! # Precision
//!
//! All time scales use split Julian Date storage (jd1, jd2) to preserve
//! nanosecond precision. When adding time, the fraction of a day goes into
//! jd2 and whole days are carried into jd1.

#[macro_use]
mod macros;

pub mod common;
pub mod conversions;
pub mod gps;
pub mod tai;
pub mod tcb;
pub mod tcg;
pub mod tdb;
pub mod tt;
pub mod ut1;
pub mod utc;
