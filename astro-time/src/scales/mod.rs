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
//! use astro_time::{JulianDate, TAI, TT, UTC};
//! use astro_time::scales::{tai_from_calendar, tt_from_calendar};
//!
//! // From Julian Date
//! let tai = TAI::from_julian_date(JulianDate::new(2451545.0, 0.0));
//!
//! // From calendar components
//! let tt = tt_from_calendar(2000, 1, 1, 12, 0, 0.0);
//! ```
//!
//! # Conversions
//!
//! Convert between scales using traits from the [`conversions`] submodule:
//!
//! ```
//! use astro_time::{JulianDate, GPS, TAI, TT};
//! use astro_time::scales::conversions::{ToTAI, ToTT, ToGPS};
//!
//! let tai = TAI::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let tt = tai.to_tt().unwrap();
//! let gps = tai.to_gps().unwrap();
//! ```
//!
//! Some conversions chain through intermediate scales internally. For example,
//! GPS to TT converts GPS -> TAI -> TT.
//!
//! # Precision
//!
//! All time scales use split Julian Date storage (jd1, jd2) to preserve
//! nanosecond precision. When adding offsets, the offset is applied to
//! the smaller-magnitude component.

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

pub use gps::gps_from_calendar;
pub use gps::GPS;
pub use tai::tai_from_calendar;
pub use tai::TAI;
pub use tcb::tcb_from_calendar;
pub use tcb::TCB;
pub use tcg::tcg_from_calendar;
pub use tcg::TCG;
pub use tdb::tdb_from_calendar;
pub use tdb::TDB;
pub use tt::tt_from_calendar;
pub use tt::TT;
pub use ut1::ut1_from_calendar;
pub use ut1::UT1;
pub use utc::utc_from_calendar;
pub use utc::UTC;

pub use conversions::{ToTAI, ToTCB, ToTCG, ToTCGFromTCB, ToTDB, ToTT, ToTTFromTDB};
