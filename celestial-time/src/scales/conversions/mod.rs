//! Time scale conversions between astronomical time systems.
//!
//! This module provides traits and implementations for converting between the eight
//! major astronomical time scales: GPS, TAI, TT, TCG, TCB, TDB, UT1, and UTC.
//!
//! # Time Scale Overview
//!
//! | Scale | Full Name | Basis | Primary Use |
//! |-------|-----------|-------|-------------|
//! | UTC | Coordinated Universal Time | Atomic + leap seconds | Civil timekeeping |
//! | TAI | International Atomic Time | Atomic clocks | Reference for other scales |
//! | TT | Terrestrial Time | TAI + 32.184s | Geocentric ephemerides |
//! | UT1 | Universal Time 1 | Earth rotation | Sidereal time, telescope pointing |
//! | GPS | GPS Time | Atomic (no leap seconds) | Satellite navigation |
//! | TCG | Geocentric Coordinate Time | Relativistic (Earth center) | Precise geocentric dynamics |
//! | TCB | Barycentric Coordinate Time | Relativistic (solar system) | Solar system dynamics |
//! | TDB | Barycentric Dynamical Time | TCB rescaled | Solar system ephemerides |
//!
//!
//! # Fixed vs Variable Offsets
//!
//! Some conversions use constant offsets:
//!
//! - **TAI <-> TT**:  Fixed 32.184 seconds
//! - **TAI <-> GPS**: Fixed 19.0 seconds
//! - **TT <-> TCG**:  Secular rate 6.969290134e-10 (IAU 2000)
//!
//! Others require external data that changes over time:
//!
//! - **UTC <-> TAI**: Leap second table (currently 37 seconds as of 2017)
//! - **UT1 <-> TAI**: IERS Earth Orientation Parameters (changes daily)
//! - **UT1 <-> TT**:  Delta-T from historical tables or predictions
//!
//! TT <-> TDB is an analytic series (~1.66ms annual oscillation) that needs the
//! observer's location rather than external data.
//!
//! # Trait Pattern
//!
//! Each target scale has a conversion trait:
//!
//! - [`ToTAI`]: Convert to International Atomic Time
//! - [`ToTT`]:  Convert to Terrestrial Time
//! - [`ToGPS`]: Convert to GPS Time
//! - [`ToUTC`]: Convert to Coordinated Universal Time
//! - [`ToTCG`]: Convert to Geocentric Coordinate Time
//!
//! Additional traits handle conversions requiring parameters:
//!
//! - [`ToUT1WithDUT1`](utc_ut1::ToUT1WithDUT1), [`ToUTCWithDUT1`](utc_ut1::ToUTCWithDUT1): UTC <-> UT1 with IERS DUT1
//! - [`ToUT1WithOffset`](ut1_tai::ToUT1WithOffset), [`ToTAIWithOffset`](ut1_tai::ToTAIWithOffset): UT1 <-> TAI with IERS offset
//! - [`ToTTWithDeltaT`](ut1_tai::ToTTWithDeltaT), [`ToUT1WithDeltaT`](ut1_tai::ToUT1WithDeltaT):  UT1 <-> TT with historical Delta-T
//! - [`ToTDB`](tt_tdb::ToTDB), [`ToTTFromTDB`](tt_tdb::ToTTFromTDB):               TT <-> TDB with observer location
//!
//! # Usage
//!
//! Simple fixed-offset conversions work directly:
//!
//! ```
//! use celestial_time::scales::gps::GPS;
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::conversions::{ToTAI, ToTT};
//! use celestial_time::julian::JulianDate;
//!
//! let tai = TAI::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let tt = tai.to_tt().unwrap();   // TAI + 32.184s
//! ```
//!
//! Conversions requiring external data take parameters:
//!
//! ```
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::ut1::UT1;
//! use celestial_time::scales::conversions::ut1_tai::{ToTAIWithOffset, ToUT1WithOffset};
//! use celestial_time::julian::JulianDate;
//!
//! // UT1-UTC from IERS Bulletin A (0.3554 s) minus TAI-UTC (32 s) at J2000
//! let ut1_tai_offset = -31.645;  // seconds
//!
//! let tai = TAI::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let ut1 = tai.to_ut1_with_offset(ut1_tai_offset).unwrap();
//! let back = ut1.to_tai_with_offset(ut1_tai_offset).unwrap();
//! ```
//!
//! # Precision Notes
//!
//! All conversions preserve precision by applying offsets to the smaller-magnitude
//! component of the two-part Julian Date. Round-trip conversions maintain
//! sub-nanosecond accuracy for fixed-offset scales.

pub mod gps_tai;
pub mod tai_tt;
pub mod tcb_tdb;
pub mod tt_tcg;
pub mod tt_tdb;
pub mod ut1_tai;
pub mod utc_tai;
pub mod utc_ut1;

use crate::scales::gps::GPS;
use crate::scales::tai::TAI;
use crate::scales::tcg::TCG;
use crate::scales::tt::TT;
use crate::scales::utc::UTC;
use crate::TimeResult;

/// Convert a time scale to GPS Time.
///
/// GPS Time runs at the same rate as TAI but with a fixed offset of -19 seconds
/// (GPS = TAI - 19s). It does not include leap seconds, making it continuous
/// since its epoch of January 6, 1980.
///
/// Implemented for: TAI
pub trait ToGPS {
    /// Convert to GPS Time.
    fn to_gps(&self) -> TimeResult<GPS>;
}

/// Convert a time scale to International Atomic Time (TAI).
///
/// TAI is the fundamental atomic time scale, maintained by a weighted average of
/// atomic clocks worldwide. It serves as the reference for most other time scales.
///
/// Implemented for: UTC, TT, GPS, TCG
///
/// For UT1 → TAI, use [`ToTAIWithOffset`](ut1_tai::ToTAIWithOffset) which requires the IERS UT1-TAI offset.
pub trait ToTAI {
    /// Convert to TAI.
    fn to_tai(&self) -> TimeResult<TAI>;
}

/// Convert a time scale to Terrestrial Time (TT).
///
/// TT is the time scale for geocentric ephemerides and Earth-based observations.
/// It differs from TAI by exactly 32.184 seconds (TT = TAI + 32.184s), a value
/// chosen for continuity with the older ET (Ephemeris Time) scale.
///
/// Implemented for: TAI, TCG
///
/// For TDB → TT, use [`ToTTFromTDB`](tt_tdb::ToTTFromTDB) which requires observer location.
/// For UT1 → TT, use [`ToTTWithDeltaT`](ut1_tai::ToTTWithDeltaT) which requires Delta-T.
pub trait ToTT {
    /// Convert to Terrestrial Time.
    fn to_tt(&self) -> TimeResult<TT>;
}

/// Convert a time scale to Geocentric Coordinate Time (TCG).
///
/// TCG is the proper time for a clock at the geocenter in the absence of
/// Earth's mass. It runs faster than TT by a rate of
/// approximately 6.969290134e-10 (about 22 ms/year).
///
/// Implemented for: TT, TAI
pub trait ToTCG {
    /// Convert to Geocentric Coordinate Time.
    fn to_tcg(&self) -> TimeResult<TCG>;
}

/// Convert a time scale to Coordinated Universal Time (UTC).
///
/// UTC is the basis for civil timekeeping. It tracks TAI but is adjusted by
/// leap seconds to stay within 0.9 seconds of UT1 (Earth rotation time).
/// Leap seconds are inserted (or theoretically removed) based on IERS
/// announcements.
///
/// Implemented for: TAI
///
/// Note: This conversion requires the leap second table in this crate,
/// which must be updated when new leap seconds are announced.
pub trait ToUTC {
    /// Convert to Coordinated Universal Time.
    fn to_utc(&self) -> TimeResult<UTC>;
}
