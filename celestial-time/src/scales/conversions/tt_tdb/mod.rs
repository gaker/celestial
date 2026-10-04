//! TT (Terrestrial Time) and TDB (Barycentric Dynamical Time) conversions.
//!
//! TDB is the independent time argument for solar system barycentric ephemerides.
//! Unlike most time scale pairs, the TT-TDB relationship is **location-dependent**
//! because it accounts for relativistic effects at the observer's position.
//!
//! # The TT-TDB Difference
//!
//! TDB is TCB linearly rescaled to TT's mean rate. TDB-TT has no secular drift,
//! only periodic variations due to:
//!
//! - Earth's orbital eccentricity (main ~1.66ms annual term)
//! - Lunar and planetary perturbations (smaller terms)
//! - Observer's position on Earth (diurnal terms, ~microsecond level)
//!
//! The difference TDB-TT oscillates with a peak-to-peak amplitude of about 3.3ms,
//! dominated by a ~1.66ms sinusoidal annual variation.
//!
//! # Algorithm
//!
//! Uses the Fairhead & Bretagnon (1990) series with 787 terms (the FAIRHD coefficients).
//! This provides sub-microsecond accuracy for dates within a few centuries of J2000.0.
//!
//! # Usage Patterns
//!
//! Two traits provide conversion:
//!
//! - [`ToTDB`]:       Convert TT to TDB (or TDB to itself)
//! - [`ToTTFromTDB`]: Convert TDB to TT with location awareness
//!
//! ```
//! use celestial_time::scales::tdb::TDB;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::conversions::tt_tdb::{ToTDB, ToTTFromTDB};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::location::Location;
//!
//! // TT → TDB: requires observer location
//! let tt = TT::from_julian_date(JulianDate::new(2451545.0, 0.5));
//! let tdb = tt.to_tdb_greenwich().unwrap();  // Uses Greenwich as reference
//!
//! // Or with explicit location
//! let tokyo = Location::from_degrees(35.6762, 139.6503, 40.0).unwrap();
//! let tdb = tt.to_tdb_with_location(&tokyo).unwrap();
//!
//! // TDB → TT: also requires location
//! let tdb = TDB::from_julian_date(JulianDate::new(2451545.0, 0.5));
//! let tt = tdb.to_tt_greenwich().unwrap();
//! ```
//!
//! # Why TDB Has No `to_tt()`
//!
//! TDB does not implement [`ToTT`](super::ToTT), because a conversion without a
//! location would silently introduce errors of up to ~2 µs. Use `to_tt_greenwich()`
//! or `to_tt_with_location()` instead.
//!
//! # UT1 Offset Parameter
//!
//! For highest precision, provide UT1-TT in seconds (that is, -ΔT). The diurnal
//! terms in the TDB-TT difference depend on the observer's local solar time,
//! which requires UT1. If omitted (UT1 = TT), the error is about 10 nanoseconds.

mod series;

use super::ut1_tai::ToUT1WithDeltaT;
use crate::julian::{finite_arg, finite_jd, JulianDate};
use crate::scales::tdb::TDB;
use crate::scales::tt::TT;
use crate::TimeResult;
use celestial_core::constants::SECONDS_PER_DAY_F64;
use celestial_core::location::Location;
use celestial_core::math::fmod;
use series::calculate_tdb_tt_difference;

/// Returns the Royal Observatory Greenwich location.
///
/// Used as the default reference point for TT-TDB conversions when no
/// observer location is specified.
fn greenwich_location() -> TimeResult<Location> {
    Ok(Location::from_degrees(51.477928, 0.0, 46.0)?)
}

const DEFAULT_UT1_MINUS_TT_SECONDS: f64 = 0.0;

// JD days start at noon, but the series' diurnal terms take the UT1 day
// fraction counted from midnight. Taking it from the two-part UT1 rounds the
// same way as ERFA's chain (eraTtut1, then the fraction), so dtr matches it to
// the last bit.
fn ut1_day_fraction(date: &JulianDate, ut1_minus_tt_seconds: f64) -> TimeResult<f64> {
    let ut1 = TT::from_julian_date(*date)
        .to_ut1_with_delta_t(-ut1_minus_tt_seconds)?
        .to_julian_date();
    let fraction = fmod(fmod(ut1.jd1(), 1.0) + fmod(ut1.jd2(), 1.0) + 0.5, 1.0);
    Ok(if fraction < 0.0 {
        fraction + 1.0
    } else {
        fraction
    })
}

/// Compute the TDB-TT offset in seconds for a given date and observer location.
///
/// This is the core calculation using Fairhead & Bretagnon (1990) coefficients.
/// The result is TDB - TT in seconds; add to TT to get TDB, subtract from TDB to get TT.
///
/// # Arguments
///
/// - `date_jd`: Julian Date (two-part for precision)
/// - `ut1_fraction`: Fraction of UT1 day (0.0 to 1.0), used for diurnal terms
/// - `location`: Observer's geographic location
///
/// # Returns
///
/// TDB - TT offset in seconds. Typical magnitude is < 0.002 seconds.
pub fn compute_tdb_tt_offset(
    date_jd: &JulianDate,
    ut1_fraction: f64,
    location: &Location,
) -> TimeResult<f64> {
    let date_jd = finite_jd(*date_jd)?;
    finite_arg("ut1_fraction", ut1_fraction)?;
    let (u, v) = location.to_geocentric_km();

    let dtr = calculate_tdb_tt_difference(
        date_jd.jd1(),
        date_jd.jd2(),
        ut1_fraction,
        location.longitude(),
        u,
        v,
    );

    Ok(dtr)
}

/// Convert a time scale to TDB (Barycentric Dynamical Time).
///
/// Implemented for TT. The conversion requires the observer's location.
///
/// # Methods
///
/// - `to_tdb_greenwich()` - Convert using Greenwich as reference location
/// - `to_tdb_with_location()` - Convert using explicit observer location
/// - `to_tdb_with_location_and_ut1_offset()` - Convert with location and UT1-TT offset
/// - `to_tdb_with_offset()` - Convert using pre-computed TDB-TT offset in seconds
pub trait ToTDB {
    /// Convert to TDB using Greenwich Observatory as reference.
    fn to_tdb_greenwich(&self) -> TimeResult<TDB>;
    /// Convert to TDB using the specified observer location.
    fn to_tdb_with_location(&self, location: &Location) -> TimeResult<TDB>;
    /// Convert to TDB with location and UT1-TT offset for maximum precision.
    fn to_tdb_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_minus_tt_seconds: f64,
    ) -> TimeResult<TDB>;
    /// Convert to TDB using a pre-computed offset (TDB-TT) in seconds.
    fn to_tdb_with_offset(&self, dtr_seconds: f64) -> TimeResult<TDB>;
}

impl ToTDB for TT {
    fn to_tdb_greenwich(&self) -> TimeResult<TDB> {
        let location = greenwich_location()?;
        self.to_tdb_with_location_and_ut1_offset(&location, DEFAULT_UT1_MINUS_TT_SECONDS)
    }

    fn to_tdb_with_location(&self, location: &Location) -> TimeResult<TDB> {
        self.to_tdb_with_location_and_ut1_offset(location, DEFAULT_UT1_MINUS_TT_SECONDS)
    }

    fn to_tdb_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_minus_tt_seconds: f64,
    ) -> TimeResult<TDB> {
        finite_arg("ut1_minus_tt_seconds", ut1_minus_tt_seconds)?;
        let tt_jd = finite_jd(self.to_julian_date())?;
        let ut1_fraction = ut1_day_fraction(&tt_jd, ut1_minus_tt_seconds)?;
        let dtr = compute_tdb_tt_offset(&tt_jd, ut1_fraction, location)?;
        self.to_tdb_with_offset(dtr)
    }

    fn to_tdb_with_offset(&self, dtr_seconds: f64) -> TimeResult<TDB> {
        let tt_jd = finite_jd(self.to_julian_date())?;

        let dtr_days = finite_arg("dtr_seconds", dtr_seconds)? / SECONDS_PER_DAY_F64;

        let tdb_jd = tt_jd.map_smaller_part(|_, small| small + dtr_days);

        Ok(TDB::from_julian_date(tdb_jd))
    }
}

/// Convert TDB to TT with location awareness.
///
/// This trait exists because the generic [`ToTT`](super::ToTT) trait cannot provide accurate
/// TDB→TT conversion without knowing the observer's location. Rather than
/// silently use a default, the design requires explicit location specification.
///
/// # Methods
///
/// - `to_tt_greenwich()` - Convert using Greenwich as reference location
/// - `to_tt_with_location()` - Convert using explicit observer location
/// - `to_tt_with_location_and_ut1_offset()` - Convert with location and UT1-TT offset
/// - `to_tt_with_offset()` - Convert using pre-computed TDB-TT offset in seconds
///
/// # Inverse Operation
///
/// The TDB-TT offset is strictly a function of TT, but this conversion evaluates
/// it at TDB in a single pass. The ~1.7 ms difference between the two moves the
/// offset by less than a picosecond.
pub trait ToTTFromTDB {
    /// Convert to TT using Greenwich Observatory as reference.
    fn to_tt_greenwich(&self) -> TimeResult<TT>;
    /// Convert to TT using the specified observer location.
    fn to_tt_with_location(&self, location: &Location) -> TimeResult<TT>;
    /// Convert to TT with location and UT1-TT offset for maximum precision.
    fn to_tt_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_minus_tt_seconds: f64,
    ) -> TimeResult<TT>;
    /// Convert to TT using a pre-computed offset (TDB-TT) in seconds.
    fn to_tt_with_offset(&self, dtr_seconds: f64) -> TimeResult<TT>;
}

impl ToTTFromTDB for TDB {
    fn to_tt_greenwich(&self) -> TimeResult<TT> {
        let location = greenwich_location()?;
        self.to_tt_with_location_and_ut1_offset(&location, DEFAULT_UT1_MINUS_TT_SECONDS)
    }

    fn to_tt_with_location(&self, location: &Location) -> TimeResult<TT> {
        self.to_tt_with_location_and_ut1_offset(location, DEFAULT_UT1_MINUS_TT_SECONDS)
    }

    fn to_tt_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_minus_tt_seconds: f64,
    ) -> TimeResult<TT> {
        finite_arg("ut1_minus_tt_seconds", ut1_minus_tt_seconds)?;
        // The offset is a function of TT, but evaluating it at TDB instead moves
        // it by under 1e-12 s. That is ERFA's chain, so we match it bit for bit.
        let tdb_jd = finite_jd(self.to_julian_date())?;
        let ut1_fraction = ut1_day_fraction(&tdb_jd, ut1_minus_tt_seconds)?;
        let dtr = compute_tdb_tt_offset(&tdb_jd, ut1_fraction, location)?;
        self.to_tt_with_offset(dtr)
    }

    fn to_tt_with_offset(&self, dtr_seconds: f64) -> TimeResult<TT> {
        let tdb_jd = finite_jd(self.to_julian_date())?;

        let dtr_days = finite_arg("dtr_seconds", dtr_seconds)? / SECONDS_PER_DAY_F64;

        let tt_jd = tdb_jd.map_smaller_part(|_, small| small - dtr_days);

        Ok(TT::from_julian_date(tt_jd))
    }
}

#[cfg(test)]
mod tests;
