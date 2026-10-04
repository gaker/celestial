//! Terrestrial Time (TT) time scale.
//!
//! TT is the modern successor to Ephemeris Time (ET) and Terrestrial Dynamical Time (TDT).
//! It provides a uniform time scale for geocentric ephemerides and is the basis for
//! planetary position calculations referenced to Earth's center.
//!
//! # Relationship to TAI
//!
//! TT differs from TAI by a fixed offset:
//!
//! ```text
//! TT = TAI + 32.184 seconds
//! ```
//!
//! The 32.184s offset was chosen to maintain continuity with ET at the 1977 epoch.
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tt::TT;
//!
//! // Create TT at J2000.0 epoch
//! let tt = TT::j2000();
//!
//! // From calendar date
//! use celestial_time::scales::tt::tt_from_calendar;
//! let tt = tt_from_calendar(2000, 1, 1, 12, 0, 0.0).unwrap();
//!
//! // Parse from ISO 8601
//! let tt: TT = "2000-01-01T12:00:00".parse().unwrap();
//!
//! // Julian centuries since J2000.0 (for precession/nutation)
//! let centuries = tt.centuries_since_j2000().unwrap();
//! ```
//!
//! # Precision
//!
//! TT uses split Julian Date storage internally, preserving microsecond accuracy
//! across the full date range. The `centuries_since_j2000()` method provides
//! the T parameter used in IAU precession and nutation models.

time_scale! {
    /// Terrestrial Time representation.
    ///
    /// Wraps a split Julian Date for high-precision time storage.
    /// TT is the primary time scale for geocentric ephemeris calculations.
    TT,
    j2000_note = "This is the J2000.0 epoch, the fundamental epoch for modern astronomical calculations."
}

uniform_day_scale!(TT, tt_from_calendar);

use crate::julian::finite_jd;
use crate::{TimeError, TimeResult};
use celestial_core::constants::MAX_CENTURIES_FROM_J2000;

impl TT {
    /// Returns the Julian year corresponding to this TT instant.
    ///
    /// Julian year = 2000.0 + (JD - J2000_JD) / 365.25
    pub fn julian_year(&self) -> f64 {
        self.0.to_julian_year()
    }

    /// Returns Julian centuries since J2000.0 (the T parameter).
    ///
    /// This is the time argument used in IAU precession and nutation series.
    /// One Julian century = 36525 days. Fails if the date is not finite or is
    /// more than 20 centuries from J2000.0, outside those models' range.
    pub fn centuries_since_j2000(&self) -> TimeResult<f64> {
        let jd = finite_jd(self.0)?;
        let t = celestial_core::utils::jd_to_centuries(jd.jd1(), jd.jd2());
        if libm::fabs(t) > MAX_CENTURIES_FROM_J2000 {
            return Err(TimeError::InvalidEpoch(format!(
                "Epoch too far from J2000.0 for the IAU models: {:.1} centuries",
                t
            )));
        }
        Ok(t)
    }
}

#[cfg(test)]
mod tests;
