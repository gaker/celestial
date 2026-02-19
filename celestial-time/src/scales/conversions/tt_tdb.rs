//! TT (Terrestrial Time) and TDB (Barycentric Dynamical Time) conversions.
//!
//! TDB is the independent time argument for solar system barycentric ephemerides.
//! Unlike most time scale pairs, the TT-TDB relationship is **location-dependent**
//! because it accounts for relativistic effects at the observer's position.
//!
//! # The TT-TDB Difference
//!
//! TDB tracks time at the solar system barycenter. An observer on Earth experiences
//! periodic variations due to:
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
//! use celestial_time::scales::{TT, TDB};
//! use celestial_time::scales::conversions::{ToTDB, ToTTFromTDB};
//! use celestial_time::JulianDate;
//! use celestial_core::Location;
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
//! # Why `TDB.to_tt()` Returns an Error
//!
//! The base [`ToTT`] trait method deliberately returns an error for TDB because
//! location-independent conversion would silently introduce errors of up to ~1.7ms.
//! Use `to_tt_greenwich()` or `to_tt_with_location()` instead.
//!
//! # UT1 Offset Parameter
//!
//! For highest precision (~microsecond), provide the UT1-UTC offset. The diurnal
//! terms in the TDB-TT difference depend on the observer's local sidereal time,
//! which requires UT1. If omitted (set to 0), the error is typically < 10 microseconds.

use super::ToTT;
use crate::{
    constants::FAIRHD,
    julian::JulianDate,
    scales::{TDB, TT},
};

use crate::{TimeError, TimeResult};
use celestial_core::constants::{
    DAYS_PER_JULIAN_MILLENNIUM, DEG_TO_RAD, J2000_JD, SECONDS_PER_DAY_F64, TWOPI,
};
use celestial_core::math::fmod;
use celestial_core::Location;

/// Returns the Royal Observatory Greenwich location.
///
/// Used as the default reference point for TT-TDB conversions when no
/// observer location is specified.
fn greenwich_location() -> Location {
    Location::from_degrees(51.477928, 0.0, 46.0).expect("Greenwich coordinates should be valid")
}

const DEFAULT_UT1_OFFSET_SECONDS: f64 = 0.0;

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
    let (u, v) = location.to_geocentric_km()?;

    let dtr = calculate_tdb_tt_difference(
        date_jd.jd1(),
        date_jd.jd2(),
        ut1_fraction,
        location.longitude,
        u,
        v,
    );

    Ok(dtr)
}

/// Convert a time scale to TDB (Barycentric Dynamical Time).
///
/// Implemented for TT and TDB. TT conversion requires observer location; TDB-to-TDB
/// is an identity operation.
///
/// # Methods
///
/// - `to_tdb_greenwich()` - Convert using Greenwich as reference location
/// - `to_tdb_with_location()` - Convert using explicit observer location
/// - `to_tdb_with_location_and_ut1_offset()` - Convert with location and UT1-UTC offset
/// - `to_tdb_with_offset()` - Convert using pre-computed TDB-TT offset in seconds
pub trait ToTDB {
    /// Convert to TDB using Greenwich Observatory as reference.
    fn to_tdb_greenwich(&self) -> TimeResult<TDB>;
    /// Convert to TDB using the specified observer location.
    fn to_tdb_with_location(&self, location: &Location) -> TimeResult<TDB>;
    /// Convert to TDB with location and UT1-UTC offset for maximum precision.
    fn to_tdb_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_offset_seconds: f64,
    ) -> TimeResult<TDB>;
    /// Convert to TDB using a pre-computed offset (TDB-TT) in seconds.
    fn to_tdb_with_offset(&self, dtr_seconds: f64) -> TimeResult<TDB>;
}

impl ToTDB for TDB {
    fn to_tdb_greenwich(&self) -> TimeResult<TDB> {
        Ok(*self)
    }

    fn to_tdb_with_location(&self, _location: &Location) -> TimeResult<TDB> {
        Ok(*self)
    }

    fn to_tdb_with_location_and_ut1_offset(
        &self,
        _location: &Location,
        _ut1_offset_seconds: f64,
    ) -> TimeResult<TDB> {
        Ok(*self)
    }

    fn to_tdb_with_offset(&self, _dtr_seconds: f64) -> TimeResult<TDB> {
        Ok(*self)
    }
}

impl ToTDB for TT {
    fn to_tdb_greenwich(&self) -> TimeResult<TDB> {
        let location = greenwich_location();
        self.to_tdb_with_location_and_ut1_offset(&location, DEFAULT_UT1_OFFSET_SECONDS)
    }

    fn to_tdb_with_location(&self, location: &Location) -> TimeResult<TDB> {
        self.to_tdb_with_location_and_ut1_offset(location, DEFAULT_UT1_OFFSET_SECONDS)
    }

    fn to_tdb_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_offset_seconds: f64,
    ) -> TimeResult<TDB> {
        let tt_jd = self.to_julian_date();

        let tt_f64 = tt_jd.to_f64();
        let ut1_fraction = ((tt_f64 - libm::trunc(tt_f64)) * SECONDS_PER_DAY_F64
            + ut1_offset_seconds)
            / SECONDS_PER_DAY_F64;
        let ut1_fraction = ut1_fraction - libm::floor(ut1_fraction);

        let dtr = compute_tdb_tt_offset(&tt_jd, ut1_fraction, location)?;
        self.to_tdb_with_offset(dtr)
    }

    fn to_tdb_with_offset(&self, dtr_seconds: f64) -> TimeResult<TDB> {
        let tt_jd = self.to_julian_date();

        let dtr_days = dtr_seconds / celestial_core::constants::SECONDS_PER_DAY_F64;

        let (tdb_jd1, tdb_jd2) = if tt_jd.jd1().abs() > tt_jd.jd2().abs() {
            (tt_jd.jd1(), tt_jd.jd2() + dtr_days)
        } else {
            (tt_jd.jd1() + dtr_days, tt_jd.jd2())
        };

        let tdb_jd = JulianDate::new(tdb_jd1, tdb_jd2);
        Ok(TDB::from_julian_date(tdb_jd))
    }
}

impl ToTT for TDB {
    fn to_tt(&self) -> TimeResult<TT> {
        Err(TimeError::ConversionError(
            "TDB→TT conversion requires observer location. \
             Use to_tt_greenwich() for Greenwich or to_tt_with_location() for other locations."
                .to_string(),
        ))
    }
}

/// Convert TDB to TT with location awareness.
///
/// This trait exists because the generic [`ToTT`] trait cannot provide accurate
/// TDB→TT conversion without knowing the observer's location. Rather than
/// silently use a default, the design requires explicit location specification.
///
/// # Methods
///
/// - `to_tt_greenwich()` - Convert using Greenwich as reference location
/// - `to_tt_with_location()` - Convert using explicit observer location
/// - `to_tt_with_location_and_ut1_offset()` - Convert with location and UT1-UTC offset
/// - `to_tt_with_offset()` - Convert using pre-computed TDB-TT offset in seconds
///
/// # Inverse Operation
///
/// The TDB→TT conversion uses iterative refinement because the offset depends on
/// TT (which we're solving for). Three iterations achieve sub-nanosecond convergence.
pub trait ToTTFromTDB {
    /// Convert to TT using Greenwich Observatory as reference.
    fn to_tt_greenwich(&self) -> TimeResult<TT>;
    /// Convert to TT using the specified observer location.
    fn to_tt_with_location(&self, location: &Location) -> TimeResult<TT>;
    /// Convert to TT with location and UT1-UTC offset for maximum precision.
    fn to_tt_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_offset_seconds: f64,
    ) -> TimeResult<TT>;
    /// Convert to TT using a pre-computed offset (TDB-TT) in seconds.
    fn to_tt_with_offset(&self, dtr_seconds: f64) -> TimeResult<TT>;
}

impl ToTTFromTDB for TDB {
    fn to_tt_greenwich(&self) -> TimeResult<TT> {
        let location = greenwich_location();
        self.to_tt_with_location_and_ut1_offset(&location, DEFAULT_UT1_OFFSET_SECONDS)
    }

    fn to_tt_with_location(&self, location: &Location) -> TimeResult<TT> {
        self.to_tt_with_location_and_ut1_offset(location, DEFAULT_UT1_OFFSET_SECONDS)
    }

    fn to_tt_with_location_and_ut1_offset(
        &self,
        location: &Location,
        ut1_offset_seconds: f64,
    ) -> TimeResult<TT> {
        let tdb_jd = self.to_julian_date();

        let tdb_f64 = tdb_jd.to_f64();
        let ut1_fraction_approx = ((tdb_f64 - libm::trunc(tdb_f64)) * SECONDS_PER_DAY_F64
            + ut1_offset_seconds)
            / SECONDS_PER_DAY_F64;
        let ut1_fraction_approx = ut1_fraction_approx - libm::floor(ut1_fraction_approx);

        let dtr_approx = compute_tdb_tt_offset(&tdb_jd, ut1_fraction_approx, location)?;

        let tt_approx = self.to_tt_with_offset(dtr_approx)?;
        let tt_jd_approx = tt_approx.to_julian_date();

        let tt_jd_approx_f64 = tt_jd_approx.jd1() + tt_jd_approx.jd2();
        let ut1_fraction = ((tt_jd_approx_f64 - libm::trunc(tt_jd_approx_f64))
            * SECONDS_PER_DAY_F64
            + ut1_offset_seconds)
            / SECONDS_PER_DAY_F64;
        let ut1_fraction = ut1_fraction - libm::floor(ut1_fraction);

        let dtr_refined = compute_tdb_tt_offset(&tt_jd_approx, ut1_fraction, location)?;

        let tt_refined = self.to_tt_with_offset(dtr_refined)?;
        let tt_jd_refined = tt_refined.to_julian_date();

        let tt_jd_refined_f64 = tt_jd_refined.jd1() + tt_jd_refined.jd2();
        let ut1_fraction_final = ((tt_jd_refined_f64 - libm::trunc(tt_jd_refined_f64))
            * SECONDS_PER_DAY_F64
            + ut1_offset_seconds)
            / SECONDS_PER_DAY_F64;
        let ut1_fraction_final = ut1_fraction_final - libm::floor(ut1_fraction_final);

        let dtr_final = compute_tdb_tt_offset(&tt_jd_refined, ut1_fraction_final, location)?;

        self.to_tt_with_offset(dtr_final)
    }

    fn to_tt_with_offset(&self, dtr_seconds: f64) -> TimeResult<TT> {
        let tdb_jd = self.to_julian_date();

        let dtr_days = dtr_seconds / celestial_core::constants::SECONDS_PER_DAY_F64;

        let (tt_jd1, tt_jd2) = if tdb_jd.jd1().abs() > tdb_jd.jd2().abs() {
            (tdb_jd.jd1(), tdb_jd.jd2() - dtr_days)
        } else {
            (tdb_jd.jd1() - dtr_days, tdb_jd.jd2())
        };

        let tt_jd = JulianDate::new(tt_jd1, tt_jd2);
        Ok(TT::from_julian_date(tt_jd))
    }
}

/// Compute TDB-TT difference using Fairhead & Bretagnon (1990) model.
///
/// This implements the full 787-term series for the TDB-TT difference, including:
/// - Fundamental arguments (Sun, Moon, planets mean longitudes)
/// - Diurnal terms from observer's geocentric position
/// - Jupiter/Saturn secular terms
///
/// # Arguments
///
/// - `date1`, `date2`: Two-part Julian Date
/// - `ut`: UT1 fraction of day (for local sidereal time calculation)
/// - `elong`: Observer's east longitude in radians
/// - `u`, `v`: Geocentric cylindrical coordinates in km (from Location::to_geocentric_km)
///
/// # Returns
///
/// TDB - TT in seconds. Range is approximately -0.00166 to +0.00166 seconds.
fn calculate_tdb_tt_difference(date1: f64, date2: f64, ut: f64, elong: f64, u: f64, v: f64) -> f64 {
    let t = ((date1 - J2000_JD) + date2) / DAYS_PER_JULIAN_MILLENNIUM;

    let tsol = fmod(ut, 1.0) * TWOPI + elong;

    let w = t / 3600.0;

    let elsun = fmod(280.46645683 + 1296027711.03429 * w, 360.0) * DEG_TO_RAD;

    let emsun = fmod(357.52910918 + 1295965810.481 * w, 360.0) * DEG_TO_RAD;

    let d = fmod(297.85019547 + 16029616012.090 * w, 360.0) * DEG_TO_RAD;

    let elj = fmod(34.35151874 + 109306899.89453 * w, 360.0) * DEG_TO_RAD;

    let els = fmod(50.07744430 + 44046398.47038 * w, 360.0) * DEG_TO_RAD;

    let wt = 0.00029e-10 * u * libm::sin(tsol + elsun - els)
        + 0.00100e-10 * u * libm::sin(tsol - 2.0 * emsun)
        + 0.00133e-10 * u * libm::sin(tsol - d)
        + 0.00133e-10 * u * libm::sin(tsol + elsun - elj)
        - 0.00229e-10 * u * libm::sin(tsol + 2.0 * elsun + emsun)
        - 0.02200e-10 * v * libm::cos(elsun + emsun)
        + 0.05312e-10 * u * libm::sin(tsol - emsun)
        - 0.13677e-10 * u * libm::sin(tsol + 2.0 * elsun)
        - 1.31840e-10 * v * libm::cos(elsun)
        + 3.17679e-10 * u * libm::sin(tsol);

    let mut w0 = 0.0;
    for j in (0..474).rev() {
        w0 += FAIRHD[j][0] * libm::sin(FAIRHD[j][1] * t + FAIRHD[j][2]);
    }

    let mut w1 = 0.0;
    for j in (474..679).rev() {
        w1 += FAIRHD[j][0] * libm::sin(FAIRHD[j][1] * t + FAIRHD[j][2]);
    }

    let mut w2 = 0.0;
    for j in (679..764).rev() {
        if FAIRHD[j][0] != 0.0 {
            w2 += FAIRHD[j][0] * libm::sin(FAIRHD[j][1] * t + FAIRHD[j][2]);
        }
    }

    let mut w3 = 0.0;
    for j in (764..784).rev() {
        if FAIRHD[j][0] != 0.0 {
            w3 += FAIRHD[j][0] * libm::sin(FAIRHD[j][1] * t + FAIRHD[j][2]);
        }
    }

    let mut w4 = 0.0;
    for j in (784..787).rev() {
        if FAIRHD[j][0] != 0.0 {
            w4 += FAIRHD[j][0] * libm::sin(FAIRHD[j][1] * t + FAIRHD[j][2]);
        }
    }

    let wf = t * (t * (t * (t * w4 + w3) + w2) + w1) + w0;

    let wj = 0.00065e-6 * libm::sin(6069.776754 * t + 4.021194)
        + 0.00033e-6 * libm::sin(213.299095 * t + 5.543132)
        - 0.00196e-6 * libm::sin(6208.294251 * t + 5.696701)
        - 0.00173e-6 * libm::sin(74.781599 * t + 2.435900)
        + 0.03638e-6 * t * t;

    wt + wf + wj
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD};

    #[test]
    fn test_identity_conversions() {
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let location = Location::from_degrees(45.0, 90.0, 100.0).unwrap();

        let via_greenwich = tdb.to_tdb_greenwich().unwrap();
        assert_eq!(
            tdb.to_julian_date().jd1(),
            via_greenwich.to_julian_date().jd1(),
            "TDB→TDB via Greenwich should preserve JD1"
        );
        assert_eq!(
            tdb.to_julian_date().jd2(),
            via_greenwich.to_julian_date().jd2(),
            "TDB→TDB via Greenwich should preserve JD2"
        );

        let via_location = tdb.to_tdb_with_location(&location).unwrap();
        assert_eq!(
            tdb.to_julian_date().jd1(),
            via_location.to_julian_date().jd1(),
            "TDB→TDB via location should preserve JD1"
        );
        assert_eq!(
            tdb.to_julian_date().jd2(),
            via_location.to_julian_date().jd2(),
            "TDB→TDB via location should preserve JD2"
        );

        let via_ut1_offset = tdb
            .to_tdb_with_location_and_ut1_offset(&location, 0.3)
            .unwrap();
        assert_eq!(
            tdb.to_julian_date().jd1(),
            via_ut1_offset.to_julian_date().jd1(),
            "TDB→TDB via UT1 offset should preserve JD1"
        );
        assert_eq!(
            tdb.to_julian_date().jd2(),
            via_ut1_offset.to_julian_date().jd2(),
            "TDB→TDB via UT1 offset should preserve JD2"
        );

        let via_offset = tdb.to_tdb_with_offset(0.001).unwrap();
        assert_eq!(
            tdb.to_julian_date().jd1(),
            via_offset.to_julian_date().jd1(),
            "TDB→TDB via offset should preserve JD1"
        );
        assert_eq!(
            tdb.to_julian_date().jd2(),
            via_offset.to_julian_date().jd2(),
            "TDB→TDB via offset should preserve JD2"
        );
    }

    #[test]
    fn test_tt_tdb_offset_verification() {
        let test_cases = [
            (J2000_JD, "J2000.0", greenwich_location()),
            (
                J2000_JD - DAYS_PER_JULIAN_CENTURY,
                "1900",
                greenwich_location(),
            ),
            (J2000_JD + 18262.5, "2050", greenwich_location()),
            (
                J2000_JD + DAYS_PER_JULIAN_CENTURY,
                "2100",
                greenwich_location(),
            ),
            (
                J2000_JD,
                "J2000 Tokyo",
                Location::from_degrees(35.6762, 139.6503, 40.0).unwrap(),
            ),
            (
                J2000_JD,
                "J2000 Sydney",
                Location::from_degrees(-33.8688, 151.2093, 58.0).unwrap(),
            ),
        ];

        for (jd, description, location) in test_cases {
            let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
            let tdb = tt.to_tdb_with_location(&location).unwrap();

            let diff_seconds = (tdb.to_julian_date().to_f64() - tt.to_julian_date().to_f64())
                * SECONDS_PER_DAY_F64;

            assert!(
                diff_seconds.abs() < 0.002,
                "{}: TT→TDB offset should be < 2ms, got {:.6} seconds",
                description,
                diff_seconds
            );

            let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
            let tt = tdb.to_tt_with_location(&location).unwrap();

            let diff_seconds = (tdb.to_julian_date().to_f64() - tt.to_julian_date().to_f64())
                * SECONDS_PER_DAY_F64;

            assert!(
                diff_seconds.abs() < 0.002,
                "{}: TDB→TT offset should be < 2ms, got {:.6} seconds",
                description,
                diff_seconds
            );
        }
    }

    #[test]
    fn test_tt_tdb_round_trip_precision() {
        const TOLERANCE_DAYS: f64 = 1e-11;

        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345];

        for jd2 in test_jd2_values {
            let original_tt = TT::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tdb = original_tt.to_tdb_greenwich().unwrap();
            let round_trip_tt = tdb.to_tt_greenwich().unwrap();

            let jd1_diff =
                (original_tt.to_julian_date().jd1() - round_trip_tt.to_julian_date().jd1()).abs();
            let jd2_diff =
                (original_tt.to_julian_date().jd2() - round_trip_tt.to_julian_date().jd2()).abs();

            assert!(
                jd1_diff < TOLERANCE_DAYS,
                "TT→TDB→TT JD1 diff {} exceeds tolerance for jd2={}",
                jd1_diff,
                jd2
            );
            assert!(
                jd2_diff < TOLERANCE_DAYS,
                "TT→TDB→TT JD2 diff {} exceeds tolerance for jd2={}",
                jd2_diff,
                jd2
            );

            let original_tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tt = original_tdb.to_tt_greenwich().unwrap();
            let round_trip_tdb = tt.to_tdb_greenwich().unwrap();

            let jd1_diff =
                (original_tdb.to_julian_date().jd1() - round_trip_tdb.to_julian_date().jd1()).abs();
            let jd2_diff =
                (original_tdb.to_julian_date().jd2() - round_trip_tdb.to_julian_date().jd2()).abs();

            assert!(
                jd1_diff < TOLERANCE_DAYS,
                "TDB→TT→TDB JD1 diff {} exceeds tolerance for jd2={}",
                jd1_diff,
                jd2
            );
            assert!(
                jd2_diff < TOLERANCE_DAYS,
                "TDB→TT→TDB JD2 diff {} exceeds tolerance for jd2={}",
                jd2_diff,
                jd2
            );
        }

        let alt_tt = TT::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tdb = alt_tt.to_tdb_greenwich().unwrap();
        let alt_round_trip = alt_tdb.to_tt_greenwich().unwrap();

        let jd1_diff =
            (alt_tt.to_julian_date().jd1() - alt_round_trip.to_julian_date().jd1()).abs();
        let jd2_diff =
            (alt_tt.to_julian_date().jd2() - alt_round_trip.to_julian_date().jd2()).abs();

        assert!(
            jd1_diff < TOLERANCE_DAYS,
            "Alternate split TT→TDB→TT JD1 diff {} exceeds tolerance",
            jd1_diff
        );
        assert!(
            jd2_diff < TOLERANCE_DAYS,
            "Alternate split TT→TDB→TT JD2 diff {} exceeds tolerance",
            jd2_diff
        );
    }

    #[test]
    fn test_api_equivalence() {
        let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.123456));

        let tdb_greenwich = tt.to_tdb_greenwich().unwrap();
        let tdb_explicit = tt.to_tdb_with_location(&greenwich_location()).unwrap();

        assert_eq!(
            tdb_greenwich.to_julian_date().jd1(),
            tdb_explicit.to_julian_date().jd1(),
            "TT: to_tdb_greenwich() should match to_tdb_with_location(greenwich) JD1"
        );
        assert_eq!(
            tdb_greenwich.to_julian_date().jd2(),
            tdb_explicit.to_julian_date().jd2(),
            "TT: to_tdb_greenwich() should match to_tdb_with_location(greenwich) JD2"
        );

        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.987654));

        let tt_greenwich = tdb.to_tt_greenwich().unwrap();
        let tt_explicit = tdb.to_tt_with_location(&greenwich_location()).unwrap();

        assert_eq!(
            tt_greenwich.to_julian_date().jd1(),
            tt_explicit.to_julian_date().jd1(),
            "TDB: to_tt_greenwich() should match to_tt_with_location(greenwich) JD1"
        );
        assert_eq!(
            tt_greenwich.to_julian_date().jd2(),
            tt_explicit.to_julian_date().jd2(),
            "TDB: to_tt_greenwich() should match to_tt_with_location(greenwich) JD2"
        );
    }

    #[test]
    fn test_tdb_to_tt_requires_location() {
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let result = tdb.to_tt();

        assert!(result.is_err(), "TDB.to_tt() should return error");

        match result {
            Err(TimeError::ConversionError(msg)) => {
                assert!(
                    msg.contains("location"),
                    "Error message should mention location requirement: {}",
                    msg
                );
                assert!(
                    msg.contains("greenwich") || msg.contains("Greenwich"),
                    "Error message should mention Greenwich option: {}",
                    msg
                );
            }
            _ => panic!("Expected ConversionError, got {:?}", result),
        }
    }
}
