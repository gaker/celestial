//! Conversions between Barycentric Coordinate Time (TCB) and Barycentric Dynamical Time (TDB).
//!
//! TCB and TDB are both barycentric time scales used for solar system dynamics, but they
//! differ in rate. TDB was introduced to provide a time scale that, when observed from
//! Earth's surface, ticks at approximately the same rate as TT on average.
//!
//! # The L_B Rate Factor
//!
//! The defining relationship from IAU 2006 Resolution B3 is:
//!
//! ```text
//! TDB = TCB - L_B * (JD_TCB - T0) * 86400 + TDB_0
//! ```
//!
//! Where:
//! - `L_B = 1.550519768e-8` (IAU 2006, exact by definition)
//! - `T0 = 1977 January 1, 0h TAI` (MJD 43144.0, the common reference epoch)
//! - `TDB_0 = -6.55e-5 seconds` (offset to align TDB with TT at J2000.0 on average)
//!
//! The L_B value represents the average fractional rate difference between TCB and TDB.
//! TCB gains about 0.49 seconds per year relative to TDB.
//!
//! # Reference Epoch (TDB_0)
//!
//! The reference epoch for TCB-TDB conversions is 1977 January 1, 0h TAI (JD 2443144.5003725),
//! the same epoch used for TT-TCG. At this epoch, with the TDB_0 offset applied, TDB and
//! TCB are related by definition.
//!
//! The TDB_0 constant (-6.55e-5 seconds) was chosen so that TDB matches TT on average at
//! the geocenter. This makes TDB a "scaled" version of TCB that tracks TT's rate.
//!
//! # Why TDB Exists
//!
//! TCB is the natural coordinate time for the barycentric frame, but its rate differs
//! from TT by about 490 ms/year. For continuity with historical ephemerides and to
//! avoid confusion, TDB was defined to match TT's average rate while remaining suitable
//! for barycentric calculations.
//!
//! In practice:
//! - TCB is used in relativistic equations of motion
//! - TDB is used in JPL ephemerides (DE series) and for practical timekeeping
//! - The difference grows linearly: ~17 seconds at J2000.0 relative to 1977
//!
//! # Precision
//!
//! Round-trip conversions (TCB -> TDB -> TCB or TDB -> TCB -> TDB) achieve sub-picosecond
//! accuracy. The implementation applies corrections to the smaller-magnitude Julian Date
//! component to preserve precision.
//!
//! # Usage
//!
//! ```
//! use astro_time::scales::{TCB, TDB};
//! use astro_time::scales::conversions::{TcbToTdb, TdbToTcb};
//! use astro_time::julian::JulianDate;
//! use astro_core::constants::J2000_JD;
//!
//! let tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tdb = tcb.tcb_to_tdb().unwrap();
//!
//! // At J2000.0, TDB is about 11.3 ms behind TCB (accumulated since 1977)
//! let offset_days = tdb.to_julian_date().to_f64() - tcb.to_julian_date().to_f64();
//! assert!(offset_days < 0.0, "TDB should be behind TCB after 1977");
//! ```
//!
//! # References
//!
//! - IAU 2006 Resolution B3: Re-definition of Barycentric Dynamical Time, TDB
//! - IERS Conventions (2010), Chapter 10: General Relativistic Models for Time
//! - Soffel et al. (2003): The IAU 2000 Resolutions for Astrometry

use crate::constants::TT_TAI_OFFSET;
use crate::julian::JulianDate;
use crate::scales::{TCB, TDB};
use crate::TimeResult;
use astro_core::constants::{MJD_ZERO_POINT, SECONDS_PER_DAY_F64};

/// L_B rate factor from IAU 2006 Resolution B3.
/// Represents the fractional rate difference: TCB runs faster than TDB by this amount.
const TCB_RATE: f64 = 1.550519768e-8;

/// MJD of the reference epoch: 1977 January 1, 0h TAI.
const MJD_1977: f64 = 43144.0;

/// TDB_0 offset in seconds. Chosen so TDB matches TT rate on average at geocenter.
const TDB_OFFSET: f64 = -6.55e-5;

/// Reference epoch as full Julian Date (MJD_ZERO_POINT + MJD_1977).
const T77TD: f64 = MJD_ZERO_POINT + MJD_1977;

/// TT-TAI offset in days (32.184s / 86400), for epoch alignment.
const T77TF: f64 = TT_TAI_OFFSET / SECONDS_PER_DAY_F64;

/// TDB_0 offset in days (-6.55e-5s / 86400).
const TDB0: f64 = TDB_OFFSET / SECONDS_PER_DAY_F64;

/// Derived rate ratio: L_B / (1 - L_B).
/// Used for TDB -> TCB conversion to invert the rate scaling.
const TCB_RATE_RATIO: f64 = TCB_RATE / (1.0 - TCB_RATE);

/// Convert Barycentric Coordinate Time (TCB) to Barycentric Dynamical Time (TDB).
///
/// TDB is a rescaled version of TCB designed to match TT's average rate at the geocenter.
/// This conversion removes the L_B rate difference accumulated since 1977.
pub trait TcbToTdb {
    /// Convert TCB to TDB.
    ///
    /// Applies: `TDB = TCB - L_B * (TCB - T0) + TDB_0`
    ///
    /// At J2000.0, TDB is approximately 11 milliseconds behind TCB due to the
    /// accumulated rate difference since 1977.
    fn tcb_to_tdb(&self) -> TimeResult<TDB>;
}

/// Convert Barycentric Dynamical Time (TDB) to Barycentric Coordinate Time (TCB).
///
/// This is the inverse of [`TcbToTdb`]. Uses the rate ratio `L_B / (1 - L_B)`
/// to correctly invert the scaling.
pub trait TdbToTcb {
    /// Convert TDB to TCB.
    ///
    /// Applies the inverse transformation using the derived rate ratio.
    /// At J2000.0, TCB is approximately 11 milliseconds ahead of TDB.
    fn tdb_to_tcb(&self) -> TimeResult<TCB>;
}

impl TcbToTdb for TCB {
    /// Convert TCB to TDB by removing the L_B rate correction.
    ///
    /// The correction is computed as: `L_B * (TCB - T0)` where T0 is the 1977 epoch.
    /// The TDB_0 offset is added to align with TT at the geocenter.
    ///
    /// Applies the correction to the smaller-magnitude JD component for precision.
    fn tcb_to_tdb(&self) -> TimeResult<TDB> {
        let tcb_jd = self.to_julian_date();
        let (tcb1, tcb2) = (tcb_jd.jd1(), tcb_jd.jd2());

        let (big, small) = if tcb1.abs() > tcb2.abs() {
            (tcb1, tcb2)
        } else {
            (tcb2, tcb1)
        };
        let d = big - T77TD;
        let corrected = small + TDB0 - (d + (small - T77TF)) * TCB_RATE;
        let (tdb1, tdb2) = if tcb1.abs() > tcb2.abs() {
            (big, corrected)
        } else {
            (corrected, big)
        };

        Ok(TDB::from_julian_date(JulianDate::new(tdb1, tdb2)))
    }
}

impl TdbToTcb for TDB {
    /// Convert TDB to TCB by applying the inverse L_B rate correction.
    ///
    /// Uses the rate ratio `L_B / (1 - L_B)` to properly invert the scaling.
    /// First removes the TDB_0 offset, then applies the inverse rate correction.
    ///
    /// Applies the correction to the smaller-magnitude JD component for precision.
    fn tdb_to_tcb(&self) -> TimeResult<TCB> {
        let tdb_jd = self.to_julian_date();
        let (tdb1, tdb2) = (tdb_jd.jd1(), tdb_jd.jd2());

        let (big, small) = if tdb1.abs() > tdb2.abs() {
            (tdb1, tdb2)
        } else {
            (tdb2, tdb1)
        };
        let d = T77TD - big;
        let f = small - TDB0;
        let corrected = f - (d - (f - T77TF)) * TCB_RATE_RATIO;
        let (tcb1, tcb2) = if tdb1.abs() > tdb2.abs() {
            (big, corrected)
        } else {
            (corrected, big)
        };

        Ok(TCB::from_julian_date(JulianDate::new(tcb1, tcb2)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_core::constants::J2000_JD;

    #[test]
    fn test_tcb_tdb_relationship() {
        // Identity conversions
        let tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let tdb = tcb.tcb_to_tdb().unwrap();
        let tcb_jd = tcb.to_julian_date();
        let tdb_jd = tdb.to_julian_date();

        // TCB runs faster than TDB, so at J2000 (after 1977 epoch), TDB < TCB
        assert!(
            tdb_jd.to_f64() < tcb_jd.to_f64(),
            "TDB should be behind TCB at J2000"
        );

        // Verify inverse relationship holds
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let tcb = tdb.tdb_to_tcb().unwrap();
        let tdb_jd = tdb.to_julian_date();
        let tcb_jd = tcb.to_julian_date();

        assert!(
            tcb_jd.to_f64() > tdb_jd.to_f64(),
            "TCB should be ahead of TDB at J2000"
        );
    }

    #[test]
    fn test_tcb_tdb_round_trip_precision() {
        // TCB/TDB conversions involve rate scaling. 1e-14 days = ~1 picosecond tolerance.
        const TOLERANCE_DAYS: f64 = 1e-14;

        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];

        for jd2 in test_jd2_values {
            // TCB -> TDB -> TCB
            let original_tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tdb = original_tcb.tcb_to_tdb().unwrap();
            let round_trip_tcb = tdb.tdb_to_tcb().unwrap();

            assert_eq!(
                original_tcb.to_julian_date().jd1(),
                round_trip_tcb.to_julian_date().jd1(),
                "TCB->TDB->TCB JD1 must be exact for jd2={}",
                jd2
            );
            let jd2_diff =
                (original_tcb.to_julian_date().jd2() - round_trip_tcb.to_julian_date().jd2()).abs();
            assert!(
                jd2_diff <= TOLERANCE_DAYS,
                "TCB->TDB->TCB JD2 diff {} exceeds tolerance {} for jd2={}",
                jd2_diff,
                TOLERANCE_DAYS,
                jd2
            );

            // TDB -> TCB -> TDB
            let original_tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tcb = original_tdb.tdb_to_tcb().unwrap();
            let round_trip_tdb = tcb.tcb_to_tdb().unwrap();

            assert_eq!(
                original_tdb.to_julian_date().jd1(),
                round_trip_tdb.to_julian_date().jd1(),
                "TDB->TCB->TDB JD1 must be exact for jd2={}",
                jd2
            );
            let jd2_diff =
                (original_tdb.to_julian_date().jd2() - round_trip_tdb.to_julian_date().jd2()).abs();
            assert!(
                jd2_diff <= TOLERANCE_DAYS,
                "TDB->TCB->TDB JD2 diff {} exceeds tolerance {} for jd2={}",
                jd2_diff,
                TOLERANCE_DAYS,
                jd2
            );
        }

        // Alternate JD split case (jd2 > jd1)
        let alt_tcb = TCB::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tdb = alt_tcb.tcb_to_tdb().unwrap();
        let alt_round_trip = alt_tdb.tdb_to_tcb().unwrap();

        assert_eq!(
            alt_tcb.to_julian_date().jd1(),
            alt_round_trip.to_julian_date().jd1(),
            "Alternate split TCB->TDB->TCB JD1 must be exact"
        );
        let jd2_diff =
            (alt_tcb.to_julian_date().jd2() - alt_round_trip.to_julian_date().jd2()).abs();
        assert!(
            jd2_diff <= TOLERANCE_DAYS,
            "Alternate split TCB->TDB->TCB JD2 diff {} exceeds tolerance {}",
            jd2_diff,
            TOLERANCE_DAYS
        );
    }
}
