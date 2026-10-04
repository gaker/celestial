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
//! - `T0 = 1977 January 1, 0h TAI` (MJD 43144.0003725 in TT/TCG/TCB, the common reference epoch)
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
//! - The difference grows linearly: ~11.25 seconds at J2000.0 relative to 1977
//!
//! # Precision
//!
//! Round-trip conversions (TCB -> TDB -> TCB or TDB -> TCB -> TDB) come back within
//! 1 ulp of jd2 (at most ~10 ps). The implementation applies corrections to the smaller-magnitude Julian Date
//! component to preserve precision.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::tcb::TCB;
//! use celestial_time::scales::tdb::TDB;
//! use celestial_time::scales::conversions::tcb_tdb::{TcbToTdb, TdbToTcb};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! let tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tdb = tcb.tcb_to_tdb().unwrap();
//!
//! // At J2000.0, TDB is about 11.25 s behind TCB (accumulated since 1977)
//! assert_eq!(tdb.to_julian_date(), JulianDate::new(J2000_JD, -0.00013025216543700573));
//! ```
//!
//! # References
//!
//! - IAU 2006 Resolution B3: Re-definition of Barycentric Dynamical Time, TDB
//! - IERS Conventions (2010), Chapter 10: General Relativistic Models for Time
//! - Soffel et al. (2003): The IAU 2000 Resolutions for Astrometry

use crate::constants::{
    MJD_1977_JAN_1, TCB_RATE_LB, TCB_RATE_RATIO, TDB_OFFSET_1977, TT_TAI_OFFSET,
};
use crate::julian::finite_jd;
use crate::scales::tcb::TCB;
use crate::scales::tdb::TDB;
use crate::TimeResult;
use celestial_core::constants::{MJD_ZERO_POINT, SECONDS_PER_DAY_F64};

/// Reference epoch as full Julian Date (MJD_ZERO_POINT + MJD_1977_JAN_1).
const T77TD: f64 = MJD_ZERO_POINT + MJD_1977_JAN_1;

/// TT-TAI offset in days (32.184s / 86400), for epoch alignment.
const T77TF: f64 = TT_TAI_OFFSET / SECONDS_PER_DAY_F64;

/// TDB_0 offset in days (-6.55e-5s / 86400). L_B matches the rate; this offset
/// keeps TDB-TT near zero at the geocenter.
const TDB0: f64 = TDB_OFFSET_1977 / SECONDS_PER_DAY_F64;

/// Convert Barycentric Coordinate Time (TCB) to Barycentric Dynamical Time (TDB).
///
/// TDB is a rescaled version of TCB designed to match TT's average rate at the geocenter.
/// This conversion removes the L_B rate difference accumulated since 1977.
pub trait TcbToTdb {
    /// Convert TCB to TDB.
    ///
    /// Applies: `TDB = TCB - L_B * (TCB - T0) + TDB_0`
    ///
    /// At J2000.0, TDB is approximately 11.25 seconds behind TCB due to the
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
    /// At J2000.0, TCB is approximately 11.25 seconds ahead of TDB.
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
        let tcb_jd = finite_jd(self.to_julian_date())?;
        let tdb_jd = tcb_jd.map_smaller_part(|big, small| {
            small + TDB0 - ((big - T77TD) + (small - T77TF)) * TCB_RATE_LB
        });
        Ok(TDB::from_julian_date(tdb_jd))
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
        let tdb_jd = finite_jd(self.to_julian_date())?;
        let tcb_jd = tdb_jd.map_smaller_part(|big, small| {
            let f = small - TDB0;
            f - ((T77TD - big) - (f - T77TF)) * TCB_RATE_RATIO
        });
        Ok(TCB::from_julian_date(tcb_jd))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use crate::TimeError;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_tcb_tdb_matches_erfa() {
        // eraTcbtdb and eraTdbtcb at J2000.0, where TCB is about 11.25 s ahead.
        let tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        assert_eq!(
            tcb.tcb_to_tdb().unwrap().to_julian_date().parts(),
            (J2000_JD, -0.00013025216543700573)
        );
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        assert_eq!(
            tdb.tdb_to_tcb().unwrap().to_julian_date().parts(),
            (J2000_JD, 0.00013025216745659132)
        );
    }

    #[test]
    fn test_tcb_tdb_round_trip_matches_erfa() {
        // eraTcbtdb then eraTdbtcb, and the reverse. Only TDB at the (JD, 0) split
        // doesn't come back exactly, by under 3e-20 d.
        let cases = [
            (J2000_JD, 0.0, 0.0, -2.710505431213761e-20),
            (J2000_JD, 0.5, 0.5, 0.5),
            (
                J2000_JD,
                0.123456789012345,
                0.123456789012345,
                0.123456789012345,
            ),
            (
                J2000_JD,
                -0.123456789012345,
                -0.123456789012345,
                -0.123456789012345,
            ),
            (J2000_JD, 0.987654321, 0.987654321, 0.987654321),
        ];
        for (jd1, jd2, tcb_back, tdb_back) in cases {
            let tcb = TCB::from_julian_date(JulianDate::new(jd1, jd2));
            let back = tcb.tcb_to_tdb().unwrap().tdb_to_tcb().unwrap();
            assert_eq!(back.to_julian_date().parts(), (jd1, tcb_back));

            let tdb = TDB::from_julian_date(JulianDate::new(jd1, jd2));
            let back = tdb.tdb_to_tcb().unwrap().tcb_to_tdb().unwrap();
            assert_eq!(back.to_julian_date().parts(), (jd1, tdb_back));
        }

        let tcb = TCB::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let back = tcb.tcb_to_tdb().unwrap().tdb_to_tcb().unwrap();
        assert_eq!(back.to_julian_date().parts(), tcb.to_julian_date().parts());
        let tdb = TDB::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let back = tdb.tdb_to_tcb().unwrap().tcb_to_tdb().unwrap();
        assert_eq!(back.to_julian_date().parts(), tdb.to_julian_date().parts());
    }

    #[test]
    fn test_non_finite_julian_date_is_rejected() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let jd = JulianDate::new(bad, 0.0);
            let expected =
                TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
            assert_eq!(
                TCB::from_julian_date(jd).tcb_to_tdb().unwrap_err(),
                expected
            );
            assert_eq!(
                TDB::from_julian_date(jd).tdb_to_tcb().unwrap_err(),
                expected
            );
        }
    }
}
