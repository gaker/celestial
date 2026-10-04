//! Conversions between TAI, TT, and TCG time scales.
//!
//! This module implements the fixed-offset and linear-rate conversions between:
//!
//! - **TAI (International Atomic Time)**: The reference atomic time scale.
//! - **TT (Terrestrial Time)**: Idealized time on the geoid. TT = TAI + 32.184s exactly.
//! - **TCG (Geocentric Coordinate Time)**: Coordinate time at the geocenter.
//!
//! # Conversion Relationships
//!
//! ```text
//! TAI <-> TT     Fixed offset: TT = TAI + 32.184 seconds
//! TAI <-> TCG   Chains through TT: TAI → TT → TCG
//! ```
//!
//! The TAI-TT offset is defined by the IAU to be exactly 32.184 seconds. This offset
//! accounts for the historical difference between atomic time and ephemeris time.
//!
//! # Precision Preservation
//!
//! All conversions add offsets to the smaller-magnitude Julian Date component to
//! preserve full f64 precision. Round-trip conversions (TAI → TT → TAI) come back
//! within 1 ulp of jd2.
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::conversions::{ToTT, ToTAI};
//! use celestial_core::constants::J2000_JD;
//!
//! let tai = TAI::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tt = tai.to_tt().unwrap();
//!
//! // TT is 32.184 seconds ahead of TAI
//! let tai_jd = tai.to_julian_date();
//! let tt_jd = tt.to_julian_date();
//! let diff_days = (tt_jd.jd1() - tai_jd.jd1()) + (tt_jd.jd2() - tai_jd.jd2());
//! let diff_seconds = diff_days * 86400.0;
//! assert_eq!(diff_seconds, 32.184);
//!
//! // Round-trip is exact
//! let back_to_tai = tt.to_tai().unwrap();
//! assert_eq!(tai.to_julian_date().jd1(), back_to_tai.to_julian_date().jd1());
//! assert_eq!(tai.to_julian_date().jd2(), back_to_tai.to_julian_date().jd2());
//! ```

use super::{ToTAI, ToTCG, ToTT};
use crate::constants::TT_TAI_OFFSET;
use crate::julian::finite_jd;
use crate::scales::tai::TAI;
use crate::scales::tcg::TCG;
use crate::scales::tt::TT;
use crate::TimeResult;
use celestial_core::constants::SECONDS_PER_DAY_F64;

/// Convert TAI to TT by adding the fixed 32.184 second offset.
///
/// The offset is added to whichever Julian Date component has smaller magnitude
/// to preserve maximum precision in the two-part representation.
impl ToTT for TAI {
    fn to_tt(&self) -> TimeResult<TT> {
        let tai_jd = finite_jd(self.to_julian_date())?;
        let dtat = TT_TAI_OFFSET / SECONDS_PER_DAY_F64;

        let tt_jd = tai_jd.map_smaller_part(|_, small| small + dtat);

        Ok(TT::from_julian_date(tt_jd))
    }
}

/// Convert TT to TAI by subtracting the fixed 32.184 second offset.
///
/// The offset is subtracted from whichever Julian Date component has smaller magnitude
/// to preserve maximum precision in the two-part representation.
impl ToTAI for TT {
    fn to_tai(&self) -> TimeResult<TAI> {
        let tt_jd = finite_jd(self.to_julian_date())?;
        let dtat = TT_TAI_OFFSET / SECONDS_PER_DAY_F64;

        let tai_jd = tt_jd.map_smaller_part(|_, small| small - dtat);

        Ok(TAI::from_julian_date(tai_jd))
    }
}

/// Convert TAI to TCG by chaining through TT.
///
/// TAI has no direct conversion to TCG. This chains: TAI → TT → TCG.
impl ToTCG for TAI {
    fn to_tcg(&self) -> TimeResult<TCG> {
        self.to_tt()?.to_tcg()
    }
}

/// Convert TCG to TAI by chaining through TT.
///
/// TCG has no direct conversion to TAI. This chains: TCG → TT → TAI.
impl ToTAI for TCG {
    fn to_tai(&self) -> TimeResult<TAI> {
        self.to_tt()?.to_tai()
    }
}

#[cfg(test)]
mod tests;
