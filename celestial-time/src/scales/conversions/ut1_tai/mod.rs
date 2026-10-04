//! Conversions between UT1, TAI, and TT time scales.
//!
//! UT1 (Universal Time 1) is tied to Earth's actual rotation. Unlike atomic time scales
//! (TAI, TT), UT1 drifts unpredictably as Earth's rotation varies due to tidal friction,
//! core-mantle coupling, and atmospheric effects.
//!
//! # Why External Offsets Are Required
//!
//! The relationship between UT1 and atomic scales cannot be computed from first principles.
//! It must be measured by the IERS (International Earth Rotation and Reference Systems
//! Service) and published as:
//!
//! - **UT1-TAI**:          Direct offset, typically around -37 seconds (as of 2024)
//! - **Delta-T (TT-UT1)**: Historical parameter, ~63.8 seconds at J2000.0
//!
//! These values change continuously. IERS Bulletin A provides predictions; Bulletin B
//! provides final values after the fact. The offset changes by roughly 1-2 ms/day.
//!
//! # Conversion Paths
//!
//! ```text
//! UT1 <-(UT1-TAI offset)-> TAI
//! UT1 <-----(Delta-T)-----> TT
//! ```
//!
//! Both require externally-supplied offset values. This module provides the traits;
//! you provide the offset from EOP (Earth Orientation Parameters) data.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::ut1::UT1;
//! use celestial_time::scales::conversions::ut1_tai::{ToTAIWithOffset, ToUT1WithOffset};
//! use celestial_time::scales::conversions::ut1_tai::{ToTTWithDeltaT, ToUT1WithDeltaT};
//! use celestial_time::julian::JulianDate;
//!
//! // UT1-UTC from IERS Bulletin A (0.3554 s) minus TAI-UTC (32 s) at J2000
//! let ut1_tai_offset = -31.645;
//!
//! let ut1 = UT1::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let tai = ut1.to_tai_with_offset(ut1_tai_offset).unwrap();
//! let back = tai.to_ut1_with_offset(ut1_tai_offset).unwrap();
//!
//! // Delta-T from historical tables or prediction models
//! let delta_t = 63.8;  // seconds at J2000.0
//!
//! let tt = ut1.to_tt_with_delta_t(delta_t).unwrap();
//! let back = tt.to_ut1_with_delta_t(delta_t).unwrap();
//! ```
//!
//! # Precision Notes
//!
//! Offsets are applied to the smaller-magnitude Julian Date component to preserve
//! precision. Round-trip conversions maintain sub-nanosecond accuracy.

use crate::julian::{finite_arg, finite_jd};
use crate::scales::tai::TAI;
use crate::scales::tt::TT;
use crate::scales::ut1::UT1;
use crate::TimeResult;
use celestial_core::constants::SECONDS_PER_DAY_F64;

/// Convert TAI to UT1 using a supplied UT1-TAI offset.
///
/// The offset comes from IERS Earth Orientation Parameters. Typical values
/// are around -37 seconds (as of 2024), becoming more negative over time
/// as leap seconds accumulate.
///
/// Note: The offset is UT1-TAI, so it's negative when UT1 is behind TAI.
pub trait ToUT1WithOffset {
    /// Convert to UT1 using the given UT1-TAI offset in seconds.
    ///
    /// The offset should be UT1-TAI (typically negative). To find UT1:
    /// `UT1 = TAI + (UT1-TAI)`
    fn to_ut1_with_offset(&self, ut1_tai_offset_seconds: f64) -> TimeResult<UT1>;
}

/// Convert UT1 to TAI using a supplied UT1-TAI offset.
///
/// The offset comes from IERS Earth Orientation Parameters. This is the
/// inverse operation of [`ToUT1WithOffset`].
pub trait ToTAIWithOffset {
    /// Convert to TAI using the given UT1-TAI offset in seconds.
    ///
    /// The offset should be UT1-TAI (typically negative). To find TAI:
    /// `TAI = UT1 - (UT1-TAI)`
    fn to_tai_with_offset(&self, ut1_tai_offset_seconds: f64) -> TimeResult<TAI>;
}

impl ToTAIWithOffset for UT1 {
    fn to_tai_with_offset(&self, ut1_tai_offset_seconds: f64) -> TimeResult<TAI> {
        let ut1_jd = finite_jd(self.to_julian_date())?;
        let offset_days =
            finite_arg("ut1_tai_offset_seconds", ut1_tai_offset_seconds)? / SECONDS_PER_DAY_F64;

        let tai_jd = ut1_jd.map_smaller_part(|_, small| small - offset_days);

        Ok(TAI::from_julian_date(tai_jd))
    }
}

impl ToUT1WithOffset for TAI {
    fn to_ut1_with_offset(&self, ut1_tai_offset_seconds: f64) -> TimeResult<UT1> {
        let tai_jd = finite_jd(self.to_julian_date())?;
        let offset_days =
            finite_arg("ut1_tai_offset_seconds", ut1_tai_offset_seconds)? / SECONDS_PER_DAY_F64;

        let ut1_jd = tai_jd.map_smaller_part(|_, small| small + offset_days);

        Ok(UT1::from_julian_date(ut1_jd))
    }
}

/// Convert UT1 to TT using Delta-T.
///
/// Delta-T is defined as TT - UT1. Unlike the fixed TAI-TT offset (32.184s),
/// Delta-T varies with Earth's rotation:
///
/// - At J2000.0: ~63.8 seconds
/// - In 2024: ~69 seconds
/// - Historical values go back centuries (reconstructed from eclipse records)
///
/// Delta-T combines two effects:
/// - The fixed TT-TAI offset (32.184s)
/// - The variable TAI-UT1 difference (leap seconds + sub-second drift)
///
/// Use this for direct UT1 <-> TT conversion when you have Delta-T from
/// historical tables or prediction models. For modern dates with EOP data,
/// chaining through TAI may be more accurate.
pub trait ToTTWithDeltaT {
    /// Convert to TT using the given Delta-T in seconds.
    ///
    /// Delta-T = TT - UT1, so: `TT = UT1 + Delta-T`
    fn to_tt_with_delta_t(&self, delta_t_seconds: f64) -> TimeResult<TT>;
}

/// Convert TT to UT1 using Delta-T.
///
/// This is the inverse of [`ToTTWithDeltaT`]. See that trait for Delta-T details.
pub trait ToUT1WithDeltaT {
    /// Convert to UT1 using the given Delta-T in seconds.
    ///
    /// Delta-T = TT - UT1, so: `UT1 = TT - Delta-T`
    fn to_ut1_with_delta_t(&self, delta_t_seconds: f64) -> TimeResult<UT1>;
}

impl ToTTWithDeltaT for UT1 {
    fn to_tt_with_delta_t(&self, delta_t_seconds: f64) -> TimeResult<TT> {
        let ut1_jd = finite_jd(self.to_julian_date())?;
        let delta_t_days = finite_arg("delta_t_seconds", delta_t_seconds)? / SECONDS_PER_DAY_F64;

        let tt_jd = ut1_jd.map_smaller_part(|_, small| small + delta_t_days);

        Ok(TT::from_julian_date(tt_jd))
    }
}

impl ToUT1WithDeltaT for TT {
    fn to_ut1_with_delta_t(&self, delta_t_seconds: f64) -> TimeResult<UT1> {
        let tt_jd = finite_jd(self.to_julian_date())?;
        let delta_t_days = finite_arg("delta_t_seconds", delta_t_seconds)? / SECONDS_PER_DAY_F64;

        let ut1_jd = tt_jd.map_smaller_part(|_, small| small - delta_t_days);

        Ok(UT1::from_julian_date(ut1_jd))
    }
}

#[cfg(test)]
mod tests;
