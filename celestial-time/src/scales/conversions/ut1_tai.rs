//! Conversions between UT1, TAI, and TT time scales.
//!
//! UT1 (Universal Time 1) is tied to Earth's actual rotation. Unlike atomic time scales
//! (TAI, TT), UT1 drifts unpredictably as Earth's rotation varies due to tidal friction,
//! core-mantle coupling, and atmospheric effects.
//!
//! # Why External Offsets Are Required
//!
//! The relationship between UT1 and atomic scales cannot be computed from first principles.
//! It must be measured by IERS (International Earth Rotation Service) and published as:
//!
//! - **UT1-TAI**:          Direct offset, typically around -37 seconds (as of 2024)
//! - **Delta-T (TT-UT1)**: Historical parameter, ~69 seconds at J2000.0
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
//! use celestial_time::scales::{TAI, TT, UT1};
//! use celestial_time::scales::conversions::{ToUT1WithOffset, ToTAIWithOffset};
//! use celestial_time::scales::conversions::{ToUT1WithDeltaT, ToTTWithDeltaT};
//! use celestial_time::julian::JulianDate;
//!
//! // UT1-TAI offset from IERS Bulletin A (example: -37.0 seconds)
//! let ut1_tai_offset = -37.0;
//!
//! let ut1 = UT1::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let tai = ut1.to_tai_with_offset(ut1_tai_offset).unwrap();
//! let back = tai.to_ut1_with_offset(ut1_tai_offset).unwrap();
//!
//! // Delta-T from historical tables or prediction models
//! let delta_t = 69.0;  // seconds at J2000.0
//!
//! let tt = ut1.to_tt_with_delta_t(delta_t).unwrap();
//! let back = tt.to_ut1_with_delta_t(delta_t).unwrap();
//! ```
//!
//! # Precision Notes
//!
//! Offsets are applied to the smaller-magnitude Julian Date component to preserve
//! precision. Round-trip conversions maintain sub-nanosecond accuracy.

use super::ToUT1;
use crate::julian::JulianDate;
use crate::scales::{TAI, TT, UT1};
use crate::TimeResult;
use celestial_core::constants::SECONDS_PER_DAY_F64;

impl ToUT1 for UT1 {
    fn to_ut1(&self) -> TimeResult<UT1> {
        Ok(*self)
    }
}

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
        let ut1_jd = self.to_julian_date();
        let offset_days = ut1_tai_offset_seconds / SECONDS_PER_DAY_F64;

        // TAI = UT1 - (UT1-TAI), so subtract the offset.
        // Apply to smaller-magnitude component for precision.
        let (tai_jd1, tai_jd2) = if ut1_jd.jd1().abs() > ut1_jd.jd2().abs() {
            (ut1_jd.jd1(), ut1_jd.jd2() - offset_days)
        } else {
            (ut1_jd.jd1() - offset_days, ut1_jd.jd2())
        };

        Ok(TAI::from_julian_date(JulianDate::new(tai_jd1, tai_jd2)))
    }
}

impl ToUT1WithOffset for TAI {
    fn to_ut1_with_offset(&self, ut1_tai_offset_seconds: f64) -> TimeResult<UT1> {
        let tai_jd = self.to_julian_date();
        let offset_days = ut1_tai_offset_seconds / SECONDS_PER_DAY_F64;

        // UT1 = TAI + (UT1-TAI), so add the offset.
        // Apply to smaller-magnitude component for precision.
        let (ut1_jd1, ut1_jd2) = if tai_jd.jd1().abs() > tai_jd.jd2().abs() {
            (tai_jd.jd1(), tai_jd.jd2() + offset_days)
        } else {
            (tai_jd.jd1() + offset_days, tai_jd.jd2())
        };

        Ok(UT1::from_julian_date(JulianDate::new(ut1_jd1, ut1_jd2)))
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
        let ut1_jd = self.to_julian_date();
        let delta_t_days = delta_t_seconds / SECONDS_PER_DAY_F64;

        // TT = UT1 + Delta-T, so add.
        // Apply to smaller-magnitude component for precision.
        let (tt_jd1, tt_jd2) = if ut1_jd.jd1().abs() > ut1_jd.jd2().abs() {
            (ut1_jd.jd1(), ut1_jd.jd2() + delta_t_days)
        } else {
            (ut1_jd.jd1() + delta_t_days, ut1_jd.jd2())
        };

        Ok(TT::from_julian_date(JulianDate::new(tt_jd1, tt_jd2)))
    }
}

impl ToUT1WithDeltaT for TT {
    fn to_ut1_with_delta_t(&self, delta_t_seconds: f64) -> TimeResult<UT1> {
        let tt_jd = self.to_julian_date();
        let delta_t_days = delta_t_seconds / SECONDS_PER_DAY_F64;

        // UT1 = TT - Delta-T, so subtract.
        // Apply to smaller-magnitude component for precision.
        let (ut1_jd1, ut1_jd2) = if tt_jd.jd1().abs() > tt_jd.jd2().abs() {
            (tt_jd.jd1(), tt_jd.jd2() - delta_t_days)
        } else {
            (tt_jd.jd1() - delta_t_days, tt_jd.jd2())
        };

        Ok(UT1::from_julian_date(JulianDate::new(ut1_jd1, ut1_jd2)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_ut1_identity_conversion() {
        let ut1 = UT1::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let identity_ut1 = ut1.to_ut1().unwrap();

        assert_eq!(
            ut1.to_julian_date().jd1(),
            identity_ut1.to_julian_date().jd1(),
            "UT1 identity conversion should preserve JD1 exactly"
        );
        assert_eq!(
            ut1.to_julian_date().jd2(),
            identity_ut1.to_julian_date().jd2(),
            "UT1 identity conversion should preserve JD2 exactly"
        );
    }

    #[test]
    fn test_ut1_tai_offset_applied_correctly() {
        let test_dates = [
            (J2000_JD, "J2000.0"),
            (2455197.5, "2010-01-01"),
            (2459580.5, "2022-01-01"),
        ];
        let ut1_tai_offset = -32.3;
        let delta_t = 69.0;

        for (jd, description) in test_dates {
            // UT1 -> TAI: TAI = UT1 - (UT1-TAI), so TAI should be ahead by 32.3s
            let ut1 = UT1::from_julian_date(JulianDate::new(jd, 0.0));
            let tai = ut1.to_tai_with_offset(ut1_tai_offset).unwrap();

            let ut1_jd = ut1.to_julian_date();
            let tai_jd = tai.to_julian_date();

            let offset_days = (tai_jd.jd1() - ut1_jd.jd1()) + (tai_jd.jd2() - ut1_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, -ut1_tai_offset,
                "{}: UT1->TAI offset must be exactly {} seconds",
                description, -ut1_tai_offset
            );

            // TAI -> UT1: UT1 = TAI + (UT1-TAI), so UT1 should be behind by 32.3s
            let tai = TAI::from_julian_date(JulianDate::new(jd, 0.0));
            let ut1 = tai.to_ut1_with_offset(ut1_tai_offset).unwrap();

            let tai_jd = tai.to_julian_date();
            let ut1_jd = ut1.to_julian_date();

            let offset_days = (tai_jd.jd1() - ut1_jd.jd1()) + (tai_jd.jd2() - ut1_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, -ut1_tai_offset,
                "{}: TAI->UT1 means TAI is {} seconds ahead",
                description, -ut1_tai_offset
            );

            // UT1 -> TT: TT = UT1 + Delta-T, so TT should be ahead by 69s
            let ut1 = UT1::from_julian_date(JulianDate::new(jd, 0.0));
            let tt = ut1.to_tt_with_delta_t(delta_t).unwrap();

            let ut1_jd = ut1.to_julian_date();
            let tt_jd = tt.to_julian_date();

            let offset_days = (tt_jd.jd1() - ut1_jd.jd1()) + (tt_jd.jd2() - ut1_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, delta_t,
                "{}: UT1->TT offset must be exactly {} seconds",
                description, delta_t
            );

            // TT -> UT1: UT1 = TT - Delta-T, so UT1 should be behind by 69s
            let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
            let ut1 = tt.to_ut1_with_delta_t(delta_t).unwrap();

            let tt_jd = tt.to_julian_date();
            let ut1_jd = ut1.to_julian_date();

            let offset_days = (tt_jd.jd1() - ut1_jd.jd1()) + (tt_jd.jd2() - ut1_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, delta_t,
                "{}: TT->UT1 means TT is {} seconds ahead",
                description, delta_t
            );
        }
    }

    #[test]
    fn test_ut1_tai_round_trip_precision() {
        // Division by SECONDS_PER_DAY introduces ~5 picosecond rounding.
        // 1e-14 days = ~1 picosecond tolerance.
        const TOLERANCE_DAYS: f64 = 1e-14;

        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];
        let test_offsets = [-32.0, -31.8, -32.5, -33.1, -30.9];

        for jd2 in test_jd2_values {
            for &offset in &test_offsets {
                // UT1 -> TAI -> UT1
                let original_ut1 = UT1::from_julian_date(JulianDate::new(J2000_JD, jd2));
                let tai = original_ut1.to_tai_with_offset(offset).unwrap();
                let round_trip_ut1 = tai.to_ut1_with_offset(offset).unwrap();

                assert_eq!(
                    original_ut1.to_julian_date().jd1(),
                    round_trip_ut1.to_julian_date().jd1(),
                    "UT1->TAI->UT1 JD1 must be exact for jd2={}, offset={}",
                    jd2,
                    offset
                );
                let jd2_diff = (original_ut1.to_julian_date().jd2()
                    - round_trip_ut1.to_julian_date().jd2())
                .abs();
                assert!(
                    jd2_diff <= TOLERANCE_DAYS,
                    "UT1->TAI->UT1 JD2 diff {} exceeds tolerance {} for jd2={}, offset={}",
                    jd2_diff,
                    TOLERANCE_DAYS,
                    jd2,
                    offset
                );

                // TAI -> UT1 -> TAI
                let original_tai = TAI::from_julian_date(JulianDate::new(J2000_JD, jd2));
                let ut1 = original_tai.to_ut1_with_offset(offset).unwrap();
                let round_trip_tai = ut1.to_tai_with_offset(offset).unwrap();

                assert_eq!(
                    original_tai.to_julian_date().jd1(),
                    round_trip_tai.to_julian_date().jd1(),
                    "TAI->UT1->TAI JD1 must be exact for jd2={}, offset={}",
                    jd2,
                    offset
                );
                let jd2_diff = (original_tai.to_julian_date().jd2()
                    - round_trip_tai.to_julian_date().jd2())
                .abs();
                assert!(
                    jd2_diff <= TOLERANCE_DAYS,
                    "TAI->UT1->TAI JD2 diff {} exceeds tolerance {} for jd2={}, offset={}",
                    jd2_diff,
                    TOLERANCE_DAYS,
                    jd2,
                    offset
                );
            }
        }

        // Alternate JD split case (jd2 > jd1)
        let alt_ut1 = UT1::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tai = alt_ut1.to_tai_with_offset(-32.0).unwrap();
        let alt_round_trip = alt_tai.to_ut1_with_offset(-32.0).unwrap();

        assert_eq!(
            alt_ut1.to_julian_date().jd1(),
            alt_round_trip.to_julian_date().jd1(),
            "Alternate split UT1->TAI->UT1 JD1 must be exact"
        );
        let jd2_diff =
            (alt_ut1.to_julian_date().jd2() - alt_round_trip.to_julian_date().jd2()).abs();
        assert!(
            jd2_diff <= TOLERANCE_DAYS,
            "Alternate split UT1->TAI->UT1 JD2 diff {} exceeds tolerance {}",
            jd2_diff,
            TOLERANCE_DAYS
        );
    }

    #[test]
    fn test_ut1_tt_round_trip_precision() {
        // Division by SECONDS_PER_DAY introduces ~5 picosecond rounding.
        // 1e-14 days = ~1 picosecond tolerance.
        const TOLERANCE_DAYS: f64 = 1e-14;

        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];
        let test_delta_t_values = [63.8, 69.0, 70.5, 65.2];

        for jd2 in test_jd2_values {
            for &delta_t in &test_delta_t_values {
                // UT1 -> TT -> UT1
                let original_ut1 = UT1::from_julian_date(JulianDate::new(J2000_JD, jd2));
                let tt = original_ut1.to_tt_with_delta_t(delta_t).unwrap();
                let round_trip_ut1 = tt.to_ut1_with_delta_t(delta_t).unwrap();

                assert_eq!(
                    original_ut1.to_julian_date().jd1(),
                    round_trip_ut1.to_julian_date().jd1(),
                    "UT1->TT->UT1 JD1 must be exact for jd2={}, delta_t={}",
                    jd2,
                    delta_t
                );
                let jd2_diff = (original_ut1.to_julian_date().jd2()
                    - round_trip_ut1.to_julian_date().jd2())
                .abs();
                assert!(
                    jd2_diff <= TOLERANCE_DAYS,
                    "UT1->TT->UT1 JD2 diff {} exceeds tolerance {} for jd2={}, delta_t={}",
                    jd2_diff,
                    TOLERANCE_DAYS,
                    jd2,
                    delta_t
                );

                // TT -> UT1 -> TT
                let original_tt = TT::from_julian_date(JulianDate::new(J2000_JD, jd2));
                let ut1 = original_tt.to_ut1_with_delta_t(delta_t).unwrap();
                let round_trip_tt = ut1.to_tt_with_delta_t(delta_t).unwrap();

                assert_eq!(
                    original_tt.to_julian_date().jd1(),
                    round_trip_tt.to_julian_date().jd1(),
                    "TT->UT1->TT JD1 must be exact for jd2={}, delta_t={}",
                    jd2,
                    delta_t
                );
                let jd2_diff = (original_tt.to_julian_date().jd2()
                    - round_trip_tt.to_julian_date().jd2())
                .abs();
                assert!(
                    jd2_diff <= TOLERANCE_DAYS,
                    "TT->UT1->TT JD2 diff {} exceeds tolerance {} for jd2={}, delta_t={}",
                    jd2_diff,
                    TOLERANCE_DAYS,
                    jd2,
                    delta_t
                );
            }
        }

        // Alternate JD split case (jd2 > jd1)
        let alt_ut1 = UT1::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tt = alt_ut1.to_tt_with_delta_t(69.0).unwrap();
        let alt_round_trip = alt_tt.to_ut1_with_delta_t(69.0).unwrap();

        assert_eq!(
            alt_ut1.to_julian_date().jd1(),
            alt_round_trip.to_julian_date().jd1(),
            "Alternate split UT1->TT->UT1 JD1 must be exact"
        );
        let jd2_diff =
            (alt_ut1.to_julian_date().jd2() - alt_round_trip.to_julian_date().jd2()).abs();
        assert!(
            jd2_diff <= TOLERANCE_DAYS,
            "Alternate split UT1->TT->UT1 JD2 diff {} exceeds tolerance {}",
            jd2_diff,
            TOLERANCE_DAYS
        );
    }
}
