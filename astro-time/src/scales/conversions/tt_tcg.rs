//! Conversions between Terrestrial Time (TT) and Geocentric Coordinate Time (TCG).
//!
//! TT and TCG differ by a constant rate defined by the IAU. TCG runs faster than TT
//! because TT accounts for gravitational time dilation at Earth's geoid, while TCG
//! is the proper time for a clock at the geocenter (in the absence of Earth's mass).
//!
//! # The L_G Rate Factor
//!
//! The defining relationship is:
//!
//! ```text
//! TCG - TT = L_G * (JD_TT - T0) * 86400
//! ```
//!
//! Where:
//! - `L_G = 6.969290134e-10` (IAU 2000 Resolution B1.9, exact by definition)
//! - `T0 = 1977 January 1, 0h TAI` (reference epoch where TCG = TT)
//! - The factor 86400 converts days to seconds
//!
//! This means TCG gains about 22 milliseconds per year relative to TT.
//!
//! # Reference Epoch
//!
//! At the reference epoch T0 (MJD 43144.0003725 in TT), TCG and TT are equal.
//! Before T0, TCG is behind TT; after T0, TCG is ahead.
//!
//! # Precision
//!
//! Round-trip conversions (TT -> TCG -> TT or TCG -> TT -> TCG) achieve sub-picosecond
//! accuracy for dates within a few centuries of J2000.0. The implementation applies
//! corrections to the smaller-magnitude Julian Date component to preserve precision.
//!
//! # Usage
//!
//! ```
//! use astro_time::scales::{TT, TCG};
//! use astro_time::scales::conversions::{ToTT, ToTCG};
//! use astro_time::julian::JulianDate;
//! use astro_core::constants::J2000_JD;
//!
//! let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tcg = tt.to_tcg().unwrap();
//!
//! // At J2000.0, TCG is about 0.506 seconds ahead of TT
//! let offset_days = tcg.to_julian_date().jd2() - tt.to_julian_date().jd2();
//! ```

use super::{ToTCG, ToTT};
use crate::constants::{TCG_RATE_LG, TCG_RATE_RATIO, TCG_REFERENCE_EPOCH};
use crate::julian::JulianDate;
use crate::scales::{TCG, TT};
use crate::TimeResult;
use astro_core::constants::MJD_ZERO_POINT;

impl ToTCG for TCG {
    /// Identity conversion. Returns self unchanged.
    fn to_tcg(&self) -> TimeResult<TCG> {
        Ok(*self)
    }
}

impl ToTT for TCG {
    /// Convert TCG to TT by removing the L_G rate correction.
    ///
    /// Computes: `TT = TCG - L_G * (JD_TCG - T0) * 86400 / 86400`
    ///
    /// The correction is subtracted because TCG runs faster than TT.
    /// At J2000.0, this removes about 0.506 seconds.
    fn to_tt(&self) -> TimeResult<TT> {
        let tcg_jd = self.to_julian_date();

        let (tt_jd1, tt_jd2) = if tcg_jd.jd1().abs() > tcg_jd.jd2().abs() {
            let correction = ((tcg_jd.jd1() - MJD_ZERO_POINT)
                + (tcg_jd.jd2() - TCG_REFERENCE_EPOCH))
                * TCG_RATE_LG;
            (tcg_jd.jd1(), tcg_jd.jd2() - correction)
        } else {
            let correction = ((tcg_jd.jd2() - MJD_ZERO_POINT)
                + (tcg_jd.jd1() - TCG_REFERENCE_EPOCH))
                * TCG_RATE_LG;
            (tcg_jd.jd1() - correction, tcg_jd.jd2())
        };

        let tt_jd = JulianDate::new(tt_jd1, tt_jd2);
        Ok(TT::from_julian_date(tt_jd))
    }
}

impl ToTCG for TT {
    /// Convert TT to TCG by applying the L_G rate correction.
    ///
    /// Uses `L_G / (1 - L_G)` as the rate ratio for the forward transformation.
    /// This ratio accounts for the fact that we're computing TCG from TT, not vice versa.
    ///
    /// At J2000.0, this adds about 0.506 seconds.
    fn to_tcg(&self) -> TimeResult<TCG> {
        let tt_jd = self.to_julian_date();

        let (tcg_jd1, tcg_jd2) = if tt_jd.jd1().abs() > tt_jd.jd2().abs() {
            let correction = ((tt_jd.jd1() - MJD_ZERO_POINT) + (tt_jd.jd2() - TCG_REFERENCE_EPOCH))
                * TCG_RATE_RATIO;
            (tt_jd.jd1(), tt_jd.jd2() + correction)
        } else {
            let correction = ((tt_jd.jd2() - MJD_ZERO_POINT) + (tt_jd.jd1() - TCG_REFERENCE_EPOCH))
                * TCG_RATE_RATIO;
            (tt_jd.jd1() + correction, tt_jd.jd2())
        };

        let tcg_jd = JulianDate::new(tcg_jd1, tcg_jd2);
        Ok(TCG::from_julian_date(tcg_jd))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_core::constants::{J2000_JD, MJD_ZERO_POINT, SECONDS_PER_DAY_F64};

    #[test]
    fn test_tcg_identity_conversion() {
        let tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let identity_tcg = tcg.to_tcg().unwrap();

        assert_eq!(
            tcg.to_julian_date().jd1(),
            identity_tcg.to_julian_date().jd1(),
            "TCG identity conversion should preserve JD1"
        );
        assert_eq!(
            tcg.to_julian_date().jd2(),
            identity_tcg.to_julian_date().jd2(),
            "TCG identity conversion should preserve JD2"
        );
    }

    #[test]
    fn test_tt_tcg_offset() {
        let test_cases = [
            (J2000_JD, 0.5058332857, "J2000.0"),
            (2455197.5, 0.7257673560, "2010-01-01"),
            (2458849.5, 0.9456713190, "2020-01-01"),
            (2469807.5, 1.6055036373, "2050-01-01"),
        ];

        let tolerance_seconds = 1e-6;

        for (jd, expected_offset_seconds, description) in test_cases {
            let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
            let tcg = tt.to_tcg().unwrap();

            let tt_jd = tt.to_julian_date();
            let tcg_jd = tcg.to_julian_date();

            let offset_days = (tcg_jd.jd1() - tt_jd.jd1()) + (tcg_jd.jd2() - tt_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            let diff = (offset_seconds - expected_offset_seconds).abs();
            assert!(
                diff < tolerance_seconds,
                "{}: TT->TCG offset should be {:.10}s, got {:.10}s (diff: {:.2e}s)",
                description,
                expected_offset_seconds,
                offset_seconds,
                diff
            );

            let tcg = TCG::from_julian_date(JulianDate::new(jd, 0.0));
            let tt = tcg.to_tt().unwrap();

            let tcg_jd = tcg.to_julian_date();
            let tt_jd = tt.to_julian_date();

            let offset_days = (tcg_jd.jd1() - tt_jd.jd1()) + (tcg_jd.jd2() - tt_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            let diff = (offset_seconds - expected_offset_seconds).abs();
            assert!(
                diff < tolerance_seconds,
                "{}: TCG->TT means TCG is {:.10}s ahead, got {:.10}s (diff: {:.2e}s)",
                description,
                expected_offset_seconds,
                offset_seconds,
                diff
            );
        }
    }

    #[test]
    fn test_tt_tcg_at_reference_epoch() {
        let reference_epoch_jd = MJD_ZERO_POINT + TCG_REFERENCE_EPOCH;

        let tt = TT::from_julian_date(JulianDate::new(reference_epoch_jd, 0.0));
        let tcg = tt.to_tcg().unwrap();

        let tt_jd = tt.to_julian_date();
        let tcg_jd = tcg.to_julian_date();

        let offset_days = (tcg_jd.jd1() - tt_jd.jd1()) + (tcg_jd.jd2() - tt_jd.jd2());
        let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

        let tolerance_seconds = 1e-12;
        assert!(
            offset_seconds.abs() < tolerance_seconds,
            "At reference epoch T0, TCG-TT should be 0, got {:.2e}s",
            offset_seconds
        );
    }

    #[test]
    fn test_tt_tcg_round_trip_precision() {
        // TCG conversions involve multiplicative scaling (LG rate). Precision loss
        // varies by jd2 magnitude: ~220 attoseconds for jd2~0, up to ~2 picoseconds
        // for jd2~+/-0.25 due to f64 magnitude mismatch when adding small corrections.
        const TOLERANCE_DAYS: f64 = 1e-14; // ~1 picosecond

        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];

        for jd2 in test_jd2_values {
            let original_tt = TT::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tcg = original_tt.to_tcg().unwrap();
            let round_trip_tt = tcg.to_tt().unwrap();

            assert_eq!(
                original_tt.to_julian_date().jd1(),
                round_trip_tt.to_julian_date().jd1(),
                "TT->TCG->TT JD1 must be exact for jd2={}",
                jd2
            );
            let jd2_diff =
                (original_tt.to_julian_date().jd2() - round_trip_tt.to_julian_date().jd2()).abs();
            assert!(
                jd2_diff <= TOLERANCE_DAYS,
                "TT->TCG->TT JD2 difference {} exceeds tolerance {} for jd2={}",
                jd2_diff,
                TOLERANCE_DAYS,
                jd2
            );

            let original_tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tt = original_tcg.to_tt().unwrap();
            let round_trip_tcg = tt.to_tcg().unwrap();

            assert_eq!(
                original_tcg.to_julian_date().jd1(),
                round_trip_tcg.to_julian_date().jd1(),
                "TCG->TT->TCG JD1 must be exact for jd2={}",
                jd2
            );
            let jd2_diff =
                (original_tcg.to_julian_date().jd2() - round_trip_tcg.to_julian_date().jd2()).abs();
            assert!(
                jd2_diff <= TOLERANCE_DAYS,
                "TCG->TT->TCG JD2 difference {} exceeds tolerance {} for jd2={}",
                jd2_diff,
                TOLERANCE_DAYS,
                jd2
            );
        }

        let alt_tt = TT::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tcg = alt_tt.to_tcg().unwrap();
        let alt_round_trip = alt_tcg.to_tt().unwrap();

        assert_eq!(
            alt_tt.to_julian_date().jd1(),
            alt_round_trip.to_julian_date().jd1(),
            "Alternate split TT->TCG->TT JD1 must be exact"
        );
        let jd2_diff =
            (alt_tt.to_julian_date().jd2() - alt_round_trip.to_julian_date().jd2()).abs();
        assert!(
            jd2_diff <= TOLERANCE_DAYS,
            "Alternate split TT->TCG->TT JD2 difference {} exceeds tolerance {}",
            jd2_diff,
            TOLERANCE_DAYS
        );
    }
}
