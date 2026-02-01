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
//! preserve full f64 precision. Round-trip conversions (TAI → TT → TAI) are exact
//! to the bit level.
//!
//! # Usage
//!
//! ```
//! use astro_time::{JulianDate, TAI, TT};
//! use astro_time::scales::conversions::{ToTT, ToTAI};
//! use astro_core::constants::J2000_JD;
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
use crate::scales::{TAI, TCG, TT};
use crate::TimeResult;
use astro_core::constants::SECONDS_PER_DAY_F64;

/// Identity conversion for TAI.
impl ToTAI for TAI {
    fn to_tai(&self) -> TimeResult<TAI> {
        Ok(*self)
    }
}

/// Identity conversion for TT.
impl ToTT for TT {
    fn to_tt(&self) -> TimeResult<TT> {
        Ok(*self)
    }
}

/// Convert TAI to TT by adding the fixed 32.184 second offset.
///
/// The offset is added to whichever Julian Date component has smaller magnitude
/// to preserve maximum precision in the two-part representation.
impl ToTT for TAI {
    fn to_tt(&self) -> TimeResult<TT> {
        let tai_jd = self.to_julian_date();
        let dtat = TT_TAI_OFFSET / SECONDS_PER_DAY_F64;

        let jd1_raw = tai_jd.jd1().to_bits();
        let jd2_raw = tai_jd.jd2().to_bits();
        let jd1_magnitude = jd1_raw & 0x7FFFFFFFFFFFFFFF;
        let jd2_magnitude = jd2_raw & 0x7FFFFFFFFFFFFFFF;
        let (tt_jd1, tt_jd2) = if jd1_magnitude > jd2_magnitude {
            (tai_jd.jd1(), tai_jd.jd2() + dtat)
        } else {
            (tai_jd.jd1() + dtat, tai_jd.jd2())
        };

        Ok(TT::from_julian_date_raw(tt_jd1, tt_jd2))
    }
}

/// Convert TT to TAI by subtracting the fixed 32.184 second offset.
///
/// The offset is subtracted from whichever Julian Date component has smaller magnitude
/// to preserve maximum precision in the two-part representation.
impl ToTAI for TT {
    fn to_tai(&self) -> TimeResult<TAI> {
        let tt_jd = self.to_julian_date();
        let dtat = TT_TAI_OFFSET / SECONDS_PER_DAY_F64;

        let jd1_raw = tt_jd.jd1().to_bits();
        let jd2_raw = tt_jd.jd2().to_bits();
        let jd1_magnitude = jd1_raw & 0x7FFFFFFFFFFFFFFF;
        let jd2_magnitude = jd2_raw & 0x7FFFFFFFFFFFFFFF;
        let (tai_jd1, tai_jd2) = if jd1_magnitude > jd2_magnitude {
            (tt_jd.jd1(), tt_jd.jd2() - dtat)
        } else {
            (tt_jd.jd1() - dtat, tt_jd.jd2())
        };

        Ok(TAI::from_julian_date_raw(tai_jd1, tai_jd2))
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
mod tests {
    use super::*;
    use crate::JulianDate;
    use astro_core::constants::J2000_JD;

    #[test]
    fn test_identity_conversions() {
        let tai = TAI::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let identity_tai = tai.to_tai().unwrap();

        assert_eq!(
            tai.to_julian_date().jd1(),
            identity_tai.to_julian_date().jd1(),
            "TAI identity conversion should preserve JD1"
        );
        assert_eq!(
            tai.to_julian_date().jd2(),
            identity_tai.to_julian_date().jd2(),
            "TAI identity conversion should preserve JD2"
        );

        let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let identity_tt = tt.to_tt().unwrap();

        assert_eq!(
            tt.to_julian_date().jd1(),
            identity_tt.to_julian_date().jd1(),
            "TT identity conversion should preserve JD1"
        );
        assert_eq!(
            tt.to_julian_date().jd2(),
            identity_tt.to_julian_date().jd2(),
            "TT identity conversion should preserve JD2"
        );
    }

    #[test]
    fn test_tai_tt_offset_32_184_seconds() {
        let test_dates = [
            (J2000_JD, "J2000.0"),
            (2455197.5, "2010-01-01"),
            (2459580.5, "2022-01-01"),
            (2440587.5, "1970-01-01 Unix epoch"),
        ];

        for (jd, description) in test_dates {
            let tai = TAI::from_julian_date(JulianDate::new(jd, 0.0));
            let tt = tai.to_tt().unwrap();

            let tai_jd = tai.to_julian_date();
            let tt_jd = tt.to_julian_date();

            let offset_days = (tt_jd.jd1() - tai_jd.jd1()) + (tt_jd.jd2() - tai_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, 32.184,
                "{}: TAI->TT offset must be exactly 32.184 seconds",
                description
            );

            let tt = TT::from_julian_date(JulianDate::new(jd, 0.0));
            let tai = tt.to_tai().unwrap();

            let tt_jd = tt.to_julian_date();
            let tai_jd = tai.to_julian_date();

            let offset_days = (tt_jd.jd1() - tai_jd.jd1()) + (tt_jd.jd2() - tai_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, 32.184,
                "{}: TT->TAI means TT is 32.184 seconds ahead",
                description
            );
        }
    }

    #[test]
    fn test_tai_tt_round_trip_precision() {
        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];

        for jd2 in test_jd2_values {
            let original_tai = TAI::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tt = original_tai.to_tt().unwrap();
            let round_trip_tai = tt.to_tai().unwrap();

            assert_eq!(
                original_tai.to_julian_date().jd1(),
                round_trip_tai.to_julian_date().jd1(),
                "TAI->TT->TAI JD1 must be exact for jd2={}",
                jd2
            );
            assert_eq!(
                original_tai.to_julian_date().jd2(),
                round_trip_tai.to_julian_date().jd2(),
                "TAI->TT->TAI JD2 must be exact for jd2={}",
                jd2
            );

            let original_tt = TT::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tai = original_tt.to_tai().unwrap();
            let round_trip_tt = tai.to_tt().unwrap();

            assert_eq!(
                original_tt.to_julian_date().jd1(),
                round_trip_tt.to_julian_date().jd1(),
                "TT->TAI->TT JD1 must be exact for jd2={}",
                jd2
            );
            assert_eq!(
                original_tt.to_julian_date().jd2(),
                round_trip_tt.to_julian_date().jd2(),
                "TT->TAI->TT JD2 must be exact for jd2={}",
                jd2
            );
        }

        let alt_tai = TAI::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tt = alt_tai.to_tt().unwrap();
        let alt_round_trip = alt_tt.to_tai().unwrap();

        assert_eq!(
            alt_tai.to_julian_date().jd1(),
            alt_round_trip.to_julian_date().jd1(),
            "Alternate split TAI->TT->TAI JD1 must be exact"
        );
        assert_eq!(
            alt_tai.to_julian_date().jd2(),
            alt_round_trip.to_julian_date().jd2(),
            "Alternate split TAI->TT->TAI JD2 must be exact"
        );
    }

    #[test]
    fn test_tai_tcg_chain_round_trip() {
        // TCG conversions involve multiplicative scaling (LG rate). Precision loss
        // varies by jd2 magnitude: ~220 attoseconds for jd2≈0, up to ~2 picoseconds
        // for jd2≈±0.25 due to f64 magnitude mismatch when adding small corrections.
        const TOLERANCE_DAYS: f64 = 1e-14; // ~1 picosecond

        let test_jd2_values = [0.0, 0.123456789, 0.5, -0.25];

        for jd2 in test_jd2_values {
            let original_tai = TAI::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tcg = original_tai.to_tcg().unwrap();
            let round_trip_tai = tcg.to_tai().unwrap();

            assert_eq!(
                original_tai.to_julian_date().jd1(),
                round_trip_tai.to_julian_date().jd1(),
                "TAI->TCG->TAI JD1 must be exact for jd2={}",
                jd2
            );
            let jd2_diff =
                (original_tai.to_julian_date().jd2() - round_trip_tai.to_julian_date().jd2()).abs();
            assert!(
                jd2_diff <= TOLERANCE_DAYS,
                "TAI->TCG->TAI JD2 difference {} exceeds tolerance {} for jd2={}",
                jd2_diff,
                TOLERANCE_DAYS,
                jd2
            );

            let original_tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tai = original_tcg.to_tai().unwrap();
            let round_trip_tcg = tai.to_tcg().unwrap();

            assert_eq!(
                original_tcg.to_julian_date().jd1(),
                round_trip_tcg.to_julian_date().jd1(),
                "TCG->TAI->TCG JD1 must be exact for jd2={}",
                jd2
            );
            let jd2_diff =
                (original_tcg.to_julian_date().jd2() - round_trip_tcg.to_julian_date().jd2()).abs();
            assert!(
                jd2_diff <= TOLERANCE_DAYS,
                "TCG->TAI->TCG JD2 difference {} exceeds tolerance {} for jd2={}",
                jd2_diff,
                TOLERANCE_DAYS,
                jd2
            );
        }
    }
}
