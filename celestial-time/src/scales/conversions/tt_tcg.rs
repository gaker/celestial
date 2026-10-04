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
//! TCG - TT = L_G * (JD_TCG - T0) * 86400
//! ```
//!
//! Where:
//! - `L_G = 6.969290134e-10` (IAU 2000 Resolution B1.9, exact by definition)
//! - `T0 = 1977 January 1, 0h TAI` (reference epoch where TCG = TT)
//! - The factor 86400 converts days to seconds
//!
//! Measured from JD_TT instead, the rate is L_G / (1 - L_G).
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
//! Round-trip conversions (TT -> TCG -> TT or TCG -> TT -> TCG) come back within
//! 1 ulp of jd2 (at most ~10 ps) for dates within a few centuries of J2000.0. The implementation applies
//! corrections to the smaller-magnitude Julian Date component to preserve precision.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::tcg::TCG;
//! use celestial_time::scales::tt::TT;
//! use celestial_time::scales::conversions::{ToTT, ToTCG};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tcg = tt.to_tcg().unwrap();
//!
//! // At J2000.0, TCG is about 0.506 seconds ahead of TT
//! let offset_days = tcg.to_julian_date().jd2() - tt.to_julian_date().jd2();
//! ```

use super::{ToTCG, ToTT};
use crate::constants::{TCG_RATE_LG, TCG_RATE_RATIO, TCG_REFERENCE_EPOCH};
use crate::julian::finite_jd;
use crate::scales::tcg::TCG;
use crate::scales::tt::TT;
use crate::TimeResult;
use celestial_core::constants::MJD_ZERO_POINT;

impl ToTT for TCG {
    /// Convert TCG to TT by removing the L_G rate correction.
    ///
    /// Computes: `TT = TCG - L_G * (JD_TCG - T0) * 86400 / 86400`
    ///
    /// The correction is subtracted because TCG runs faster than TT.
    /// At J2000.0, this removes about 0.506 seconds.
    fn to_tt(&self) -> TimeResult<TT> {
        let tcg_jd = finite_jd(self.to_julian_date())?;

        let tt_jd = tcg_jd.map_smaller_part(|big, small| {
            small - ((big - MJD_ZERO_POINT) + (small - TCG_REFERENCE_EPOCH)) * TCG_RATE_LG
        });
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
        let tt_jd = finite_jd(self.to_julian_date())?;

        let tcg_jd = tt_jd.map_smaller_part(|big, small| {
            small + ((big - MJD_ZERO_POINT) + (small - TCG_REFERENCE_EPOCH)) * TCG_RATE_RATIO
        });
        Ok(TCG::from_julian_date(tcg_jd))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use crate::TimeError;
    use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT};

    #[test]
    fn test_tt_tcg_matches_erfa() {
        // eraTttcg and eraTcgtt.
        let cases = [
            (J2000_JD, 5.85455192154085e-6, -5.854551917460643e-6),
            (2455197.5, 8.400085144758405e-6, -8.400085138904144e-6),
            (2458849.5, 1.0945269903469018e-5, -1.0945269895840943e-5),
            (2469807.5, 1.8582218037628628e-5, -1.8582218024678145e-5),
        ];
        for (jd, tcg_jd2, tt_jd2) in cases {
            let tcg = TT::from_julian_date(JulianDate::new(jd, 0.0))
                .to_tcg()
                .unwrap();
            assert_eq!(tcg.to_julian_date().parts(), (jd, tcg_jd2));
            let tt = TCG::from_julian_date(JulianDate::new(jd, 0.0))
                .to_tt()
                .unwrap();
            assert_eq!(tt.to_julian_date().parts(), (jd, tt_jd2));
        }

        // At the 1977 reference epoch they agree, up to the rounding of the epoch.
        let epoch = MJD_ZERO_POINT + TCG_REFERENCE_EPOCH;
        let tcg = TT::from_julian_date(JulianDate::new(epoch, 0.0))
            .to_tcg()
            .unwrap();
        assert_eq!(
            tcg.to_julian_date().parts(),
            (epoch, 1.1155817123279573e-19)
        );
    }

    #[test]
    fn test_tt_tcg_round_trip_matches_erfa() {
        // eraTttcg then eraTcgtt, and the reverse. Only the (JD, 0) split doesn't
        // come back exactly, by under 3e-21 d.
        let cases = [
            (J2000_JD, 0.0, 2.541098841762901e-21, 1.6940658945086007e-21),
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
        for (jd1, jd2, tt_back, tcg_back) in cases {
            let tt = TT::from_julian_date(JulianDate::new(jd1, jd2));
            let back = tt.to_tcg().unwrap().to_tt().unwrap();
            assert_eq!(back.to_julian_date().parts(), (jd1, tt_back));

            let tcg = TCG::from_julian_date(JulianDate::new(jd1, jd2));
            let back = tcg.to_tt().unwrap().to_tcg().unwrap();
            assert_eq!(back.to_julian_date().parts(), (jd1, tcg_back));
        }

        let tt = TT::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let back = tt.to_tcg().unwrap().to_tt().unwrap();
        assert_eq!(back.to_julian_date().parts(), tt.to_julian_date().parts());
        let tcg = TCG::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let back = tcg.to_tt().unwrap().to_tcg().unwrap();
        assert_eq!(back.to_julian_date().parts(), tcg.to_julian_date().parts());
    }

    #[test]
    fn test_non_finite_julian_date_is_rejected() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let jd = JulianDate::new(bad, 0.0);
            let expected =
                TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
            assert_eq!(TCG::from_julian_date(jd).to_tt().unwrap_err(), expected);
            assert_eq!(TT::from_julian_date(jd).to_tcg().unwrap_err(), expected);
        }
    }
}
