//! Conversions between Geocentric Coordinate Time (TCG) and Barycentric Coordinate Time (TCB).
//!
//! TCG and TCB are coordinate time scales for the geocentric and barycentric reference frames,
//! respectively. TCB runs faster than TCG because the solar system's gravitational potential
//! at Earth's orbit causes additional time dilation beyond Earth's own gravitational field.
//!
//! # The L_B Rate Factor
//!
//! The defining relationship from IAU 2006 Resolution B3 is:
//!
//! ```text
//! TCB - TCG = L_B * (JD_TCG - T0) * 86400
//! ```
//!
//! Where:
//! - `L_B = 1.550519768e-8` (IAU 2006, exact by definition)
//! - `T0 = 1977 January 1, 0h TAI` (MJD 43144.0, the common reference epoch)
//! - The factor 86400 converts days to seconds
//!
//! L_B represents the average fractional rate difference between TCB and TCG due to the
//! Sun's gravitational potential at Earth's orbit. TCB gains about 489 milliseconds per
//! year relative to TCG.
//!
//! # Reference Epoch
//!
//! At the reference epoch T0 (1977 January 1, 0h TAI, JD 2443144.5003725), TCG and TCB
//! are defined to be equal. This is the same epoch used for the TT-TCG relationship.
//!
//! Before T0, TCB is behind TCG; after T0, TCB is ahead.
//!
//! # Physical Interpretation
//!
//! The L_B rate difference arises from general relativity:
//!
//! - **TCG**: Proper time for a clock at the geocenter (removing Earth's gravitational
//!   potential but still in the Sun's potential well)
//! - **TCB**: Proper time for a clock at the solar system barycenter (outside all
//!   gravitational potentials of the solar system)
//!
//! Since Earth orbits within the Sun's gravitational potential, clocks at the geocenter
//! run slower than clocks at the barycenter. The L_B value encapsulates this difference.
//!
//! # Accumulated Offset at J2000.0
//!
//! At J2000.0 (about 23 years after the 1977 reference epoch), TCB is approximately
//! 11.25 seconds ahead of TCG:
//!
//! ```text
//! TCB - TCG at J2000.0 = L_B * 23 years * 86400 * 365.25 days/year
//!                      = 1.550519768e-8 * 7.26e8 seconds
//!                      ≈ 11.25 seconds
//! ```
//!
//! # Precision
//!
//! Round-trip conversions (TCG -> TCB -> TCG or TCB -> TCG -> TCB) achieve sub-picosecond
//! accuracy. The implementation applies corrections to the smaller-magnitude Julian Date
//! component to preserve precision.
//!
//! # Usage
//!
//! ```
//! use celestial_time::scales::{TCG, TCB};
//! use celestial_time::scales::conversions::{ToTCB, ToTCGFromTCB};
//! use celestial_time::julian::JulianDate;
//! use celestial_core::constants::J2000_JD;
//!
//! let tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, 0.0));
//! let tcb = tcg.to_tcb().unwrap();
//!
//! // At J2000.0, TCB is about 11.25 seconds ahead of TCG
//! let offset_days = tcb.to_julian_date().to_f64() - tcg.to_julian_date().to_f64();
//! assert!(offset_days > 0.0, "TCB should be ahead of TCG after 1977");
//! ```
//!
//! # References
//!
//! - IAU 2006 Resolution B3: Re-definition of Barycentric Dynamical Time, TDB
//! - IAU 2000 Resolution B1.9: Definition of TCG
//! - IERS Conventions (2010), Chapter 10: General Relativistic Models for Time
//! - Petit & Luzum (2010): IERS Technical Note 36

use crate::constants::{TCB_RATE_LB, TCB_RATE_RATIO, TCB_REFERENCE_EPOCH};
use crate::julian::JulianDate;
use crate::scales::{TCB, TCG};
use crate::TimeResult;
use celestial_core::constants::MJD_ZERO_POINT;

/// Convert Geocentric Coordinate Time (TCG) to Barycentric Coordinate Time (TCB).
///
/// TCB runs faster than TCG by the L_B rate factor due to the solar system's
/// gravitational potential at Earth's orbit. This trait applies the rate correction
/// accumulated since the 1977 reference epoch.
pub trait ToTCB {
    /// Convert to Barycentric Coordinate Time (TCB).
    ///
    /// Applies: `TCB = TCG + L_B / (1 - L_B) * (TCG - T0)`
    ///
    /// At J2000.0, this adds approximately 11.25 seconds.
    fn to_tcb(&self) -> TimeResult<TCB>;
}

/// Convert Barycentric Coordinate Time (TCB) to Geocentric Coordinate Time (TCG).
///
/// This is the inverse of [`ToTCB`]. The conversion removes the L_B rate difference
/// that accumulates between the barycentric and geocentric coordinate times.
pub trait ToTCGFromTCB {
    /// Convert to Geocentric Coordinate Time (TCG).
    ///
    /// Applies: `TCG = TCB - L_B * (TCB - T0)`
    ///
    /// At J2000.0, this subtracts approximately 11.25 seconds.
    fn to_tcg(&self) -> TimeResult<TCG>;
}

impl ToTCB for TCB {
    /// Identity conversion. Returns self unchanged.
    fn to_tcb(&self) -> TimeResult<TCB> {
        Ok(*self)
    }
}

impl ToTCB for TCG {
    /// Convert TCG to TCB by applying the L_B rate correction.
    ///
    /// Uses `L_B / (1 - L_B)` as the rate ratio for the forward transformation.
    /// This ratio accounts for the fact that we're computing TCB from TCG, not vice versa.
    ///
    /// The correction is computed relative to the 1977 reference epoch where TCG = TCB.
    /// Applies the correction to the smaller-magnitude JD component for precision.
    fn to_tcb(&self) -> TimeResult<TCB> {
        let tcg_jd = self.to_julian_date();

        let (tcb_jd1, tcb_jd2) = if tcg_jd.jd1().abs() > tcg_jd.jd2().abs() {
            let correction = ((tcg_jd.jd1() - MJD_ZERO_POINT)
                + (tcg_jd.jd2() - TCB_REFERENCE_EPOCH))
                * TCB_RATE_RATIO;
            (tcg_jd.jd1(), tcg_jd.jd2() + correction)
        } else {
            let correction = ((tcg_jd.jd2() - MJD_ZERO_POINT)
                + (tcg_jd.jd1() - TCB_REFERENCE_EPOCH))
                * TCB_RATE_RATIO;
            (tcg_jd.jd1() + correction, tcg_jd.jd2())
        };

        let tcb_jd = JulianDate::new(tcb_jd1, tcb_jd2);
        Ok(TCB::from_julian_date(tcb_jd))
    }
}

impl ToTCGFromTCB for TCB {
    /// Convert TCB to TCG by removing the L_B rate correction.
    ///
    /// Computes: `TCG = TCB - L_B * (JD_TCB - T0)`
    ///
    /// The correction is subtracted because TCB runs faster than TCG.
    /// At J2000.0, this removes about 11.25 seconds.
    ///
    /// Applies the correction to the smaller-magnitude JD component for precision.
    fn to_tcg(&self) -> TimeResult<TCG> {
        let tcb_jd = self.to_julian_date();

        let (tcg_jd1, tcg_jd2) = if tcb_jd.jd1().abs() > tcb_jd.jd2().abs() {
            let correction = ((tcb_jd.jd1() - MJD_ZERO_POINT)
                + (tcb_jd.jd2() - TCB_REFERENCE_EPOCH))
                * TCB_RATE_LB;
            (tcb_jd.jd1(), tcb_jd.jd2() - correction)
        } else {
            let correction = ((tcb_jd.jd2() - MJD_ZERO_POINT)
                + (tcb_jd.jd1() - TCB_REFERENCE_EPOCH))
                * TCB_RATE_LB;
            (tcb_jd.jd1() - correction, tcb_jd.jd2())
        };

        let tcg_jd = JulianDate::new(tcg_jd1, tcg_jd2);
        Ok(TCG::from_julian_date(tcg_jd))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::MJD_1977_JAN_1;
    use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT, SECONDS_PER_DAY_F64};

    #[test]
    fn test_identity_conversions() {
        let tcb = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let identity_tcb = tcb.to_tcb().unwrap();

        assert_eq!(
            tcb.to_julian_date().jd1(),
            identity_tcb.to_julian_date().jd1()
        );
        assert_eq!(
            tcb.to_julian_date().jd2(),
            identity_tcb.to_julian_date().jd2()
        );

        let tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, 0.999999999999999));
        let tcb_converted = tcg.to_tcb().unwrap();
        let tcg_back = tcb_converted.to_tcg().unwrap();
        let tcb_again = tcg_back.to_tcb().unwrap();

        assert_eq!(
            tcb_converted.to_julian_date().jd1(),
            tcb_again.to_julian_date().jd1()
        );
        assert_eq!(
            tcb_converted.to_julian_date().jd2(),
            tcb_again.to_julian_date().jd2()
        );
    }

    #[test]
    fn test_tcg_tcb_offset_at_j2000() {
        let tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let tcb = tcg.to_tcb().unwrap();
        let tcb_jd = tcb.to_julian_date().to_f64();

        assert!(tcb_jd > J2000_JD, "TCB should be ahead of TCG");

        let diff_seconds = (tcb_jd - J2000_JD) * SECONDS_PER_DAY_F64;
        assert!(
            diff_seconds > 11.0 && diff_seconds < 12.0,
            "TCB-TCG at J2000.0 should be ~11.25 seconds: {:.6} seconds",
            diff_seconds
        );

        let tcb_at_j2000 = TCB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let tcg_from_tcb = tcb_at_j2000.to_tcg().unwrap();
        let tcg_jd = tcg_from_tcb.to_julian_date().to_f64();

        assert!(tcg_jd < J2000_JD, "TCG should be behind TCB");

        let reverse_diff = (J2000_JD - tcg_jd) * SECONDS_PER_DAY_F64;
        assert!(
            reverse_diff > 11.0 && reverse_diff < 12.0,
            "TCG-TCB reverse difference should be ~11.25s: {:.6} seconds",
            reverse_diff
        );
    }

    #[test]
    fn test_tcg_tcb_rate_relationship() {
        assert_eq!(TCB_RATE_LB, 1.550519768e-8);

        let reference_epoch = TCG::from_julian_date(JulianDate::new(TCB_REFERENCE_EPOCH, 0.0));
        let one_day_later = TCG::from_julian_date(JulianDate::new(TCB_REFERENCE_EPOCH + 1.0, 0.0));

        let tcb_ref = reference_epoch.to_tcb().unwrap();
        let tcb_day = one_day_later.to_tcb().unwrap();

        let tcb_diff = tcb_day.to_julian_date().to_f64() - tcb_ref.to_julian_date().to_f64();
        let expected_diff = 1.0 + TCB_RATE_LB / (1.0 - TCB_RATE_LB);

        let relative_error = (tcb_diff - expected_diff).abs() / expected_diff;
        assert!(
            relative_error < 1e-12,
            "TCB rate should match expected relativistic correction: {:.2e}",
            relative_error
        );

        let ten_years_days = 3652.5;
        let tcg_j2000 = TCG::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let tcg_j2010 = TCG::from_julian_date(JulianDate::new(J2000_JD + ten_years_days, 0.0));

        let tcb_j2000 = tcg_j2000.to_tcb().unwrap();
        let tcb_j2010 = tcg_j2010.to_tcb().unwrap();

        let tcb_interval =
            tcb_j2010.to_julian_date().to_f64() - tcb_j2000.to_julian_date().to_f64();
        let expected_drift = ten_years_days * TCB_RATE_RATIO;
        let actual_drift = tcb_interval - ten_years_days;

        let drift_error = (actual_drift - expected_drift).abs() / expected_drift;
        assert!(
            drift_error < 1e-4,
            "10-year secular drift error: {:.2e}",
            drift_error
        );
    }

    #[test]
    fn test_tcg_tcb_round_trip_precision() {
        let tolerance = 1e-14;
        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345];

        for jd2 in test_jd2_values {
            let tcg = TCG::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tcb = tcg.to_tcb().unwrap();
            let back_tcg = tcb.to_tcg().unwrap();

            let total_diff =
                (tcg.to_julian_date().to_f64() - back_tcg.to_julian_date().to_f64()).abs();
            assert!(
                total_diff < tolerance,
                "TCG round trip for jd2={} exceeded tolerance: {:.2e}",
                jd2,
                total_diff
            );

            let tcb_rt = TCB::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tcg_from = tcb_rt.to_tcg().unwrap();
            let back_tcb = tcg_from.to_tcb().unwrap();

            let tcb_diff =
                (tcb_rt.to_julian_date().to_f64() - back_tcb.to_julian_date().to_f64()).abs();
            assert!(
                tcb_diff < tolerance,
                "TCB round trip for jd2={} exceeded tolerance: {:.2e}",
                jd2,
                tcb_diff
            );
        }

        let tcg_alt = TCG::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let tcb_alt = tcg_alt.to_tcb().unwrap();
        let back_alt = tcb_alt.to_tcg().unwrap();
        let alt_diff =
            (tcg_alt.to_julian_date().to_f64() - back_alt.to_julian_date().to_f64()).abs();
        assert!(
            alt_diff < tolerance,
            "Alternate JD split round trip exceeded tolerance: {:.2e}",
            alt_diff
        );

        let tcb_alt2 = TCB::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let tcg_alt2 = tcb_alt2.to_tcg().unwrap();
        let back_alt2 = tcg_alt2.to_tcb().unwrap();
        let alt2_diff =
            (tcb_alt2.to_julian_date().to_f64() - back_alt2.to_julian_date().to_f64()).abs();
        assert!(
            alt2_diff < tolerance,
            "Alternate TCB split round trip exceeded tolerance: {:.2e}",
            alt2_diff
        );
    }

    #[test]
    fn test_reference_epoch_behavior() {
        let tcg_at_ref =
            TCG::from_julian_date(JulianDate::new(MJD_ZERO_POINT + MJD_1977_JAN_1, 0.0));
        let tcb_at_ref = tcg_at_ref.to_tcb().unwrap();

        let diff_seconds = (tcb_at_ref.to_julian_date().to_f64()
            - tcg_at_ref.to_julian_date().to_f64())
            * SECONDS_PER_DAY_F64;

        assert!(
            diff_seconds.abs() < 1.0,
            "At 1977 Jan 1 reference epoch, TCB and TCG should be nearly equal: {:.6} seconds",
            diff_seconds
        );
    }
}
