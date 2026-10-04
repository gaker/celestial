//! GPS and TAI time scale conversions.
//!
//! Provides bidirectional conversion between GPS Time and International Atomic Time (TAI).
//! The relationship is a fixed offset: TAI = GPS + 19 seconds.
//!
//! # Background
//!
//! GPS Time started on January 6, 1980 (GPS epoch) when TAI-GPS was exactly 19 seconds.
//! Unlike UTC, GPS does not include leap seconds, so this offset remains constant.
//!
//! ```text
//! TAI = GPS + 19.0 seconds
//! GPS = TAI - 19.0 seconds
//! ```
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::gps::GPS;
//! use celestial_time::scales::tai::TAI;
//! use celestial_time::scales::conversions::{ToGPS, ToTAI};
//!
//! let gps = GPS::from_julian_date(JulianDate::new(2451545.0, 0.0));
//! let tai = gps.to_tai().unwrap();
//!
//! let back_to_gps = tai.to_gps().unwrap();
//! ```
//!
//! # Precision
//!
//! Only a fixed offset is added or subtracted. Round-trip GPS -> TAI -> GPS comes back
//! within 1 ulp of jd2, from rounding when the offset is added.

use super::{ToGPS, ToTAI};
use crate::constants::GPS_TO_TAI_OFFSET_SECONDS;
use crate::julian::finite_jd;
use crate::scales::gps::GPS;
use crate::scales::tai::TAI;
use crate::TimeResult;
use celestial_core::constants::SECONDS_PER_DAY_F64;

/// GPS to TAI conversion. Adds 19 seconds.
impl ToTAI for GPS {
    /// Converts GPS time to TAI by adding the fixed 19-second offset.
    ///
    /// The offset is added to the smaller-magnitude JD component to preserve precision.
    fn to_tai(&self) -> TimeResult<TAI> {
        let gps_jd = finite_jd(self.to_julian_date())?;
        let offset_days = GPS_TO_TAI_OFFSET_SECONDS / SECONDS_PER_DAY_F64;

        let tai_jd = gps_jd.map_smaller_part(|_, small| small + offset_days);

        Ok(TAI::from_julian_date(tai_jd))
    }
}

/// TAI to GPS conversion. Subtracts 19 seconds.
impl ToGPS for TAI {
    /// Converts TAI to GPS time by subtracting the fixed 19-second offset.
    ///
    /// The offset is subtracted from the smaller-magnitude JD component to preserve precision.
    fn to_gps(&self) -> TimeResult<GPS> {
        let tai_jd = finite_jd(self.to_julian_date())?;
        let offset_days = GPS_TO_TAI_OFFSET_SECONDS / SECONDS_PER_DAY_F64;

        let gps_jd = tai_jd.map_smaller_part(|_, small| small - offset_days);

        Ok(GPS::from_julian_date(gps_jd))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use crate::TimeError;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_gps_tai_offset_19_seconds() {
        let test_dates = [
            (2444244.5, "GPS epoch 1980-01-06"),
            (J2000_JD, "J2000.0"),
            (2455197.5, "2010-01-01"),
            (2459580.5, "2022-01-01"),
        ];

        for (jd, description) in test_dates {
            let gps = GPS::from_julian_date(JulianDate::new(jd, 0.0));
            let tai = gps.to_tai().unwrap();

            let gps_jd = gps.to_julian_date();
            let tai_jd = tai.to_julian_date();

            let offset_days = (tai_jd.jd1() - gps_jd.jd1()) + (tai_jd.jd2() - gps_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, 19.0,
                "{}: GPS->TAI offset must be exactly 19 seconds",
                description
            );

            let tai = TAI::from_julian_date(JulianDate::new(jd, 0.0));
            let gps = tai.to_gps().unwrap();

            let tai_jd = tai.to_julian_date();
            let gps_jd = gps.to_julian_date();

            let offset_days = (tai_jd.jd1() - gps_jd.jd1()) + (tai_jd.jd2() - gps_jd.jd2());
            let offset_seconds = offset_days * SECONDS_PER_DAY_F64;

            assert_eq!(
                offset_seconds, 19.0,
                "{}: TAI->GPS means TAI is 19 seconds ahead",
                description
            );
        }
    }

    #[test]
    fn test_gps_tai_round_trip_precision() {
        let test_jd2_values = [0.0, 0.5, 0.123456789012345, -0.123456789012345, 0.987654321];

        for jd2 in test_jd2_values {
            let original_gps = GPS::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let tai = original_gps.to_tai().unwrap();
            let round_trip_gps = tai.to_gps().unwrap();

            assert_eq!(
                original_gps.to_julian_date().jd1(),
                round_trip_gps.to_julian_date().jd1(),
                "GPS->TAI->GPS JD1 must be exact for jd2={}",
                jd2
            );
            assert_eq!(
                original_gps.to_julian_date().jd2(),
                round_trip_gps.to_julian_date().jd2(),
                "GPS->TAI->GPS JD2 must be exact for jd2={}",
                jd2
            );

            let original_tai = TAI::from_julian_date(JulianDate::new(J2000_JD, jd2));
            let gps = original_tai.to_gps().unwrap();
            let round_trip_tai = gps.to_tai().unwrap();

            assert_eq!(
                original_tai.to_julian_date().jd1(),
                round_trip_tai.to_julian_date().jd1(),
                "TAI->GPS->TAI JD1 must be exact for jd2={}",
                jd2
            );
            assert_eq!(
                original_tai.to_julian_date().jd2(),
                round_trip_tai.to_julian_date().jd2(),
                "TAI->GPS->TAI JD2 must be exact for jd2={}",
                jd2
            );
        }

        let alt_gps = GPS::from_julian_date(JulianDate::new(0.5, J2000_JD));
        let alt_tai = alt_gps.to_tai().unwrap();
        let alt_round_trip = alt_tai.to_gps().unwrap();

        assert_eq!(
            alt_gps.to_julian_date().jd1(),
            alt_round_trip.to_julian_date().jd1(),
            "Alternate split GPS->TAI->GPS JD1 must be exact"
        );
        assert_eq!(
            alt_gps.to_julian_date().jd2(),
            alt_round_trip.to_julian_date().jd2(),
            "Alternate split GPS->TAI->GPS JD2 must be exact"
        );
    }

    #[test]
    fn test_non_finite_julian_date_is_rejected() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let jd = JulianDate::new(bad, 0.0);
            let expected =
                TimeError::InvalidEpoch(format!("Julian Date ({}, 0) is not finite", bad));
            assert_eq!(GPS::from_julian_date(jd).to_tai().unwrap_err(), expected);
            assert_eq!(TAI::from_julian_date(jd).to_gps().unwrap_err(), expected);
        }
    }
}
