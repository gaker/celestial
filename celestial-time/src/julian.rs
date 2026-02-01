use crate::constants::{SECONDS_TO_DAYS, UNIX_EPOCH_JD};
use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT, SECONDS_PER_DAY_F64};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct JulianDate {
    pub jd1: f64,
    pub jd2: f64,
}

impl JulianDate {
    pub fn new(jd1: f64, jd2: f64) -> Self {
        Self { jd1, jd2 }
    }

    pub fn from_f64(jd: f64) -> Self {
        Self::new(jd, 0.0)
    }

    pub fn j2000() -> Self {
        Self::new(J2000_JD, 0.0)
    }

    pub fn unix_epoch() -> Self {
        Self::new(UNIX_EPOCH_JD, 0.0)
    }

    pub fn jd1(&self) -> f64 {
        self.jd1
    }

    pub fn jd2(&self) -> f64 {
        self.jd2
    }

    pub fn to_f64(&self) -> f64 {
        self.jd1 + self.jd2
    }

    pub fn add_days(&self, days: f64) -> Self {
        Self::new(self.jd1, self.jd2 + days)
    }

    pub fn add_seconds(&self, seconds: f64) -> Self {
        self.add_days(seconds * SECONDS_TO_DAYS)
    }

    pub fn from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> Self {
        // Algorithm matches ERFA's eraCal2jd + eraDtf2d convention:
        // jd1 = full Julian Date at midnight (integer-ish)
        // jd2 = fraction of day
        let my = (month as i32 - 14) / 12;
        let iypmy = year + my;

        // Compute MJD for 0h of the given day (same algorithm as ERFA eraCal2jd)
        let mjd = ((1461 * (iypmy + 4800)) / 4 + (367 * (month as i32 - 2 - 12 * my)) / 12
            - (3 * ((iypmy + 4900) / 100)) / 4
            + day as i32
            - 2432076) as f64;

        // Full Julian Date at midnight = MJD epoch + MJD
        let jd1 = MJD_ZERO_POINT + mjd;

        // Day fraction from time components
        let jd2 = (60.0 * (60 * hour as i32 + minute as i32) as f64 + second) / SECONDS_PER_DAY_F64;

        Self::new(jd1, jd2)
    }
}

impl fmt::Display for JulianDate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "JD {:.9}", self.to_f64())
    }
}

impl From<f64> for JulianDate {
    fn from(jd: f64) -> Self {
        Self::from_f64(jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_julian_date_creation() {
        let jd = JulianDate::new(J2000_JD, 0.5);
        assert_eq!(jd.jd1(), J2000_JD);
        assert_eq!(jd.jd2(), 0.5);
        assert_eq!(jd.to_f64(), 2451545.5);
    }

    #[test]
    fn test_j2000_epoch() {
        let j2000 = JulianDate::j2000();
        assert_eq!(j2000.to_f64(), J2000_JD);
    }

    #[test]
    fn test_unix_epoch() {
        let unix = JulianDate::unix_epoch();
        assert_eq!(unix.to_f64(), crate::constants::UNIX_EPOCH_JD);
    }

    #[test]
    fn test_arithmetic() {
        let jd = JulianDate::new(J2000_JD, 0.0);
        let jd_plus_day = jd.add_days(1.0);
        assert_eq!(jd_plus_day.to_f64(), 2451546.0);

        let jd_plus_hour = jd.add_seconds(3600.0);
        assert!((jd_plus_hour.to_f64() - 2_451_545.041_666_666_5).abs() < 1e-15);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serde_round_trip() {
        let test_cases = [
            JulianDate::new(J2000_JD, 0.0),          // J2000.0
            JulianDate::new(2451545.5, 0.123456789), // J2000.0 + 12h + fraction
            JulianDate::new(2440587.5, 0.0),         // Unix epoch
            JulianDate::new(J2000_JD, 0.999999999),  // High precision
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: JulianDate = serde_json::from_str(&json).unwrap();

            assert_eq!(
                original.jd1(),
                deserialized.jd1(),
                "JD1 precision lost in serde round-trip"
            );
            assert_eq!(
                original.jd2(),
                deserialized.jd2(),
                "JD2 precision lost in serde round-trip"
            );
            assert_eq!(
                original, deserialized,
                "JulianDate equality lost in serde round-trip"
            );
        }
    }
}
