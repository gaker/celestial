//! Barycentric Dynamical Time (TDB) scale.
//!
//! TDB is the independent time argument for barycentric ephemerides of the solar system.
//! It differs from TT by small periodic terms (max ~1.7 ms) due to relativistic effects
//! from Earth's orbital motion around the solar system barycenter.
//!
//! # Background
//!
//! TDB runs at a different rate than TT due to gravitational time dilation and velocity
//! effects. The difference TDB-TT is dominated by a ~1.7 ms amplitude term with a period
//! of one year, plus smaller terms. For most applications, TDB ≈ TT to within 2 ms.
//!
//! The IAU recommends using TCB (Barycentric Coordinate Time) for rigorous relativistic
//! work. TDB is defined as a linear transformation of TCB that keeps TDB-TT bounded.
//!
//! # Usage
//!
//! ```
//! use celestial_time::{JulianDate, TDB};
//!
//! let tdb = TDB::j2000();
//! let tdb_plus_day = tdb.add_days(1.0);
//!
//! let tdb_from_cal = celestial_time::scales::tdb::tdb_from_calendar(2000, 1, 1, 12, 0, 0.0);
//! ```
//!
//! # Precision
//!
//! Internally stores time as a split Julian Date for sub-microsecond precision.
//! Arithmetic operations preserve precision by operating on the underlying JulianDate.

use crate::constants::UNIX_EPOCH_JD;
use crate::julian::JulianDate;
use crate::parsing::parse_iso8601;
use crate::{TimeError, TimeResult};
use celestial_core::constants::SECONDS_PER_DAY_F64;
use std::fmt;
use std::str::FromStr;

/// Barycentric Dynamical Time.
///
/// A time scale for solar system barycentric ephemerides. Wraps a split Julian Date
/// for high-precision arithmetic. TDB tracks TT to within ~2 ms over centuries.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TDB(JulianDate);

impl TDB {
    /// Creates TDB from Unix timestamp components.
    ///
    /// Converts seconds and nanoseconds since 1970-01-01 to TDB Julian Date.
    pub fn new(seconds: i64, nanos: u32) -> Self {
        let total_seconds =
            seconds as f64 + nanos as f64 / celestial_core::constants::NANOSECONDS_PER_SECOND_F64;
        let jd = JulianDate::from_f64(UNIX_EPOCH_JD + total_seconds / SECONDS_PER_DAY_F64);
        Self(jd)
    }

    /// Creates TDB from a Julian Date.
    pub fn from_julian_date(jd: JulianDate) -> Self {
        Self(jd)
    }

    /// Returns TDB at the J2000.0 epoch (2000-01-01T12:00:00 TDB).
    pub fn j2000() -> Self {
        Self(JulianDate::j2000())
    }

    /// Returns the underlying Julian Date.
    pub fn to_julian_date(&self) -> JulianDate {
        self.0
    }

    /// Adds seconds to this TDB instant, returning a new TDB.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self(self.0.add_seconds(seconds))
    }

    /// Adds days to this TDB instant, returning a new TDB.
    pub fn add_days(&self, days: f64) -> Self {
        Self(self.0.add_days(days))
    }
}

impl fmt::Display for TDB {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TDB {}", self.0)
    }
}

/// Conversion from JulianDate. No transformation applied.
impl From<JulianDate> for TDB {
    fn from(jd: JulianDate) -> Self {
        Self::from_julian_date(jd)
    }
}

/// Parses ISO 8601 string as TDB.
///
/// The string is interpreted directly as TDB without any scale conversion.
impl FromStr for TDB {
    type Err = TimeError;

    fn from_str(s: &str) -> TimeResult<Self> {
        let parsed = parse_iso8601(s)?;
        Ok(Self::from_julian_date(parsed.to_julian_date()))
    }
}

/// Creates TDB from calendar components.
///
/// Interprets the calendar date directly as TDB. For high-precision work,
/// prefer constructing from a Julian Date directly.
pub fn tdb_from_calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> TDB {
    let jd = JulianDate::from_calendar(year, month, day, hour, minute, second);
    TDB::from_julian_date(jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::UNIX_EPOCH_JD;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_tdb_constructors() {
        assert_eq!(TDB::new(0, 0).to_julian_date().to_f64(), UNIX_EPOCH_JD);
        assert_eq!(TDB::j2000().to_julian_date().to_f64(), J2000_JD);
        assert_eq!(
            tdb_from_calendar(2000, 1, 1, 12, 0, 0.0).to_julian_date().to_f64(),
            J2000_JD
        );
    }

    #[test]
    fn test_tdb_arithmetic() {
        let tdb = TDB::j2000();
        assert_eq!(tdb.add_days(1.0).to_julian_date().to_f64(), J2000_JD + 1.0);
        assert_eq!(
            tdb.add_seconds(3600.0).to_julian_date().to_f64(),
            J2000_JD + 1.0 / 24.0
        );
    }

    #[test]
    fn test_tdb_from_julian_date_trait() {
        let jd = JulianDate::new(J2000_JD, 0.123456789);
        let tdb_direct = TDB::from_julian_date(jd);
        let tdb_from_trait: TDB = jd.into();

        assert_eq!(
            tdb_direct.to_julian_date().jd1(),
            tdb_from_trait.to_julian_date().jd1()
        );
        assert_eq!(
            tdb_direct.to_julian_date().jd2(),
            tdb_from_trait.to_julian_date().jd2()
        );
    }

    #[test]
    fn test_tdb_display() {
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let display_str = format!("{}", tdb);

        assert!(display_str.starts_with("TDB"));
        assert!(display_str.contains("2451545"));
    }

    #[test]
    fn test_tdb_string_parsing() {
        assert_eq!(
            TDB::from_str("2000-01-01T12:00:00").unwrap().to_julian_date().to_f64(),
            TDB::j2000().to_julian_date().to_f64()
        );
        assert!(TDB::from_str("invalid-date").is_err());

        let result = TDB::from_str("2000-01-01T12:00:00.123").unwrap();
        let expected_jd = J2000_JD + 0.123 / SECONDS_PER_DAY_F64;
        let diff = (result.to_julian_date().to_f64() - expected_jd).abs();
        assert!(diff < 1e-14, "fractional seconds diff: {:.2e}", diff);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_tdb_serde_round_trip() {
        let test_cases = [
            TDB::j2000(),
            TDB::new(0, 0),
            tdb_from_calendar(2024, 6, 15, 14, 30, 45.123),
            tdb_from_calendar(1990, 12, 31, 23, 59, 59.999999999),
        ];

        for original in test_cases {
            let json = serde_json::to_string(&original).unwrap();
            let deserialized: TDB = serde_json::from_str(&json).unwrap();

            let jd1_diff =
                (original.to_julian_date().jd1() - deserialized.to_julian_date().jd1()).abs();
            let jd2_diff =
                (original.to_julian_date().jd2() - deserialized.to_julian_date().jd2()).abs();
            let total_diff =
                (original.to_julian_date().to_f64() - deserialized.to_julian_date().to_f64()).abs();

            assert!(jd1_diff < 1e-14, "jd1 diff: {:.2e}", jd1_diff);
            assert!(jd2_diff < 1e-14, "jd2 diff: {:.2e}", jd2_diff);
            assert!(total_diff < 1e-14, "total diff: {:.2e}", total_diff);
        }
    }
}
