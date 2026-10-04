use crate::constants::UNIX_EPOCH_JD;
use crate::scales::conversions::utc_tai::calendar_to_julian;
use crate::{TimeError, TimeResult};
use celestial_core::constants::{
    DAYS_PER_JULIAN_YEAR, J2000_JD, MJD_ZERO_POINT, NANOSECONDS_PER_SECOND_F64, SECONDS_PER_DAY,
    SECONDS_PER_DAY_F64,
};
use std::fmt;
use std::ops::Sub;

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct JulianDate {
    jd1: f64,
    jd2: f64,
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

    // Whole days are carried into jd1. Left in jd2 they would grow it until
    // small steps (a microsecond after a century) round away.
    pub fn add_days(&self, days: f64) -> Self {
        let whole_days = libm::trunc(days);
        let jd2 = self.jd2 + (days - whole_days);
        let carry = libm::trunc(jd2);
        Self::new(self.jd1 + (whole_days + carry), jd2 - carry)
    }

    pub fn add_seconds(&self, seconds: f64) -> Self {
        self.add_days(seconds / SECONDS_PER_DAY_F64)
    }

    // Whole days and the time of day are kept apart so nanoseconds survive.
    pub(crate) fn from_unix_time(seconds: i64, nanos: u32) -> Self {
        let jd1 = UNIX_EPOCH_JD + (seconds / SECONDS_PER_DAY) as f64;
        let second_of_day =
            (seconds % SECONDS_PER_DAY) as f64 + f64::from(nanos) / NANOSECONDS_PER_SECOND_F64;
        Self::new(jd1, second_of_day / SECONDS_PER_DAY_F64)
    }

    pub fn from_calendar(
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: f64,
    ) -> TimeResult<Self> {
        // Same split as ERFA's dtf2d: midnight's JD in jd1, the time of day in jd2.
        let (mjd_zero, mjd) = calendar_to_julian(year, month.into(), day.into())?;
        let jd2 = day_fraction(hour, minute, second, 0.0)?;
        Ok(Self::new(mjd_zero + mjd, jd2))
    }

    // Subtracting J2000 from jd1 before adding jd2 keeps jd2's low bits, as
    // eraEpj does.
    pub fn to_julian_year(&self) -> f64 {
        2000.0 + ((self.jd1 - J2000_JD) + self.jd2) / DAYS_PER_JULIAN_YEAR
    }

    // eraEpj2jd's split: the MJD zero point in jd1 and the MJD in jd2.
    pub fn from_julian_year(year: f64) -> Self {
        let mjd = (J2000_JD - MJD_ZERO_POINT) + (year - 2000.0) * DAYS_PER_JULIAN_YEAR;
        Self::new(MJD_ZERO_POINT, mjd)
    }

    // Offsets and rate corrections rewrite only the smaller-magnitude part, so the
    // larger part stays exact. On a tie jd1 counts as the smaller, as in ERFA.
    pub(crate) fn map_smaller_part(self, f: impl FnOnce(f64, f64) -> f64) -> Self {
        if libm::fabs(self.jd1) > libm::fabs(self.jd2) {
            Self::new(self.jd1, f(self.jd1, self.jd2))
        } else {
            Self::new(f(self.jd2, self.jd1), self.jd2)
        }
    }
}

// `leap` is how much longer than 86,400 s the day is: the leap second at the
// end of a UTC leap-second day, and zero for every other scale. The extra time
// goes into the day and into the seconds of its last minute.
pub(crate) fn day_fraction(hour: u8, minute: u8, second: f64, leap: f64) -> TimeResult<f64> {
    let second_limit = if hour == 23 && minute == 59 {
        60.0 + leap
    } else {
        60.0
    };
    if hour > 23 || minute > 59 || !(0.0..second_limit).contains(&second) {
        return Err(TimeError::InvalidDate(format!(
            "time {:02}:{:02}:{} is out of range",
            hour, minute, second
        )));
    }
    let seconds_into_day = 60.0 * f64::from(60 * u32::from(hour) + u32::from(minute)) + second;
    Ok(seconds_into_day / (SECONDS_PER_DAY_F64 + leap))
}

pub(crate) fn finite_jd(jd: JulianDate) -> TimeResult<JulianDate> {
    if jd.jd1.is_finite() && jd.jd2.is_finite() {
        return Ok(jd);
    }
    Err(TimeError::InvalidEpoch(format!(
        "Julian Date ({}, {}) is not finite",
        jd.jd1, jd.jd2
    )))
}

pub(crate) fn finite_arg(name: &str, value: f64) -> TimeResult<f64> {
    if value.is_finite() {
        return Ok(value);
    }
    Err(TimeError::ConversionError(format!(
        "{} must be finite, got {}",
        name, value
    )))
}

impl fmt::Display for JulianDate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "JD {} + {}", self.jd1, self.jd2)
    }
}

// Dates are equal when jd1 + jd2 is the same number, however it is split. The
// rounded sum and its exact rounding error together pin that number, where the
// rounded sum alone would merge dates less than an ulp apart.
impl PartialEq for JulianDate {
    fn eq(&self, other: &Self) -> bool {
        let (sum, error) = exact_sum(self.jd1, self.jd2);
        let (other_sum, other_error) = exact_sum(other.jd1, other.jd2);
        if !sum.is_finite() || !other_sum.is_finite() {
            // The rounding error is NaN here, so only identical parts can match.
            return (self.jd1, self.jd2) == (other.jd1, other.jd2);
        }
        sum == other_sum && error == other_error
    }
}

// Knuth's two-sum: a + b == sum + error exactly, for any finite sum.
fn exact_sum(a: f64, b: f64) -> (f64, f64) {
    let sum = a + b;
    let b_rounded = sum - a;
    let error = (a - (sum - b_rounded)) + (b - b_rounded);
    (sum, error)
}

// Tests that pin ERFA's split compare the parts, since `==` ignores the split.
#[cfg(test)]
impl JulianDate {
    pub(crate) fn parts(&self) -> (f64, f64) {
        (self.jd1, self.jd2)
    }
}

// The difference in days. Subtracting part by part keeps the precision of a
// short interval between two large dates.
impl Sub<Self> for JulianDate {
    type Output = f64;

    fn sub(self, other: Self) -> f64 {
        (self.jd1 - other.jd1) + (self.jd2 - other.jd2)
    }
}

#[cfg(test)]
mod tests;
