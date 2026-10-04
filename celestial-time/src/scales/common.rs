use super::conversions::utc_tai::calendar_to_julian;
use crate::constants::{PRE_LEAP_SECOND_ENTRIES, TAI_UTC_OFFSETS, UTC_DRIFT_CORRECTIONS};
use crate::{TimeError, TimeResult};

// The calendar arithmetic needs iypmy + 4800 >= 0; ERFA's cal2jd has the same limit.
const MIN_CALENDAR_YEAR: i32 = -4799;

pub fn get_tai_utc_offset(year: i32, month: i32, day: i32, fraction: f64) -> TimeResult<f64> {
    if !(0.0..=1.0).contains(&fraction) {
        return Err(TimeError::ConversionError(format!(
            "Day fraction {} is outside [0, 1]",
            fraction
        )));
    }
    let (_, modified_jd) = calendar_to_julian(year, month, day)?;
    let Some(i) = offset_table_index(year, month) else {
        return Ok(0.0);
    };

    let mut tai_minus_utc = TAI_UTC_OFFSETS[i].2;
    if i < PRE_LEAP_SECOND_ENTRIES {
        let (drift_mjd, drift_rate) = UTC_DRIFT_CORRECTIONS[i];
        tai_minus_utc += (modified_jd + fraction - drift_mjd) * drift_rate;
    }
    Ok(tai_minus_utc)
}

fn offset_table_index(year: i32, month: i32) -> Option<usize> {
    match TAI_UTC_OFFSETS.binary_search_by(|&(entry_year, entry_month, _)| {
        (entry_year, entry_month).cmp(&(year, month))
    }) {
        Ok(i) => Some(i),
        Err(i) => i.checked_sub(1),
    }
}

pub(crate) fn next_calendar_day(year: i32, month: i32, day: i32) -> TimeResult<(i32, i32, i32)> {
    validate_calendar_date(year, month, day)?;
    if Some(day) != days_in_month(year, month) {
        Ok((year, month, day + 1))
    } else if month < 12 {
        Ok((year, month + 1, 1))
    } else {
        let next_year = year
            .checked_add(1)
            .ok_or_else(|| invalid_date(year, month, day, "the next day is out of range"))?;
        Ok((next_year, 1, 1))
    }
}

pub(crate) fn validate_calendar_date(year: i32, month: i32, day: i32) -> TimeResult<()> {
    if year < MIN_CALENDAR_YEAR {
        return Err(invalid_date(year, month, day, "year is before -4799"));
    }
    match days_in_month(year, month) {
        None => Err(invalid_date(year, month, day, "month is out of range")),
        Some(last) if !(1..=last).contains(&day) => {
            Err(invalid_date(year, month, day, "day is out of range"))
        }
        Some(_) => Ok(()),
    }
}

fn days_in_month(year: i32, month: i32) -> Option<i32> {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => Some(31),
        4 | 6 | 9 | 11 => Some(30),
        2 if is_leap_year(year) => Some(29),
        2 => Some(28),
        _ => None,
    }
}

fn invalid_date(year: i32, month: i32, day: i32, reason: &str) -> TimeError {
    TimeError::InvalidDate(format!("{:04}-{:02}-{:02}: {}", year, month, day, reason))
}

fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0)
}
