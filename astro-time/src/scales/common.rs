use crate::{TimeError, TimeResult};

pub fn get_tai_utc_offset(year: i32, month: i32, day: i32, fraction: f64) -> f64 {
    use crate::constants::{PRE_LEAP_SECOND_ENTRIES, TAI_UTC_OFFSETS, UTC_DRIFT_CORRECTIONS};

    if !(0.0..=1.0).contains(&fraction) {
        return 0.0;
    }

    let my = (month - 14) / 12;
    let iypmy = year + my;
    let modified_jd = ((1461 * (iypmy + 4800)) / 4 + (367 * (month - 2 - 12 * my)) / 12
        - (3 * ((iypmy + 4900) / 100)) / 4
        + day
        - 2432076) as f64;

    if year < TAI_UTC_OFFSETS[0].0 {
        return 0.0;
    }

    let m = 12 * year + month;

    let i = match TAI_UTC_OFFSETS
        .binary_search_by(|&(entry_year, entry_month, _)| (12 * entry_year + entry_month).cmp(&m))
    {
        Ok(idx) => idx, // Exact match found
        Err(idx) => {
            if idx == 0 {
                return 0.0; // Before the first entry
            }
            idx - 1 // Use the entry just before the insertion point
        }
    };

    let mut tai_minus_utc = TAI_UTC_OFFSETS[i].2;

    if i < PRE_LEAP_SECOND_ENTRIES {
        let (drift_mjd, drift_rate) = UTC_DRIFT_CORRECTIONS[i];
        tai_minus_utc += (modified_jd + fraction - drift_mjd) * drift_rate;
    }

    tai_minus_utc
}

pub fn next_calendar_day(year: i32, month: i32, day: i32) -> TimeResult<(i32, i32, i32)> {
    let days_in_month = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if is_leap_year(year) {
                29
            } else {
                28
            }
        }
        _ => {
            return Err(TimeError::ConversionError(format!(
                "Invalid month: {}",
                month
            )))
        }
    };

    if day < days_in_month {
        Ok((year, month, day + 1))
    } else if month < 12 {
        Ok((year, month + 1, 1))
    } else {
        Ok((year + 1, 1, 1))
    }
}

pub fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0)
}
