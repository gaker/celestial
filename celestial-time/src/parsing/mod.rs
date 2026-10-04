use crate::julian::JulianDate;
use crate::scales::common::validate_calendar_date;
use crate::{TimeError, TimeResult};
use std::fmt;

#[derive(Debug, Clone)]
pub(crate) struct ParsedDateTime {
    pub(crate) year: i32,
    pub(crate) month: u8,
    pub(crate) day: u8,
    pub(crate) hour: u8,
    pub(crate) minute: u8,
    pub(crate) second: f64,
}

impl ParsedDateTime {
    pub(crate) fn to_julian_date(&self) -> TimeResult<JulianDate> {
        JulianDate::from_calendar(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
        )
    }
}

// Z designates UTC, so the UTC parser strips it before calling this and every
// other scale rejects it here.
pub(crate) fn parse_iso8601(s: &str) -> TimeResult<ParsedDateTime> {
    let s = s.trim();
    if s.ends_with('Z') {
        return Err(TimeError::ParseError(format!(
            "Only UTC takes a Z suffix: '{}'",
            s
        )));
    }
    let (date, time) = s.split_once(['T', ' ']).ok_or_else(|| {
        TimeError::ParseError(format!(
            "Invalid datetime format: '{}'. Expected YYYY-MM-DDTHH:MM:SS",
            s
        ))
    })?;
    let [year, month, day] = split_fields(date, '-', "date", "YYYY-MM-DD")?;
    let [hour, minute, second] = split_fields(time, ':', "time", "HH:MM:SS")?;
    let parsed = ParsedDateTime {
        year: i32::from(digits(year, 4, "year")?),
        month: digits(month, 2, "month")? as u8,
        day: digits(day, 2, "day")? as u8,
        hour: digits(hour, 2, "hour")? as u8,
        minute: digits(minute, 2, "minute")? as u8,
        second: parse_seconds(second)?,
    };
    check_ranges(&parsed)?;
    Ok(parsed)
}

fn split_fields<'a>(
    s: &'a str,
    separator: char,
    name: &str,
    shape: &str,
) -> TimeResult<[&'a str; 3]> {
    let mut fields = s.split(separator);
    match (fields.next(), fields.next(), fields.next(), fields.next()) {
        (Some(a), Some(b), Some(c), None) => Ok([a, b, c]),
        _ => Err(TimeError::ParseError(format!(
            "Invalid {} format: '{}'. Expected {}",
            name, s, shape
        ))),
    }
}

// ISO 8601 fields have a fixed width, so `2000-1-1T0:0:0` is not a date.
fn digits(field: &str, width: usize, name: &str) -> TimeResult<u16> {
    if field.len() != width || !field.bytes().all(|b| b.is_ascii_digit()) {
        return Err(TimeError::ParseError(format!(
            "Invalid {}: '{}'",
            name, field
        )));
    }
    Ok(field.bytes().fold(0, |n, b| n * 10 + u16::from(b - b'0')))
}

// A UTC minute can hold a leap second, so 60.x is syntactically valid; each
// scale then checks whether the day it names has that second.
fn check_ranges(parsed: &ParsedDateTime) -> TimeResult<()> {
    validate_calendar_date(parsed.year, parsed.month.into(), parsed.day.into())?;
    if parsed.hour > 23 {
        return Err(out_of_range("Hour", parsed.hour));
    }
    if parsed.minute > 59 {
        return Err(out_of_range("Minute", parsed.minute));
    }
    if parsed.second >= 61.0 {
        return Err(out_of_range("Second", parsed.second));
    }
    Ok(())
}

fn out_of_range(name: &str, value: impl fmt::Display) -> TimeError {
    TimeError::ParseError(format!("{} out of range: {}", name, value))
}

// f64's FromStr also takes signs, exponents, "NaN" and "inf", none of which
// belong in an ISO 8601 time, so check the shape before converting. The
// fraction may have any number of digits.
fn parse_seconds(s: &str) -> TimeResult<f64> {
    let is_digits = |part: &str| !part.is_empty() && part.bytes().all(|b| b.is_ascii_digit());
    let well_formed = match s.split_once('.') {
        Some((whole, fraction)) => whole.len() == 2 && is_digits(whole) && is_digits(fraction),
        None => s.len() == 2 && is_digits(s),
    };
    let invalid = || TimeError::ParseError(format!("Invalid second: '{}'", s));
    if !well_formed {
        return Err(invalid());
    }
    s.parse::<f64>().map_err(|_| invalid())
}

#[cfg(test)]
mod tests;
