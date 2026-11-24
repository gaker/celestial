// Times in the report, truncated to the second.

use std::str::FromStr;

use celestial_time::scales::utc::UTC;
use celestial_time::TimeResult;

// `utc` as YYYY-MM-DD HH:MM:SS, truncated to the second.
pub(super) fn to_the_second(utc: &UTC) -> String {
    match truncated_iso8601(utc) {
        Ok(text) => text.replacen('T', " ", 1),
        Err(e) => format!("(no calendar date: {e})"),
    }
}

// The ISO 8601 form of `utc` is rounded to the millisecond, so from 0.5 ms before a
// second it names that second. When the second it names starts after `utc`, the text is
// that of the instant half a second earlier, which rounding to the millisecond does not
// carry into the next second. celestial_time then handles leap seconds and day ends.
fn truncated_iso8601(utc: &UTC) -> TimeResult<String> {
    let text = whole_seconds(utc)?;
    let start = UTC::from_str(&text)?;
    let (s, u) = (start.to_julian_date(), utc.to_julian_date());
    if (s.jd1() - u.jd1()) + (s.jd2() - u.jd2()) > 0.0 {
        return whole_seconds(&start.add_seconds(-0.5));
    }
    Ok(text)
}

// The ISO 8601 form of `utc` without the fraction of a second.
fn whole_seconds(utc: &UTC) -> TimeResult<String> {
    let iso = utc.to_iso8601()?;
    Ok(iso.split('.').next().unwrap_or_default().to_string())
}
