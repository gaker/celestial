use chrono::{Datelike, NaiveDate, NaiveDateTime, NaiveTime, Timelike};

use celestial_images::formats::Image;
use celestial_time::{utc_from_calendar, JulianDate};

use super::error::MetadataError;

/// Reads the observation epoch from an image's headers as a UTC Julian Date.
///
/// Tries headers in this order:
///
/// 1. `MJD-OBS` — numeric Modified Julian Date (highest priority, highest precision)
/// 2. `DATE-OBS` — ISO 8601 timestamp (e.g. `"2024-06-15T03:22:45.123"`) or a date-only
///    string paired with `TIME-OBS` for the time component
///
/// # Errors
///
/// Returns [`MetadataError::MissingHeader`] if neither header exists,
/// [`MetadataError::MissingTimeOfDay`] if `DATE-OBS` is date-only without `TIME-OBS`,
/// or [`MetadataError::DateParse`] / [`MetadataError::TimeParse`] for malformed values.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::jd_from_image;
///
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// let jd = jd_from_image(&img)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn jd_from_image(img: &Image) -> Result<JulianDate, MetadataError> {
    if let Some(kw) = img.get_keyword("MJD-OBS") {
        let mjd = kw
            .value
            .as_ref()
            .and_then(|v| v.as_real())
            .ok_or(MetadataError::InvalidHeaderType {
                keyword: "MJD-OBS",
                kind: "numeric",
            })?;
        return Ok(JulianDate::new(2_400_000.5, mjd));
    }

    let date_obs = img
        .get_keyword("DATE-OBS")
        .and_then(|kw| kw.value.as_ref()?.as_string().map(|s| s.to_owned()))
        .ok_or(MetadataError::MissingHeader("DATE-OBS or MJD-OBS"))?;

    let dt = parse_date_obs(img, &date_obs)?;
    let sec = dt.second() as f64 + dt.nanosecond() as f64 / 1_000_000_000.0;
    let utc = utc_from_calendar(
        dt.year(),
        dt.month() as u8,
        dt.day() as u8,
        dt.hour() as u8,
        dt.minute() as u8,
        sec,
    );
    Ok(utc.to_julian_date())
}

fn parse_date_obs(img: &Image, date_obs: &str) -> Result<NaiveDateTime, MetadataError> {
    if date_obs.contains('T') {
        return NaiveDateTime::parse_from_str(date_obs, "%Y-%m-%dT%H:%M:%S%.f").map_err(|e| {
            MetadataError::DateParse {
                value: date_obs.to_string(),
                source: e,
            }
        });
    }

    let date =
        NaiveDate::parse_from_str(date_obs, "%Y-%m-%d").map_err(|e| MetadataError::DateParse {
            value: date_obs.to_string(),
            source: e,
        })?;
    let time = img
        .get_keyword("TIME-OBS")
        .and_then(|kw| kw.value.as_ref()?.as_string().map(|s| s.to_owned()))
        .ok_or(MetadataError::MissingTimeOfDay)?;
    let t = NaiveTime::parse_from_str(&time, "%H:%M:%S%.f").map_err(|e| {
        MetadataError::TimeParse {
            value: time,
            source: e,
        }
    })?;
    Ok(date.and_time(t))
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::Keyword;
    use celestial_images::formats::PixelData;

    fn stub_image() -> Image {
        Image::new(PixelData::F32(vec![0.0; 4]), [2, 2])
    }

    #[test]
    fn mjd_obs_takes_priority_and_is_returned_verbatim() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("MJD-OBS", 60476.140798));
        img.set_keyword(Keyword::string("DATE-OBS", "1999-01-01T00:00:00.000"));
        let jd = jd_from_image(&img).unwrap();
        assert_eq!(jd.to_f64(), 2_400_000.5 + 60476.140798);
    }

    #[test]
    fn mjd_obs_non_numeric_is_rejected() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("MJD-OBS", "bogus"));
        let err = jd_from_image(&img).unwrap_err();
        assert!(
            matches!(err, MetadataError::InvalidHeaderType { keyword: "MJD-OBS", .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn iso_datetime_is_parsed_to_matching_jd() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15T03:22:45.123"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 6, 15, 3, 22, 45.123).to_julian_date();
        assert_eq!(jd.to_f64(), expected.to_f64());
    }

    #[test]
    fn date_plus_separate_time_obs_yields_same_jd_as_iso() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15"));
        img.set_keyword(Keyword::string("TIME-OBS", "03:22:45.123"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 6, 15, 3, 22, 45.123).to_julian_date();
        assert_eq!(jd.to_f64(), expected.to_f64());
    }

    #[test]
    fn missing_all_headers_returns_missing_header() {
        let img = stub_image();
        let err = jd_from_image(&img).unwrap_err();
        assert!(matches!(
            err,
            MetadataError::MissingHeader("DATE-OBS or MJD-OBS")
        ));
    }

    #[test]
    fn date_only_without_time_obs_returns_missing_time_of_day() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15"));
        let err = jd_from_image(&img).unwrap_err();
        assert!(matches!(err, MetadataError::MissingTimeOfDay));
    }

    #[test]
    fn malformed_iso_date_returns_date_parse_error() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-13-40T99:99:99"));
        let err = jd_from_image(&img).unwrap_err();
        assert!(
            matches!(&err, MetadataError::DateParse { value, .. } if value == "2024-13-40T99:99:99"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn malformed_date_without_t_returns_date_parse_error() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "not-a-date"));
        let err = jd_from_image(&img).unwrap_err();
        assert!(matches!(err, MetadataError::DateParse { .. }));
    }

    #[test]
    fn malformed_time_obs_returns_time_parse_error() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15"));
        img.set_keyword(Keyword::string("TIME-OBS", "not-a-time"));
        let err = jd_from_image(&img).unwrap_err();
        assert!(
            matches!(&err, MetadataError::TimeParse { value, .. } if value == "not-a-time"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn fractional_seconds_are_preserved_in_jd() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15T12:00:00.5"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 6, 15, 12, 0, 0.5).to_julian_date();
        assert_eq!(jd.to_f64(), expected.to_f64());
    }
}
