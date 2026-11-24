use std::str::FromStr;

use celestial_images::formats::image::Image;
use celestial_time::julian::JulianDate;
use celestial_time::scales::utc::UTC;
use celestial_time::TimeError;

use super::error::MetadataError;

/// Reads the observation epoch from an image's headers as a UTC Julian Date.
///
/// Tries headers in this order:
///
/// 1. `MJD-OBS`, a numeric Modified Julian Date (highest priority, highest precision)
/// 2. `DATE-OBS`, an ISO 8601 timestamp such as `"2024-06-15T03:22:45.123"`, with an
///    optional trailing `Z`, or a date-only string paired with `TIME-OBS` for the time
///    of day
///
/// A key that repeats is read from its last card.
///
/// # Errors
///
/// Returns [`MetadataError::MissingHeader`] if neither header exists,
/// [`MetadataError::NonFiniteValue`] if `MJD-OBS` is NaN or infinite,
/// [`MetadataError::MissingTimeOfDay`] if `DATE-OBS` is date-only without `TIME-OBS`,
/// [`MetadataError::DateParse`] / [`MetadataError::TimeParse`] for malformed values,
/// or [`MetadataError::InvalidTime`] for a time that UTC doesn't have, such as
/// `23:59:60` on a day without a leap second.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::time::jd_from_image;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// let jd = jd_from_image(&img)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn jd_from_image(img: &Image) -> Result<JulianDate, MetadataError> {
    if img.get_last("MJD-OBS").exists() {
        return jd_from_mjd_obs(img);
    }
    let utc = utc_from_date_obs(img)?.ok_or(MetadataError::MissingHeader("DATE-OBS or MJD-OBS"))?;
    Ok(utc.to_julian_date())
}

fn jd_from_mjd_obs(img: &Image) -> Result<JulianDate, MetadataError> {
    let mjd = img
        .get_last("MJD-OBS")
        .as_f64()
        .ok_or(MetadataError::InvalidHeaderType {
            keyword: "MJD-OBS",
            kind: "numeric",
        })?;
    if !mjd.is_finite() {
        return Err(MetadataError::NonFiniteValue {
            field: "MJD-OBS",
            value: mjd,
        });
    }
    Ok(JulianDate::new(2_400_000.5, mjd))
}

// The instant DATE-OBS names, taking the time of day from TIME-OBS when DATE-OBS is a
// date alone. None when there is no DATE-OBS string. celestial_time parses every field.
// The date is parsed on its own first, so that an error names the card that is wrong.
fn utc_from_date_obs(img: &Image) -> Result<Option<UTC>, MetadataError> {
    let Some(date_obs) = img.get_last("DATE-OBS").as_str() else {
        return Ok(None);
    };
    let date_obs = date_obs.trim();
    let date_parse = |source| MetadataError::DateParse {
        value: date_obs.to_string(),
        source,
    };
    // T or a space separates the date from the time, as celestial_time reads them.
    let date = date_obs
        .split_once(['T', ' '])
        .map_or(date_obs, |(date, _)| date);
    UTC::from_str(&format!("{date}T00:00:00")).map_err(date_parse)?;
    if date == date_obs {
        return utc_from_time_obs(img, date).map(Some);
    }
    match UTC::from_str(date_obs) {
        Ok(utc) => Ok(Some(utc)),
        Err(e @ TimeError::ParseError(_)) => Err(date_parse(e)),
        Err(e) => Err(MetadataError::InvalidTime(e)),
    }
}

fn utc_from_time_obs(img: &Image, date: &str) -> Result<UTC, MetadataError> {
    let time = img
        .get_last("TIME-OBS")
        .as_str()
        .ok_or(MetadataError::MissingTimeOfDay)?
        .trim();
    match UTC::from_str(&format!("{date}T{time}")) {
        Ok(utc) => Ok(utc),
        Err(e @ TimeError::ParseError(_)) => Err(MetadataError::TimeParse {
            value: time.to_string(),
            source: e,
        }),
        Err(e) => Err(MetadataError::InvalidTime(e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::keywords::Keyword;
    use celestial_images::formats::pixel_data::PixelData;
    use celestial_time::scales::utc::utc_from_calendar;

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
            matches!(
                err,
                MetadataError::InvalidHeaderType {
                    keyword: "MJD-OBS",
                    ..
                }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn iso_datetime_is_parsed_to_matching_jd() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15T03:22:45.123"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 6, 15, 3, 22, 45.123)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd.to_f64(), expected.to_f64());
    }

    #[test]
    fn date_plus_separate_time_obs_yields_same_jd_as_iso() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15"));
        img.set_keyword(Keyword::string("TIME-OBS", "03:22:45.123"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 6, 15, 3, 22, 45.123)
            .unwrap()
            .to_julian_date();
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
        let expected = utc_from_calendar(2024, 6, 15, 12, 0, 0.5)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd.to_f64(), expected.to_f64());
    }

    #[test]
    fn leap_second_is_accepted_on_a_leap_second_day() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2016-12-31T23:59:60.5"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2016, 12, 31, 23, 59, 60.5)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd, expected);
    }

    #[test]
    fn leap_second_on_an_ordinary_day_returns_invalid_time() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2016-12-30T23:59:60"));
        let err = jd_from_image(&img).unwrap_err();
        assert!(
            matches!(err, MetadataError::InvalidTime(_)),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn non_finite_mjd_obs_is_rejected() {
        let cases = [
            (f64::NAN, "MJD-OBS must be finite, got NaN"),
            (f64::INFINITY, "MJD-OBS must be finite, got inf"),
        ];
        for (mjd, expected) in cases {
            let mut img = stub_image();
            img.set_keyword(Keyword::real("MJD-OBS", mjd));
            let err = jd_from_image(&img).err();
            assert_eq!(err.map(|e| e.to_string()).as_deref(), Some(expected));
        }
    }

    #[test]
    fn date_obs_with_a_z_suffix_is_read_as_utc() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-01-15T03:04:05.123Z"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 1, 15, 3, 4, 5.123)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd, expected);
    }

    #[test]
    fn date_only_with_time_obs_on_a_leap_second_day_keeps_the_leap_second() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2016-12-31"));
        img.set_keyword(Keyword::string("TIME-OBS", "23:59:60.5"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2016, 12, 31, 23, 59, 60.5)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd, expected);
    }

    #[test]
    fn repeated_date_obs_reads_the_last_card() {
        let mut img = stub_image();
        img.keywords
            .push(Keyword::string("DATE-OBS", "1999-01-01T00:00:00"));
        img.keywords
            .push(Keyword::string("DATE-OBS", "2024-06-15T03:22:45.123"));
        let jd = jd_from_image(&img).unwrap();
        let expected = utc_from_calendar(2024, 6, 15, 3, 22, 45.123)
            .unwrap()
            .to_julian_date();
        assert_eq!(jd, expected);
    }
}
