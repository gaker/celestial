use celestial_images::formats::Image;

/// Reads the raw `DATE-OBS` string from image headers.
///
/// Returns the header value verbatim — no parsing, no normalization. Returns `None` if
/// `DATE-OBS` is missing or non-string. Use [`crate::metadata::jd_from_image`] when
/// you need a parsed [`celestial_time::JulianDate`].
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::observation::date_obs_from_image;
///
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// let label = date_obs_from_image(&img).unwrap_or_default();
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn date_obs_from_image(img: &Image) -> Option<String> {
    img.get_keyword("DATE-OBS")
        .and_then(|kw| kw.value.as_ref()?.as_string().map(|s| s.to_owned()))
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
    fn returns_string_value_verbatim() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15T03:22:45.123"));
        assert_eq!(
            date_obs_from_image(&img).as_deref(),
            Some("2024-06-15T03:22:45.123")
        );
    }

    #[test]
    fn missing_header_returns_none() {
        let img = stub_image();
        assert!(date_obs_from_image(&img).is_none());
    }

    #[test]
    fn non_string_value_returns_none() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("DATE-OBS", 60476.14));
        assert!(date_obs_from_image(&img).is_none());
    }
}
