use celestial_images::formats::Image;

/// Observer location parsed from FITS `OBSGEO-*` headers.
///
/// Longitude and latitude are decimal degrees. Altitude is meters above the ellipsoid
/// and is optional — `OBSGEO-H` isn't always present.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::observation::geodetic_from_image;
///
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// if let Some(geo) = geodetic_from_image(&img) {
///     println!("observer: {:.4}° {:.4}°", geo.lat_deg, geo.lon_deg);
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Geodetic {
    /// Longitude in decimal degrees, east positive (from `OBSGEO-L`).
    pub lon_deg: f64,
    /// Latitude in decimal degrees (from `OBSGEO-B`).
    pub lat_deg: f64,
    /// Altitude in meters above ellipsoid (from `OBSGEO-H`). Optional.
    pub alt_m: Option<f64>,
}

/// Reads observer location from `OBSGEO-B` (lat), `OBSGEO-L` (lon), `OBSGEO-H` (alt).
///
/// Returns `None` if either of the required headers (`OBSGEO-B`, `OBSGEO-L`) is missing
/// or non-numeric. `OBSGEO-H` is optional — if missing, `alt_m` is `None` in the result.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::observation::geodetic_from_image;
///
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// let maybe_geo = geodetic_from_image(&img);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn geodetic_from_image(img: &Image) -> Option<Geodetic> {
    let lat_deg = img
        .get_keyword("OBSGEO-B")
        .and_then(|kw| kw.value.as_ref()?.as_real())?;
    let lon_deg = img
        .get_keyword("OBSGEO-L")
        .and_then(|kw| kw.value.as_ref()?.as_real())?;
    let alt_m = img
        .get_keyword("OBSGEO-H")
        .and_then(|kw| kw.value.as_ref()?.as_real());
    Some(Geodetic {
        lon_deg,
        lat_deg,
        alt_m,
    })
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
    fn all_three_headers_populate_all_fields() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBSGEO-B", 35.75));
        img.set_keyword(Keyword::real("OBSGEO-L", -95.25));
        img.set_keyword(Keyword::real("OBSGEO-H", 350.0));
        let geo = geodetic_from_image(&img).unwrap();
        assert_eq!(geo.lat_deg, 35.75);
        assert_eq!(geo.lon_deg, -95.25);
        assert_eq!(geo.alt_m, Some(350.0));
    }

    #[test]
    fn altitude_is_optional() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBSGEO-B", 0.0));
        img.set_keyword(Keyword::real("OBSGEO-L", 0.0));
        let geo = geodetic_from_image(&img).unwrap();
        assert_eq!(geo.alt_m, None);
    }

    #[test]
    fn missing_latitude_returns_none() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBSGEO-L", 0.0));
        img.set_keyword(Keyword::real("OBSGEO-H", 0.0));
        assert!(geodetic_from_image(&img).is_none());
    }

    #[test]
    fn missing_longitude_returns_none() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBSGEO-B", 0.0));
        img.set_keyword(Keyword::real("OBSGEO-H", 0.0));
        assert!(geodetic_from_image(&img).is_none());
    }

    #[test]
    fn no_geodetic_headers_returns_none() {
        let img = stub_image();
        assert!(geodetic_from_image(&img).is_none());
    }

    #[test]
    fn non_numeric_header_value_is_treated_as_missing() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("OBSGEO-B", "bogus"));
        img.set_keyword(Keyword::real("OBSGEO-L", 0.0));
        assert!(geodetic_from_image(&img).is_none());
    }
}
