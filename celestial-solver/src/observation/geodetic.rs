use celestial_core::angle::parse_dms;
use celestial_images::formats::image::Image;

/// Observer location parsed from FITS `OBSGEO-*` or `SITE*` headers.
///
/// Longitude and latitude are decimal degrees, longitude east positive. `SITELONG` is
/// read east positive, but some capture programs write it west positive, so `OBSGEO-L`
/// wins when both are present. Altitude is meters above the ellipsoid or sea level and
/// is optional; `OBSGEO-H` and `SITEELEV` aren't always present.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::observation::geodetic::geodetic_from_image;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// if let Some(geo) = geodetic_from_image(&img) {
///     println!("observer: {:.4}° {:.4}°", geo.lat_deg, geo.lon_deg);
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Geodetic {
    /// Longitude in decimal degrees, east positive (from `OBSGEO-L` or `SITELONG`).
    pub lon_deg: f64,
    /// Latitude in decimal degrees (from `OBSGEO-B` or `SITELAT`).
    pub lat_deg: f64,
    /// Altitude in meters (from `OBSGEO-H` or `SITEELEV`). Optional.
    pub alt_m: Option<f64>,
}

/// Reads observer location from `OBSGEO-B` (lat), `OBSGEO-L` (lon), `OBSGEO-H` (alt),
/// or from `SITELAT`, `SITELONG` and `SITEELEV`.
///
/// `OBSGEO-B` and `OBSGEO-L` take precedence when both are present; otherwise
/// `SITELAT` and `SITELONG` are read, each a number of degrees or a string holding a
/// sexagesimal angle such as `'+39 00 25.56'` or a decimal number. `SITELONG` is read
/// east positive; some capture programs write it west positive, which is why
/// `OBSGEO-L` wins when both are present, and a warning is logged when `SITELONG` is
/// `OBSGEO-L` with the sign flipped. Returns `None` if neither pair is present with
/// finite values. The altitude is `OBSGEO-H`, else `SITEELEV`, else `None`. A key that
/// repeats is read from its last card.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::observation::geodetic::geodetic_from_image;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// let maybe_geo = geodetic_from_image(&img);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn geodetic_from_image(img: &Image) -> Option<Geodetic> {
    let obsgeo = || Some((number(img, "OBSGEO-B")?, number(img, "OBSGEO-L")?));
    let site = || Some((degrees(img, "SITELAT")?, degrees(img, "SITELONG")?));
    let (lat_deg, lon_deg) = match obsgeo() {
        Some(position) => {
            warn_if_sitelong_is_west_positive(img, position.1);
            position
        }
        None => site()?,
    };
    let alt_m = number(img, "OBSGEO-H").or_else(|| number(img, "SITEELEV"));
    Some(Geodetic {
        lon_deg,
        lat_deg,
        alt_m,
    })
}

// The finite numeric value of `key`.
fn number(img: &Image, key: &str) -> Option<f64> {
    finite(key, img.get_last(key).as_f64()?)
}

// Two longitudes within 1 arcmin are taken as the same. A SITELONG and an OBSGEO-L
// written from one site setting differ by the rounding of the written values, at most
// 0.004 arcsec in the files at hand. An east-positive SITELONG is about OBSGEO-L, so it
// differs from -OBSGEO-L by about 2 |lon|; within 1 arcmin of longitude 0 that is too
// small to tell the signs apart, and the check is skipped.
const SAME_LONGITUDE_DEG: f64 = 1.0 / 60.0;

fn warn_if_sitelong_is_west_positive(img: &Image, obsgeo_lon: f64) {
    let Some(site_lon) = degrees(img, "SITELONG") else {
        return;
    };
    if sign_is_flipped(site_lon, obsgeo_lon) {
        log::warn!(
            "SITELONG = {site_lon} is OBSGEO-L = {obsgeo_lon} with the sign flipped; \
             SITELONG is read east positive, so OBSGEO-L is used"
        );
    }
}

fn sign_is_flipped(site_lon: f64, obsgeo_lon: f64) -> bool {
    libm::fabs(site_lon + obsgeo_lon) < SAME_LONGITUDE_DEG
        && libm::fabs(obsgeo_lon) >= SAME_LONGITUDE_DEG
}

// `key` in degrees, from a number or a string. PixInsight reads SITELAT and SITELONG
// the same way, with SITELONG east positive.
fn degrees(img: &Image, key: &str) -> Option<f64> {
    let field = img.get_last(key);
    if let Some(value) = field.as_f64() {
        return finite(key, value);
    }
    let text = field.as_str()?;
    let value = parse_dms(text)
        .map(|angle| angle.degrees())
        .ok()
        .or_else(|| text.trim().parse::<f64>().ok());
    match value {
        Some(value) => finite(key, value),
        None => {
            log::warn!("{key} = {text:?} is not an angle; the observer location skips it");
            None
        }
    }
}

fn finite(key: &str, value: f64) -> Option<f64> {
    if value.is_finite() {
        return Some(value);
    }
    log::warn!("{key} = {value} is not finite; the observer location skips it");
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::keywords::Keyword;
    use celestial_images::formats::pixel_data::PixelData;

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

    // Capture programs write SITELAT and SITELONG as sexagesimal strings or as numbers,
    // with longitude east positive.
    #[test]
    fn site_keywords_give_the_location_without_obsgeo() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("SITELAT", "+39 00 25.56"));
        img.set_keyword(Keyword::string("SITELONG", "-92:20:04.56"));
        img.set_keyword(Keyword::real("SITEELEV", 228.0));
        let geo = geodetic_from_image(&img).unwrap();
        let lat = parse_dms("+39 00 25.56").unwrap().degrees();
        let lon = parse_dms("-92:20:04.56").unwrap().degrees();
        assert_eq!(
            (geo.lat_deg, geo.lon_deg, geo.alt_m),
            (lat, lon, Some(228.0))
        );

        let mut img = stub_image();
        img.set_keyword(Keyword::real("SITELAT", 39.0071));
        img.set_keyword(Keyword::string("SITELONG", "-92.3346"));
        let geo = geodetic_from_image(&img).unwrap();
        assert_eq!(
            (geo.lat_deg, geo.lon_deg, geo.alt_m),
            (39.0071, -92.3346, None)
        );
    }

    #[test]
    fn obsgeo_takes_precedence_over_site_keywords() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBSGEO-B", 39.0071));
        img.set_keyword(Keyword::real("OBSGEO-L", -92.3346));
        img.set_keyword(Keyword::real("OBSGEO-H", 228.0));
        img.set_keyword(Keyword::string("SITELAT", "+39 00 25.56"));
        img.set_keyword(Keyword::string("SITELONG", "+92 20 04.56"));
        img.set_keyword(Keyword::real("SITEELEV", 300.0));
        let geo = geodetic_from_image(&img).unwrap();
        assert_eq!(
            (geo.lat_deg, geo.lon_deg, geo.alt_m),
            (39.0071, -92.3346, Some(228.0))
        );
    }

    #[test]
    fn site_values_that_are_not_angles_are_missing() {
        for (lat, lon) in [("north", "-92.3346"), ("+39 00 25.56", "NaN")] {
            let mut img = stub_image();
            img.set_keyword(Keyword::string("SITELAT", lat));
            img.set_keyword(Keyword::string("SITELONG", lon));
            assert!(geodetic_from_image(&img).is_none(), "{lat} {lon}");
        }
    }

    #[test]
    fn repeated_obsgeo_reads_the_last_card() {
        let mut img = stub_image();
        img.keywords.push(Keyword::real("OBSGEO-B", 10.0));
        img.keywords.push(Keyword::real("OBSGEO-B", 39.0071));
        img.keywords.push(Keyword::real("OBSGEO-L", -92.3346));
        let geo = geodetic_from_image(&img).unwrap();
        assert_eq!((geo.lat_deg, geo.lon_deg), (39.0071, -92.3346));
    }

    // SITELONG beside OBSGEO-L as two files in the repository have them, each with the
    // sign flipped, then a SITELONG that agrees with OBSGEO-L, one 2 arcmin away, and
    // one within 1 arcmin of longitude 0, where the sign cannot be told.
    #[test]
    fn sign_is_flipped_matches_sitelong_against_minus_obsgeo_l() {
        let rosette = parse_dms("+92 20 04.42").unwrap().degrees();
        let bode = parse_dms("+92 20 04.56").unwrap().degrees();
        assert!(sign_is_flipped(rosette, -92.33456));
        assert!(sign_is_flipped(bode, -92.3346));
        assert!(!sign_is_flipped(-92.3346, -92.3346));
        assert!(!sign_is_flipped(92.3346 + 2.0 / 60.0, -92.3346));
        assert!(!sign_is_flipped(0.005, -0.005));
    }
}
