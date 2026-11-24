use std::cell::Cell;

use celestial_images::fits::header::parser::Field;
use celestial_images::formats::image::Image;
use celestial_wcs::builder::WcsBuilder;
use celestial_wcs::header::KeywordProvider;

use crate::pixel_frame::detection_to_sky;

// What an image's CTYPE cards and the WCS they describe say about the position hint.
pub(super) enum HeaderWcs {
    // No CTYPE1 or CTYPE2 card.
    Absent,
    // CTYPE names axes other than right ascension and declination, so CRVAL is not a
    // sky position.
    OtherAxes,
    // The sky position, (ra, dec) in degrees, of the image centre.
    Centre(f64, f64),
    // CTYPE names right ascension and declination, but the WCS has no sky position for
    // the image centre, for the reason given.
    Unusable(String),
}

// Reads the header WCS and evaluates it at the image centre. Every card is the last one
// of its name, as cfitsio and WCSLIB read a repeated key.
pub(super) fn header_wcs(img: &Image) -> HeaderWcs {
    let ctypes = ["CTYPE1", "CTYPE2"].map(|key| img.get_last(key));
    if ctypes.iter().all(|ctype| !ctype.exists()) {
        return HeaderWcs::Absent;
    }
    if !ctypes.iter().all(names_equatorial_axis) {
        return HeaderWcs::OtherAxes;
    }
    match centre(img) {
        Ok((ra, dec)) => HeaderWcs::Centre(ra, dec),
        Err(reason) => HeaderWcs::Unusable(reason),
    }
}

// An absent card names no axis, so it does not rule out an equatorial pair; the WCS
// builder then reports the missing card.
fn names_equatorial_axis(ctype: &Field<'_>) -> bool {
    match ctype.as_str() {
        Some(text) => {
            let axis = text.trim().split('-').next().unwrap_or_default();
            axis == "RA" || axis == "DEC"
        }
        None => !ctype.exists(),
    }
}

// The builder takes an absent CRPIX or CRVAL as 0 and an absent CDELT as 1 deg, which
// would put the hint anywhere, so the cards that place the image must be present.
fn centre(img: &Image) -> Result<(f64, f64), String> {
    if let Some(key) = missing_placement_card(img) {
        return Err(format!("no {key}"));
    }
    let keywords = ImageKeywords::new(img);
    let wcs = WcsBuilder::from_header(&keywords)
        .and_then(WcsBuilder::build)
        .map_err(|e| e.to_string())?;
    if keywords.mistyped.get() {
        return Err("a WCS card has a value of the wrong type".to_string());
    }
    // The image centre in detection pixels, FITS pixel ((w + 1) / 2, (h + 1) / 2).
    let (x, y) = (
        (img.width() as f64 - 1.0) / 2.0,
        (img.height() as f64 - 1.0) / 2.0,
    );
    let (ra, dec) = detection_to_sky(&wcs, x, y).map_err(|e| e.to_string())?;
    if !(ra.is_finite() && dec.is_finite()) {
        return Err(format!("the image centre maps to ({ra}, {dec})"));
    }
    Ok((ra, dec))
}

fn missing_placement_card(img: &Image) -> Option<&'static str> {
    let present = |key: &str| img.get_last(key).exists();
    let reference = ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2"];
    if let Some(key) = reference.into_iter().find(|key| !present(key)) {
        return Some(key);
    }
    let cd = ["CD1_1", "CD1_2", "CD2_1", "CD2_2"]
        .into_iter()
        .any(present);
    if !cd && !(present("CDELT1") && present("CDELT2")) {
        return Some("CD matrix or CDELT1/2");
    }
    None
}

// The image's cards for the WCS builder, which treats a card with a value of the wrong
// type as absent; `mistyped` records that one was read.
struct ImageKeywords<'a> {
    img: &'a Image,
    mistyped: Cell<bool>,
}

impl<'a> ImageKeywords<'a> {
    fn new(img: &'a Image) -> Self {
        Self {
            img,
            mistyped: Cell::new(false),
        }
    }

    fn typed<T>(&self, key: &str, read: impl Fn(&Field<'a>) -> Option<T>) -> Option<T> {
        let field = self.img.get_last(key);
        let value = read(&field);
        if field.exists() && value.is_none() {
            self.mistyped.set(true);
        }
        value
    }
}

impl KeywordProvider for ImageKeywords<'_> {
    fn get_string(&self, key: &str) -> Option<String> {
        self.typed(key, |field| field.as_str().map(str::to_string))
    }

    fn get_float(&self, key: &str) -> Option<f64> {
        self.typed(key, Field::as_f64)
    }

    fn get_int(&self, key: &str) -> Option<i64> {
        self.typed(key, Field::as_i64)
    }
}
