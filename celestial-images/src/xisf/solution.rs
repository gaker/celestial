use crate::fits::header::keywords::{find_last_keyword, Keyword, KeywordValue};
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::property::{XisfProperty, XisfPropertyValue};
use celestial_wcs::builder::{Wcs, WcsBuilder};
use celestial_wcs::header::KeywordProvider;
use celestial_wcs::spherical::Projection;
use std::cell::Cell;

// These ids follow the PCL property names. They have not been checked
// against files written by PixInsight.
pub(crate) const PROJECTION_SYSTEM: &str = "PCL:AstrometricSolution:ProjectionSystem";
pub(crate) const REFERENCE_CELESTIAL: &str =
    "PCL:AstrometricSolution:ReferenceCelestialCoordinates";
pub(crate) const REFERENCE_IMAGE: &str = "PCL:AstrometricSolution:ReferenceImageCoordinates";
pub(crate) const LINEAR_MATRIX: &str = "PCL:AstrometricSolution:LinearTransformationMatrix";
pub(crate) const REFERENCE_NATIVE: &str = "PCL:AstrometricSolution:ReferenceNativeCoordinates";
pub(crate) const CELESTIAL_POLE_NATIVE: &str =
    "PCL:AstrometricSolution:CelestialPoleNativeCoordinates";

pub(crate) const SOLUTION_IDS: [&str; 6] = [
    PROJECTION_SYSTEM,
    REFERENCE_CELESTIAL,
    REFERENCE_IMAGE,
    LINEAR_MATRIX,
    REFERENCE_NATIVE,
    CELESTIAL_POLE_NATIVE,
];

// XISF 1.0 Table 16 names these seven projections and no others.
const PROJECTION_SYSTEMS: [(&str, &str); 7] = [
    ("TAN", "Gnomonic"),
    ("STG", "Stereographic"),
    ("ZEA", "ZenithalEqualArea"),
    ("SIN", "Orthographic"),
    ("CAR", "PlateCarree"),
    ("MER", "Mercator"),
    ("AIT", "HammerAitoff"),
];

pub(crate) fn projection_system_name(code: &str) -> Option<&'static str> {
    PROJECTION_SYSTEMS
        .iter()
        .find(|(c, _)| *c == code)
        .map(|(_, name)| *name)
}

fn projection_code(name: &str) -> Option<&'static str> {
    PROJECTION_SYSTEMS
        .iter()
        .find(|(_, n)| *n == name)
        .map(|(code, _)| *code)
}

// Each Table 16 projection without PV2 parameters, as Paper II defines it.
pub(crate) fn plain_projection(code: &str) -> Option<Projection> {
    match code {
        "TAN" => Some(Projection::tan()),
        "STG" => Some(Projection::stg()),
        "ZEA" => Some(Projection::zea()),
        "SIN" => Some(Projection::sin()),
        "CAR" => Some(Projection::car()),
        "MER" => Some(Projection::mer()),
        "AIT" => Some(Projection::ait()),
        _ => None,
    }
}

// Reads keywords for the WCS builder and records a keyword it finds with a
// value of the wrong type, which the builder alone would treat as absent.
struct KeywordList<'a> {
    keywords: &'a [Keyword],
    mistyped: Cell<bool>,
}

impl<'a> KeywordList<'a> {
    fn new(keywords: &'a [Keyword]) -> Self {
        Self {
            keywords,
            mistyped: Cell::new(false),
        }
    }

    fn typed<T>(&self, key: &str, read: impl Fn(&KeywordValue) -> Option<T>) -> Option<T> {
        let keyword = self.keywords.get(find_last_keyword(self.keywords, key)?)?;
        let value = keyword.value.as_ref().and_then(read);
        if value.is_none() {
            self.mistyped.set(true);
        }
        value
    }
}

impl KeywordProvider for KeywordList<'_> {
    fn get_string(&self, key: &str) -> Option<String> {
        self.typed(key, |v| v.as_string().map(str::to_string))
    }

    // FITS readers type a card such as CRPIX1 = 512 as an integer.
    fn get_float(&self, key: &str) -> Option<f64> {
        self.typed(key, KeywordValue::as_real)
    }

    fn get_int(&self, key: &str) -> Option<i64> {
        self.typed(key, KeywordValue::as_integer)
    }
}

// A celestial axis type has a four-character code after a dash, as in
// "RA---TAN". Linear axes such as "LINEAR" or "PIXEL" have none.
fn has_celestial_ctype(keywords: &[Keyword]) -> bool {
    let list = KeywordList::new(keywords);
    ["CTYPE1", "CTYPE2"].iter().any(|key| {
        list.get_string(key)
            .is_some_and(|ctype| ctype.find('-').is_some_and(|i| i > 0))
    })
}

// Returns None when the keywords describe no celestial WCS, and an error
// when they describe one that does not build. The builder takes a missing
// CRVAL as 0 and a value of the wrong type as missing, so both are errors.
pub(crate) fn celestial_wcs(keywords: &[Keyword]) -> Result<Option<Wcs>> {
    if !has_celestial_ctype(keywords) {
        return Ok(None);
    }
    let list = KeywordList::new(keywords);
    let wcs = WcsBuilder::from_header(&list)
        .and_then(|builder| builder.build())
        .map_err(|e| XisfError::AstrometricSolution(e.to_string()))?;
    if list.mistyped.get() {
        return Err(unusable("a WCS keyword has a value of the wrong type"));
    }
    if ["CRVAL1", "CRVAL2"]
        .iter()
        .any(|key| list.get_float(key).is_none())
    {
        return Err(unusable("CRVAL1 and CRVAL2 are required"));
    }
    Ok(Some(wcs))
}

fn unusable(reason: &str) -> XisfError {
    XisfError::AstrometricSolution(reason.to_string())
}

// Returns the FITS WCS keywords the solution properties describe, or None
// when they describe no usable solution. XISF 1.0 section 11.5.3.7 makes an
// unknown projection name or a missing required property leave the whole
// solution unavailable. Pixel i has its centre at i + 0.5 in XISF image
// coordinates and at i + 1 in FITS, so CRPIX is the reference plus 0.5.
pub(crate) fn solution_keywords(properties: &[XisfProperty]) -> Option<Vec<Keyword>> {
    let code = projection_code(string_property(properties, PROJECTION_SYSTEM)?)?;
    let [alpha_0, delta_0] = pair(properties, REFERENCE_CELESTIAL)?;
    let [x_ref, y_ref] = pair(properties, REFERENCE_IMAGE)?;
    let cd = matrix_2x2(properties, LINEAR_MATRIX)?;
    let mut keywords = vec![
        Keyword::string("CTYPE1", format!("RA---{code}").as_str()),
        Keyword::string("CTYPE2", format!("DEC--{code}").as_str()),
        Keyword::real("CRPIX1", x_ref + 0.5),
        Keyword::real("CRPIX2", y_ref + 0.5),
        Keyword::real("CRVAL1", alpha_0),
        Keyword::real("CRVAL2", delta_0),
        Keyword::real("CD1_1", cd[0]),
        Keyword::real("CD1_2", cd[1]),
        Keyword::real("CD2_1", cd[2]),
        Keyword::real("CD2_2", cd[3]),
    ];
    keywords.extend(native_keywords(properties, code)?);
    celestial_wcs(&keywords).ok().flatten()?;
    Some(keywords)
}

// PV1_1 and PV1_2 are written only when the reference point is not the
// projection's default, and LONPOLE and LATPOLE only when the pole is given.
// A property that is present but malformed makes the solution unusable.
fn native_keywords(properties: &[XisfProperty], code: &str) -> Option<Vec<Keyword>> {
    let mut keywords = Vec::new();
    if let Some(property) = find(properties, REFERENCE_NATIVE) {
        let [phi_0, theta_0] = as_pair(property)?;
        if (phi_0, theta_0) != plain_projection(code)?.native_reference() {
            keywords.push(Keyword::real("PV1_1", phi_0));
            keywords.push(Keyword::real("PV1_2", theta_0));
        }
    }
    if let Some(property) = find(properties, CELESTIAL_POLE_NATIVE) {
        let [phi_p, theta_p] = as_pair(property)?;
        keywords.push(Keyword::real("LONPOLE", phi_p));
        keywords.push(Keyword::real("LATPOLE", theta_p));
    }
    Some(keywords)
}

fn find<'a>(properties: &'a [XisfProperty], id: &str) -> Option<&'a XisfProperty> {
    properties.iter().find(|p| p.id == id)
}

fn string_property<'a>(properties: &'a [XisfProperty], id: &str) -> Option<&'a str> {
    match &find(properties, id)?.value {
        XisfPropertyValue::String(text) => Some(text.trim()),
        _ => None,
    }
}

fn pair(properties: &[XisfProperty], id: &str) -> Option<[f64; 2]> {
    as_pair(find(properties, id)?)
}

fn as_pair(property: &XisfProperty) -> Option<[f64; 2]> {
    match &property.value {
        XisfPropertyValue::F64Vector(values) => values.as_slice().try_into().ok(),
        _ => None,
    }
}

fn matrix_2x2(properties: &[XisfProperty], id: &str) -> Option<[f64; 4]> {
    match &find(properties, id)?.value {
        XisfPropertyValue::F64Matrix {
            rows: 2,
            cols: 2,
            data,
        } => data.as_slice().try_into().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gnomonic() -> Vec<XisfProperty> {
        vec![
            XisfProperty::string(PROJECTION_SYSTEM, "Gnomonic"),
            XisfProperty::f64_vector(REFERENCE_CELESTIAL, vec![37.5, 30.0]),
            XisfProperty::f64_vector(REFERENCE_IMAGE, vec![511.5, 255.25]),
            XisfProperty::f64_matrix(LINEAR_MATRIX, 2, 2, vec![-1e-3, 2e-6, 3e-6, 1e-3]),
        ]
    }

    fn names(keywords: &[Keyword]) -> Vec<&str> {
        keywords.iter().map(|k| k.name.as_str()).collect()
    }

    fn replace(properties: &mut Vec<XisfProperty>, property: XisfProperty) {
        properties.retain(|p| p.id != property.id);
        properties.push(property);
    }

    #[test]
    fn required_properties_give_the_linear_keywords() {
        let keywords = solution_keywords(&gnomonic()).unwrap();
        assert_eq!(
            names(&keywords),
            [
                "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2",
                "CD2_1", "CD2_2"
            ]
        );
        assert_eq!(keywords[2].value, Some(KeywordValue::Real(512.0)));
        assert_eq!(keywords[3].value, Some(KeywordValue::Real(255.75)));
    }

    #[test]
    fn default_native_reference_adds_no_pv_and_the_pole_adds_lonpole_latpole() {
        let mut properties = gnomonic();
        properties.push(XisfProperty::f64_vector(REFERENCE_NATIVE, vec![0.0, 90.0]));
        properties.push(XisfProperty::f64_vector(
            CELESTIAL_POLE_NATIVE,
            vec![180.0, 30.0],
        ));
        let keywords = solution_keywords(&properties).unwrap();
        assert_eq!(&names(&keywords)[10..], ["LONPOLE", "LATPOLE"]);
        assert_eq!(keywords[10].value, Some(KeywordValue::Real(180.0)));
        assert_eq!(keywords[11].value, Some(KeywordValue::Real(30.0)));
    }

    #[test]
    fn other_native_reference_adds_pv1_1_and_pv1_2() {
        let mut properties = gnomonic();
        replace(
            &mut properties,
            XisfProperty::string(PROJECTION_SYSTEM, "PlateCarree"),
        );
        properties.push(XisfProperty::f64_vector(REFERENCE_NATIVE, vec![5.0, 20.0]));
        let keywords = solution_keywords(&properties).unwrap();
        assert_eq!(
            keywords[0].value,
            Some(KeywordValue::String("RA---CAR".into()))
        );
        assert_eq!(&names(&keywords)[10..], ["PV1_1", "PV1_2"]);
        assert_eq!(keywords[11].value, Some(KeywordValue::Real(20.0)));
    }

    #[test]
    fn unusable_solutions_give_no_keywords() {
        let unknown_name = XisfProperty::string(PROJECTION_SYSTEM, "ZenithalEquidistant");
        let short_vector = XisfProperty::f64_vector(REFERENCE_CELESTIAL, vec![37.5]);
        let wrong_shape = XisfProperty::f64_matrix(LINEAR_MATRIX, 1, 4, vec![1.0, 0.0, 0.0, 1.0]);
        let singular = XisfProperty::f64_matrix(LINEAR_MATRIX, 2, 2, vec![0.0; 4]);
        let bad_pole = XisfProperty::f64_vector(CELESTIAL_POLE_NATIVE, vec![1.0, 2.0, 3.0]);
        let bad_type = XisfProperty::float64(REFERENCE_NATIVE, 90.0);
        for property in [
            unknown_name,
            short_vector,
            wrong_shape,
            singular,
            bad_pole,
            bad_type,
        ] {
            let mut properties = gnomonic();
            replace(&mut properties, property.clone());
            assert_eq!(solution_keywords(&properties), None, "{property:?}");
        }
        for missing in [
            PROJECTION_SYSTEM,
            REFERENCE_CELESTIAL,
            REFERENCE_IMAGE,
            LINEAR_MATRIX,
        ] {
            let mut properties = gnomonic();
            properties.retain(|p| p.id != missing);
            assert_eq!(solution_keywords(&properties), None, "{missing}");
        }
    }
}
