use crate::xisf::header::XisfProperty;

pub fn wcs_to_xisf_properties(wcs_keywords: &[celestial_wcs::WcsKeyword]) -> Vec<XisfProperty> {
    let mut properties = Vec::new();
    let mut crval1: Option<f64> = None;
    let mut crval2: Option<f64> = None;
    let mut crpix1: Option<f64> = None;
    let mut crpix2: Option<f64> = None;
    let mut cd1_1: Option<f64> = None;
    let mut cd1_2: Option<f64> = None;
    let mut cd2_1: Option<f64> = None;
    let mut cd2_2: Option<f64> = None;
    let mut proj_code: Option<String> = None;
    let mut lonpole: Option<f64> = None;
    let mut latpole: Option<f64> = None;

    for kw in wcs_keywords {
        match kw.name.as_str() {
            "CTYPE1" => {
                if let celestial_wcs::WcsKeywordValue::String(s) = &kw.value {
                    proj_code = extract_projection_code(s);
                }
            }
            "CRVAL1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crval1 = Some(*v);
                }
            }
            "CRVAL2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crval2 = Some(*v);
                }
            }
            "CRPIX1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crpix1 = Some(*v);
                }
            }
            "CRPIX2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crpix2 = Some(*v);
                }
            }
            "CD1_1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd1_1 = Some(*v);
                }
            }
            "CD1_2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd1_2 = Some(*v);
                }
            }
            "CD2_1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd2_1 = Some(*v);
                }
            }
            "CD2_2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd2_2 = Some(*v);
                }
            }
            "LONPOLE" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    lonpole = Some(*v);
                }
            }
            "LATPOLE" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    latpole = Some(*v);
                }
            }
            _ => {}
        }
    }

    if let Some(code) = proj_code {
        properties.push(XisfProperty::string(
            "PCL:AstrometricSolution:ProjectionSystem",
            code,
        ));
    }

    if let (Some(ra), Some(dec)) = (crval1, crval2) {
        properties.push(XisfProperty::f64_vector(
            "PCL:AstrometricSolution:ReferenceCelestialCoordinates",
            vec![ra, dec],
        ));
    }

    if let (Some(x), Some(y)) = (crpix1, crpix2) {
        properties.push(XisfProperty::f64_vector(
            "PCL:AstrometricSolution:ReferenceImageCoordinates",
            vec![x, y],
        ));
    }

    if let (Some(c11), Some(c12), Some(c21), Some(c22)) = (cd1_1, cd1_2, cd2_1, cd2_2) {
        properties.push(XisfProperty::f64_matrix(
            "PCL:AstrometricSolution:LinearTransformationMatrix",
            2,
            2,
            vec![c11, c12, c21, c22],
        ));
    }

    properties.push(XisfProperty::f64_vector(
        "PCL:AstrometricSolution:ReferenceNativeCoordinates",
        vec![0.0, 90.0],
    ));

    let phi_p = lonpole.unwrap_or_else(|| {
        if crval2.unwrap_or(0.0) >= 90.0 {
            0.0
        } else {
            180.0
        }
    });
    let theta_p = latpole.unwrap_or_else(|| crval2.unwrap_or(0.0));
    properties.push(XisfProperty::f64_vector(
        "PCL:AstrometricSolution:CelestialPoleNativeCoordinates",
        vec![phi_p, theta_p],
    ));

    properties
}

fn extract_projection_code(ctype: &str) -> Option<String> {
    let trimmed = ctype.trim();
    if trimmed.len() < 8 {
        return None;
    }
    let dashes_pos = trimmed.find('-')?;
    let after_dashes = &trimmed[dashes_pos..];
    let code_start = after_dashes.trim_start_matches('-');
    if code_start.is_empty() {
        return None;
    }
    Some(projection_system_name(code_start))
}

fn projection_system_name(code: &str) -> String {
    match code {
        "TAN" => "Gnomonic",
        "SIN" => "Orthographic",
        "STG" => "Stereographic",
        "ARC" => "ZenithalEquidistant",
        "ZEA" => "ZenithalEqualArea",
        "AZP" => "ZenithalPerspective",
        "SZP" => "SlantZenithalPerspective",
        "ZPN" => "ZenithalPolynomial",
        "AIR" => "Airy",
        "CYP" => "CylindricalPerspective",
        "CEA" => "CylindricalEqualArea",
        "CAR" => "PlateCarree",
        "MER" => "Mercator",
        "SFL" => "SansonFlamsteed",
        "PAR" => "Parabolic",
        "MOL" => "Mollweide",
        "AIT" => "HammerAitoff",
        "COP" => "ConicPerspective",
        "COE" => "ConicEqualArea",
        "COD" => "ConicEquidistant",
        "COO" => "ConicOrthomorphic",
        other => return other.to_string(),
    }
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xisf::header::XisfPropertyValue;
    use celestial_wcs::{WcsKeyword, WcsKeywordValue};

    fn find_prop<'a>(props: &'a [XisfProperty], id: &str) -> Option<&'a XisfProperty> {
        props.iter().find(|p| p.id == id)
    }

    #[test]
    fn extract_projection_code_rejects_short_ctype() {
        assert_eq!(extract_projection_code(""), None);
        assert_eq!(extract_projection_code("RA--TAN"), None);
        assert_eq!(extract_projection_code("short"), None);
    }

    #[test]
    fn extract_projection_code_rejects_ctype_without_dashes() {
        assert_eq!(extract_projection_code("RA_____"), None);
        assert_eq!(extract_projection_code("12345678"), None);
    }

    #[test]
    fn extract_projection_code_rejects_dashes_only_trailing() {
        assert_eq!(extract_projection_code("RA------"), None);
    }

    #[test]
    fn extract_projection_code_finds_known_code_after_dashes() {
        assert_eq!(extract_projection_code("RA---TAN"), Some("Gnomonic".to_string()));
        assert_eq!(extract_projection_code("DEC--SIN"), Some("Orthographic".to_string()));
    }

    #[test]
    fn extract_projection_code_returns_unknown_code_unchanged() {
        assert_eq!(extract_projection_code("FOO---XYZ"), Some("XYZ".to_string()));
    }

    #[test]
    fn extract_projection_code_trims_whitespace() {
        assert_eq!(
            extract_projection_code("  RA---TAN  "),
            Some("Gnomonic".to_string())
        );
    }

    #[test]
    fn projection_system_name_maps_every_known_code() {
        let cases = [
            ("TAN", "Gnomonic"),
            ("SIN", "Orthographic"),
            ("STG", "Stereographic"),
            ("ARC", "ZenithalEquidistant"),
            ("ZEA", "ZenithalEqualArea"),
            ("AZP", "ZenithalPerspective"),
            ("SZP", "SlantZenithalPerspective"),
            ("ZPN", "ZenithalPolynomial"),
            ("AIR", "Airy"),
            ("CYP", "CylindricalPerspective"),
            ("CEA", "CylindricalEqualArea"),
            ("CAR", "PlateCarree"),
            ("MER", "Mercator"),
            ("SFL", "SansonFlamsteed"),
            ("PAR", "Parabolic"),
            ("MOL", "Mollweide"),
            ("AIT", "HammerAitoff"),
            ("COP", "ConicPerspective"),
            ("COE", "ConicEqualArea"),
            ("COD", "ConicEquidistant"),
            ("COO", "ConicOrthomorphic"),
        ];
        for (code, expected) in cases {
            assert_eq!(projection_system_name(code), expected);
        }
    }

    #[test]
    fn projection_system_name_echoes_unknown_code() {
        assert_eq!(projection_system_name("XYZ"), "XYZ");
        assert_eq!(projection_system_name(""), "");
        assert_eq!(projection_system_name("BON"), "BON");
    }

    #[test]
    fn empty_wcs_input_emits_only_constant_properties() {
        let props = wcs_to_xisf_properties(&[]);
        assert_eq!(props.len(), 2);

        let native = find_prop(&props, "PCL:AstrometricSolution:ReferenceNativeCoordinates").unwrap();
        assert_eq!(
            native.value,
            XisfPropertyValue::F64Vector(vec![0.0, 90.0])
        );

        let pole = find_prop(&props, "PCL:AstrometricSolution:CelestialPoleNativeCoordinates")
            .unwrap();
        assert_eq!(pole.value, XisfPropertyValue::F64Vector(vec![180.0, 0.0]));
    }

    #[test]
    fn partial_input_emits_only_available_sections() {
        let kws = vec![
            WcsKeyword::real("CRVAL1", 150.0),
            WcsKeyword::real("CRVAL2", 30.0),
        ];
        let props = wcs_to_xisf_properties(&kws);
        let celestial =
            find_prop(&props, "PCL:AstrometricSolution:ReferenceCelestialCoordinates").unwrap();
        assert_eq!(
            celestial.value,
            XisfPropertyValue::F64Vector(vec![150.0, 30.0])
        );
        assert!(find_prop(&props, "PCL:AstrometricSolution:ReferenceImageCoordinates").is_none());
        assert!(find_prop(&props, "PCL:AstrometricSolution:LinearTransformationMatrix").is_none());
        assert!(find_prop(&props, "PCL:AstrometricSolution:ProjectionSystem").is_none());
    }

    #[test]
    fn full_input_emits_every_section() {
        let kws = vec![
            WcsKeyword {
                name: "CTYPE1".to_string(),
                value: WcsKeywordValue::String("RA---TAN".to_string()),
            },
            WcsKeyword::real("CRVAL1", 10.0),
            WcsKeyword::real("CRVAL2", 20.0),
            WcsKeyword::real("CRPIX1", 512.0),
            WcsKeyword::real("CRPIX2", 512.0),
            WcsKeyword::real("CD1_1", 0.001),
            WcsKeyword::real("CD1_2", 0.0),
            WcsKeyword::real("CD2_1", 0.0),
            WcsKeyword::real("CD2_2", 0.001),
        ];
        let props = wcs_to_xisf_properties(&kws);

        let proj = find_prop(&props, "PCL:AstrometricSolution:ProjectionSystem").unwrap();
        assert_eq!(proj.value, XisfPropertyValue::String("Gnomonic".to_string()));

        let celestial =
            find_prop(&props, "PCL:AstrometricSolution:ReferenceCelestialCoordinates").unwrap();
        assert_eq!(
            celestial.value,
            XisfPropertyValue::F64Vector(vec![10.0, 20.0])
        );

        let image =
            find_prop(&props, "PCL:AstrometricSolution:ReferenceImageCoordinates").unwrap();
        assert_eq!(
            image.value,
            XisfPropertyValue::F64Vector(vec![512.0, 512.0])
        );

        let matrix = find_prop(&props, "PCL:AstrometricSolution:LinearTransformationMatrix")
            .unwrap();
        match &matrix.value {
            XisfPropertyValue::F64Matrix { rows, cols, data } => {
                assert_eq!(*rows, 2);
                assert_eq!(*cols, 2);
                assert_eq!(data, &vec![0.001, 0.0, 0.0, 0.001]);
            }
            _ => panic!("expected matrix"),
        }

        let pole = find_prop(&props, "PCL:AstrometricSolution:CelestialPoleNativeCoordinates")
            .unwrap();
        // CRVAL2=20 < 90 → phi_p=180, theta_p defaults to CRVAL2=20
        assert_eq!(
            pole.value,
            XisfPropertyValue::F64Vector(vec![180.0, 20.0])
        );
    }

    #[test]
    fn phi_p_defaults_to_zero_when_crval2_is_at_pole() {
        let kws = vec![
            WcsKeyword::real("CRVAL1", 0.0),
            WcsKeyword::real("CRVAL2", 90.0),
        ];
        let props = wcs_to_xisf_properties(&kws);
        let pole = find_prop(&props, "PCL:AstrometricSolution:CelestialPoleNativeCoordinates")
            .unwrap();
        assert_eq!(pole.value, XisfPropertyValue::F64Vector(vec![0.0, 90.0]));
    }

    #[test]
    fn phi_p_defaults_to_180_below_pole() {
        let kws = vec![WcsKeyword::real("CRVAL2", -30.0)];
        let props = wcs_to_xisf_properties(&kws);
        let pole = find_prop(&props, "PCL:AstrometricSolution:CelestialPoleNativeCoordinates")
            .unwrap();
        assert_eq!(
            pole.value,
            XisfPropertyValue::F64Vector(vec![180.0, -30.0])
        );
    }

    #[test]
    fn explicit_lonpole_and_latpole_override_defaults() {
        let kws = vec![
            WcsKeyword::real("CRVAL2", 30.0),
            WcsKeyword::real("LONPOLE", 45.0),
            WcsKeyword::real("LATPOLE", 75.0),
        ];
        let props = wcs_to_xisf_properties(&kws);
        let pole = find_prop(&props, "PCL:AstrometricSolution:CelestialPoleNativeCoordinates")
            .unwrap();
        assert_eq!(pole.value, XisfPropertyValue::F64Vector(vec![45.0, 75.0]));
    }

    #[test]
    fn non_real_valued_crval_is_ignored() {
        let kws = vec![
            WcsKeyword {
                name: "CRVAL1".to_string(),
                value: WcsKeywordValue::Integer(5),
            },
            WcsKeyword::real("CRVAL2", 10.0),
        ];
        let props = wcs_to_xisf_properties(&kws);
        assert!(
            find_prop(&props, "PCL:AstrometricSolution:ReferenceCelestialCoordinates").is_none()
        );
    }

    #[test]
    fn ctype1_with_non_string_value_is_ignored() {
        let kws = vec![WcsKeyword {
            name: "CTYPE1".to_string(),
            value: WcsKeywordValue::Real(1.0),
        }];
        let props = wcs_to_xisf_properties(&kws);
        assert!(find_prop(&props, "PCL:AstrometricSolution:ProjectionSystem").is_none());
    }

    #[test]
    fn unknown_wcs_keywords_are_dropped_silently() {
        let kws = vec![
            WcsKeyword::real("UNKNOWN", 42.0),
            WcsKeyword::real("CRVAL1", 10.0),
            WcsKeyword::real("CRVAL2", 20.0),
        ];
        let props = wcs_to_xisf_properties(&kws);
        let celestial =
            find_prop(&props, "PCL:AstrometricSolution:ReferenceCelestialCoordinates").unwrap();
        assert_eq!(
            celestial.value,
            XisfPropertyValue::F64Vector(vec![10.0, 20.0])
        );
    }
}
