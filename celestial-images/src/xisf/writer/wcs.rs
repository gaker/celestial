use crate::fits::header::keywords::{find_last_keyword, Keyword};
use crate::xisf::header::property::XisfProperty;
use crate::xisf::solution::{
    celestial_wcs, plain_projection, projection_system_name, CELESTIAL_POLE_NATIVE, LINEAR_MATRIX,
    PROJECTION_SYSTEM, REFERENCE_CELESTIAL, REFERENCE_IMAGE, REFERENCE_NATIVE,
};
use celestial_wcs::builder::{CoordType, Wcs};

// Returns the astrometric solution properties for the celestial WCS the
// keywords describe. A WCS that does not build, or that the properties
// cannot express, gets none, and the FITSKeyword elements still hold it.
pub(crate) fn solution_properties(keywords: &[Keyword]) -> Vec<XisfProperty> {
    let Ok(Some(wcs)) = celestial_wcs(keywords) else {
        return Vec::new();
    };
    if !is_icrs_compatible(keywords) {
        return Vec::new();
    }
    wcs_to_xisf_properties(&wcs).unwrap_or_default()
}

// XISF 1.0 section 11.5.3.7.2 refers a solution without a reference system
// property to the ICRS. FK5 at J2000 differs from it by tens of
// milliarcseconds and is accepted; other systems would be mislabelled.
// Paper II reads a header without RADESYS as FK4 when EQUINOX is below
// 1984 and FK5 otherwise. EPOCH is the older name for EQUINOX.
fn is_icrs_compatible(keywords: &[Keyword]) -> bool {
    let value = |name: &str| {
        keywords
            .get(find_last_keyword(keywords, name)?)?
            .value
            .as_ref()
    };
    let equinox = match value("EQUINOX").or_else(|| value("EPOCH")) {
        None => None,
        Some(v) => match v.as_real() {
            Some(year) => Some(year),
            None => return false,
        },
    };
    let system = match value("RADESYS") {
        None => None,
        Some(v) => match v.as_string() {
            Some(text) => Some(text.trim()),
            None => return false,
        },
    };
    match (system, equinox) {
        (Some("ICRS"), _) | (Some("FK5") | None, None) => true,
        (Some("FK5") | None, Some(year)) => year == 2000.0,
        _ => false,
    }
}

// ReferenceImageCoordinates is CRPIX - 0.5 because XISF puts the centre of
// pixel i at i + 0.5 and FITS puts it at i + 1. A SIP or TPV distortion is
// left out: the linear part is the solution's first layer, and the FITS
// keywords keep the distortion terms.
fn wcs_to_xisf_properties(wcs: &Wcs) -> Option<Vec<XisfProperty>> {
    let name = projection_system(wcs)?;
    let [crpix1, crpix2] = wcs.crpix();
    let cd = wcs.cd_matrix();
    let (phi_0, theta_0) = wcs.native_reference();
    let pole = wcs.rotation();
    Some(vec![
        XisfProperty::string(PROJECTION_SYSTEM, name),
        XisfProperty::f64_vector(REFERENCE_CELESTIAL, wcs.crval().to_vec()),
        XisfProperty::f64_vector(REFERENCE_IMAGE, vec![crpix1 - 0.5, crpix2 - 0.5]),
        XisfProperty::f64_matrix(LINEAR_MATRIX, 2, 2, cd.as_flattened().to_vec()),
        XisfProperty::f64_vector(REFERENCE_NATIVE, vec![phi_0, theta_0]),
        XisfProperty::f64_vector(
            CELESTIAL_POLE_NATIVE,
            vec![pole.phi_p_deg(), pole.delta_p_deg()],
        ),
    ])
}

// The solution properties hold equatorial coordinates, one of the Table 16
// projections without PV2 parameters, and no PV1_0 offset of the fiducial
// point. Any other WCS has no projection system name.
fn projection_system(wcs: &Wcs) -> Option<&'static str> {
    let code = wcs.projection_code();
    let equatorial = wcs.coord_type() == CoordType::Equatorial;
    let offset = wcs.to_keywords().iter().any(|k| k.name == "PV1_0");
    let plain = plain_projection(code).as_ref() == Some(wcs.projection());
    if !equatorial || offset || !plain {
        return None;
    }
    projection_system_name(code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xisf::header::property::XisfPropertyValue;
    use celestial_wcs::builder::WcsBuilder;
    use celestial_wcs::spherical::Projection;

    fn find_prop<'a>(props: &'a [XisfProperty], id: &str) -> Option<&'a XisfProperty> {
        props.iter().find(|p| p.id == id)
    }

    fn vector(props: &[XisfProperty], id: &str) -> Option<Vec<f64>> {
        match &find_prop(props, id)?.value {
            XisfPropertyValue::F64Vector(v) => Some(v.clone()),
            _ => None,
        }
    }

    fn built(projection: Projection) -> Wcs {
        WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(37.5, 30.0)
            .cd_matrix([[-1e-3, 0.0], [0.0, 1e-3]])
            .projection(projection)
            .build()
            .unwrap()
    }

    // A TAN header with CRVAL (10, 20), CRPIX (512, 512), CD 0.001 I.
    fn tan_keywords() -> Vec<Keyword> {
        vec![
            Keyword::string("CTYPE1", "RA---TAN"),
            Keyword::string("CTYPE2", "DEC--TAN"),
            Keyword::real("CRVAL1", 10.0),
            Keyword::real("CRVAL2", 20.0),
            Keyword::real("CRPIX1", 512.0),
            Keyword::real("CRPIX2", 512.0),
            Keyword::real("CD1_1", 0.001),
            Keyword::real("CD1_2", 0.0),
            Keyword::real("CD2_1", 0.0),
            Keyword::real("CD2_2", 0.001),
        ]
    }

    fn set(keywords: &mut Vec<Keyword>, keyword: Keyword) {
        keywords.retain(|k| k.name != keyword.name);
        keywords.push(keyword);
    }

    // Paper II: CAR's fiducial point is (0, 0); with δ0 = 30 < θ0 = 0 is false,
    // so φp = φ0 = 0, and δp = 90 - 30 = 60.
    #[test]
    fn car_solution_gets_cylindrical_native_reference_and_pole() {
        let props = wcs_to_xisf_properties(&built(Projection::car())).unwrap();
        assert_eq!(vector(&props, REFERENCE_NATIVE), Some(vec![0.0, 0.0]));
        assert_eq!(vector(&props, CELESTIAL_POLE_NATIVE), Some(vec![0.0, 60.0]));
    }

    // XISF 1.0 Table 16 names seven projection systems; COP is not one.
    #[test]
    fn cop_solution_writes_no_projection_name_outside_table_16() {
        assert_eq!(
            wcs_to_xisf_properties(&built(Projection::cop(45.0, 25.0))),
            None
        );
    }

    #[test]
    fn projections_outside_table_16_give_no_properties() {
        for projection in [
            Projection::arc(),
            Projection::cea(),
            Projection::mol(),
            Projection::coe(45.0, 25.0),
            Projection::tsc(),
            Projection::sin_with_params(0.1, 0.2),
        ] {
            let result = wcs_to_xisf_properties(&built(projection.clone()));
            assert_eq!(result, None, "{projection:?}");
        }
    }

    #[test]
    fn projection_system_name_maps_the_table_16_codes() {
        let cases = [
            ("TAN", "Gnomonic"),
            ("SIN", "Orthographic"),
            ("STG", "Stereographic"),
            ("ZEA", "ZenithalEqualArea"),
            ("CAR", "PlateCarree"),
            ("MER", "Mercator"),
            ("AIT", "HammerAitoff"),
        ];
        for (code, expected) in cases {
            assert_eq!(projection_system_name(code), Some(expected));
            let props = wcs_to_xisf_properties(&built(plain_projection(code).unwrap())).unwrap();
            assert_eq!(
                find_prop(&props, PROJECTION_SYSTEM).map(|p| &p.value),
                Some(&XisfPropertyValue::String(expected.to_string()))
            );
        }
    }

    #[test]
    fn projection_system_name_has_no_name_for_other_codes() {
        for code in ["ARC", "COP", "CEA", "MOL", "TPV", "XYZ", ""] {
            assert_eq!(projection_system_name(code), None, "{code}");
        }
    }

    #[test]
    fn keywords_without_a_celestial_wcs_give_no_properties() {
        let linear = vec![
            Keyword::string("CTYPE1", "LINEAR"),
            Keyword::real("CRVAL1", 150.0),
        ];
        let numeric_ctype = vec![Keyword::real("CTYPE1", 1.0), Keyword::real("CRVAL2", 30.0)];
        for keywords in [vec![], linear, numeric_ctype] {
            assert_eq!(solution_properties(&keywords), vec![]);
        }
    }

    #[test]
    fn celestial_keywords_that_do_not_build_give_no_properties() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::real("CD1_1", 0.0));
        set(&mut keywords, Keyword::real("CD2_2", 0.0));
        assert_eq!(solution_properties(&keywords), vec![]);
    }

    // D7: XISF puts the centre of pixel i at i + 0.5, FITS at i + 1.
    #[test]
    fn full_input_emits_every_section() {
        let props = solution_properties(&tan_keywords());
        assert_eq!(
            props.iter().map(|p| p.id.as_str()).collect::<Vec<_>>(),
            crate::xisf::solution::SOLUTION_IDS
        );
        assert_eq!(vector(&props, REFERENCE_CELESTIAL), Some(vec![10.0, 20.0]));
        assert_eq!(vector(&props, REFERENCE_IMAGE), Some(vec![511.5, 511.5]));
        assert_eq!(
            find_prop(&props, LINEAR_MATRIX).map(|p| &p.value),
            Some(&XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 2,
                data: vec![0.001, 0.0, 0.0, 0.001]
            })
        );
        assert_eq!(vector(&props, REFERENCE_NATIVE), Some(vec![0.0, 90.0]));
        // Zenithal: δ0 = 20 < θ0 = 90, so φp = 180, and δp = δ0.
        assert_eq!(
            vector(&props, CELESTIAL_POLE_NATIVE),
            Some(vec![180.0, 20.0])
        );
    }

    #[test]
    fn repeated_wcs_keywords_use_the_last_card() {
        let mut keywords = tan_keywords();
        keywords.extend([
            Keyword::real("CRVAL1", 30.0),
            Keyword::string("RADESYS", "FK4"),
            Keyword::string("RADESYS", "ICRS"),
        ]);
        let props = solution_properties(&keywords);
        assert_eq!(vector(&props, REFERENCE_CELESTIAL), Some(vec![30.0, 20.0]));
    }

    #[test]
    fn integer_valued_keywords_are_read_as_reals() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::integer("CRVAL1", 5));
        set(&mut keywords, Keyword::integer("CRPIX2", 100));
        let props = solution_properties(&keywords);
        assert_eq!(vector(&props, REFERENCE_CELESTIAL), Some(vec![5.0, 20.0]));
        assert_eq!(vector(&props, REFERENCE_IMAGE), Some(vec![511.5, 99.5]));
    }

    #[test]
    fn phi_p_defaults_to_zero_when_crval2_is_at_pole() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::real("CRVAL2", 90.0));
        let props = solution_properties(&keywords);
        assert_eq!(vector(&props, CELESTIAL_POLE_NATIVE), Some(vec![0.0, 90.0]));
    }

    #[test]
    fn phi_p_defaults_to_180_below_pole() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::real("CRVAL2", -30.0));
        let props = solution_properties(&keywords);
        assert_eq!(
            vector(&props, CELESTIAL_POLE_NATIVE),
            Some(vec![180.0, -30.0])
        );
    }

    // For a zenithal projection δp is δ0, so LATPOLE does not change it.
    #[test]
    fn explicit_lonpole_sets_phi_p_and_latpole_does_not_move_a_zenithal_pole() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::real("CRVAL2", 30.0));
        set(&mut keywords, Keyword::real("LONPOLE", 45.0));
        set(&mut keywords, Keyword::real("LATPOLE", 75.0));
        let props = solution_properties(&keywords);
        assert_eq!(
            vector(&props, CELESTIAL_POLE_NATIVE),
            Some(vec![45.0, 30.0])
        );
    }

    #[test]
    fn unknown_wcs_keywords_are_ignored() {
        let mut keywords = tan_keywords();
        keywords.push(Keyword::real("UNKNOWN", 42.0));
        let props = solution_properties(&keywords);
        assert_eq!(vector(&props, REFERENCE_CELESTIAL), Some(vec![10.0, 20.0]));
    }

    #[test]
    fn galactic_axes_give_no_properties() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::string("CTYPE1", "GLON-TAN"));
        set(&mut keywords, Keyword::string("CTYPE2", "GLAT-TAN"));
        assert_eq!(solution_properties(&keywords), vec![]);
    }

    #[test]
    fn fiducial_offset_gives_no_properties() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::string("CTYPE1", "RA---CAR"));
        set(&mut keywords, Keyword::string("CTYPE2", "DEC--CAR"));
        keywords.push(Keyword::real("PV1_0", 1.0));
        keywords.push(Keyword::real("PV1_2", 10.0));
        assert_eq!(solution_properties(&keywords), vec![]);
    }

    #[test]
    fn reference_systems_other_than_icrs_and_fk5_j2000_give_no_properties() {
        let accepted = [
            vec![Keyword::string("RADESYS", "ICRS")],
            vec![Keyword::string("RADESYS", "FK5")],
            vec![
                Keyword::string("RADESYS", "FK5"),
                Keyword::integer("EQUINOX", 2000),
            ],
        ];
        for extra in accepted {
            let mut keywords = tan_keywords();
            keywords.extend(extra);
            assert_eq!(solution_properties(&keywords).len(), 6);
        }
        let refused = [
            vec![Keyword::string("RADESYS", "FK4")],
            vec![
                Keyword::string("RADESYS", "FK5"),
                Keyword::real("EQUINOX", 1950.0),
            ],
            vec![Keyword::string("RADESYS", "GAPPT")],
        ];
        for extra in refused {
            let mut keywords = tan_keywords();
            keywords.extend(extra);
            assert_eq!(solution_properties(&keywords), vec![]);
        }
    }

    // SIP terms stay in the FITS keywords; the properties hold the linear part.
    #[test]
    fn sip_keywords_give_the_linear_solution() {
        let mut keywords = tan_keywords();
        set(&mut keywords, Keyword::string("CTYPE1", "RA---TAN-SIP"));
        set(&mut keywords, Keyword::string("CTYPE2", "DEC--TAN-SIP"));
        keywords.extend([
            Keyword::integer("A_ORDER", 2),
            Keyword::integer("B_ORDER", 2),
            Keyword::real("A_2_0", 1e-6),
            Keyword::real("B_0_2", -2e-6),
        ]);
        assert_eq!(
            solution_properties(&keywords),
            solution_properties(&tan_keywords())
        );
    }
}
