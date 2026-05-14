//! Writing solved WCS into image headers and files.
//!
//! Methods on [`SolveResult`] for producing FITS keywords, XISF properties, and for
//! saving a solved image to disk.

use std::path::Path;

use celestial_catalog::query::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;
use celestial_images::fits::header::Keyword;
use celestial_images::formats::Image;
use celestial_images::xisf::{XisfProperty, XisfPropertyValue};
use celestial_images::Result as ImageResult;

use crate::fit_wcs::{SipSolution, WcsSolution};
use crate::match_field::StarPair;
use crate::solve::SolveResult;

const SPLINE_PREFIX: &str = "PCL:AstrometricSolution:SplineWorldTransformation:";

impl SolveResult {
    /// Builds the FITS WCS keywords (plus SIP coefficients, if SIP was fit).
    ///
    /// Returns the 12 standard TAN WCS cards (`CTYPE1/2`, `CUNIT1/2`, `CRPIX1/2`,
    /// `CRVAL1/2`, `CD1_1/CD1_2/CD2_1/CD2_2`) and, if SIP was fit, the `A_ORDER`,
    /// `B_ORDER`, and `A_p_q` / `B_p_q` coefficient cards.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// for kw in result.to_fits_keywords() {
    ///     println!("{}", kw.name);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn to_fits_keywords(&self) -> Vec<Keyword> {
        let mut out = wcs_keywords(&self.wcs);
        if let Some(sip) = &self.sip {
            out.extend(sip_keywords(sip));
        }
        out
    }

    /// Builds the PixInsight spline-transformation XISF properties.
    ///
    /// Returns an empty vector when there are fewer than 5 matched pairs (the minimum
    /// for a spline fit). Otherwise returns the 11 properties under
    /// `PCL:AstrometricSolution:SplineWorldTransformation:*` that PixInsight reads for
    /// its WCS overlay.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// let props = result.to_xisf_properties();
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn to_xisf_properties(&self) -> Vec<XisfProperty> {
        if self.pairs.len() >= 5 {
            pcl_spline_properties(&self.wcs, &self.pairs)
        } else {
            Vec::new()
        }
    }

    /// Applies WCS/SIP keywords and XISF properties to `img` in place.
    ///
    /// The XISF properties are added unconditionally — if `img` is later saved as FITS,
    /// the FITS writer ignores them. Existing keywords with the same names are
    /// overwritten.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    ///
    /// let mut with_wcs = img.clone();
    /// result.write_into(&mut with_wcs);
    /// with_wcs.save("m31_solved.fits")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_into(&self, img: &mut Image) {
        for kw in self.to_fits_keywords() {
            img.set_keyword(kw);
        }
        img.add_xisf_properties(self.to_xisf_properties());
    }

    /// Clones `img`, applies the solved WCS, and writes it to `path`.
    ///
    /// The input `img` is not modified. The output format is dispatched on the path
    /// extension — `.fits` / `.fit` / `.fts` for FITS, `.xisf` for XISF. Missing
    /// parent directories are created.
    ///
    /// # Errors
    ///
    /// Returns [`celestial_images::ImageError`] for filesystem errors, unsupported
    /// path extensions, or write failures.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// result.save_with(&img, "out/solved.fits")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn save_with<P: AsRef<Path>>(&self, img: &Image, path: P) -> ImageResult<()> {
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let mut out = img.clone();
        self.write_into(&mut out);
        out.save(path)
    }
}

fn wcs_keywords(wcs: &WcsSolution) -> Vec<Keyword> {
    vec![
        Keyword::string("CTYPE1", "RA---TAN"),
        Keyword::string("CTYPE2", "DEC--TAN"),
        Keyword::string("CUNIT1", "deg"),
        Keyword::string("CUNIT2", "deg"),
        Keyword::real("CRPIX1", wcs.crpix1),
        Keyword::real("CRPIX2", wcs.crpix2),
        Keyword::real("CRVAL1", wcs.crval1),
        Keyword::real("CRVAL2", wcs.crval2),
        Keyword::real("CD1_1", wcs.cd1_1),
        Keyword::real("CD1_2", wcs.cd1_2),
        Keyword::real("CD2_1", wcs.cd2_1),
        Keyword::real("CD2_2", wcs.cd2_2),
    ]
}

fn sip_keywords(sip: &SipSolution) -> Vec<Keyword> {
    let mut out = vec![
        Keyword::integer("A_ORDER", sip.a_order as i64),
        Keyword::integer("B_ORDER", sip.b_order as i64),
    ];
    let mut a_keys: Vec<_> = sip.a_coeffs.iter().collect();
    a_keys.sort_by_key(|((p, q), _)| (*p, *q));
    for ((p, q), c) in a_keys {
        out.push(Keyword::real(format!("A_{}_{}", p, q), *c));
    }
    let mut b_keys: Vec<_> = sip.b_coeffs.iter().collect();
    b_keys.sort_by_key(|((p, q), _)| (*p, *q));
    for ((p, q), c) in b_keys {
        out.push(Keyword::real(format!("B_{}_{}", p, q), *c));
    }
    out
}

fn pcl_spline_properties(wcs: &WcsSolution, pairs: &[StarPair]) -> Vec<XisfProperty> {
    let (world, image) = build_control_points(wcs, pairs);
    let linear_approx = linear_world_to_image(wcs);

    vec![
        XisfProperty::string(format!("{SPLINE_PREFIX}Version"), "2.0"),
        XisfProperty::string(format!("{SPLINE_PREFIX}RBFType"), "ThinPlateSpline"),
        XisfProperty::int32(format!("{SPLINE_PREFIX}SplineOrder"), 2),
        XisfProperty::float64(format!("{SPLINE_PREFIX}SplineSmoothness"), 0.025),
        XisfProperty::int32(format!("{SPLINE_PREFIX}MaxSplinePoints"), 2100),
        XisfProperty::new(
            format!("{SPLINE_PREFIX}UseSimplifiers"),
            XisfPropertyValue::Boolean(true),
        ),
        XisfProperty::float64(format!("{SPLINE_PREFIX}SimplifierRejectFraction"), 0.10),
        XisfProperty::new(
            format!("{SPLINE_PREFIX}Truncated"),
            XisfPropertyValue::Boolean(false),
        ),
        XisfProperty::f64_vector(format!("{SPLINE_PREFIX}ControlPoints:World"), world),
        XisfProperty::f64_vector(format!("{SPLINE_PREFIX}ControlPoints:Image"), image),
        XisfProperty::f64_matrix(
            format!("{SPLINE_PREFIX}LinearApproximation"),
            2,
            3,
            linear_approx,
        ),
    ]
}

fn build_control_points(wcs: &WcsSolution, pairs: &[StarPair]) -> (Vec<f64>, Vec<f64>) {
    let mut world = Vec::with_capacity(pairs.len() * 2);
    let mut image = Vec::with_capacity(pairs.len() * 2);
    for p in pairs {
        if let Some((xi_rad, eta_rad)) =
            tan_project_star(p.ra_deg, p.dec_deg, wcs.crval1, wcs.crval2)
        {
            world.push(xi_rad * RAD_TO_DEG);
            world.push(eta_rad * RAD_TO_DEG);
            image.push(p.px_x);
            image.push(p.px_y);
        }
    }
    (world, image)
}

fn linear_world_to_image(wcs: &WcsSolution) -> Vec<f64> {
    let det = wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1;
    let a = wcs.cd2_2 / det;
    let b = -wcs.cd1_2 / det;
    let c = -wcs.cd2_1 / det;
    let d = wcs.cd1_1 / det;
    vec![a, b, wcs.crpix1, c, d, wcs.crpix2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::DetectedStar;
    use crate::match_field::FieldMatch;
    use celestial_images::fits::header::KeywordValue;
    use std::collections::HashMap;

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            crpix1: 1024.0,
            crpix2: 768.0,
            crval1: 180.0,
            crval2: 30.0,
            cd1_1: 0.001,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            width: 2048,
            height: 1536,
            focal_mm: Some(1000.0),
            pixel_um: Some(3.76),
            n_stars: 50,
            rms_px: 0.3,
            weighted_rms_px: 0.25,
            residuals: Vec::new(),
        }
    }

    fn sample_sip(order: u32) -> SipSolution {
        let mut a_coeffs = HashMap::new();
        let mut b_coeffs = HashMap::new();
        a_coeffs.insert((2, 0), 1e-6);
        a_coeffs.insert((0, 2), 2e-6);
        a_coeffs.insert((1, 1), 3e-6);
        b_coeffs.insert((0, 2), 4e-6);
        b_coeffs.insert((2, 0), 5e-6);

        SipSolution {
            crpix1: 1024.0,
            crpix2: 768.0,
            crval1: 180.0,
            crval2: 30.0,
            cd1_1: 0.001,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            width: 2048,
            height: 1536,
            a_order: order,
            b_order: order,
            a_coeffs,
            b_coeffs,
            n_stars: 25,
            linear_rms_px: 0.4,
            rms_px: 0.2,
            weighted_rms_px: 0.18,
        }
    }

    fn pair(px_x: f64, px_y: f64, ra_deg: f64, dec_deg: f64) -> StarPair {
        StarPair {
            px_x,
            px_y,
            ra_deg,
            dec_deg,
            votes: 4,
            snr: 20.0,
        }
    }

    fn build_result(wcs: WcsSolution, sip: Option<SipSolution>, pairs: Vec<StarPair>) -> SolveResult {
        let field_match = FieldMatch {
            image_stars: Vec::new(),
            catalog_stars: Vec::new(),
            matches: Vec::new(),
            pairs: pairs.clone(),
        };
        SolveResult {
            stars: Vec::<DetectedStar>::new(),
            field_match,
            initial_wcs: wcs.clone(),
            wcs,
            sip,
            refined: true,
            pairs,
        }
    }

    fn find(keywords: &[Keyword], name: &str) -> Option<KeywordValue> {
        keywords
            .iter()
            .find(|k| k.name == name)
            .and_then(|k| k.value.clone())
    }

    #[test]
    fn fits_keywords_without_sip_emit_only_wcs_cards() {
        let result = build_result(sample_wcs(), None, Vec::new());
        let kws = result.to_fits_keywords();
        assert_eq!(kws.len(), 12);
        let names: Vec<&str> = kws.iter().map(|k| k.name.as_str()).collect();
        assert_eq!(
            names,
            vec![
                "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                "CD1_1", "CD1_2", "CD2_1", "CD2_2",
            ]
        );
        assert_eq!(
            find(&kws, "CTYPE1"),
            Some(KeywordValue::String("RA---TAN".to_string()))
        );
        assert_eq!(
            find(&kws, "CTYPE2"),
            Some(KeywordValue::String("DEC--TAN".to_string()))
        );
        assert_eq!(find(&kws, "CRVAL1"), Some(KeywordValue::Real(180.0)));
        assert_eq!(find(&kws, "CRVAL2"), Some(KeywordValue::Real(30.0)));
        assert_eq!(find(&kws, "CRPIX1"), Some(KeywordValue::Real(1024.0)));
        assert_eq!(find(&kws, "CRPIX2"), Some(KeywordValue::Real(768.0)));
        assert_eq!(find(&kws, "CD1_1"), Some(KeywordValue::Real(0.001)));
        assert_eq!(find(&kws, "CD2_2"), Some(KeywordValue::Real(0.001)));
    }

    #[test]
    fn fits_keywords_with_sip_appends_sorted_a_and_b_coefficients() {
        let result = build_result(sample_wcs(), Some(sample_sip(3)), Vec::new());
        let kws = result.to_fits_keywords();
        assert_eq!(kws.len(), 12 + 2 + 3 + 2);

        assert_eq!(find(&kws, "A_ORDER"), Some(KeywordValue::Integer(3)));
        assert_eq!(find(&kws, "B_ORDER"), Some(KeywordValue::Integer(3)));

        let a_order: Vec<&str> = kws
            .iter()
            .filter(|k| k.name.starts_with("A_") && k.name != "A_ORDER")
            .map(|k| k.name.as_str())
            .collect();
        assert_eq!(a_order, vec!["A_0_2", "A_1_1", "A_2_0"]);

        let b_order: Vec<&str> = kws
            .iter()
            .filter(|k| k.name.starts_with("B_") && k.name != "B_ORDER")
            .map(|k| k.name.as_str())
            .collect();
        assert_eq!(b_order, vec!["B_0_2", "B_2_0"]);

        assert_eq!(find(&kws, "A_1_1"), Some(KeywordValue::Real(3e-6)));
        assert_eq!(find(&kws, "B_0_2"), Some(KeywordValue::Real(4e-6)));
    }

    #[test]
    fn xisf_properties_empty_when_pairs_below_threshold() {
        let pairs = vec![pair(100.0, 100.0, 180.01, 30.01); 4];
        let result = build_result(sample_wcs(), None, pairs);
        assert!(result.to_xisf_properties().is_empty());
    }

    #[test]
    fn xisf_properties_emitted_when_pairs_meet_threshold() {
        let pairs: Vec<StarPair> = (0..5)
            .map(|i| {
                let off = i as f64 * 0.01;
                pair(100.0 + i as f64, 200.0 + i as f64, 180.0 + off, 30.0 + off)
            })
            .collect();
        let result = build_result(sample_wcs(), None, pairs);
        let props = result.to_xisf_properties();

        assert_eq!(props.len(), 11);

        let ids: Vec<&str> = props.iter().map(|p| p.id.as_str()).collect();
        let expect_ids = [
            "Version",
            "RBFType",
            "SplineOrder",
            "SplineSmoothness",
            "MaxSplinePoints",
            "UseSimplifiers",
            "SimplifierRejectFraction",
            "Truncated",
            "ControlPoints:World",
            "ControlPoints:Image",
            "LinearApproximation",
        ];
        for (id, suffix) in ids.iter().zip(expect_ids.iter()) {
            assert!(
                id.ends_with(suffix),
                "expected id to end with {suffix:?}, got {id:?}"
            );
        }

        let get = |suffix: &str| props.iter().find(|p| p.id.ends_with(suffix)).unwrap();

        assert!(matches!(&get("Version").value, XisfPropertyValue::String(s) if s == "2.0"));
        assert!(
            matches!(&get("RBFType").value, XisfPropertyValue::String(s) if s == "ThinPlateSpline")
        );
        assert!(matches!(get("SplineOrder").value, XisfPropertyValue::Int32(2)));
        assert!(matches!(get("MaxSplinePoints").value, XisfPropertyValue::Int32(2100)));
        assert!(matches!(get("UseSimplifiers").value, XisfPropertyValue::Boolean(true)));
        assert!(matches!(get("Truncated").value, XisfPropertyValue::Boolean(false)));
        assert!(matches!(get("SplineSmoothness").value, XisfPropertyValue::Float64(s) if s == 0.025));
        assert!(matches!(
            get("SimplifierRejectFraction").value,
            XisfPropertyValue::Float64(s) if s == 0.10
        ));

        if let XisfPropertyValue::F64Vector(v) = &get("ControlPoints:World").value {
            assert_eq!(v.len(), 10);
        } else {
            panic!("ControlPoints:World should be F64Vector");
        }
        if let XisfPropertyValue::F64Vector(v) = &get("ControlPoints:Image").value {
            assert_eq!(v.len(), 10);
            assert_eq!(v[0], 100.0);
            assert_eq!(v[1], 200.0);
        } else {
            panic!("ControlPoints:Image should be F64Vector");
        }

        if let XisfPropertyValue::F64Matrix { rows, cols, data } = &get("LinearApproximation").value
        {
            assert_eq!(*rows, 2);
            assert_eq!(*cols, 3);
            assert_eq!(data.len(), 6);
            // CD is diag(0.001, 0.001): inverse is diag(1000, 1000)
            assert!((data[0] - 1000.0).abs() < 1e-9);
            assert!((data[1] - 0.0).abs() < 1e-9);
            assert_eq!(data[2], 1024.0);
            assert!((data[3] - 0.0).abs() < 1e-9);
            assert!((data[4] - 1000.0).abs() < 1e-9);
            assert_eq!(data[5], 768.0);
        } else {
            panic!("LinearApproximation should be F64Matrix");
        }
    }

    #[test]
    fn xisf_properties_skip_pairs_that_fail_tan_projection() {
        let mut pairs: Vec<StarPair> = (0..5)
            .map(|i| {
                let off = i as f64 * 0.01;
                pair(100.0 + i as f64, 200.0 + i as f64, 180.0 + off, 30.0 + off)
            })
            .collect();
        pairs.push(pair(500.0, 500.0, 0.0, -30.0));
        let result = build_result(sample_wcs(), None, pairs);
        let props = result.to_xisf_properties();

        let world = props.iter().find(|p| p.id.ends_with("ControlPoints:World")).unwrap();
        if let XisfPropertyValue::F64Vector(v) = &world.value {
            assert_eq!(v.len(), 10, "anti-center pair should be dropped");
        } else {
            panic!("expected F64Vector");
        }
    }

    #[test]
    fn linear_approximation_reflects_cd_inverse_for_rotated_matrix() {
        let mut wcs = sample_wcs();
        wcs.cd1_1 = 0.0;
        wcs.cd1_2 = 0.001;
        wcs.cd2_1 = -0.001;
        wcs.cd2_2 = 0.0;
        let pairs: Vec<StarPair> = (0..5)
            .map(|i| {
                let off = i as f64 * 0.01;
                pair(100.0 + i as f64, 200.0 + i as f64, 180.0 + off, 30.0 + off)
            })
            .collect();
        let result = build_result(wcs, None, pairs);
        let props = result.to_xisf_properties();
        let linear = props
            .iter()
            .find(|p| p.id.ends_with("LinearApproximation"))
            .unwrap();
        if let XisfPropertyValue::F64Matrix { data, .. } = &linear.value {
            // det = 0*0 - 0.001*(-0.001) = 1e-6
            // a = d22/det = 0/1e-6 = 0
            // b = -d12/det = -0.001/1e-6 = -1000
            // c = -d21/det = 0.001/1e-6 = 1000
            // d = d11/det = 0/1e-6 = 0
            assert!((data[0] - 0.0).abs() < 1e-9);
            assert!((data[1] - (-1000.0)).abs() < 1e-9);
            assert!((data[3] - 1000.0).abs() < 1e-9);
            assert!((data[4] - 0.0).abs() < 1e-9);
        } else {
            panic!("expected F64Matrix");
        }
    }

    fn blank_image() -> Image {
        use celestial_images::formats::PixelData;
        Image::new(PixelData::F32(vec![0.0; 16]), [4, 4])
    }

    #[test]
    fn write_into_applies_wcs_keywords_to_image() {
        let mut img = blank_image();
        let result = build_result(sample_wcs(), None, Vec::new());
        result.write_into(&mut img);

        assert!(img.get_keyword("CRVAL1").is_some());
        assert!(img.get_keyword("CRVAL2").is_some());
        assert!(img.get_keyword("CD1_1").is_some());
        assert!(img.get_keyword("CTYPE1").is_some());
    }

    #[test]
    fn write_into_includes_sip_keywords_when_present() {
        let mut img = blank_image();
        let result = build_result(sample_wcs(), Some(sample_sip(3)), Vec::new());
        result.write_into(&mut img);

        assert_eq!(
            img.get_keyword("A_ORDER").unwrap().value.clone().unwrap(),
            KeywordValue::Integer(3)
        );
        assert!(img.get_keyword("B_1_1").is_none(), "unset coeff should not appear");
        assert!(img.get_keyword("A_1_1").is_some());
    }

    #[test]
    fn save_with_writes_fits_and_round_trips_keywords() {
        let tmp = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
        let img = blank_image();
        let result = build_result(sample_wcs(), None, Vec::new());
        result.save_with(&img, tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert!(restored.get_keyword("CRVAL1").is_some());
        assert!(restored.get_keyword("CD1_1").is_some());
    }

    #[test]
    fn save_with_does_not_mutate_caller_image() {
        let tmp = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
        let img = blank_image();
        let result = build_result(sample_wcs(), None, Vec::new());
        result.save_with(&img, tmp.path()).unwrap();

        assert!(img.get_keyword("CRVAL1").is_none(), "caller's image must not be modified");
    }

    #[test]
    fn save_with_creates_missing_parent_directories() {
        let tmp = tempfile::tempdir().unwrap();
        let nested = tmp.path().join("a/b/c/out.fits");
        let img = blank_image();
        let result = build_result(sample_wcs(), None, Vec::new());
        result.save_with(&img, &nested).unwrap();

        assert!(nested.exists());
    }

    #[test]
    fn save_with_unsupported_extension_errors() {
        let tmp = tempfile::Builder::new().suffix(".bogus").tempfile().unwrap();
        let img = blank_image();
        let result = build_result(sample_wcs(), None, Vec::new());
        assert!(result.save_with(&img, tmp.path()).is_err());
    }
}
