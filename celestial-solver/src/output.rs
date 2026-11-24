//! Writing solved WCS into image headers and files.
//!
//! Methods on [`SolveResult`] for producing FITS keywords, XISF properties, and for
//! saving a solved image to disk.

use std::path::Path;

use celestial_catalog::query::quad::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;
use celestial_images::core::errors::Result as ImageResult;
use celestial_images::fits::header::keywords::Keyword;
use celestial_images::formats::image::Image;
use celestial_images::xisf::header::property::{XisfProperty, XisfPropertyValue};

use crate::fit_wcs::solution::WcsSolution;
use crate::match_field::StarPair;
use crate::pixel_frame::{detection_to_pixinsight, pixinsight_linear_transformation};
use crate::solve::SolveResult;

const SOLUTION_PREFIX: &str = "PCL:AstrometricSolution:";
const SPLINE_PREFIX: &str = "PCL:AstrometricSolution:SplineWorldTransformation:";

impl SolveResult {
    /// Builds the FITS WCS keywords of [`SolveResult::fits_wcs`].
    ///
    /// Returns `CTYPE1/2`, `CUNIT1/2`, `CRPIX1/2` (1-based FITS pixel coordinates),
    /// `CRVAL1/2`, `CD1_1/CD1_2/CD2_1/CD2_2` and `LATPOLE`, which equals `CRVAL2` for the
    /// TAN projection. If SIP was fit, the CTYPEs end in `-SIP` and the `A_ORDER`,
    /// `B_ORDER`, `A_p_q`, `B_p_q` cards follow, with the inverse `AP_ORDER`, `BP_ORDER`,
    /// `AP_p_q` and `BP_p_q` cards. `RADESYS = 'ICRS'` comes last.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// for kw in result.to_fits_keywords() {
    ///     println!("{}", kw.name);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn to_fits_keywords(&self) -> Vec<Keyword> {
        let keywords = self.fits_wcs.to_keywords().into_iter().map(Keyword::from);
        // The catalog positions are ICRS.
        keywords
            .chain([Keyword::string("RADESYS", "ICRS")])
            .collect()
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
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
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
    /// Every earlier WCS keyword of axes 1 and 2 (`CTYPEi`, `CRPIXi`, `CRVALi`, `CUNITi`,
    /// `CDi_j`, `PCi_j` and their `CDiiijjj` and `PCiiijjj` forms, `CDELTi`, `CROTAi`,
    /// `PVi_m`, SIP `A_*`, `B_*`, `AP_*`, `BP_*`, the lookup-table cards `CPDISi`,
    /// `CPERRi`, `CQDISi`, `CQERRi`, `DPi`, `DQi`, `D2IMDISi`, `D2IMERRi`, and `WCSAXES`,
    /// `LONPOLE`, `LATPOLE`, `RADESYS`, `RADECSYS`) and every earlier
    /// `PCL:AstrometricSolution:*` property is removed first, so the image holds one
    /// solution. `EQUINOX`, `EPOCH` and the keywords of other axes are kept. The XISF properties are added unconditionally; if
    /// `img` is later saved as FITS, the FITS writer ignores them.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    ///
    /// let mut with_wcs = img.clone();
    /// result.write_into(&mut with_wcs);
    /// with_wcs.save("m31_solved.fits")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_into(&self, img: &mut Image) {
        img.keywords.retain(|k| !is_wcs_keyword(&k.name));
        for kw in self.to_fits_keywords() {
            img.set_keyword(kw);
        }
        img.xisf_properties
            .retain(|p| !p.id.starts_with(SOLUTION_PREFIX));
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
    /// Returns [`celestial_images::core::errors::ImageError`] for filesystem errors, unsupported
    /// path extensions, or write failures.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
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

// A keyword of the celestial WCS on image axes 1 and 2, in any form a reader would
// combine with the keywords `to_fits_keywords` writes: the current and pre-1996 matrix
// forms, the SIP and distortion-paper lookup-table cards, and the reference frame.
// EQUINOX and EPOCH are not: a capture program may give them for its own coordinates,
// and the RADESYS = 'ICRS' that `to_fits_keywords` writes makes the solution
// independent of them.
fn is_wcs_keyword(name: &str) -> bool {
    const NAMES: [&str; 11] = [
        "WCSAXES", "LONPOLE", "LATPOLE", "RADESYS", "RADECSYS", "A_ORDER", "B_ORDER", "AP_ORDER",
        "BP_ORDER", "A_DMAX", "B_DMAX",
    ];
    const PER_AXIS: [&str; 14] = [
        "CTYPE", "CRPIX", "CRVAL", "CUNIT", "CDELT", "CROTA", "CPDIS", "CPERR", "CQDIS", "CQERR",
        "DP", "DQ", "D2IMDIS", "D2IMERR",
    ];
    const SIP: [&str; 4] = ["A_", "B_", "AP_", "BP_"];
    let has = |prefixes: &[&str], rest: &dyn Fn(&str) -> bool| {
        prefixes
            .iter()
            .any(|p| name.strip_prefix(p).is_some_and(rest))
    };
    NAMES.contains(&name)
        || has(&PER_AXIS, &is_axis)
        || has(&["CD", "PC"], &is_axis_pair)
        || has(&["PV"], &|s| {
            s.split_once('_')
                .is_some_and(|(i, m)| is_axis(i) && is_index(m))
        })
        || has(&SIP, &|s| {
            s.split_once('_')
                .is_some_and(|(p, q)| is_index(p) && is_index(q))
        })
}

fn is_index(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

fn is_axis(s: &str) -> bool {
    s == "1" || s == "2"
}

// `i_j`, or the pre-1996 `iiijjj`, with i and j each axis 1 or 2.
fn is_axis_pair(s: &str) -> bool {
    if let Some((i, j)) = s.split_once('_') {
        return is_axis(i) && is_axis(j);
    }
    let old_axis = |a: &str| a == "001" || a == "002";
    s.len() == 6 && s.is_char_boundary(3) && old_axis(&s[..3]) && old_axis(&s[3..])
}

fn pcl_spline_properties(wcs: &WcsSolution, pairs: &[StarPair]) -> Vec<XisfProperty> {
    let (world, image) = build_control_points(wcs, pairs);
    let linear_approx = pixinsight_linear_transformation(wcs).to_vec();

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

// The pairs as PixInsight control points: tangent-plane positions about CRVAL in
// degrees, and pixels in PixInsight's convention.
fn build_control_points(wcs: &WcsSolution, pairs: &[StarPair]) -> (Vec<f64>, Vec<f64>) {
    let mut world = Vec::with_capacity(pairs.len() * 2);
    let mut image = Vec::with_capacity(pairs.len() * 2);
    for p in pairs {
        if let Some((xi_rad, eta_rad)) =
            tan_project_star(p.ra_deg, p.dec_deg, wcs.crval1, wcs.crval2)
        {
            world.push(xi_rad * RAD_TO_DEG);
            world.push(eta_rad * RAD_TO_DEG);
            let (x, y) = detection_to_pixinsight((p.px_x, p.px_y));
            image.push(x);
            image.push(y);
        }
    }
    (world, image)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fit_wcs::sip::SipSolution;
    use celestial_images::fits::header::keywords::KeywordValue;
    use std::collections::BTreeMap;

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            focal_mm: Some(1000.0),
            pixel_um: Some(3.76),
            n_stars: 50,
            rms_px: 0.3,
            weighted_rms_px: 0.25,
            ..WcsSolution::diagonal((180.0, 30.0), (0.001, 0.001), (2048, 1536))
        }
    }

    fn sample_sip(order: u32) -> SipSolution {
        let mut a_coeffs = BTreeMap::new();
        let mut b_coeffs = BTreeMap::new();
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

    fn build_result(
        wcs: WcsSolution,
        sip: Option<SipSolution>,
        pairs: Vec<StarPair>,
    ) -> SolveResult {
        SolveResult::of_fit(wcs, sip, pairs)
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
        assert_eq!(kws.len(), 14);
        let names: Vec<&str> = kws.iter().map(|k| k.name.as_str()).collect();
        assert_eq!(
            names,
            vec![
                "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                "CD1_1", "CD1_2", "CD2_1", "CD2_2", "LATPOLE", "RADESYS",
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
        assert_eq!(find(&kws, "CRPIX1"), Some(KeywordValue::Real(1025.0)));
        assert_eq!(find(&kws, "CRPIX2"), Some(KeywordValue::Real(769.0)));
        assert_eq!(find(&kws, "CD1_1"), Some(KeywordValue::Real(0.001)));
        assert_eq!(find(&kws, "CD2_2"), Some(KeywordValue::Real(0.001)));
    }

    #[test]
    fn fits_keywords_with_sip_appends_sorted_a_and_b_coefficients() {
        let result = build_result(sample_wcs(), Some(sample_sip(3)), Vec::new());
        let kws = result.to_fits_keywords();
        // 13 TAN cards, A_ORDER, B_ORDER, 5 A and B terms, AP_ORDER, BP_ORDER, the 15
        // terms of each order-4 inverse polynomial, and RADESYS.
        assert_eq!(kws.len(), 13 + 2 + 5 + 2 + 30 + 1);

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
        assert!(matches!(
            get("SplineOrder").value,
            XisfPropertyValue::Int32(2)
        ));
        assert!(matches!(
            get("MaxSplinePoints").value,
            XisfPropertyValue::Int32(2100)
        ));
        assert!(matches!(
            get("UseSimplifiers").value,
            XisfPropertyValue::Boolean(true)
        ));
        assert!(matches!(
            get("Truncated").value,
            XisfPropertyValue::Boolean(false)
        ));
        assert!(
            matches!(get("SplineSmoothness").value, XisfPropertyValue::Float64(s) if s == 0.025)
        );
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
            assert_eq!(v[0], 100.5);
            assert_eq!(v[1], 200.5);
        } else {
            panic!("ControlPoints:Image should be F64Vector");
        }

        if let XisfPropertyValue::F64Matrix { rows, cols, data } = &get("LinearApproximation").value
        {
            assert_eq!(*rows, 2);
            assert_eq!(*cols, 3);
            assert_eq!(data.len(), 6);
            // CD is diag(0.001, 0.001), and the PixInsight reference pixel is CRPIX + 0.5.
            assert_eq!(data[0], 0.001);
            assert_eq!(data[1], 0.0);
            assert_eq!(data[2], -(0.001 * 1024.5 + 0.0 * 768.5));
            assert_eq!(data[3], 0.0);
            assert_eq!(data[4], 0.001);
            assert_eq!(data[5], -(0.0 * 1024.5 + 0.001 * 768.5));
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

        let world = props
            .iter()
            .find(|p| p.id.ends_with("ControlPoints:World"))
            .unwrap();
        if let XisfPropertyValue::F64Vector(v) = &world.value {
            assert_eq!(v.len(), 10, "anti-center pair should be dropped");
        } else {
            panic!("expected F64Vector");
        }
    }

    #[test]
    fn linear_approximation_reflects_cd_for_rotated_matrix() {
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
            // The 2x2 part maps PixInsight pixels to the tangent plane: it is CD.
            assert_eq!(data[0], 0.0);
            assert_eq!(data[1], 0.001);
            assert_eq!(data[3], -0.001);
            assert_eq!(data[4], 0.0);
        } else {
            panic!("expected F64Matrix");
        }
    }

    fn blank_image() -> Image {
        use celestial_images::formats::pixel_data::PixelData;
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
        assert!(
            img.get_keyword("B_1_1").is_none(),
            "unset coeff should not appear"
        );
        assert!(img.get_keyword("A_1_1").is_some());
    }

    fn five_pairs() -> Vec<StarPair> {
        let pairs = (0..5).map(|i| {
            let off = i as f64 * 0.01;
            pair(100.0 + i as f64, 200.0 + i as f64, 180.0 + off, 30.0 + off)
        });
        pairs.collect()
    }

    // Detection pixel centres are integers from 0; FITS pixel centres are integers from 1.
    #[test]
    fn fits_keywords_write_crpix_in_the_one_based_fits_frame() {
        let result = build_result(sample_wcs(), None, Vec::new());
        let kws = result.to_fits_keywords();
        assert_eq!(find(&kws, "CRPIX1"), Some(KeywordValue::Real(1025.0)));
        assert_eq!(find(&kws, "CRPIX2"), Some(KeywordValue::Real(769.0)));
    }

    #[test]
    fn fits_keywords_with_sip_mark_ctype_sip() {
        let result = build_result(sample_wcs(), Some(sample_sip(3)), Vec::new());
        let kws = result.to_fits_keywords();
        let sip = |s: &str| Some(KeywordValue::String(s.to_string()));
        assert_eq!(find(&kws, "CTYPE1"), sip("RA---TAN-SIP"));
        assert_eq!(find(&kws, "CTYPE2"), sip("DEC--TAN-SIP"));
    }

    // One order above the forward polynomial, as astrometry.net's
    // sip_ensure_inverse_polynomials chooses.
    #[test]
    fn fits_keywords_with_sip_include_the_inverse_polynomial() {
        let result = build_result(sample_wcs(), Some(sample_sip(3)), Vec::new());
        let kws = result.to_fits_keywords();
        assert_eq!(find(&kws, "AP_ORDER"), Some(KeywordValue::Integer(4)));
        assert_eq!(find(&kws, "BP_ORDER"), Some(KeywordValue::Integer(4)));
    }

    fn image_with_stale_wcs() -> Image {
        let mut img = blank_image();
        let stale = [("CTYPE3", "FREQ"), ("CUNIT3", "Hz"), ("RADESYS", "FK5")];
        for (name, value) in stale {
            img.keywords.push(Keyword::string(name, value));
        }
        let stale = [
            "CRPIX3", "CRVAL3", "PC1_1", "PC1_2", "PC2_1", "PC2_2", "CDELT1", "CDELT2", "CROTA1",
            "CROTA2", "PV2_1", "PV1_3", "LONPOLE", "LATPOLE", "EQUINOX", "A_2_0", "B_0_2",
            "AP_1_0", "BP_0_1", "A_DMAX", "B_DMAX",
        ];
        for name in stale {
            img.keywords.push(Keyword::real(name, 1.5));
        }
        for name in ["WCSAXES", "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER"] {
            img.keywords.push(Keyword::integer(name, 2));
        }
        img.keywords.push(Keyword::string("OBJECT", "M42"));
        img.keywords.push(Keyword::real("EXPTIME", 120.0));
        img
    }

    #[test]
    fn write_into_removes_stale_wcs_cards() {
        let mut img = image_with_stale_wcs();
        let result = build_result(sample_wcs(), None, Vec::new());
        result.write_into(&mut img);
        let mut names: Vec<&str> = img.keywords.iter().map(|k| k.name.as_str()).collect();
        names.sort_unstable();
        assert_eq!(
            names,
            [
                "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CRPIX1", "CRPIX2", "CRPIX3", "CRVAL1",
                "CRVAL2", "CRVAL3", "CTYPE1", "CTYPE2", "CTYPE3", "CUNIT1", "CUNIT2", "CUNIT3",
                "EQUINOX", "EXPTIME", "LATPOLE", "OBJECT", "RADESYS",
            ]
        );
    }

    fn names_after_write_into(keywords: Vec<Keyword>) -> Vec<String> {
        let mut img = blank_image();
        img.keywords.extend(keywords);
        build_result(sample_wcs(), None, Vec::new()).write_into(&mut img);
        let mut names: Vec<String> = img.keywords.into_iter().map(|k| k.name).collect();
        names.sort_unstable();
        names
    }

    fn written_names() -> Vec<String> {
        names_after_write_into(Vec::new())
    }

    // The pre-1996 PC and CD forms, `PCiiijjj` and `CDiiijjj`, which readers still apply.
    #[test]
    fn write_into_removes_old_style_pc_and_cd_cards() {
        let stale = [
            "PC001001", "PC001002", "PC002001", "PC002002", "CD001001", "CD002002",
        ];
        let keywords = stale.map(|name| Keyword::real(name, 0.5)).to_vec();
        assert_eq!(names_after_write_into(keywords), written_names());
    }

    // The distortion-paper lookup tables (`CPDISj`, `CPERRj`, the record-valued `DPj`,
    // and their `CQ`/`DQ` forms) and the detector-to-image tables (`D2IMDISj`, `D2IMERRj`).
    #[test]
    fn write_into_removes_distortion_lookup_cards() {
        let mut keywords = Vec::new();
        for name in ["CPDIS1", "CPDIS2", "CQDIS1", "D2IMDIS1", "D2IMDIS2"] {
            keywords.push(Keyword::string(name, "LOOKUP"));
        }
        for name in ["DP1", "DP2", "DQ1"] {
            keywords.push(Keyword::string(name, "NAXES: 2"));
        }
        for name in ["CPERR1", "CPERR2", "CQERR2", "D2IMERR1", "D2IMERR2"] {
            keywords.push(Keyword::real(name, 0.1));
        }
        assert_eq!(names_after_write_into(keywords), written_names());
    }

    // The solution describes axes 1 and 2. A third axis keeps its cards, and EQUINOX, which
    // a capture program may give for its own coordinates, stays: RADESYS = 'ICRS' makes the
    // solution independent of it.
    #[test]
    fn write_into_keeps_a_third_axis_and_equinox() {
        let kept = [
            "CTYPE3", "CRPIX3", "CRVAL3", "CDELT3", "PC3_3", "PC1_3", "EQUINOX",
        ];
        let keywords = kept.map(|name| Keyword::real(name, 3.0)).to_vec();
        let mut expected = written_names();
        expected.extend(kept.map(String::from));
        expected.sort_unstable();
        assert_eq!(names_after_write_into(keywords), expected);
        let result = build_result(sample_wcs(), None, Vec::new());
        assert_eq!(
            find(&result.to_fits_keywords(), "RADESYS"),
            Some(KeywordValue::String("ICRS".to_string()))
        );
    }

    #[test]
    fn write_into_replaces_earlier_astrometric_properties() {
        let mut img = blank_image();
        img.add_xisf_properties([
            XisfProperty::string("PCL:AstrometricSolution:ProjectionSystem", "Orthographic"),
            XisfProperty::f64_vector(format!("{SPLINE_PREFIX}ControlPoints:Image"), vec![1.0]),
            XisfProperty::string("Observation:Object:Name", "M42"),
        ]);
        let result = build_result(sample_wcs(), None, five_pairs());
        result.write_into(&mut img);
        result.write_into(&mut img);
        let mut ids: Vec<&str> = img.xisf_properties.iter().map(|p| p.id.as_str()).collect();
        ids.sort_unstable();
        let mut expected: Vec<String> = result
            .to_xisf_properties()
            .into_iter()
            .map(|p| p.id)
            .chain(["Observation:Object:Name".to_string()])
            .collect();
        expected.sort_unstable();
        assert_eq!(ids, expected);
    }

    // PixInsight's control points are pixel coordinates whose pixel centres are at
    // integers plus one half.
    #[test]
    fn xisf_control_points_use_the_pixinsight_pixel_frame() {
        let result = build_result(sample_wcs(), None, five_pairs());
        let props = result.to_xisf_properties();
        let image = props
            .iter()
            .find(|p| p.id.ends_with("ControlPoints:Image"))
            .unwrap();
        let XisfPropertyValue::F64Vector(v) = &image.value else {
            panic!("ControlPoints:Image should be F64Vector");
        };
        assert_eq!(v[..4], [100.5, 200.5, 101.5, 201.5]);
    }

    // PixInsight's LinearApproximation maps a PixInsight pixel to the tangent plane,
    // world = M * (pixel - reference), as [M | -M * reference].
    #[test]
    fn xisf_linear_approximation_maps_image_to_world() {
        let mut wcs = sample_wcs();
        wcs.cd1_2 = 2e-5;
        wcs.cd2_1 = -3e-5;
        let result = build_result(wcs, None, five_pairs());
        let props = result.to_xisf_properties();
        let linear = props
            .iter()
            .find(|p| p.id.ends_with("LinearApproximation"))
            .unwrap();
        let XisfPropertyValue::F64Matrix { rows, cols, data } = &linear.value else {
            panic!("LinearApproximation should be F64Matrix");
        };
        let (rx, ry) = (1024.5, 768.5);
        let expected = [
            0.001,
            2e-5,
            -(0.001 * rx + 2e-5 * ry),
            -3e-5,
            0.001,
            -(-3e-5 * rx + 0.001 * ry),
        ];
        assert_eq!((*rows, *cols), (2, 3));
        assert_eq!(data[..], expected);
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

        assert!(
            img.get_keyword("CRVAL1").is_none(),
            "caller's image must not be modified"
        );
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
        let tmp = tempfile::Builder::new()
            .suffix(".bogus")
            .tempfile()
            .unwrap();
        let img = blank_image();
        let result = build_result(sample_wcs(), None, Vec::new());
        assert!(result.save_with(&img, tmp.path()).is_err());
    }
}
