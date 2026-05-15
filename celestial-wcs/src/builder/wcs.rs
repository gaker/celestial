
use celestial_core::Angle;

use crate::coordinate::{CelestialCoord, IntermediateCoord, PixelCoord};
use crate::distortion::DistortionModel;
use crate::error::WcsResult;
use crate::linear::LinearTransform;
use crate::spherical::{Projection, SphericalRotation};

use super::keyword::WcsKeyword;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoordType {
    #[default]
    Equatorial,
    Galactic,
    Ecliptic,
    Helioecliptic,
    Supergalactic,
    Generic,
}

impl CoordType {
    pub fn from_ctype_prefix(prefix: &str) -> Self {
        match prefix {
            "RA" | "DEC" => Self::Equatorial,
            "GLON" | "GLAT" => Self::Galactic,
            "ELON" | "ELAT" => Self::Ecliptic,
            "HLON" | "HLAT" => Self::Helioecliptic,
            "SLON" | "SLAT" => Self::Supergalactic,
            _ => Self::Generic,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Wcs {
    linear: LinearTransform,
    projection: Projection,
    rotation: SphericalRotation,
    coord_type: CoordType,
    proj_code: String,
    crval_deg: (f64, f64),
    distortion: Option<DistortionModel>,
}

impl Wcs {
    pub fn new(
        linear: LinearTransform,
        projection: Projection,
        rotation: SphericalRotation,
        coord_type: CoordType,
        proj_code: String,
        crval_deg: (f64, f64),
        distortion: Option<DistortionModel>,
    ) -> Self {
        Self {
            linear,
            projection,
            rotation,
            coord_type,
            proj_code,
            crval_deg,
            distortion,
        }
    }

    pub fn pixel_to_celestial(&self, pixel: PixelCoord) -> WcsResult<CelestialCoord> {
        let pixel = match &self.distortion {
            Some(DistortionModel::Sip(sip)) => {
                let (x, y) = sip.apply(pixel.x(), pixel.y());
                PixelCoord::new(x, y)
            }
            _ => pixel,
        };
        let intermediate = self.linear.pixel_to_intermediate(pixel);
        let intermediate = match &self.distortion {
            Some(DistortionModel::Tpv(tpv)) => {
                let (x, y) = tpv
                    .as_ref()
                    .apply(intermediate.x_deg(), intermediate.y_deg());
                IntermediateCoord::new(x, y)
            }
            Some(DistortionModel::Tnx(tnx)) => {
                let (x, y) = tnx.apply(intermediate.x_deg(), intermediate.y_deg());
                IntermediateCoord::new(x, y)
            }
            _ => intermediate,
        };
        let native = self.projection.deproject(intermediate)?;
        self.rotation.native_to_celestial(native)
    }

    pub fn celestial_to_pixel(&self, celestial: CelestialCoord) -> WcsResult<PixelCoord> {
        let native = self.rotation.celestial_to_native(celestial)?;
        let intermediate = self.projection.project(native)?;
        let intermediate = match &self.distortion {
            Some(DistortionModel::Tpv(tpv)) => {
                let (x, y) = tpv
                    .as_ref()
                    .apply_inverse(intermediate.x_deg(), intermediate.y_deg())?;
                IntermediateCoord::new(x, y)
            }
            Some(DistortionModel::Tnx(tnx)) => {
                let (x, y) = tnx.apply_inverse(intermediate.x_deg(), intermediate.y_deg())?;
                IntermediateCoord::new(x, y)
            }
            _ => intermediate,
        };
        let pixel = self.linear.intermediate_to_pixel(intermediate);
        match &self.distortion {
            Some(DistortionModel::Sip(sip)) => {
                let (x, y) = sip.apply_inverse(pixel.x(), pixel.y())?;
                Ok(PixelCoord::new(x, y))
            }
            _ => Ok(pixel),
        }
    }

    pub fn pix2world(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        let pixel = PixelCoord::new(x, y);
        let celestial = self.pixel_to_celestial(pixel)?;
        Ok((celestial.alpha().degrees(), celestial.delta().degrees()))
    }

    pub fn world2pix(&self, lon: f64, lat: f64) -> WcsResult<(f64, f64)> {
        let celestial = CelestialCoord::new(Angle::from_degrees(lon), Angle::from_degrees(lat));
        let pixel = self.celestial_to_pixel(celestial)?;
        Ok((pixel.x(), pixel.y()))
    }

    #[inline]
    pub fn projection_code(&self) -> &str {
        &self.proj_code
    }

    #[inline]
    pub fn coord_type(&self) -> CoordType {
        self.coord_type
    }

    #[inline]
    pub fn crpix(&self) -> [f64; 2] {
        self.linear.crpix()
    }

    #[inline]
    pub fn crval(&self) -> (f64, f64) {
        self.crval_deg
    }

    #[inline]
    pub fn pixel_scale(&self) -> f64 {
        self.linear.pixel_scale()
    }

    pub fn to_keywords(&self) -> Vec<WcsKeyword> {
        let mut keywords = Vec::new();

        keywords.extend(self.ctype_keywords());
        keywords.extend(self.crpix_keywords());
        keywords.extend(self.crval_keywords());
        keywords.extend(self.cd_keywords());
        keywords.extend(self.pole_keywords());
        keywords.extend(self.pv_keywords());

        keywords
    }

    fn ctype_keywords(&self) -> Vec<WcsKeyword> {
        let (prefix1, prefix2) = ctype_prefixes(&self.coord_type);
        let code = &self.proj_code;

        vec![
            WcsKeyword::string("CTYPE1", format_ctype(prefix1, code)),
            WcsKeyword::string("CTYPE2", format_ctype(prefix2, code)),
        ]
    }

    fn crpix_keywords(&self) -> Vec<WcsKeyword> {
        let crpix = self.linear.crpix();
        vec![
            WcsKeyword::real("CRPIX1", crpix[0]),
            WcsKeyword::real("CRPIX2", crpix[1]),
        ]
    }

    fn crval_keywords(&self) -> Vec<WcsKeyword> {
        vec![
            WcsKeyword::real("CRVAL1", self.crval_deg.0),
            WcsKeyword::real("CRVAL2", self.crval_deg.1),
        ]
    }

    fn cd_keywords(&self) -> Vec<WcsKeyword> {
        let cd = self.linear.cd_matrix();
        vec![
            WcsKeyword::real("CD1_1", cd[0][0]),
            WcsKeyword::real("CD1_2", cd[0][1]),
            WcsKeyword::real("CD2_1", cd[1][0]),
            WcsKeyword::real("CD2_2", cd[1][1]),
        ]
    }

    fn pole_keywords(&self) -> Vec<WcsKeyword> {
        let mut keywords = Vec::new();
        let phi_p = self.rotation.phi_p_degrees();
        let delta_p = self.rotation.delta_p_degrees();

        let default_phi_p = default_lonpole(&self.coord_type, self.crval_deg.1, &self.projection);
        if (phi_p - default_phi_p).abs() > 1e-10 {
            keywords.push(WcsKeyword::real("LONPOLE", phi_p));
        }

        if (delta_p - 90.0).abs() > 1e-10 {
            keywords.push(WcsKeyword::real("LATPOLE", delta_p));
        }

        keywords
    }

    fn pv_keywords(&self) -> Vec<WcsKeyword> {
        projection_pv_keywords(&self.projection)
    }
}

pub(super) fn ctype_prefixes(coord_type: &CoordType) -> (&'static str, &'static str) {
    match coord_type {
        CoordType::Equatorial => ("RA", "DEC"),
        CoordType::Galactic => ("GLON", "GLAT"),
        CoordType::Ecliptic => ("ELON", "ELAT"),
        CoordType::Helioecliptic => ("HLON", "HLAT"),
        CoordType::Supergalactic => ("SLON", "SLAT"),
        CoordType::Generic => ("XLON", "XLAT"),
    }
}

pub(super) fn format_ctype(prefix: &str, proj_code: &str) -> String {
    let padding_len = 4 - prefix.len();
    let dashes = "-".repeat(padding_len + 1);
    format!("{}{}{}", prefix, dashes, proj_code)
}

pub(super) fn default_lonpole(_coord_type: &CoordType, crval_lat: f64, projection: &Projection) -> f64 {
    let (_, theta_0) = projection.native_reference();
    if crval_lat >= theta_0 {
        0.0
    } else {
        180.0
    }
}

pub(super) fn projection_pv_keywords(projection: &Projection) -> Vec<WcsKeyword> {
    match projection {
        Projection::Sin { xi, eta } if *xi != 0.0 || *eta != 0.0 => {
            vec![
                WcsKeyword::real("PV2_1", *xi),
                WcsKeyword::real("PV2_2", *eta),
            ]
        }
        Projection::Azp { mu, gamma } => {
            vec![
                WcsKeyword::real("PV2_1", *mu),
                WcsKeyword::real("PV2_2", *gamma),
            ]
        }
        Projection::Szp { mu, phi_c, theta_c } => {
            vec![
                WcsKeyword::real("PV2_1", *mu),
                WcsKeyword::real("PV2_2", *phi_c),
                WcsKeyword::real("PV2_3", *theta_c),
            ]
        }
        Projection::Zpn { coeffs } => coeffs
            .iter()
            .enumerate()
            .filter(|(_, &v)| v != 0.0)
            .map(|(i, &v)| WcsKeyword::real(format!("PV2_{}", i), v))
            .collect(),
        Projection::Air { theta_b } => {
            vec![WcsKeyword::real("PV2_1", *theta_b)]
        }
        Projection::Cea { lambda } if *lambda != 1.0 => {
            vec![WcsKeyword::real("PV2_1", *lambda)]
        }
        Projection::Cyp { mu, lambda } => {
            vec![
                WcsKeyword::real("PV2_1", *mu),
                WcsKeyword::real("PV2_2", *lambda),
            ]
        }
        Projection::Cop { theta_a, eta }
        | Projection::Coe { theta_a, eta }
        | Projection::Cod { theta_a, eta }
        | Projection::Coo { theta_a, eta } => {
            vec![
                WcsKeyword::real("PV2_1", *theta_a),
                WcsKeyword::real("PV2_2", *eta),
            ]
        }
        Projection::Bon { theta_1 } => {
            vec![WcsKeyword::real("PV2_1", *theta_1)]
        }
        _ => Vec::new(),
    }
}

pub(super) fn projection_to_code(proj: &Projection) -> String {
    match proj {
        Projection::Tan => "TAN",
        Projection::Sin { .. } => "SIN",
        Projection::Arc => "ARC",
        Projection::Stg => "STG",
        Projection::Zea => "ZEA",
        Projection::Azp { .. } => "AZP",
        Projection::Szp { .. } => "SZP",
        Projection::Zpn { .. } => "ZPN",
        Projection::Air { .. } => "AIR",
        Projection::Car => "CAR",
        Projection::Mer => "MER",
        Projection::Cea { .. } => "CEA",
        Projection::Cyp { .. } => "CYP",
        Projection::Sfl => "SFL",
        Projection::Par => "PAR",
        Projection::Mol => "MOL",
        Projection::Ait => "AIT",
        Projection::Cop { .. } => "COP",
        Projection::Coe { .. } => "COE",
        Projection::Cod { .. } => "COD",
        Projection::Coo { .. } => "COO",
        Projection::Bon { .. } => "BON",
        Projection::Pco => "PCO",
        Projection::Tsc => "TSC",
        Projection::Csc => "CSC",
        Projection::Qsc => "QSC",
    }
    .to_string()
}


#[cfg(test)]
mod tests {
    use super::*;

    fn create_simple_tan_wcs() -> crate::error::WcsResult<Wcs> {
        let crpix = [512.0, 512.0];
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let linear = LinearTransform::from_cd(crpix, cd)?;

        let projection = Projection::tan();
        let (_, theta_0) = projection.native_reference();
        let crval_lon = 180.0;
        let crval_lat = 45.0;

        let rotation = SphericalRotation::from_crval(
            Angle::from_degrees(crval_lon),
            Angle::from_degrees(crval_lat),
            Angle::from_degrees(theta_0),
            None,
            None,
        )?;

        Ok(Wcs::new(
            linear,
            projection,
            rotation,
            CoordType::Equatorial,
            "TAN".to_string(),
            (crval_lon, crval_lat),
            None,
        ))
    }

    #[test]
    fn test_wcs_simple_tan_pipeline() {
        use crate::PixelCoord;
        use celestial_core::assert_ulp_lt;

        let wcs = create_simple_tan_wcs().unwrap();

        // CRPIX maps to CRVAL by definition.
        let celestial =
            wcs.pixel_to_celestial(PixelCoord::new(512.0, 512.0)).unwrap();
        assert_ulp_lt!(celestial.alpha().degrees(), 180.0, 10);
        assert_ulp_lt!(celestial.delta().degrees(), 45.0, 10);

        // Round trip via the typed API (PixelCoord <-> CelestialCoord), with a
        // tight ULP bound at CRPIX and an absolute floor off-centre.
        let at_crpix = PixelCoord::new(512.0, 512.0);
        let recovered = wcs
            .celestial_to_pixel(wcs.pixel_to_celestial(at_crpix).unwrap())
            .unwrap();
        assert_ulp_lt!(at_crpix.x(), recovered.x(), 10);
        assert_ulp_lt!(at_crpix.y(), recovered.y(), 10);

        let off_center = PixelCoord::new(256.0, 768.0);
        let recovered = wcs
            .celestial_to_pixel(wcs.pixel_to_celestial(off_center).unwrap())
            .unwrap();
        assert!((off_center.x() - recovered.x()).abs() < 1e-9);
        assert!((off_center.y() - recovered.y()).abs() < 1e-9);

        // The f64 convenience API (pix2world/world2pix) round-trips as well.
        let (x, y) = (300.0, 700.0);
        let (ra, dec) = wcs.pix2world(x, y).unwrap();
        let (xr, yr) = wcs.world2pix(ra, dec).unwrap();
        // Tolerance accounts for ARM vs x86 FPU differences in trig functions.
        assert!((x - xr).abs() < 1e-8);
        assert!((y - yr).abs() < 1e-8);
    }

    #[test]
    fn test_wcs_accessors() {
        let wcs = create_simple_tan_wcs().unwrap();
        assert_eq!(wcs.projection_code(), "TAN");
        assert_eq!(wcs.coord_type(), CoordType::Equatorial);
        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
        assert_eq!(wcs.pixel_scale(), 0.001);
    }

    #[test]
    fn test_coord_type_from_ctype_prefix() {
        assert_eq!(CoordType::from_ctype_prefix("RA"), CoordType::Equatorial);
        assert_eq!(CoordType::from_ctype_prefix("DEC"), CoordType::Equatorial);
        assert_eq!(CoordType::from_ctype_prefix("GLON"), CoordType::Galactic);
        assert_eq!(CoordType::from_ctype_prefix("GLAT"), CoordType::Galactic);
        assert_eq!(CoordType::from_ctype_prefix("ELON"), CoordType::Ecliptic);
        assert_eq!(CoordType::from_ctype_prefix("ELAT"), CoordType::Ecliptic);
        assert_eq!(
            CoordType::from_ctype_prefix("HLON"),
            CoordType::Helioecliptic
        );
        assert_eq!(
            CoordType::from_ctype_prefix("HLAT"),
            CoordType::Helioecliptic
        );
        assert_eq!(
            CoordType::from_ctype_prefix("SLON"),
            CoordType::Supergalactic
        );
        assert_eq!(
            CoordType::from_ctype_prefix("SLAT"),
            CoordType::Supergalactic
        );
        assert_eq!(CoordType::from_ctype_prefix("UNKNOWN"), CoordType::Generic);
    }

    #[test]
    fn test_wcs_roundtrip_across_projections_and_cds() {
        use crate::PixelCoord;

        // Three independent configurations cover:
        //   - rotated CD matrix with TAN at southern declination,
        //   - diagonal CD with ARC near the pole,
        //   - diagonal CD with STG very close to the pole (where rotation
        //     conditioning is most fragile).
        let scale = 0.0005;
        let (sin_th, cos_th) = libm::sincos(celestial_core::constants::PI / 6.0);

        struct Case {
            crpix: [f64; 2],
            cd: [[f64; 2]; 2],
            projection: Projection,
            code: &'static str,
            crval: (f64, f64),
            pixel: (f64, f64),
        }
        let cases = [
            Case {
                crpix: [512.0, 512.0],
                cd: [
                    [scale * cos_th, -scale * sin_th],
                    [scale * sin_th, scale * cos_th],
                ],
                projection: Projection::tan(),
                code: "TAN",
                crval: (120.0, -30.0),
                pixel: (400.0, 600.0),
            },
            Case {
                crpix: [256.0, 256.0],
                cd: [[0.002, 0.0], [0.0, 0.002]],
                projection: Projection::arc(),
                code: "ARC",
                crval: (90.0, 60.0),
                pixel: (200.0, 300.0),
            },
            Case {
                crpix: [128.0, 128.0],
                cd: [[0.005, 0.0], [0.0, 0.005]],
                projection: Projection::stg(),
                code: "STG",
                crval: (0.0, 85.0),
                pixel: (100.0, 150.0),
            },
        ];

        for case in cases {
            let linear = LinearTransform::from_cd(case.crpix, case.cd).unwrap();
            let (_, theta_0) = case.projection.native_reference();
            let rotation = SphericalRotation::from_crval(
                Angle::from_degrees(case.crval.0),
                Angle::from_degrees(case.crval.1),
                Angle::from_degrees(theta_0),
                None,
                None,
            )
            .unwrap();
            let wcs = Wcs::new(
                linear,
                case.projection,
                rotation,
                CoordType::Equatorial,
                case.code.to_string(),
                case.crval,
                None,
            );

            let original = PixelCoord::new(case.pixel.0, case.pixel.1);
            let recovered = wcs
                .celestial_to_pixel(wcs.pixel_to_celestial(original).unwrap())
                .unwrap();
            assert!(
                (original.x() - recovered.x()).abs() < 1e-8,
                "{} x: {} vs {}", case.code, original.x(), recovered.x(),
            );
            assert!(
                (original.y() - recovered.y()).abs() < 1e-8,
                "{} y: {} vs {}", case.code, original.y(), recovered.y(),
            );
        }
    }

    #[test]
    fn test_projection_to_code_all() {
        let cases: &[(Projection, &str)] = &[
            (Projection::tan(), "TAN"),
            (Projection::sin(), "SIN"),
            (Projection::sin_with_params(0.1, 0.2), "SIN"),
            (Projection::arc(), "ARC"),
            (Projection::stg(), "STG"),
            (Projection::zea(), "ZEA"),
            (Projection::azp(2.0, 30.0), "AZP"),
            (Projection::szp(2.0, 45.0, 60.0), "SZP"),
            (Projection::zpn(vec![0.0, 1.0]), "ZPN"),
            (Projection::air(45.0), "AIR"),
            (Projection::car(), "CAR"),
            (Projection::mer(), "MER"),
            (Projection::cea_with_lambda(0.5), "CEA"),
            (Projection::cyp(1.0, 2.0), "CYP"),
            (Projection::sfl(), "SFL"),
            (Projection::par(), "PAR"),
            (Projection::mol(), "MOL"),
            (Projection::ait(), "AIT"),
            (Projection::cop(45.0, 0.0), "COP"),
            (Projection::coe(45.0, 0.0), "COE"),
            (Projection::cod(45.0, 0.0), "COD"),
            (Projection::coo(45.0, 0.0), "COO"),
            (Projection::bon(45.0), "BON"),
            (Projection::pco(), "PCO"),
            (Projection::tsc(), "TSC"),
            (Projection::csc(), "CSC"),
            (Projection::qsc(), "QSC"),
        ];
        for (proj, expected) in cases {
            assert_eq!(projection_to_code(proj), *expected);
        }
    }

    #[test]
    fn test_format_ctype_all() {
        let cases = [
            ("RA", "TAN", "RA---TAN"),
            ("DEC", "TAN", "DEC--TAN"),
            ("GLON", "SIN", "GLON-SIN"),
            ("GLAT", "SIN", "GLAT-SIN"),
            ("ELON", "ARC", "ELON-ARC"),
            ("ELAT", "ARC", "ELAT-ARC"),
        ];
        for (prefix, code, expected) in cases {
            let s = format_ctype(prefix, code);
            assert_eq!(s, expected);
            assert_eq!(s.len(), 8);
        }
    }

    fn tan_wcs_with(distortion: DistortionModel) -> Wcs {
        let linear = LinearTransform::from_cd([512.0, 512.0], [[0.001, 0.0], [0.0, 0.001]]).unwrap();
        let projection = Projection::tan();
        let (_, theta_0) = projection.native_reference();
        let rotation = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(theta_0),
            None,
            None,
        )
        .unwrap();
        Wcs::new(
            linear,
            projection,
            rotation,
            CoordType::Equatorial,
            "TAN".to_string(),
            (180.0, 45.0),
            Some(distortion),
        )
    }

    #[test]
    fn test_pipeline_roundtrips_through_each_distortion_kind() {
        use crate::distortion::{
            CrossTerms, SipDistortion, SurfaceType, TnxDistortion, TnxSurface, TpvDistortion,
        };

        // SIP — pixel-space distortion. A nonzero higher-order coefficient
        // exercises the apply / apply_inverse branches; identity-on-pixels would
        // hide both paths because the iterative inverse starts at the right answer.
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6);
        sip.set_b(0, 2, 1e-6);
        let sip_wcs = tan_wcs_with(DistortionModel::Sip(sip));
        let (ra, dec) = sip_wcs.pix2world(500.0, 500.0).unwrap();
        let (x, y) = sip_wcs.world2pix(ra, dec).unwrap();
        assert!((x - 500.0).abs() < 1e-6, "SIP roundtrip x: got {}", x);
        assert!((y - 500.0).abs() < 1e-6, "SIP roundtrip y: got {}", y);

        // TPV — intermediate-space distortion. PV1_1 = PV2_1 = 1.0 is the
        // linear identity; add small higher-order term to drive the iterative inverse.
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(4, 1e-4);
        tpv.set_pv2(6, 1e-4);
        let tpv_wcs = tan_wcs_with(DistortionModel::Tpv(Box::new(tpv)));
        let (ra, dec) = tpv_wcs.pix2world(500.0, 500.0).unwrap();
        let (x, y) = tpv_wcs.world2pix(ra, dec).unwrap();
        assert!((x - 500.0).abs() < 1e-6, "TPV roundtrip x: got {}", x);
        assert!((y - 500.0).abs() < 1e-6, "TPV roundtrip y: got {}", y);

        // TNX — intermediate-space distortion. Identity-by-coefficients (only the
        // x^1 and y^1 terms set to 1.0 with no_cross) exercises the dispatch.
        let lng = TnxSurface::new(
            SurfaceType::Polynomial,
            2,
            2,
            CrossTerms::None,
            (-1.0, 1.0),
            (-1.0, 1.0),
            vec![0.0, 1.0, 0.0, 1.0],
        )
        .unwrap();
        let lat = TnxSurface::new(
            SurfaceType::Polynomial,
            2,
            2,
            CrossTerms::None,
            (-1.0, 1.0),
            (-1.0, 1.0),
            vec![0.0, 1.0, 0.0, 1.0],
        )
        .unwrap();
        let tnx_wcs = tan_wcs_with(DistortionModel::Tnx(TnxDistortion::new(lng, lat)));
        let (ra, dec) = tnx_wcs.pix2world(500.0, 500.0).unwrap();
        let (x, y) = tnx_wcs.world2pix(ra, dec).unwrap();
        assert!((x - 500.0).abs() < 1e-6, "TNX roundtrip x: got {}", x);
        assert!((y - 500.0).abs() < 1e-6, "TNX roundtrip y: got {}", y);
    }

    #[test]
    fn test_projection_pv_keywords_smoke() {
        // Each match arm in projection_pv_keywords gets exercised. Doesn't
        // verify FITS values exhaustively — just that the right keyword names
        // are produced for each parameterized projection.
        let names = |kws: Vec<WcsKeyword>| kws.iter().map(|k| k.name.clone()).collect::<Vec<_>>();

        assert_eq!(names(projection_pv_keywords(&Projection::tan())), Vec::<String>::new());

        assert_eq!(
            names(projection_pv_keywords(&Projection::sin_with_params(0.1, -0.2))),
            vec!["PV2_1", "PV2_2"],
        );
        // SIN with both params zero takes the default arm (no keywords).
        assert_eq!(
            names(projection_pv_keywords(&Projection::sin_with_params(0.0, 0.0))),
            Vec::<String>::new(),
        );

        assert_eq!(
            names(projection_pv_keywords(&Projection::azp(2.0, 30.0))),
            vec!["PV2_1", "PV2_2"],
        );
        assert_eq!(
            names(projection_pv_keywords(&Projection::szp(2.0, 0.0, 90.0))),
            vec!["PV2_1", "PV2_2", "PV2_3"],
        );

        // ZPN: only nonzero coefficients are emitted, named PV2_<index>.
        let zpn_kws = projection_pv_keywords(&Projection::zpn(vec![0.0, 1.0, 0.0, 0.01]));
        assert_eq!(names(zpn_kws), vec!["PV2_1", "PV2_3"]);

        assert_eq!(
            names(projection_pv_keywords(&Projection::air(45.0))),
            vec!["PV2_1"],
        );
        // CEA with lambda = 1.0 (the spec default) emits nothing.
        assert_eq!(
            names(projection_pv_keywords(&Projection::cea_with_lambda(1.0))),
            Vec::<String>::new(),
        );
        // CEA with non-default lambda emits PV2_1.
        assert_eq!(
            names(projection_pv_keywords(&Projection::cea_with_lambda(0.5))),
            vec!["PV2_1"],
        );
        assert_eq!(
            names(projection_pv_keywords(&Projection::cyp(1.0, 0.5))),
            vec!["PV2_1", "PV2_2"],
        );

        // The four conics share a single match arm — checking one is enough.
        assert_eq!(
            names(projection_pv_keywords(&Projection::cop(45.0, 15.0))),
            vec!["PV2_1", "PV2_2"],
        );

        assert_eq!(
            names(projection_pv_keywords(&Projection::bon(30.0))),
            vec!["PV2_1"],
        );
    }
}
