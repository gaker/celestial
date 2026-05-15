use std::collections::HashMap;

use celestial_core::Angle;

use crate::distortion::DistortionModel;
use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;
use crate::linear::LinearTransform;
use crate::spherical::{Projection, SphericalRotation};

use super::parse::{create_projection_from_code, parse_ctype, parse_matrix, parse_pv_params};
use super::wcs::{projection_to_code, CoordType, Wcs};

#[derive(Debug, Clone, PartialEq, Default)]
pub(super) enum MatrixSpec {
    #[default]
    None,
    Cd([[f64; 2]; 2]),
    PcCdelt {
        pc: [[f64; 2]; 2],
        cdelt: [f64; 2],
    },
}

#[derive(Debug, Clone, Default)]
pub struct WcsBuilder {
    crpix: Option<[f64; 2]>,
    crval: Option<[f64; 2]>,
    matrix: MatrixSpec,
    projection: Option<Projection>,
    lonpole: Option<f64>,
    latpole: Option<f64>,
    pv_params: HashMap<(u8, u8), f64>,
    coord_type: Option<CoordType>,
    proj_code: Option<String>,
    distortion: Option<DistortionModel>,
}

impl WcsBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn crpix(mut self, x: f64, y: f64) -> Self {
        self.crpix = Some([x, y]);
        self
    }

    pub fn crval(mut self, lon: f64, lat: f64) -> Self {
        self.crval = Some([lon, lat]);
        self
    }

    pub fn cd_matrix(mut self, cd: [[f64; 2]; 2]) -> Self {
        self.matrix = MatrixSpec::Cd(cd);
        self
    }

    pub fn pc_cdelt(mut self, pc: [[f64; 2]; 2], cdelt: [f64; 2]) -> Self {
        self.matrix = MatrixSpec::PcCdelt { pc, cdelt };
        self
    }

    pub fn projection(mut self, proj: Projection) -> Self {
        self.projection = Some(proj);
        self
    }

    pub fn lonpole(mut self, lonpole: f64) -> Self {
        self.lonpole = Some(lonpole);
        self
    }

    pub fn latpole(mut self, latpole: f64) -> Self {
        self.latpole = Some(latpole);
        self
    }

    pub fn pv(mut self, axis: u8, index: u8, value: f64) -> Self {
        self.pv_params.insert((axis, index), value);
        self
    }

    pub fn coord_type(mut self, coord_type: CoordType) -> Self {
        self.coord_type = Some(coord_type);
        self
    }

    pub fn proj_code(mut self, code: impl Into<String>) -> Self {
        self.proj_code = Some(code.into());
        self
    }

    pub fn distortion(mut self, distortion: DistortionModel) -> Self {
        self.distortion = Some(distortion);
        self
    }

    pub fn from_header(header: &impl KeywordProvider) -> WcsResult<Self> {
        let ctype1 = header.require_string("CTYPE1")?;
        let ctype2 = header.require_string("CTYPE2")?;

        let (prefix1, proj_code1) = parse_ctype(&ctype1)?;
        let (_prefix2, proj_code2) = parse_ctype(&ctype2)?;

        if proj_code1 != proj_code2 {
            return Err(WcsError::invalid_keyword(
                "CTYPE1/CTYPE2",
                format!(
                    "Mismatched projection codes: '{}' vs '{}'",
                    proj_code1, proj_code2
                ),
            ));
        }

        let coord_type = CoordType::from_ctype_prefix(prefix1);
        let proj_code = proj_code1.to_string();

        let crpix1 = header.require_float("CRPIX1")?;
        let crpix2 = header.require_float("CRPIX2")?;

        let crval1 = header.require_float("CRVAL1")?;
        let crval2 = header.require_float("CRVAL2")?;

        let matrix = parse_matrix(header)?;

        let lonpole = header.get_float("LONPOLE");
        let latpole = header.get_float("LATPOLE");

        let pv_params = parse_pv_params(header);

        let mut builder = Self::new()
            .crpix(crpix1, crpix2)
            .crval(crval1, crval2)
            .coord_type(coord_type)
            .proj_code(proj_code);

        builder.matrix = matrix;

        if let Some(lp) = lonpole {
            builder = builder.lonpole(lp);
        }
        if let Some(lp) = latpole {
            builder = builder.latpole(lp);
        }

        builder.pv_params = pv_params;

        Ok(builder)
    }

    pub fn validate(&self) -> WcsResult<()> {
        if self.crpix.is_none() {
            return Err(WcsError::missing_keyword("Missing CRPIX"));
        }
        if self.crval.is_none() {
            return Err(WcsError::missing_keyword("Missing CRVAL"));
        }
        if self.matrix == MatrixSpec::None {
            return Err(WcsError::missing_keyword(
                "Missing transformation matrix (CD or PC+CDELT)",
            ));
        }
        if self.projection.is_none() && self.proj_code.is_none() {
            return Err(WcsError::missing_keyword(
                "Missing projection (set projection or proj_code)",
            ));
        }
        Ok(())
    }

    pub fn build(self) -> WcsResult<Wcs> {
        self.validate()?;

        let crpix = self.crpix.unwrap();
        let crval = self.crval.unwrap();

        let linear = match &self.matrix {
            MatrixSpec::Cd(cd) => LinearTransform::from_cd(crpix, *cd)?,
            MatrixSpec::PcCdelt { pc, cdelt } => {
                LinearTransform::from_pc_cdelt(crpix, *pc, *cdelt)?
            }
            MatrixSpec::None => unreachable!("validate() ensures matrix is set"),
        };

        let projection = match self.projection {
            Some(proj) => proj,
            None => {
                let code = self.proj_code.as_ref().unwrap();
                create_projection_from_code(code, &self.pv_params)?
            }
        };

        let (_, theta_0) = projection.native_reference();

        let lonpole = self.lonpole.map(Angle::from_degrees);
        let latpole = self.latpole.map(Angle::from_degrees);

        let rotation = SphericalRotation::from_crval(
            Angle::from_degrees(crval[0]),
            Angle::from_degrees(crval[1]),
            Angle::from_degrees(theta_0),
            lonpole,
            latpole,
        )?;

        let coord_type = self.coord_type.unwrap_or_default();
        let proj_code = self
            .proj_code
            .unwrap_or_else(|| projection_to_code(&projection));

        Ok(Wcs::new(
            linear,
            projection,
            rotation,
            coord_type,
            proj_code,
            (crval[0], crval[1]),
            self.distortion,
        ))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::keyword::WcsKeywordValue;
    use crate::header::KeywordMap;

    fn minimal_builder() -> WcsBuilder {
        WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
    }

    fn tan_header() -> KeywordMap {
        let mut h = KeywordMap::new();
        h.set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", 0.001)
            .set_float("CD2_2", 0.001);
        h
    }

    #[test]
    fn test_builder_new_is_empty_and_chaining_sets_every_field() {
        // Default builder has no state set.
        let empty = WcsBuilder::new();
        assert!(empty.crpix.is_none());
        assert!(empty.crval.is_none());
        assert_eq!(empty.matrix, MatrixSpec::None);
        assert!(empty.projection.is_none());
        assert!(empty.lonpole.is_none());
        assert!(empty.latpole.is_none());
        assert!(empty.pv_params.is_empty());
        assert!(empty.coord_type.is_none());
        assert!(empty.proj_code.is_none());

        // Every setter populates its slot via the fluent API.  The string-form
        // proj_code path also exercises the `impl Into<String>` arm.
        let builder = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .lonpole(180.0)
            .latpole(90.0)
            .pv(2, 1, 0.5)
            .coord_type(CoordType::Equatorial)
            .proj_code(String::from("TAN"));

        assert_eq!(builder.crpix, Some([512.0, 512.0]));
        assert_eq!(builder.crval, Some([180.0, 45.0]));
        assert_eq!(builder.matrix, MatrixSpec::Cd([[0.001, 0.0], [0.0, 0.001]]));
        assert_eq!(builder.projection, Some(Projection::Tan));
        assert_eq!(builder.lonpole, Some(180.0));
        assert_eq!(builder.latpole, Some(90.0));
        assert_eq!(builder.pv_params.get(&(2, 1)), Some(&0.5));
        assert_eq!(builder.coord_type, Some(CoordType::Equatorial));
        assert_eq!(builder.proj_code, Some("TAN".to_string()));
    }

    #[test]
    fn test_matrix_setters_overwrite_each_other() {
        let pc = [[1.0, 0.0], [0.0, 1.0]];
        let cdelt = [0.001, 0.001];
        let cd = [[0.002, 0.0], [0.0, 0.002]];

        let cd_wins = WcsBuilder::new().pc_cdelt(pc, cdelt).cd_matrix(cd);
        assert_eq!(cd_wins.matrix, MatrixSpec::Cd(cd));

        let pc_wins = WcsBuilder::new().cd_matrix(cd).pc_cdelt(pc, cdelt);
        assert_eq!(pc_wins.matrix, MatrixSpec::PcCdelt { pc, cdelt });
    }

    #[test]
    fn test_from_header_cd_matrix() {
        let mut header = tan_header();
        header.set_float("CD1_1", -0.001).set_float("CD2_2", 0.001);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.crpix, Some([512.0, 512.0]));
        assert_eq!(builder.crval, Some([180.0, 45.0]));
        assert_eq!(builder.coord_type, Some(CoordType::Equatorial));
        assert_eq!(builder.proj_code, Some("TAN".to_string()));
        assert_eq!(
            builder.matrix,
            MatrixSpec::Cd([[-0.001, 0.0], [0.0, 0.001]])
        );
    }

    #[test]
    fn test_from_header_pc_cdelt() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "GLON-ARC")
            .set_string("CTYPE2", "GLAT-ARC")
            .set_float("CRPIX1", 256.0)
            .set_float("CRPIX2", 256.0)
            .set_float("CRVAL1", 90.0)
            .set_float("CRVAL2", 30.0)
            .set_float("CDELT1", -0.002)
            .set_float("CDELT2", 0.002)
            .set_float("PC1_1", 0.866)
            .set_float("PC1_2", -0.5)
            .set_float("PC2_1", 0.5)
            .set_float("PC2_2", 0.866);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.crpix, Some([256.0, 256.0]));
        assert_eq!(builder.crval, Some([90.0, 30.0]));
        assert_eq!(builder.coord_type, Some(CoordType::Galactic));
        assert_eq!(builder.proj_code, Some("ARC".to_string()));
        assert_eq!(
            builder.matrix,
            MatrixSpec::PcCdelt {
                pc: [[0.866, -0.5], [0.5, 0.866]],
                cdelt: [-0.002, 0.002]
            }
        );
    }

    #[test]
    fn test_from_header_pc_defaults_to_identity() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---SIN")
            .set_string("CTYPE2", "DEC--SIN")
            .set_float("CRPIX1", 100.0)
            .set_float("CRPIX2", 100.0)
            .set_float("CRVAL1", 0.0)
            .set_float("CRVAL2", 0.0)
            .set_float("CDELT1", 0.001)
            .set_float("CDELT2", 0.001);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(
            builder.matrix,
            MatrixSpec::PcCdelt {
                pc: [[1.0, 0.0], [0.0, 1.0]],
                cdelt: [0.001, 0.001]
            }
        );
    }

    #[test]
    fn test_from_header_missing_crpix() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", 0.001);

        let result = WcsBuilder::from_header(&header);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("CRPIX1"));
    }

    #[test]
    fn test_from_header_missing_matrix() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0);

        let result = WcsBuilder::from_header(&header);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("no transformation matrix"));
    }

    #[test]
    fn test_from_header_pv_params() {
        let mut header = tan_header();
        header
            .set_string("CTYPE1", "RA---SIN")
            .set_string("CTYPE2", "DEC--SIN")
            .set_float("PV2_1", 0.5)
            .set_float("PV2_2", 0.25);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.pv_params.get(&(2, 1)), Some(&0.5));
        assert_eq!(builder.pv_params.get(&(2, 2)), Some(&0.25));
    }

    #[test]
    fn test_from_header_lonpole_latpole() {
        let mut header = tan_header();
        header.set_float("LONPOLE", 180.0).set_float("LATPOLE", 45.0);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.lonpole, Some(180.0));
        assert_eq!(builder.latpole, Some(45.0));
    }

    #[test]
    fn test_from_header_mismatched_projection() {
        let mut header = tan_header();
        header.set_string("CTYPE2", "DEC--SIN");

        let err = WcsBuilder::from_header(&header).unwrap_err().to_string();
        assert!(err.contains("Mismatched"), "got: {}", err);
    }

    #[test]
    fn test_from_header_ecliptic_coords() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "ELON-CAR")
            .set_string("CTYPE2", "ELAT-CAR")
            .set_float("CRPIX1", 180.0)
            .set_float("CRPIX2", 90.0)
            .set_float("CRVAL1", 0.0)
            .set_float("CRVAL2", 0.0)
            .set_float("CD1_1", 1.0)
            .set_float("CD2_2", 1.0);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.coord_type, Some(CoordType::Ecliptic));
        assert_eq!(builder.proj_code, Some("CAR".to_string()));
    }

    #[test]
    fn test_from_header_partial_cd_matrix_defaults_missing_to_zero() {
        // tan_header() only sets CD1_1 and CD2_2; CD1_2 / CD2_1 must default to 0.
        let header = tan_header();
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(builder.matrix, MatrixSpec::Cd([[0.001, 0.0], [0.0, 0.001]]));
    }

    #[test]
    fn test_build_succeeds_with_minimal_valid_config() {
        // Also covers test_validate_returns_ok_for_valid_builder: the same
        // builder must both validate and build.
        let builder = minimal_builder().projection(Projection::tan());
        assert!(builder.clone().validate().is_ok());
        let wcs = builder.build().unwrap();
        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
        assert_eq!(wcs.projection_code(), "TAN");
    }

    #[test]
    fn test_build_fails_on_each_missing_required_field() {
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let cases: &[(WcsBuilder, &str)] = &[
            (
                WcsBuilder::new()
                    .crval(180.0, 45.0)
                    .cd_matrix(cd)
                    .projection(Projection::tan()),
                "Missing CRPIX",
            ),
            (
                WcsBuilder::new()
                    .crpix(512.0, 512.0)
                    .cd_matrix(cd)
                    .projection(Projection::tan()),
                "Missing CRVAL",
            ),
            (
                WcsBuilder::new()
                    .crpix(512.0, 512.0)
                    .crval(180.0, 45.0)
                    .projection(Projection::tan()),
                "Missing transformation matrix",
            ),
            (
                WcsBuilder::new()
                    .crpix(512.0, 512.0)
                    .crval(180.0, 45.0)
                    .cd_matrix(cd),
                "Missing projection",
            ),
        ];
        for (builder, needle) in cases {
            let err = builder.clone().build().unwrap_err().to_string();
            assert!(err.contains(needle), "expected {:?}, got: {}", needle, err);
        }
    }

    #[test]
    fn test_build_unsupported_projection() {
        // Unknown projection codes must surface the code in the error.
        let err = minimal_builder().proj_code("XYZ").build().unwrap_err().to_string();
        assert!(err.contains("XYZ"), "got: {}", err);
    }

    #[test]
    fn test_build_conic_pv_requirement() {
        // The four conics plus BON all require PV2_1 (theta_a / theta_1); the
        // build must report the missing keyword.  Providing the parameter then
        // makes the same builder succeed.
        for code in &["COP", "COE", "COD", "COO", "BON"] {
            let err = minimal_builder()
                .proj_code(*code)
                .build()
                .unwrap_err()
                .to_string();
            assert!(err.contains("PV2_1"), "{} missing PV2_1: {}", code, err);

            let wcs = minimal_builder()
                .proj_code(*code)
                .pv(2, 1, 45.0)
                .build()
                .unwrap();
            assert_eq!(wcs.projection_code(), *code);
        }
    }

    #[test]
    fn test_builder_full_roundtrip() {
        use crate::PixelCoord;

        let wcs = minimal_builder()
            .proj_code("TAN")
            .coord_type(CoordType::Equatorial)
            .build()
            .unwrap();

        let original = PixelCoord::new(300.0, 700.0);
        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        let tol = 1e-8;
        assert!((original.x() - recovered.x()).abs() < tol);
        assert!((original.y() - recovered.y()).abs() < tol);
    }

    #[test]
    fn test_builder_from_header_then_build() {
        let wcs = WcsBuilder::from_header(&tan_header()).unwrap().build().unwrap();

        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
        assert_eq!(wcs.projection_code(), "TAN");
        assert_eq!(wcs.coord_type(), CoordType::Equatorial);
    }

    #[test]
    fn test_build_with_pc_cdelt() {
        let wcs = WcsBuilder::new()
            .crpix(256.0, 256.0)
            .crval(90.0, 30.0)
            .pc_cdelt([[1.0, 0.0], [0.0, 1.0]], [0.002, 0.002])
            .proj_code("ARC")
            .build()
            .unwrap();

        assert_eq!(wcs.crpix(), [256.0, 256.0]);
        assert_eq!(wcs.projection_code(), "ARC");
    }

    #[test]
    fn test_build_with_lonpole_latpole() {
        let wcs = minimal_builder()
            .projection(Projection::tan())
            .lonpole(180.0)
            .latpole(45.0)
            .build()
            .unwrap();
        assert_eq!(wcs.crval(), (180.0, 45.0));
    }

    #[test]
    fn test_build_with_proj_code_round_trips_all() {
        // Sweep over every projection that takes (or omits) PV parameters,
        // plus exercise the Projection-enum path via STG.  The enum case
        // confirms `projection()` can substitute for `proj_code()`.
        type PvTriple = (u8, u8, f64);
        let cases: &[(&str, &[PvTriple])] = &[
            ("SIN", &[]),
            ("ARC", &[]),
            ("AZP", &[(2, 1, 2.0), (2, 2, 30.0)]),
            ("SZP", &[(2, 1, 2.0), (2, 2, 45.0), (2, 3, 60.0)]),
            ("ZPN", &[(2, 0, 0.0), (2, 1, 1.0), (2, 5, 0.001)]),
            ("AIR", &[(2, 1, 45.0)]),
            ("CEA", &[(2, 1, 0.5)]),
            ("COE", &[(2, 1, 45.0)]),
            ("COD", &[(2, 1, 30.0)]),
            ("COO", &[(2, 1, 60.0)]),
            ("BON", &[(2, 1, 45.0)]),
        ];
        for (code, pvs) in cases {
            let mut builder = minimal_builder().proj_code(*code);
            for (axis, idx, value) in *pvs {
                builder = builder.pv(*axis, *idx, *value);
            }
            let wcs = builder.build()
                .unwrap_or_else(|e| panic!("{} failed to build: {}", code, e));
            assert_eq!(wcs.projection_code(), *code);
        }

        // Projection::stg() (no proj_code) also infers the code from the enum.
        let wcs = minimal_builder().projection(Projection::stg()).build().unwrap();
        assert_eq!(wcs.projection_code(), "STG");
    }

    #[test]
    fn test_ctype_keywords_emit_correct_prefixes() {
        let cases: &[(Projection, CoordType, &str, &str)] = &[
            (Projection::tan(), CoordType::Equatorial, "RA---TAN", "DEC--TAN"),
            (Projection::sin(), CoordType::Galactic, "GLON-SIN", "GLAT-SIN"),
        ];
        for (proj, coord_type, expected1, expected2) in cases {
            let wcs = minimal_builder()
                .projection(proj.clone())
                .coord_type(*coord_type)
                .build()
                .unwrap();

            let keywords = wcs.to_keywords();
            let find = |name: &str| {
                keywords
                    .iter()
                    .find(|k| k.name == name)
                    .unwrap_or_else(|| panic!("{} not found", name))
                    .value
                    .clone()
            };
            assert_eq!(find("CTYPE1"), WcsKeywordValue::String((*expected1).to_string()));
            assert_eq!(find("CTYPE2"), WcsKeywordValue::String((*expected2).to_string()));
        }
    }

}
