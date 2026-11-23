use std::collections::HashMap;

use crate::distortion::DistortionModel;
use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;
use crate::linear::LinearTransform;
use crate::spherical::{Projection, SphericalRotation};

use super::distortion_keywords::{check_distortion, header_distortion};
use super::fiducial::Fiducial;
use super::parse::{axis_pair, parse_ctype_pair, parse_matrix, parse_pv_params};
use super::units::{celestial_units, in_degrees};
use super::wcs::{CoordType, Wcs};
use super::world_axes::LongitudeFirst;

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

// One slot, so that the projection and its code cannot disagree.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum ProjectionSpec {
    Code(String),
    Given(Projection),
}

#[derive(Debug, Clone, Default)]
pub struct WcsBuilder {
    crpix: Option<[f64; 2]>,
    crval: Option<[f64; 2]>,
    matrix: MatrixSpec,
    projection: Option<ProjectionSpec>,
    lonpole: Option<f64>,
    latpole: Option<f64>,
    pv_params: HashMap<(u8, u8), f64>,
    coord_type: Option<CoordType>,
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
        self.projection = Some(ProjectionSpec::Given(proj));
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
        self.projection = Some(ProjectionSpec::Code(code.into()));
        self
    }

    pub fn distortion(mut self, distortion: DistortionModel) -> Self {
        self.distortion = Some(distortion);
        self
    }

    pub fn from_header(header: &impl KeywordProvider) -> WcsResult<Self> {
        let header = &LongitudeFirst::new(header)?;
        let ctype1 = header.require_string("CTYPE1")?;
        let ctype2 = header.require_string("CTYPE2")?;
        let (prefix, code, suffix) = parse_ctype_pair(&ctype1, &ctype2)?;
        let units = celestial_units(header)?;
        let crpix = axis_pair(header, "CRPIX", 0.0);
        let crval = in_degrees(axis_pair(header, "CRVAL", 0.0), units);
        let distortion = header_distortion(header, code, suffix, crpix)?;
        let tpv = matches!(distortion, Some(DistortionModel::Tpv(_)));
        let mut builder = Self::new()
            .crpix(crpix[0], crpix[1])
            .crval(crval[0], crval[1])
            .coord_type(CoordType::from_ctype_prefix(prefix))
            .proj_code(if tpv { "TAN" } else { code });
        builder.matrix = parse_matrix(header, units);
        builder.distortion = distortion;
        Ok(builder.with_pole_and_pv(header, tpv))
    }

    fn with_pole_and_pv(mut self, header: &impl KeywordProvider, tpv: bool) -> Self {
        self.lonpole = header.get_float("LONPOLE");
        self.latpole = header.get_float("LATPOLE");
        // TPV's PVi_m are the terms of its polynomial, already read.
        if !tpv {
            self.pv_params = parse_pv_params(header);
        }
        self
    }

    pub fn validate(&self) -> WcsResult<()> {
        self.reference_point()?;
        if self.matrix == MatrixSpec::None {
            return Err(missing_matrix());
        }
        if self.projection.is_none() {
            return Err(missing_projection());
        }
        Ok(())
    }

    pub fn build(self) -> WcsResult<Wcs> {
        self.validate()?;
        let (crpix, crval) = self.reference_point()?;
        let linear = self.linear_transform(crpix)?;
        let projection = self.resolved_projection()?;
        let fiducial = Fiducial::from_pv(&projection, &self.pv_params)?;
        check_distortion(self.distortion.as_ref(), crpix, &projection, &fiducial)?;
        let rotation = self.rotation(crval, &fiducial)?;
        let coord_type = self.coord_type.unwrap_or_default();
        let crval = (crval[0], crval[1]);
        let wcs = Wcs::new(
            linear,
            projection,
            rotation,
            coord_type,
            crval,
            self.distortion,
        );
        Ok(wcs.with_fiducial(fiducial))
    }

    fn reference_point(&self) -> WcsResult<([f64; 2], [f64; 2])> {
        let crpix = self
            .crpix
            .ok_or_else(|| WcsError::missing_keyword("Missing CRPIX"))?;
        let crval = self
            .crval
            .ok_or_else(|| WcsError::missing_keyword("Missing CRVAL"))?;
        Ok((crpix, crval))
    }

    fn linear_transform(&self, crpix: [f64; 2]) -> WcsResult<LinearTransform> {
        match &self.matrix {
            MatrixSpec::Cd(cd) => LinearTransform::from_cd(crpix, *cd),
            MatrixSpec::PcCdelt { pc, cdelt } => LinearTransform::from_pc_cdelt(crpix, *pc, *cdelt),
            MatrixSpec::None => Err(missing_matrix()),
        }
    }

    fn rotation(&self, crval: [f64; 2], fiducial: &Fiducial) -> WcsResult<SphericalRotation> {
        let lonpole = self.pole_parameter(3, self.lonpole);
        let latpole = self.pole_parameter(4, self.latpole);
        let (phi_0, theta_0) = (fiducial.phi_0(), fiducial.theta_0());
        SphericalRotation::from_crval(crval[0], crval[1], phi_0, theta_0, lonpole, latpole)
    }

    // Paper II §2.4: PV2_m are the projection's parameters, which an explicit
    // Projection already holds; beside one they would go unread.
    fn resolved_projection(&self) -> WcsResult<Projection> {
        match &self.projection {
            Some(ProjectionSpec::Code(code)) => {
                Projection::from_code(code, |m| self.pv_params.get(&(2, m)).copied())
            }
            Some(ProjectionSpec::Given(_)) if self.pv_params.keys().any(|&(axis, _)| axis == 2) => {
                Err(WcsError::invalid_parameter(
                    "PV2_m set beside an explicit Projection, which holds its own parameters",
                ))
            }
            Some(ProjectionSpec::Given(projection)) => Ok(projection.clone()),
            None => Err(missing_projection()),
        }
    }

    // Paper II §2.6: PV1_3 and PV1_4 encapsulate LONPOLE and LATPOLE, and take
    // precedence over them.
    fn pole_parameter(&self, index: u8, keyword: Option<f64>) -> Option<f64> {
        let encapsulated = self.pv_params.get(&(1, index)).copied();
        encapsulated.or(keyword)
    }
}

fn missing_matrix() -> WcsError {
    WcsError::missing_keyword("Missing transformation matrix (CD or PC+CDELT)")
}

fn missing_projection() -> WcsError {
    WcsError::missing_keyword("Missing projection (set projection or proj_code)")
}

#[cfg(test)]
mod tests {
    use super::super::keyword::{keyword_map, WcsKeywordValue};
    use super::*;
    use crate::header::KeywordMap;
    use celestial_core::constants::PI;

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
        assert_eq!(builder.lonpole, Some(180.0));
        assert_eq!(builder.latpole, Some(90.0));
        assert_eq!(builder.pv_params.get(&(2, 1)), Some(&0.5));
        assert_eq!(builder.coord_type, Some(CoordType::Equatorial));
        let tan = ProjectionSpec::Code("TAN".to_string());
        assert_eq!(builder.projection, Some(tan));
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
        let code = ProjectionSpec::Code("TAN".to_string());
        assert_eq!(builder.projection, Some(code));
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
        let code = ProjectionSpec::Code("ARC".to_string());
        assert_eq!(builder.projection, Some(code));
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
    fn test_from_header_crota2_rotates_cdelt_axes() {
        // Paper II eq. 189 with ρ = CROTA2 = 30°, whose sine rounds to 0.5 and
        // cosine to 0.8660254037844386.
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 100.0)
            .set_float("CRPIX2", 100.0)
            .set_float("CRVAL1", 0.0)
            .set_float("CRVAL2", 0.0)
            .set_float("CDELT1", -0.001)
            .set_float("CDELT2", 0.002)
            .set_float("CROTA2", 30.0);

        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(
            builder.matrix,
            MatrixSpec::Cd([
                [-0.0008660254037844386, -0.001],
                [-0.0005, 0.0017320508075688772]
            ])
        );

        // At ρ = 90° the cosine is exactly 0, so the diagonal vanishes.
        header.set_float("CROTA2", 90.0);
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(
            builder.matrix,
            MatrixSpec::Cd([[0.0, -0.002], [-0.001, 0.0]])
        );

        // Paper II §6.2 bars CROTA2 beside PC, so a PC keyword marks the
        // header as the current formalism and CROTA2 is ignored.
        header.set_float("PC1_2", 0.25);
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(
            builder.matrix,
            MatrixSpec::PcCdelt {
                pc: [[1.0, 0.25], [0.0, 1.0]],
                cdelt: [-0.001, 0.002]
            }
        );
    }

    // Paper I §2.4: CRVAL and CRPIX default to 0, CDELT to 1 and PC to δij.
    #[test]
    fn test_from_header_takes_paper_i_defaults() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN");
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(builder.crpix, Some([0.0, 0.0]));
        assert_eq!(builder.crval, Some([0.0, 0.0]));
        assert_eq!(
            builder.matrix,
            MatrixSpec::PcCdelt {
                pc: [[1.0, 0.0], [0.0, 1.0]],
                cdelt: [1.0, 1.0]
            }
        );
        builder.build().unwrap();

        header.set_float("CDELT1", -1e-3).set_float("PC1_2", 0.25);
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(
            builder.matrix,
            MatrixSpec::PcCdelt {
                pc: [[1.0, 0.25], [0.0, 1.0]],
                cdelt: [-1e-3, 1.0]
            }
        );
    }

    // Paper I Table 6 defines arcmin, arcsec and mas in degrees, and Paper II
    // §2 evaluates celestial axes in degrees.
    #[test]
    fn test_from_header_scales_angular_units_to_degrees() {
        let mut header = tan_header();
        header
            .set_string("CUNIT1", "arcsec  ")
            .set_string("CUNIT2", "arcmin")
            .set_float("CRVAL1", 445500.0)
            .set_float("CRVAL2", 2250.0)
            .set_float("CD1_1", -4.5)
            .set_float("CD2_2", 0.75);
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(builder.crval, Some([123.75, 37.5]));
        assert_eq!(
            builder.matrix,
            MatrixSpec::Cd([[-0.00125, 0.0], [0.0, 0.0125]])
        );

        header
            .set_string("CUNIT1", "rad")
            .set_string("CUNIT2", "mas")
            .set_float("CRVAL1", PI)
            .set_float("CRVAL2", 135_000_000.0);
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(builder.crval, Some([180.0, 37.5]));

        header.set_string("CUNIT2", "m");
        let result = WcsBuilder::from_header(&header);
        assert!(
            matches!(&result, Err(WcsError::InvalidKeyword { keyword, .. }) if keyword == "CUNIT2"),
            "{result:?}"
        );
    }

    // CDELT carries the units; PC is a pure number.
    #[test]
    fn test_from_header_scales_cdelt_but_not_pc() {
        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_string("CUNIT1", "arcsec")
            .set_string("CUNIT2", "arcsec")
            .set_float("CDELT1", -4.5)
            .set_float("CDELT2", 4.5)
            .set_float("PC1_2", 0.25)
            .set_float("CROTA2", 30.0);
        let builder = WcsBuilder::from_header(&header).unwrap();
        assert_eq!(
            builder.matrix,
            MatrixSpec::PcCdelt {
                pc: [[1.0, 0.25], [0.0, 1.0]],
                cdelt: [-0.00125, 0.00125]
            }
        );
    }

    // A latitude-first header and its longitude-first twin with the pixel
    // axes transposed describe the same sky: CRPIX and the matrix columns
    // trade places and nothing else changes. Paper I §2.1.2 defaults PCi_j
    // to δij, and Paper II §6.1 rotates by the latitude axis's CROTA.
    #[test]
    fn test_from_header_reads_latitude_first_axes() {
        let mut lat_first = KeywordMap::new();
        lat_first
            .set_string("CTYPE1", "DEC--SIN")
            .set_string("CTYPE2", "RA---SIN")
            .set_string("CUNIT1", "arcsec")
            .set_float("CRPIX1", 100.0)
            .set_float("CRPIX2", 200.0)
            .set_float("CRVAL1", 135000.0)
            .set_float("CRVAL2", 123.5)
            .set_float("CDELT1", 3.6)
            .set_float("CDELT2", -0.001)
            .set_float("CROTA1", 30.0)
            .set_float("PV1_1", 0.25)
            .set_float("PV2_3", 170.0);
        let mut twin = KeywordMap::new();
        twin.set_string("CTYPE1", "RA---SIN")
            .set_string("CTYPE2", "DEC--SIN")
            .set_string("CUNIT2", "arcsec")
            .set_float("CRPIX1", 200.0)
            .set_float("CRPIX2", 100.0)
            .set_float("CRVAL1", 123.5)
            .set_float("CRVAL2", 135000.0)
            .set_float("CDELT1", -0.001)
            .set_float("CDELT2", 3.6)
            .set_float("CROTA2", 30.0)
            .set_float("PV1_3", 170.0)
            .set_float("PV2_1", 0.25);
        assert_transposed_twins(&lat_first, &twin);

        lat_first.set_float("PC1_2", 0.25);
        twin.set_float("PC2_1", 0.25);
        assert_transposed_twins(&lat_first, &twin);
    }

    #[track_caller]
    fn assert_transposed_twins(header: &KeywordMap, twin: &KeywordMap) {
        let builder = WcsBuilder::from_header(header).unwrap();
        let twin = WcsBuilder::from_header(twin).unwrap();
        let swap = |m: [[f64; 2]; 2]| m.map(|row| [row[1], row[0]]);
        let expected = match twin.matrix {
            MatrixSpec::Cd(cd) => MatrixSpec::Cd(swap(cd)),
            MatrixSpec::PcCdelt { pc, cdelt } => MatrixSpec::PcCdelt {
                pc: swap(pc),
                cdelt,
            },
            MatrixSpec::None => MatrixSpec::None,
        };
        assert_eq!(builder.matrix, expected);
        assert_eq!(builder.crpix, twin.crpix.map(|[x, y]| [y, x]));
        assert_eq!(builder.crval, twin.crval);
        assert_eq!(builder.pv_params, twin.pv_params);
        assert_eq!(builder.coord_type, twin.coord_type);
        assert_eq!(builder.projection, twin.projection);
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
        header
            .set_float("LONPOLE", 180.0)
            .set_float("LATPOLE", 45.0);

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

    // Paper II §3: RA/DEC, xLON/xLAT and yzLN/yzLT pairs, x or yz alike.
    #[test]
    fn test_from_header_requires_a_celestial_pair() {
        let mut header = tan_header();
        let pairs = [
            ("GLON-TAN", "GLAT-TAN", true),
            ("HPLN-TAN", "HPLT-TAN", true),
            ("DEC--TAN", "RA---TAN", true),
            ("RA---TAN", "GLAT-TAN", false),
            ("GLON-TAN", "ELAT-TAN", false),
            ("HPLN-TAN", "HGLT-TAN", false),
            ("GLON-TAN", "GLON-TAN", false),
            ("FOO--TAN", "BAR--TAN", false),
        ];
        for (ctype1, ctype2, is_pair) in pairs {
            header
                .set_string("CTYPE1", ctype1)
                .set_string("CTYPE2", ctype2);
            let result = WcsBuilder::from_header(&header);
            let refused = matches!(&result, Err(WcsError::InvalidKeyword { .. }));
            assert_eq!(result.is_ok(), is_pair, "{ctype1}, {ctype2}: {result:?}");
            assert_eq!(refused, !is_pair, "{ctype1}, {ctype2}: {result:?}");
        }
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
        let code = ProjectionSpec::Code("CAR".to_string());
        assert_eq!(builder.projection, Some(code));
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
        // The same builder must both validate and build.
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
        let err = minimal_builder()
            .proj_code("XYZ")
            .build()
            .unwrap_err()
            .to_string();
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

    // The TAN WCS of wcs::tests::test_wcs_simple_tan_pipeline, built here
    // from its code; it gives that test's pinned values.
    #[test]
    fn test_builder_full_roundtrip() {
        let wcs = minimal_builder()
            .proj_code("TAN")
            .coord_type(CoordType::Equatorial)
            .build()
            .unwrap();
        let sky = (179.69920249673768, 45.187604542264104);
        assert_eq!(wcs.pix2world(300.0, 700.0).unwrap(), sky);
        let back = (300.0000000000147, 699.9999999999836);
        assert_eq!(wcs.world2pix(sky.0, sky.1).unwrap(), back);
    }

    #[test]
    fn test_builder_from_header_then_build() {
        let wcs = WcsBuilder::from_header(&tan_header())
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
        assert_eq!(wcs.projection_code(), "TAN");
        assert_eq!(wcs.coord_type(), CoordType::Equatorial);
    }

    // Paper I eq. 1 scales row i of PC by CDELTi. Every product here is a
    // power-of-two scaling, so PC and CDELT give exactly this CD, and column
    // scaling would give [[0.001, 0.001], [−0.0005, 0.002]].
    #[test]
    fn test_build_with_pc_cdelt() {
        let base = || WcsBuilder::new().crpix(256.0, 256.0).crval(90.0, 30.0);
        let from_pc = base()
            .pc_cdelt([[0.5, 0.25], [-0.25, 0.5]], [0.002, 0.004])
            .proj_code("ARC")
            .build()
            .unwrap();
        let from_cd = base()
            .cd_matrix([[0.001, 0.0005], [-0.001, 0.002]])
            .proj_code("ARC")
            .build()
            .unwrap();

        assert_eq!(from_pc.crpix(), [256.0, 256.0]);
        assert_eq!(from_pc.projection_code(), "ARC");
        let sky = from_cd.pix2world(200.0, 300.0).unwrap();
        assert_eq!(from_pc.pix2world(200.0, 300.0).unwrap(), sky);
        assert_eq!(
            from_pc.world2pix(sky.0, sky.1).unwrap(),
            from_cd.world2pix(sky.0, sky.1).unwrap()
        );
    }

    // LONPOLE 90 turns the TAN sky a quarter turn about the reference point,
    // and LATPOLE −90 puts the CAR native pole at δp = −45° rather than +45°
    // (Paper II eqs. 8 to 10), so neither row is where the defaults put it.
    // Against 60 digits pix2world is within 1.48 ulp and world2pix within
    // 39 ulp, at most 0.46 of the summed half-ulp rounding of the steps; most
    // of that is θ's rounding scaled by 1/CDELT.
    #[test]
    fn test_build_with_lonpole_latpole() {
        let tan = minimal_builder()
            .projection(Projection::tan())
            .lonpole(90.0)
            .build()
            .unwrap();
        let car = minimal_builder()
            .projection(Projection::car())
            .latpole(-90.0)
            .build()
            .unwrap();
        let rows = [
            (
                tan,
                (179.92916572878087, 45.099978005641276),
                (611.9999999999875, 561.9999999999869),
            ),
            (
                car,
                (179.85870209276425, 44.949912809713744),
                (612.0000000000041, 561.999999999994),
            ),
        ];
        for (wcs, sky, back) in rows {
            assert_eq!(wcs.pix2world(612.0, 562.0).unwrap(), sky);
            assert_eq!(wcs.world2pix(sky.0, sky.1).unwrap(), back);
        }
    }

    fn car_header(keywords: &[(&str, f64)]) -> KeywordMap {
        let mut h = KeywordMap::new();
        h.set_string("CTYPE1", "RA---CAR")
            .set_string("CTYPE2", "DEC--CAR")
            .set_float("CRPIX1", 100.0)
            .set_float("CRPIX2", 100.0)
            .set_float("CRVAL1", 150.0)
            .set_float("CRVAL2", 40.0)
            .set_float("CDELT1", 0.05)
            .set_float("CDELT2", 0.05);
        for &(key, value) in keywords {
            h.set_float(key, value);
        }
        h
    }

    fn header_wcs(header: &KeywordMap) -> Wcs {
        WcsBuilder::from_header(header).unwrap().build().unwrap()
    }

    #[test]
    fn test_pv1_3_and_pv1_4_take_precedence_over_lonpole_and_latpole() {
        // Paper II §2.6. With LONPOLE 0 the δp of eq. 8 has the roots ±50°,
        // between which LATPOLE chooses.
        let cases = [
            (
                vec![("LONPOLE", 0.0), ("PV1_3", 30.0)],
                vec![("LONPOLE", 30.0)],
            ),
            (
                vec![("LONPOLE", 0.0), ("LATPOLE", 90.0), ("PV1_4", -90.0)],
                vec![("LONPOLE", 0.0), ("LATPOLE", -90.0)],
            ),
        ];
        for (encapsulated, plain) in cases {
            let got = header_wcs(&car_header(&encapsulated)).pix2world(130.0, 80.0);
            let expected = header_wcs(&car_header(&plain)).pix2world(130.0, 80.0);
            assert_eq!(got.unwrap(), expected.unwrap(), "{encapsulated:?}");
        }
    }

    #[test]
    fn test_pv1_1_and_pv1_2_set_the_native_point_of_crval() {
        // Paper II §2.5. CAR puts (φ0, θ0) = (30°, 22.5°) at (x, y) = (φ0, θ0),
        // 600 and 450 pixels of 0.05° from CRPIX.
        let wcs = header_wcs(&car_header(&[("PV1_1", 30.0), ("PV1_2", 22.5)]));
        assert_eq!(wcs.pix2world(700.0, 550.0).unwrap(), (150.0, 40.0));
    }

    #[test]
    fn test_pv1_0_offsets_the_fiducial_point_onto_crpix() {
        // Paper II §2.5: a nonzero PV1_0 takes the projection of (φ0, θ0) off
        // (x, y). For CAR and (0, 22.5°) that moves every pixel 450 rows and
        // puts CRVAL at CRPIX.
        let fiducial = [("PV1_1", 0.0), ("PV1_2", 22.5)];
        let offset = header_wcs(&car_header(&[("PV1_0", 1.0), fiducial[0], fiducial[1]]));
        let plain = header_wcs(&car_header(&fiducial));
        let shifted = plain.pix2world(130.0, 530.0).unwrap();
        assert_eq!(offset.pix2world(130.0, 80.0).unwrap(), shifted);
        assert_eq!(offset.pix2world(100.0, 100.0).unwrap(), (150.0, 40.0));
        let (x, y) = plain.world2pix(150.0, 40.0).unwrap();
        assert_eq!(offset.world2pix(150.0, 40.0).unwrap(), (x, y - 450.0));
    }

    #[test]
    fn test_to_keywords_keeps_the_fiducial_point() {
        let fiducial = [("PV1_0", 1.0), ("PV1_1", 30.0), ("PV1_2", 22.5)];
        let wcs = header_wcs(&car_header(&fiducial));
        let header = keyword_map(wcs.to_keywords());
        for (name, value) in fiducial {
            assert_eq!(header.get_float(name), Some(value), "{name}");
        }
        let reread = header_wcs(&header);
        let expected = wcs.pix2world(130.0, 80.0).unwrap();
        assert_eq!(reread.pix2world(130.0, 80.0).unwrap(), expected);
    }

    #[test]
    fn test_projection_and_its_code_cannot_disagree() {
        let sin = minimal_builder()
            .proj_code("TAN")
            .projection(Projection::sin())
            .build()
            .unwrap();
        assert_eq!(sin.projection_code(), "SIN");
        let tan = minimal_builder()
            .projection(Projection::sin())
            .proj_code("TAN")
            .build()
            .unwrap();
        let plain_tan = minimal_builder().proj_code("TAN").build().unwrap();
        assert_eq!(
            tan.pix2world(100.0, 900.0).unwrap(),
            plain_tan.pix2world(100.0, 900.0).unwrap()
        );

        // An explicit Projection carries its own parameters.
        let beside = minimal_builder()
            .projection(Projection::sin())
            .pv(2, 1, 0.5)
            .build();
        assert!(beside.is_err());
    }

    #[test]
    fn test_build_with_proj_code_round_trips_all() {
        // Sweep over every projection that takes (or omits) PV parameters:
        // the code and its PV2_m must resolve to the projection built
        // directly, each parameter in its place. The enum case confirms
        // `projection()` can substitute for `proj_code()`.
        type PvTriple = (u8, u8, f64);
        let zpn = Projection::zpn(vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.001]);
        let cases: &[(&str, &[PvTriple], Projection)] = &[
            ("SIN", &[], Projection::sin()),
            ("ARC", &[], Projection::arc()),
            (
                "AZP",
                &[(2, 1, 2.0), (2, 2, 30.0)],
                Projection::azp(2.0, 30.0),
            ),
            (
                "SZP",
                &[(2, 1, 2.0), (2, 2, 45.0), (2, 3, 60.0)],
                Projection::szp(2.0, 45.0, 60.0),
            ),
            ("ZPN", &[(2, 0, 0.0), (2, 1, 1.0), (2, 5, 0.001)], zpn),
            ("AIR", &[(2, 1, 45.0)], Projection::air(45.0)),
            ("CEA", &[(2, 1, 0.5)], Projection::cea_with_lambda(0.5)),
            ("COE", &[(2, 1, 45.0)], Projection::coe(45.0, 0.0)),
            (
                "COD",
                &[(2, 1, 30.0), (2, 2, 5.0)],
                Projection::cod(30.0, 5.0),
            ),
            ("COO", &[(2, 1, 60.0)], Projection::coo(60.0, 0.0)),
            ("BON", &[(2, 1, 45.0)], Projection::bon(45.0)),
        ];
        for (code, pvs, expected) in cases {
            let mut builder = minimal_builder().proj_code(*code);
            for (axis, idx, value) in *pvs {
                builder = builder.pv(*axis, *idx, *value);
            }
            assert_eq!(builder.resolved_projection().unwrap(), *expected);
            let wcs = builder
                .build()
                .unwrap_or_else(|e| panic!("{} failed to build: {}", code, e));
            assert_eq!(wcs.projection_code(), *code);
        }

        // Projection::stg() (no proj_code) also infers the code from the enum.
        let wcs = minimal_builder()
            .projection(Projection::stg())
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "STG");
    }

    #[test]
    fn test_ctype_keywords_emit_correct_prefixes() {
        let cases: &[(Projection, CoordType, &str, &str)] = &[
            (
                Projection::tan(),
                CoordType::Equatorial,
                "RA---TAN",
                "DEC--TAN",
            ),
            (
                Projection::sin(),
                CoordType::Galactic,
                "GLON-SIN",
                "GLAT-SIN",
            ),
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
            assert_eq!(
                find("CTYPE1"),
                WcsKeywordValue::String((*expected1).to_string())
            );
            assert_eq!(
                find("CTYPE2"),
                WcsKeywordValue::String((*expected2).to_string())
            );
        }
    }
}
