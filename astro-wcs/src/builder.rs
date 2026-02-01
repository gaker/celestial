use std::collections::HashMap;

use astro_core::Angle;

use crate::coordinate::{CelestialCoord, IntermediateCoord, PixelCoord};
use crate::distortion::DistortionModel;
use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;
use crate::linear::LinearTransform;
use crate::spherical::{Projection, SphericalRotation};

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
        let pixel = self.apply_sip_forward(pixel);
        let intermediate = self.linear.pixel_to_intermediate(pixel);
        let intermediate = self.apply_tpv_tnx_forward(intermediate);
        let native = self.projection.deproject(intermediate)?;
        self.rotation.native_to_celestial(native)
    }

    pub fn celestial_to_pixel(&self, celestial: CelestialCoord) -> WcsResult<PixelCoord> {
        let native = self.rotation.celestial_to_native(celestial)?;
        let intermediate = self.projection.project(native)?;
        let intermediate = self.apply_tpv_tnx_inverse(intermediate)?;
        let pixel = self.linear.intermediate_to_pixel(intermediate);
        self.apply_sip_inverse(pixel)
    }

    fn apply_sip_forward(&self, pixel: PixelCoord) -> PixelCoord {
        match &self.distortion {
            Some(DistortionModel::Sip(sip)) => {
                let (x, y) = sip.apply(pixel.x(), pixel.y());
                PixelCoord::new(x, y)
            }
            _ => pixel,
        }
    }

    fn apply_tpv_tnx_forward(&self, intermediate: IntermediateCoord) -> IntermediateCoord {
        match &self.distortion {
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
        }
    }

    fn apply_tpv_tnx_inverse(
        &self,
        intermediate: IntermediateCoord,
    ) -> WcsResult<IntermediateCoord> {
        match &self.distortion {
            Some(DistortionModel::Tpv(tpv)) => {
                let (x, y) = tpv
                    .as_ref()
                    .apply_inverse(intermediate.x_deg(), intermediate.y_deg())?;
                Ok(IntermediateCoord::new(x, y))
            }
            Some(DistortionModel::Tnx(tnx)) => {
                let (x, y) = tnx.apply_inverse(intermediate.x_deg(), intermediate.y_deg())?;
                Ok(IntermediateCoord::new(x, y))
            }
            _ => Ok(intermediate),
        }
    }

    fn apply_sip_inverse(&self, pixel: PixelCoord) -> WcsResult<PixelCoord> {
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

    #[inline]
    pub fn projection(&self) -> &Projection {
        &self.projection
    }

    #[inline]
    pub fn rotation(&self) -> &SphericalRotation {
        &self.rotation
    }

    #[inline]
    pub fn linear(&self) -> &LinearTransform {
        &self.linear
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

fn ctype_prefixes(coord_type: &CoordType) -> (&'static str, &'static str) {
    match coord_type {
        CoordType::Equatorial => ("RA", "DEC"),
        CoordType::Galactic => ("GLON", "GLAT"),
        CoordType::Ecliptic => ("ELON", "ELAT"),
        CoordType::Helioecliptic => ("HLON", "HLAT"),
        CoordType::Supergalactic => ("SLON", "SLAT"),
        CoordType::Generic => ("XLON", "XLAT"),
    }
}

fn format_ctype(prefix: &str, proj_code: &str) -> String {
    let padding_len = 4 - prefix.len();
    let dashes = "-".repeat(padding_len + 1);
    format!("{}{}{}", prefix, dashes, proj_code)
}

fn default_lonpole(_coord_type: &CoordType, crval_lat: f64, projection: &Projection) -> f64 {
    let (_, theta_0) = projection.native_reference();
    if crval_lat >= theta_0 {
        0.0
    } else {
        180.0
    }
}

fn projection_pv_keywords(projection: &Projection) -> Vec<WcsKeyword> {
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
        Projection::Cop { theta_a }
        | Projection::Coe { theta_a }
        | Projection::Cod { theta_a }
        | Projection::Coo { theta_a } => {
            vec![WcsKeyword::real("PV2_1", *theta_a)]
        }
        Projection::Bon { theta_1 } => {
            vec![WcsKeyword::real("PV2_1", *theta_1)]
        }
        _ => Vec::new(),
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum WcsKeywordValue {
    Real(f64),
    Integer(i64),
    String(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct WcsKeyword {
    pub name: String,
    pub value: WcsKeywordValue,
}

impl WcsKeyword {
    pub fn real(name: impl Into<String>, value: f64) -> Self {
        Self {
            name: name.into(),
            value: WcsKeywordValue::Real(value),
        }
    }

    pub fn integer(name: impl Into<String>, value: i64) -> Self {
        Self {
            name: name.into(),
            value: WcsKeywordValue::Integer(value),
        }
    }

    pub fn string(name: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            value: WcsKeywordValue::String(value.into()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
enum MatrixSpec {
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

fn create_projection_from_code(
    code: &str,
    pv_params: &HashMap<(u8, u8), f64>,
) -> WcsResult<Projection> {
    match code {
        "TAN" => Ok(Projection::tan()),
        "SIN" => {
            let xi = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let eta = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            if xi == 0.0 && eta == 0.0 {
                Ok(Projection::sin())
            } else {
                Ok(Projection::sin_with_params(xi, eta))
            }
        }
        "ARC" => Ok(Projection::arc()),
        "STG" => Ok(Projection::stg()),
        "ZEA" => Ok(Projection::zea()),
        "AZP" => {
            let mu = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let gamma = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            Ok(Projection::azp(mu, gamma))
        }
        "SZP" => {
            let mu = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let phi_c = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            let theta_c = pv_params.get(&(2, 3)).copied().unwrap_or(90.0);
            Ok(Projection::szp(mu, phi_c, theta_c))
        }
        "ZPN" => {
            let mut coeffs = Vec::new();
            for i in 0..=20 {
                if let Some(&val) = pv_params.get(&(2, i)) {
                    while coeffs.len() < i as usize {
                        coeffs.push(0.0);
                    }
                    coeffs.push(val);
                }
            }
            if coeffs.is_empty() {
                coeffs.push(0.0);
                coeffs.push(1.0);
            }
            Ok(Projection::zpn(coeffs))
        }
        "AIR" => {
            let theta_b = pv_params.get(&(2, 1)).copied().unwrap_or(90.0);
            Ok(Projection::air(theta_b))
        }
        "CAR" => Ok(Projection::car()),
        "MER" => Ok(Projection::mer()),
        "CEA" => {
            let lambda = pv_params.get(&(2, 1)).copied().unwrap_or(1.0);
            Ok(Projection::cea_with_lambda(lambda))
        }
        "CYP" => {
            let mu = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let lambda = pv_params.get(&(2, 2)).copied().unwrap_or(1.0);
            Ok(Projection::cyp(mu, lambda))
        }
        "SFL" => Ok(Projection::sfl()),
        "PAR" => Ok(Projection::par()),
        "MOL" => Ok(Projection::mol()),
        "AIT" => Ok(Projection::ait()),
        "COP" => {
            let theta_a = pv_params.get(&(2, 1)).ok_or_else(|| {
                WcsError::missing_keyword("COP projection requires PV2_1 (theta_a)")
            })?;
            Ok(Projection::cop(*theta_a))
        }
        "COE" => {
            let theta_a = pv_params.get(&(2, 1)).ok_or_else(|| {
                WcsError::missing_keyword("COE projection requires PV2_1 (theta_a)")
            })?;
            Ok(Projection::coe(*theta_a))
        }
        "COD" => {
            let theta_a = pv_params.get(&(2, 1)).ok_or_else(|| {
                WcsError::missing_keyword("COD projection requires PV2_1 (theta_a)")
            })?;
            Ok(Projection::cod(*theta_a))
        }
        "COO" => {
            let theta_a = pv_params.get(&(2, 1)).ok_or_else(|| {
                WcsError::missing_keyword("COO projection requires PV2_1 (theta_a)")
            })?;
            Ok(Projection::coo(*theta_a))
        }
        "BON" => {
            let theta_1 = pv_params.get(&(2, 1)).ok_or_else(|| {
                WcsError::missing_keyword("BON projection requires PV2_1 (theta_1)")
            })?;
            Ok(Projection::bon(*theta_1))
        }
        "PCO" => Ok(Projection::pco()),
        "TSC" => Ok(Projection::tsc()),
        "CSC" => Ok(Projection::csc()),
        "QSC" => Ok(Projection::qsc()),
        _ => Err(WcsError::unsupported_projection(code)),
    }
}

fn projection_to_code(proj: &Projection) -> String {
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

fn parse_ctype(ctype: &str) -> WcsResult<(&str, &str)> {
    let trimmed = ctype.trim();

    if let Some(dash_pos) = trimmed.rfind('-') {
        if dash_pos == 0 {
            return Err(WcsError::invalid_keyword(
                "CTYPE",
                format!("Invalid CTYPE format: '{}'", ctype),
            ));
        }

        let prefix_part = &trimmed[..dash_pos];
        let proj_part = &trimmed[dash_pos + 1..];

        let prefix = prefix_part.trim_end_matches('-');

        if proj_part.is_empty() {
            return Err(WcsError::invalid_keyword(
                "CTYPE",
                format!("Missing projection code in CTYPE: '{}'", ctype),
            ));
        }

        Ok((prefix, proj_part))
    } else {
        Err(WcsError::invalid_keyword(
            "CTYPE",
            format!("Invalid CTYPE format (no dash separator): '{}'", ctype),
        ))
    }
}

fn parse_matrix(header: &impl KeywordProvider) -> WcsResult<MatrixSpec> {
    let cd11 = header.get_float("CD1_1");
    let cd12 = header.get_float("CD1_2");
    let cd21 = header.get_float("CD2_1");
    let cd22 = header.get_float("CD2_2");

    if cd11.is_some() || cd12.is_some() || cd21.is_some() || cd22.is_some() {
        let cd = [
            [cd11.unwrap_or(0.0), cd12.unwrap_or(0.0)],
            [cd21.unwrap_or(0.0), cd22.unwrap_or(0.0)],
        ];
        return Ok(MatrixSpec::Cd(cd));
    }

    let cdelt1 = header.get_float("CDELT1");
    let cdelt2 = header.get_float("CDELT2");

    if let (Some(c1), Some(c2)) = (cdelt1, cdelt2) {
        let pc11 = header.get_float("PC1_1").unwrap_or(1.0);
        let pc12 = header.get_float("PC1_2").unwrap_or(0.0);
        let pc21 = header.get_float("PC2_1").unwrap_or(0.0);
        let pc22 = header.get_float("PC2_2").unwrap_or(1.0);

        let pc = [[pc11, pc12], [pc21, pc22]];
        let cdelt = [c1, c2];

        return Ok(MatrixSpec::PcCdelt { pc, cdelt });
    }

    Err(WcsError::missing_keyword(
        "CD1_1 or CDELT1 (no transformation matrix found)",
    ))
}

fn parse_pv_params(header: &impl KeywordProvider) -> HashMap<(u8, u8), f64> {
    let mut pv_params = HashMap::new();

    for axis in 1..=2u8 {
        for index in 0..=20u8 {
            let key = format!("PV{}_{}", axis, index);
            if let Some(value) = header.get_float(&key) {
                pv_params.insert((axis, index), value);
            }
        }
    }

    pv_params
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builder_new() {
        let builder = WcsBuilder::new();
        assert!(builder.crpix.is_none());
        assert!(builder.crval.is_none());
        assert_eq!(builder.matrix, MatrixSpec::None);
        assert!(builder.projection.is_none());
        assert!(builder.lonpole.is_none());
        assert!(builder.latpole.is_none());
        assert!(builder.pv_params.is_empty());
        assert!(builder.coord_type.is_none());
        assert!(builder.proj_code.is_none());
    }

    #[test]
    fn test_builder_crpix() {
        let builder = WcsBuilder::new().crpix(512.0, 512.0);
        assert_eq!(builder.crpix, Some([512.0, 512.0]));
    }

    #[test]
    fn test_builder_crval() {
        let builder = WcsBuilder::new().crval(180.0, 45.0);
        assert_eq!(builder.crval, Some([180.0, 45.0]));
    }

    #[test]
    fn test_builder_cd_matrix() {
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let builder = WcsBuilder::new().cd_matrix(cd);
        assert_eq!(builder.matrix, MatrixSpec::Cd(cd));
    }

    #[test]
    fn test_builder_pc_cdelt() {
        let pc = [[1.0, 0.0], [0.0, 1.0]];
        let cdelt = [0.001, 0.001];
        let builder = WcsBuilder::new().pc_cdelt(pc, cdelt);
        assert_eq!(builder.matrix, MatrixSpec::PcCdelt { pc, cdelt });
    }

    #[test]
    fn test_builder_projection() {
        let builder = WcsBuilder::new().projection(Projection::tan());
        assert_eq!(builder.projection, Some(Projection::Tan));
    }

    #[test]
    fn test_builder_lonpole() {
        let builder = WcsBuilder::new().lonpole(180.0);
        assert_eq!(builder.lonpole, Some(180.0));
    }

    #[test]
    fn test_builder_latpole() {
        let builder = WcsBuilder::new().latpole(90.0);
        assert_eq!(builder.latpole, Some(90.0));
    }

    #[test]
    fn test_builder_pv() {
        let builder = WcsBuilder::new().pv(2, 1, 0.5).pv(2, 2, 0.25);
        assert_eq!(builder.pv_params.get(&(2, 1)), Some(&0.5));
        assert_eq!(builder.pv_params.get(&(2, 2)), Some(&0.25));
    }

    #[test]
    fn test_builder_coord_type() {
        let builder = WcsBuilder::new().coord_type(CoordType::Galactic);
        assert_eq!(builder.coord_type, Some(CoordType::Galactic));
    }

    #[test]
    fn test_builder_proj_code() {
        let builder = WcsBuilder::new().proj_code("TAN");
        assert_eq!(builder.proj_code, Some("TAN".to_string()));
    }

    #[test]
    fn test_builder_chaining() {
        let builder = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .lonpole(180.0)
            .latpole(90.0)
            .pv(2, 1, 0.5)
            .coord_type(CoordType::Equatorial)
            .proj_code("TAN");

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
    fn test_cd_matrix_overwrites_pc_cdelt() {
        let pc = [[1.0, 0.0], [0.0, 1.0]];
        let cdelt = [0.001, 0.001];
        let cd = [[0.002, 0.0], [0.0, 0.002]];

        let builder = WcsBuilder::new().pc_cdelt(pc, cdelt).cd_matrix(cd);

        assert_eq!(builder.matrix, MatrixSpec::Cd(cd));
    }

    #[test]
    fn test_pc_cdelt_overwrites_cd_matrix() {
        let cd = [[0.002, 0.0], [0.0, 0.002]];
        let pc = [[1.0, 0.0], [0.0, 1.0]];
        let cdelt = [0.001, 0.001];

        let builder = WcsBuilder::new().cd_matrix(cd).pc_cdelt(pc, cdelt);

        assert_eq!(builder.matrix, MatrixSpec::PcCdelt { pc, cdelt });
    }

    #[test]
    fn test_coord_type_default() {
        assert_eq!(CoordType::default(), CoordType::Equatorial);
    }

    #[test]
    fn test_builder_with_string_proj_code() {
        let builder = WcsBuilder::new().proj_code(String::from("SIN"));
        assert_eq!(builder.proj_code, Some("SIN".to_string()));
    }

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
    fn test_wcs_pixel_to_celestial_at_crpix() {
        use crate::PixelCoord;
        use astro_core::assert_ulp_lt;

        let wcs = create_simple_tan_wcs().unwrap();
        let pixel = PixelCoord::new(512.0, 512.0);
        let celestial = wcs.pixel_to_celestial(pixel).unwrap();

        assert_ulp_lt!(celestial.alpha().degrees(), 180.0, 10);
        assert_ulp_lt!(celestial.delta().degrees(), 45.0, 10);
    }

    #[test]
    fn test_wcs_roundtrip_at_crpix() {
        use crate::PixelCoord;
        use astro_core::assert_ulp_lt;

        let wcs = create_simple_tan_wcs().unwrap();
        let original = PixelCoord::new(512.0, 512.0);

        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        assert_ulp_lt!(original.x(), recovered.x(), 10);
        assert_ulp_lt!(original.y(), recovered.y(), 10);
    }

    #[test]
    fn test_wcs_roundtrip_off_center() {
        use crate::PixelCoord;

        let wcs = create_simple_tan_wcs().unwrap();
        let original = PixelCoord::new(256.0, 768.0);

        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        let tol = 1e-9;
        assert!((original.x() - recovered.x()).abs() < tol);
        assert!((original.y() - recovered.y()).abs() < tol);
    }

    #[test]
    fn test_wcs_pix2world_world2pix_roundtrip() {
        let wcs = create_simple_tan_wcs().unwrap();
        let (x, y) = (300.0, 700.0);

        let (ra, dec) = wcs.pix2world(x, y).unwrap();
        let (x_recovered, y_recovered) = wcs.world2pix(ra, dec).unwrap();

        // Tolerance accounts for ARM vs x86 FPU differences in trig functions
        let tol = 1e-8;
        assert!((x - x_recovered).abs() < tol);
        assert!((y - y_recovered).abs() < tol);
    }

    #[test]
    fn test_wcs_projection_code() {
        let wcs = create_simple_tan_wcs().unwrap();
        assert_eq!(wcs.projection_code(), "TAN");
    }

    #[test]
    fn test_wcs_coord_type() {
        let wcs = create_simple_tan_wcs().unwrap();
        assert_eq!(wcs.coord_type(), CoordType::Equatorial);
    }

    #[test]
    fn test_wcs_crpix() {
        let wcs = create_simple_tan_wcs().unwrap();
        assert_eq!(wcs.crpix(), [512.0, 512.0]);
    }

    #[test]
    fn test_wcs_crval() {
        let wcs = create_simple_tan_wcs().unwrap();
        assert_eq!(wcs.crval(), (180.0, 45.0));
    }

    #[test]
    fn test_wcs_pixel_scale() {
        let wcs = create_simple_tan_wcs().unwrap();
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
    fn test_wcs_roundtrip_with_rotated_cd_matrix() {
        use crate::PixelCoord;

        let crpix = [512.0, 512.0];
        let angle = std::f64::consts::PI / 6.0;
        let scale = 0.0005;
        let (angle_s, angle_c) = angle.sin_cos();
        let cd = [
            [scale * angle_c, -scale * angle_s],
            [scale * angle_s, scale * angle_c],
        ];
        let linear = LinearTransform::from_cd(crpix, cd).unwrap();

        let projection = Projection::tan();
        let (_, theta_0) = projection.native_reference();
        let crval_lon = 120.0;
        let crval_lat = -30.0;

        let rotation = SphericalRotation::from_crval(
            Angle::from_degrees(crval_lon),
            Angle::from_degrees(crval_lat),
            Angle::from_degrees(theta_0),
            None,
            None,
        )
        .unwrap();

        let wcs = Wcs::new(
            linear,
            projection,
            rotation,
            CoordType::Equatorial,
            "TAN".to_string(),
            (crval_lon, crval_lat),
            None,
        );

        let original = PixelCoord::new(400.0, 600.0);
        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        let tol = 1e-8;
        assert!((original.x() - recovered.x()).abs() < tol);
        assert!((original.y() - recovered.y()).abs() < tol);
    }

    #[test]
    fn test_wcs_roundtrip_arc_projection() {
        use crate::PixelCoord;

        let crpix = [256.0, 256.0];
        let cd = [[0.002, 0.0], [0.0, 0.002]];
        let linear = LinearTransform::from_cd(crpix, cd).unwrap();

        let projection = Projection::arc();
        let (_, theta_0) = projection.native_reference();
        let crval_lon = 90.0;
        let crval_lat = 60.0;

        let rotation = SphericalRotation::from_crval(
            Angle::from_degrees(crval_lon),
            Angle::from_degrees(crval_lat),
            Angle::from_degrees(theta_0),
            None,
            None,
        )
        .unwrap();

        let wcs = Wcs::new(
            linear,
            projection,
            rotation,
            CoordType::Equatorial,
            "ARC".to_string(),
            (crval_lon, crval_lat),
            None,
        );

        let original = PixelCoord::new(200.0, 300.0);
        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        let tol = 1e-9;
        assert!((original.x() - recovered.x()).abs() < tol);
        assert!((original.y() - recovered.y()).abs() < tol);
    }

    #[test]
    fn test_wcs_roundtrip_stg_projection() {
        use crate::PixelCoord;

        let crpix = [128.0, 128.0];
        let cd = [[0.005, 0.0], [0.0, 0.005]];
        let linear = LinearTransform::from_cd(crpix, cd).unwrap();

        let projection = Projection::stg();
        let (_, theta_0) = projection.native_reference();
        let crval_lon = 0.0;
        let crval_lat = 85.0;

        let rotation = SphericalRotation::from_crval(
            Angle::from_degrees(crval_lon),
            Angle::from_degrees(crval_lat),
            Angle::from_degrees(theta_0),
            None,
            None,
        )
        .unwrap();

        let wcs = Wcs::new(
            linear,
            projection,
            rotation,
            CoordType::Equatorial,
            "STG".to_string(),
            (crval_lon, crval_lat),
            None,
        );

        let original = PixelCoord::new(100.0, 150.0);
        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        let tol = 1e-9;
        assert!((original.x() - recovered.x()).abs() < tol);
        assert!((original.y() - recovered.y()).abs() < tol);
    }

    #[test]
    fn test_parse_ctype_ra_tan() {
        let (prefix, proj) = parse_ctype("RA---TAN").unwrap();
        assert_eq!(prefix, "RA");
        assert_eq!(proj, "TAN");
    }

    #[test]
    fn test_parse_ctype_dec_tan() {
        let (prefix, proj) = parse_ctype("DEC--TAN").unwrap();
        assert_eq!(prefix, "DEC");
        assert_eq!(proj, "TAN");
    }

    #[test]
    fn test_parse_ctype_glon_sin() {
        let (prefix, proj) = parse_ctype("GLON-SIN").unwrap();
        assert_eq!(prefix, "GLON");
        assert_eq!(proj, "SIN");
    }

    #[test]
    fn test_parse_ctype_glat_sin() {
        let (prefix, proj) = parse_ctype("GLAT-SIN").unwrap();
        assert_eq!(prefix, "GLAT");
        assert_eq!(proj, "SIN");
    }

    #[test]
    fn test_parse_ctype_with_whitespace() {
        let (prefix, proj) = parse_ctype("  RA---TAN  ").unwrap();
        assert_eq!(prefix, "RA");
        assert_eq!(proj, "TAN");
    }

    #[test]
    fn test_parse_ctype_invalid_no_dash() {
        let result = parse_ctype("RATAN");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("no dash separator"));
    }

    #[test]
    fn test_parse_ctype_invalid_empty_proj() {
        let result = parse_ctype("RA---");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Missing projection"));
    }

    #[test]
    fn test_from_header_cd_matrix() {
        use crate::header::KeywordMap;

        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", -0.001)
            .set_float("CD1_2", 0.0)
            .set_float("CD2_1", 0.0)
            .set_float("CD2_2", 0.001);

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
        use crate::header::KeywordMap;

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
        use crate::header::KeywordMap;

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
        use crate::header::KeywordMap;

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
        use crate::header::KeywordMap;

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
        use crate::header::KeywordMap;

        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---SIN")
            .set_string("CTYPE2", "DEC--SIN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", 0.001)
            .set_float("CD2_2", 0.001)
            .set_float("PV2_1", 0.5)
            .set_float("PV2_2", 0.25);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.pv_params.get(&(2, 1)), Some(&0.5));
        assert_eq!(builder.pv_params.get(&(2, 2)), Some(&0.25));
    }

    #[test]
    fn test_from_header_lonpole_latpole() {
        use crate::header::KeywordMap;

        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", 0.001)
            .set_float("CD2_2", 0.001)
            .set_float("LONPOLE", 180.0)
            .set_float("LATPOLE", 45.0);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.lonpole, Some(180.0));
        assert_eq!(builder.latpole, Some(45.0));
    }

    #[test]
    fn test_from_header_mismatched_projection() {
        use crate::header::KeywordMap;

        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--SIN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", 0.001);

        let result = WcsBuilder::from_header(&header);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Mismatched"));
    }

    #[test]
    fn test_from_header_ecliptic_coords() {
        use crate::header::KeywordMap;

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
    fn test_from_header_partial_cd_matrix() {
        use crate::header::KeywordMap;

        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", 0.001)
            .set_float("CD2_2", 0.001);

        let builder = WcsBuilder::from_header(&header).unwrap();

        assert_eq!(builder.matrix, MatrixSpec::Cd([[0.001, 0.0], [0.0, 0.001]]));
    }

    #[test]
    fn test_build_succeeds_with_minimal_valid_config() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build()
            .unwrap();

        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
        assert_eq!(wcs.projection_code(), "TAN");
    }

    #[test]
    fn test_build_fails_on_missing_crpix() {
        let result = WcsBuilder::new()
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build();

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Missing CRPIX"));
    }

    #[test]
    fn test_build_fails_on_missing_crval() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build();

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Missing CRVAL"));
    }

    #[test]
    fn test_build_fails_on_missing_matrix() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .projection(Projection::tan())
            .build();

        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Missing transformation matrix"));
    }

    #[test]
    fn test_build_fails_on_missing_projection() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .build();

        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Missing projection"));
    }

    #[test]
    fn test_build_creates_correct_projection_from_code() {
        let wcs = WcsBuilder::new()
            .crpix(256.0, 256.0)
            .crval(90.0, 60.0)
            .cd_matrix([[0.002, 0.0], [0.0, 0.002]])
            .proj_code("ARC")
            .build()
            .unwrap();

        assert_eq!(wcs.projection_code(), "ARC");
    }

    #[test]
    fn test_build_with_sin_params() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("SIN")
            .pv(2, 1, 0.5)
            .pv(2, 2, 0.25)
            .build()
            .unwrap();

        assert_eq!(wcs.projection_code(), "SIN");
    }

    #[test]
    fn test_build_unsupported_projection() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("XYZ")
            .build();

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("XYZ"));
    }

    #[test]
    fn test_build_conic_missing_param() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COP")
            .build();

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_build_conic_with_param() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COP")
            .pv(2, 1, 45.0)
            .build()
            .unwrap();

        assert_eq!(wcs.projection_code(), "COP");
    }

    #[test]
    fn test_builder_full_roundtrip() {
        use crate::PixelCoord;

        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("TAN")
            .coord_type(CoordType::Equatorial)
            .build()
            .unwrap();

        let original = PixelCoord::new(300.0, 700.0);
        let celestial = wcs.pixel_to_celestial(original).unwrap();
        let recovered = wcs.celestial_to_pixel(celestial).unwrap();

        // Tolerance accounts for ARM vs x86 FPU differences in trig functions
        let tol = 1e-8;
        assert!((original.x() - recovered.x()).abs() < tol);
        assert!((original.y() - recovered.y()).abs() < tol);
    }

    #[test]
    fn test_builder_from_header_then_build() {
        use crate::header::KeywordMap;

        let mut header = KeywordMap::new();
        header
            .set_string("CTYPE1", "RA---TAN")
            .set_string("CTYPE2", "DEC--TAN")
            .set_float("CRPIX1", 512.0)
            .set_float("CRPIX2", 512.0)
            .set_float("CRVAL1", 180.0)
            .set_float("CRVAL2", 45.0)
            .set_float("CD1_1", -0.001)
            .set_float("CD1_2", 0.0)
            .set_float("CD2_1", 0.0)
            .set_float("CD2_2", 0.001);

        let wcs = WcsBuilder::from_header(&header).unwrap().build().unwrap();

        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
        assert_eq!(wcs.projection_code(), "TAN");
        assert_eq!(wcs.coord_type(), CoordType::Equatorial);
    }

    #[test]
    fn test_validate_returns_ok_for_valid_builder() {
        let builder = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan());

        assert!(builder.validate().is_ok());
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
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .lonpole(180.0)
            .latpole(45.0)
            .build()
            .unwrap();

        assert_eq!(wcs.crval(), (180.0, 45.0));
    }

    #[test]
    fn test_projection_inferred_from_enum() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::stg())
            .build()
            .unwrap();

        assert_eq!(wcs.projection_code(), "STG");
    }

    // ==================== Tests for create_projection_from_code ====================

    #[test]
    fn test_create_sin_projection_default_params() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("SIN", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "SIN");
    }

    #[test]
    fn test_create_arc_projection() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("ARC", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "ARC");
    }

    #[test]
    fn test_create_stg_projection() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("STG", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "STG");
    }

    #[test]
    fn test_create_zea_projection() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("ZEA", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "ZEA");
    }

    #[test]
    fn test_create_azp_projection_with_params() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 2.0); // mu
        pv_params.insert((2, 2), 30.0); // gamma
        let proj = create_projection_from_code("AZP", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "AZP");
    }

    #[test]
    fn test_create_azp_projection_default_params() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("AZP", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "AZP");
    }

    #[test]
    fn test_create_szp_projection_with_params() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 2.0); // mu
        pv_params.insert((2, 2), 45.0); // phi_c
        pv_params.insert((2, 3), 60.0); // theta_c
        let proj = create_projection_from_code("SZP", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "SZP");
    }

    #[test]
    fn test_create_szp_projection_default_params() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("SZP", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "SZP");
    }

    #[test]
    fn test_create_zpn_projection_with_coefficients() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 0), 0.0);
        pv_params.insert((2, 1), 1.0);
        pv_params.insert((2, 3), 0.1); // sparse: index 3 with gap at index 2
        let proj = create_projection_from_code("ZPN", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "ZPN");
    }

    #[test]
    fn test_create_zpn_projection_empty_coefficients_uses_defaults() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("ZPN", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "ZPN");
    }

    #[test]
    fn test_create_air_projection_with_theta_b() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 45.0); // theta_b
        let proj = create_projection_from_code("AIR", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "AIR");
    }

    #[test]
    fn test_create_air_projection_default_theta_b() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("AIR", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "AIR");
    }

    #[test]
    fn test_create_cea_projection_with_lambda() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 0.5); // lambda
        let proj = create_projection_from_code("CEA", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "CEA");
    }

    #[test]
    fn test_create_cea_projection_default_lambda() {
        let pv_params = HashMap::new();
        let proj = create_projection_from_code("CEA", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "CEA");
    }

    #[test]
    fn test_create_cyp_projection_with_params() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 1.0); // mu
        pv_params.insert((2, 2), 2.0); // lambda
        let proj = create_projection_from_code("CYP", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "CYP");
    }

    #[test]
    fn test_create_coe_projection_with_theta_a() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 45.0); // theta_a
        let proj = create_projection_from_code("COE", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "COE");
    }

    #[test]
    fn test_create_coe_projection_missing_theta_a() {
        let pv_params = HashMap::new();
        let result = create_projection_from_code("COE", &pv_params);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_create_cod_projection_with_theta_a() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 30.0); // theta_a
        let proj = create_projection_from_code("COD", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "COD");
    }

    #[test]
    fn test_create_cod_projection_missing_theta_a() {
        let pv_params = HashMap::new();
        let result = create_projection_from_code("COD", &pv_params);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_create_coo_projection_with_theta_a() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 60.0); // theta_a
        let proj = create_projection_from_code("COO", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "COO");
    }

    #[test]
    fn test_create_coo_projection_missing_theta_a() {
        let pv_params = HashMap::new();
        let result = create_projection_from_code("COO", &pv_params);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_create_bon_projection_with_theta_1() {
        let mut pv_params = HashMap::new();
        pv_params.insert((2, 1), 45.0); // theta_1
        let proj = create_projection_from_code("BON", &pv_params).unwrap();
        assert_eq!(projection_to_code(&proj), "BON");
    }

    #[test]
    fn test_create_bon_projection_missing_theta_1() {
        let pv_params = HashMap::new();
        let result = create_projection_from_code("BON", &pv_params);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    // ==================== Tests for projection_to_code ====================

    #[test]
    fn test_projection_to_code_tan() {
        assert_eq!(projection_to_code(&Projection::tan()), "TAN");
    }

    #[test]
    fn test_projection_to_code_sin() {
        assert_eq!(projection_to_code(&Projection::sin()), "SIN");
    }

    #[test]
    fn test_projection_to_code_sin_with_params() {
        assert_eq!(
            projection_to_code(&Projection::sin_with_params(0.1, 0.2)),
            "SIN"
        );
    }

    #[test]
    fn test_projection_to_code_arc() {
        assert_eq!(projection_to_code(&Projection::arc()), "ARC");
    }

    #[test]
    fn test_projection_to_code_stg() {
        assert_eq!(projection_to_code(&Projection::stg()), "STG");
    }

    #[test]
    fn test_projection_to_code_zea() {
        assert_eq!(projection_to_code(&Projection::zea()), "ZEA");
    }

    #[test]
    fn test_projection_to_code_azp() {
        assert_eq!(projection_to_code(&Projection::azp(2.0, 30.0)), "AZP");
    }

    #[test]
    fn test_projection_to_code_szp() {
        assert_eq!(projection_to_code(&Projection::szp(2.0, 45.0, 60.0)), "SZP");
    }

    #[test]
    fn test_projection_to_code_zpn() {
        assert_eq!(projection_to_code(&Projection::zpn(vec![0.0, 1.0])), "ZPN");
    }

    #[test]
    fn test_projection_to_code_air() {
        assert_eq!(projection_to_code(&Projection::air(45.0)), "AIR");
    }

    #[test]
    fn test_projection_to_code_car() {
        assert_eq!(projection_to_code(&Projection::car()), "CAR");
    }

    #[test]
    fn test_projection_to_code_mer() {
        assert_eq!(projection_to_code(&Projection::mer()), "MER");
    }

    #[test]
    fn test_projection_to_code_cea() {
        assert_eq!(projection_to_code(&Projection::cea_with_lambda(0.5)), "CEA");
    }

    #[test]
    fn test_projection_to_code_cyp() {
        assert_eq!(projection_to_code(&Projection::cyp(1.0, 2.0)), "CYP");
    }

    #[test]
    fn test_projection_to_code_sfl() {
        assert_eq!(projection_to_code(&Projection::sfl()), "SFL");
    }

    #[test]
    fn test_projection_to_code_par() {
        assert_eq!(projection_to_code(&Projection::par()), "PAR");
    }

    #[test]
    fn test_projection_to_code_mol() {
        assert_eq!(projection_to_code(&Projection::mol()), "MOL");
    }

    #[test]
    fn test_projection_to_code_ait() {
        assert_eq!(projection_to_code(&Projection::ait()), "AIT");
    }

    #[test]
    fn test_projection_to_code_cop() {
        assert_eq!(projection_to_code(&Projection::cop(45.0)), "COP");
    }

    #[test]
    fn test_projection_to_code_coe() {
        assert_eq!(projection_to_code(&Projection::coe(45.0)), "COE");
    }

    #[test]
    fn test_projection_to_code_cod() {
        assert_eq!(projection_to_code(&Projection::cod(45.0)), "COD");
    }

    #[test]
    fn test_projection_to_code_coo() {
        assert_eq!(projection_to_code(&Projection::coo(45.0)), "COO");
    }

    #[test]
    fn test_projection_to_code_bon() {
        assert_eq!(projection_to_code(&Projection::bon(45.0)), "BON");
    }

    #[test]
    fn test_projection_to_code_pco() {
        assert_eq!(projection_to_code(&Projection::pco()), "PCO");
    }

    #[test]
    fn test_projection_to_code_tsc() {
        assert_eq!(projection_to_code(&Projection::tsc()), "TSC");
    }

    #[test]
    fn test_projection_to_code_csc() {
        assert_eq!(projection_to_code(&Projection::csc()), "CSC");
    }

    #[test]
    fn test_projection_to_code_qsc() {
        assert_eq!(projection_to_code(&Projection::qsc()), "QSC");
    }

    // ==================== Tests for parse_ctype edge cases ====================

    #[test]
    fn test_parse_ctype_leading_dash_invalid() {
        let result = parse_ctype("-TAN");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid CTYPE"));
    }

    #[test]
    fn test_parse_ctype_single_dash_separator() {
        let (prefix, proj) = parse_ctype("GLON-TAN").unwrap();
        assert_eq!(prefix, "GLON");
        assert_eq!(proj, "TAN");
    }

    #[test]
    fn test_parse_ctype_multiple_dashes() {
        let (prefix, proj) = parse_ctype("RA---ZEA").unwrap();
        assert_eq!(prefix, "RA");
        assert_eq!(proj, "ZEA");
    }

    // ==================== WcsBuilder integration tests for projection codes ====================

    #[test]
    fn test_build_with_sin_default_params() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("SIN")
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "SIN");
    }

    #[test]
    fn test_build_with_azp_params() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("AZP")
            .pv(2, 1, 2.0)
            .pv(2, 2, 30.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "AZP");
    }

    #[test]
    fn test_build_with_szp_params() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("SZP")
            .pv(2, 1, 2.0)
            .pv(2, 2, 45.0)
            .pv(2, 3, 60.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "SZP");
    }

    #[test]
    fn test_build_with_zpn_sparse_coefficients() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("ZPN")
            .pv(2, 0, 0.0)
            .pv(2, 1, 1.0)
            .pv(2, 5, 0.001) // sparse coefficient
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "ZPN");
    }

    #[test]
    fn test_build_with_air_projection() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("AIR")
            .pv(2, 1, 45.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "AIR");
    }

    #[test]
    fn test_build_with_cea_lambda() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("CEA")
            .pv(2, 1, 0.5)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "CEA");
    }

    #[test]
    fn test_build_with_coe_theta_a() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COE")
            .pv(2, 1, 45.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "COE");
    }

    #[test]
    fn test_build_with_cod_theta_a() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COD")
            .pv(2, 1, 30.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "COD");
    }

    #[test]
    fn test_build_with_coo_theta_a() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COO")
            .pv(2, 1, 60.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "COO");
    }

    #[test]
    fn test_build_with_bon_theta_1() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("BON")
            .pv(2, 1, 45.0)
            .build()
            .unwrap();
        assert_eq!(wcs.projection_code(), "BON");
    }

    #[test]
    fn test_build_coe_missing_required_param() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COE")
            .build();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_build_cod_missing_required_param() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COD")
            .build();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_build_coo_missing_required_param() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("COO")
            .build();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_build_bon_missing_required_param() {
        let result = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .proj_code("BON")
            .build();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("PV2_1"));
    }

    #[test]
    fn test_format_ctype_equatorial() {
        assert_eq!(format_ctype("RA", "TAN"), "RA---TAN");
        assert_eq!(format_ctype("DEC", "TAN"), "DEC--TAN");
    }

    #[test]
    fn test_format_ctype_galactic() {
        assert_eq!(format_ctype("GLON", "SIN"), "GLON-SIN");
        assert_eq!(format_ctype("GLAT", "SIN"), "GLAT-SIN");
    }

    #[test]
    fn test_format_ctype_ecliptic() {
        assert_eq!(format_ctype("ELON", "ARC"), "ELON-ARC");
        assert_eq!(format_ctype("ELAT", "ARC"), "ELAT-ARC");
    }

    #[test]
    fn test_format_ctype_length() {
        assert_eq!(format_ctype("RA", "TAN").len(), 8);
        assert_eq!(format_ctype("DEC", "TAN").len(), 8);
        assert_eq!(format_ctype("GLON", "TAN").len(), 8);
        assert_eq!(format_ctype("GLAT", "TAN").len(), 8);
    }

    #[test]
    fn test_ctype_keywords_equatorial() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .coord_type(CoordType::Equatorial)
            .build()
            .unwrap();

        let keywords = wcs.to_keywords();
        let ctype1 = keywords
            .iter()
            .find(|k| k.name == "CTYPE1")
            .expect("CTYPE1 not found");
        let ctype2 = keywords
            .iter()
            .find(|k| k.name == "CTYPE2")
            .expect("CTYPE2 not found");

        assert_eq!(
            ctype1.value,
            WcsKeywordValue::String("RA---TAN".to_string())
        );
        assert_eq!(
            ctype2.value,
            WcsKeywordValue::String("DEC--TAN".to_string())
        );
    }

    #[test]
    fn test_ctype_keywords_galactic() {
        let wcs = WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::sin())
            .coord_type(CoordType::Galactic)
            .build()
            .unwrap();

        let keywords = wcs.to_keywords();
        let ctype1 = keywords
            .iter()
            .find(|k| k.name == "CTYPE1")
            .expect("CTYPE1 not found");
        let ctype2 = keywords
            .iter()
            .find(|k| k.name == "CTYPE2")
            .expect("CTYPE2 not found");

        assert_eq!(
            ctype1.value,
            WcsKeywordValue::String("GLON-SIN".to_string())
        );
        assert_eq!(
            ctype2.value,
            WcsKeywordValue::String("GLAT-SIN".to_string())
        );
    }
}
