use crate::coordinate::{CelestialCoord, IntermediateCoord, PixelCoord};
use crate::distortion::DistortionModel;
use crate::error::WcsResult;
use crate::linear::LinearTransform;
use crate::spherical::{Projection, SphericalRotation};

use super::distortion_keywords::{ctype_code, distortion_keywords};
use super::fiducial::Fiducial;
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
    pub(crate) fn from_ctype_prefix(prefix: &str) -> Self {
        let named = [
            Self::Equatorial,
            Self::Galactic,
            Self::Ecliptic,
            Self::Helioecliptic,
            Self::Supergalactic,
        ];
        named
            .into_iter()
            .find(|coord_type| {
                let (lon, lat) = coord_type.ctype_prefixes();
                prefix == lon || prefix == lat
            })
            .unwrap_or(Self::Generic)
    }

    fn ctype_prefixes(self) -> (&'static str, &'static str) {
        match self {
            Self::Equatorial => ("RA", "DEC"),
            Self::Galactic => ("GLON", "GLAT"),
            Self::Ecliptic => ("ELON", "ELAT"),
            Self::Helioecliptic => ("HLON", "HLAT"),
            Self::Supergalactic => ("SLON", "SLAT"),
            Self::Generic => ("XLON", "XLAT"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Wcs {
    linear: LinearTransform,
    projection: Projection,
    rotation: SphericalRotation,
    coord_type: CoordType,
    crval_deg: (f64, f64),
    distortion: Option<DistortionModel>,
    fiducial: Fiducial,
}

impl Wcs {
    pub(super) fn new(
        linear: LinearTransform,
        projection: Projection,
        rotation: SphericalRotation,
        coord_type: CoordType,
        crval_deg: (f64, f64),
        distortion: Option<DistortionModel>,
    ) -> Self {
        Self {
            linear,
            fiducial: Fiducial::of(&projection),
            projection,
            rotation,
            coord_type,
            crval_deg,
            distortion,
        }
    }

    pub(super) fn with_fiducial(self, fiducial: Fiducial) -> Self {
        Self { fiducial, ..self }
    }

    pub fn pixel_to_celestial(&self, pixel: PixelCoord) -> WcsResult<CelestialCoord> {
        let pixel = self.pixel_distortion(pixel);
        let intermediate = self.plane_distortion(self.linear.pixel_to_intermediate(pixel));
        let native = self
            .projection
            .deproject(self.fiducial.add_offset(intermediate))?;
        self.rotation.native_to_celestial(native)
    }

    pub fn celestial_to_pixel(&self, celestial: CelestialCoord) -> WcsResult<PixelCoord> {
        let native = self.rotation.celestial_to_native(celestial)?;
        let plane = self.projection.project(native)?;
        let intermediate = self.plane_distortion_inverse(self.fiducial.remove_offset(plane))?;
        self.pixel_distortion_inverse(self.linear.intermediate_to_pixel(intermediate))
    }

    // SIP corrects pixel offsets ahead of the linear transform; TPV corrects
    // intermediate coordinates after it.
    fn pixel_distortion(&self, pixel: PixelCoord) -> PixelCoord {
        match &self.distortion {
            Some(DistortionModel::Sip(sip)) => {
                let (x, y) = sip.apply(pixel.x(), pixel.y());
                PixelCoord::new(x, y)
            }
            _ => pixel,
        }
    }

    fn pixel_distortion_inverse(&self, pixel: PixelCoord) -> WcsResult<PixelCoord> {
        match &self.distortion {
            Some(DistortionModel::Sip(sip)) => {
                let (x, y) = sip.apply_inverse(pixel.x(), pixel.y())?;
                Ok(PixelCoord::new(x, y))
            }
            _ => Ok(pixel),
        }
    }

    fn plane_distortion(&self, plane: IntermediateCoord) -> IntermediateCoord {
        let (x, y) = (plane.x_deg(), plane.y_deg());
        let (x, y) = match &self.distortion {
            Some(DistortionModel::Tpv(tpv)) => tpv.apply(x, y),
            _ => return plane,
        };
        IntermediateCoord::new(x, y)
    }

    fn plane_distortion_inverse(&self, plane: IntermediateCoord) -> WcsResult<IntermediateCoord> {
        let (x, y) = (plane.x_deg(), plane.y_deg());
        let (x, y) = match &self.distortion {
            Some(DistortionModel::Tpv(tpv)) => tpv.apply_inverse(x, y)?,
            _ => return Ok(plane),
        };
        Ok(IntermediateCoord::new(x, y))
    }

    pub fn pix2world(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        let pixel = PixelCoord::new(x, y);
        let celestial = self.pixel_to_celestial(pixel)?;
        Ok((celestial.alpha_deg(), celestial.delta_deg()))
    }

    pub fn world2pix(&self, lon: f64, lat: f64) -> WcsResult<(f64, f64)> {
        let celestial = CelestialCoord::new(lon, lat);
        let pixel = self.celestial_to_pixel(celestial)?;
        Ok((pixel.x(), pixel.y()))
    }

    #[inline]
    pub fn projection_code(&self) -> &'static str {
        self.projection.code()
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
        // Paper II §2 requires 'deg' for celestial axes, given explicitly or not.
        keywords.push(WcsKeyword::string("CUNIT1", "deg"));
        keywords.push(WcsKeyword::string("CUNIT2", "deg"));
        keywords.extend(self.crpix_keywords());
        keywords.extend(self.crval_keywords());
        keywords.extend(self.cd_keywords());
        keywords.extend(self.pole_keywords());
        keywords.extend(self.pv_keywords());
        keywords.extend(distortion_keywords(self.distortion.as_ref()));

        keywords
    }

    fn ctype_keywords(&self) -> Vec<WcsKeyword> {
        let (prefix1, prefix2) = self.coord_type.ctype_prefixes();
        let code = &ctype_code(self.projection_code(), self.distortion.as_ref());

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

        let default_phi_p = SphericalRotation::default_lonpole(
            self.crval_deg.1,
            self.fiducial.phi_0(),
            self.fiducial.theta_0(),
        );
        if (phi_p - default_phi_p).abs() > 1e-10 {
            keywords.push(WcsKeyword::real("LONPOLE", phi_p));
        }

        if (delta_p - 90.0).abs() > 1e-10 {
            keywords.push(WcsKeyword::real("LATPOLE", delta_p));
        }

        keywords
    }

    fn pv_keywords(&self) -> Vec<WcsKeyword> {
        let mut keywords = self.fiducial.keywords(&self.projection);
        let pv = self.projection.pv_values().into_iter();
        keywords.extend(pv.map(|(m, value)| WcsKeyword::real(format!("PV2_{m}"), value)));
        keywords
    }
}

pub(super) fn format_ctype(prefix: &str, proj_code: &str) -> String {
    let padding_len = 4 - prefix.len();
    let dashes = "-".repeat(padding_len + 1);
    format!("{}{}{}", prefix, dashes, proj_code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::WcsBuilder;
    use crate::distortion::{SipDistortion, TpvDistortion};

    fn tan_builder() -> WcsBuilder {
        WcsBuilder::new()
            .crpix(512.0, 512.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
    }

    fn create_simple_tan_wcs() -> crate::error::WcsResult<Wcs> {
        tan_builder().build()
    }

    type PixelRow = ((f64, f64), (f64, f64), (f64, f64));

    // (pixel, sky, pixel): pix2world of the first pair and world2pix of the
    // second, each pinned one way.
    fn assert_pixel_pins(wcs: &Wcs, rows: &[PixelRow]) {
        for &(pixel, sky, back) in rows {
            let got = wcs.pix2world(pixel.0, pixel.1).unwrap();
            assert_eq!(got, sky, "pix2world {pixel:?}");
            assert_eq!(
                wcs.world2pix(sky.0, sky.1).unwrap(),
                back,
                "world2pix {sky:?}"
            );
        }
    }

    // CRPIX maps to CRVAL and back exactly. Against 60 digits the other rows
    // are within 2.08 ulp on the sky and 102 ulp in pixels, at most 0.57 of
    // the summed half-ulp rounding of the steps. The pixel carries the
    // rounding of θ near 90°: half an ulp of θ, 7.1e-15°, is 7.1e-12 px at
    // 0.001°/px.
    #[test]
    fn test_wcs_simple_tan_pipeline() {
        let wcs = create_simple_tan_wcs().unwrap();
        let rows = [
            ((512.0, 512.0), (180.0, 45.0), (512.0, 512.0)),
            (
                (256.0, 768.0),
                (179.63634134715633, 45.25542127650822),
                (256.00000000000557, 767.9999999999987),
            ),
            (
                (300.0, 700.0),
                (179.69920249673768, 45.187604542264104),
                (300.0000000000147, 699.9999999999836),
            ),
        ];
        assert_pixel_pins(&wcs, &rows);
    }

    #[test]
    fn test_world2pix_rejects_latitude_beyond_a_pole() {
        let wcs = create_simple_tan_wcs().unwrap();
        for lat in [100.0, -100.0] {
            let result = wcs.world2pix(180.0, lat);
            assert!(result.is_err(), "{lat}: {result:?}");
        }
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

    // A rotated CD matrix with TAN at southern declination, then ARC and STG
    // nearer the pole. Against 60 digits, sky to pixel is within 182 ulp, at
    // most 0.71 of the summed half-ulp rounding of the steps, and carries the
    // rounding of θ near 90° as above. Pixel to sky is within 1.51 ulp except
    // α at δ = 85°, 116 ulp of a small longitude and 0.36 of the budget: θ's
    // half ulp moves it 302 ulp. On the sky that offset is 2.2e-15°.
    #[test]
    fn test_wcs_pinned_across_projections_and_cds() {
        let rotated = [
            [4.330127018922194e-4, -2.5e-4],
            [2.5e-4, 4.330127018922194e-4],
        ];
        let diagonal = |scale| [[scale, 0.0], [0.0, scale]];
        let cases = [
            (
                Projection::tan(),
                (120.0, -30.0),
                512.0,
                rotated,
                (400.0, 600.0),
            ),
            (
                Projection::arc(),
                (90.0, 60.0),
                256.0,
                diagonal(0.002),
                (200.0, 300.0),
            ),
            (
                Projection::stg(),
                (0.0, 85.0),
                128.0,
                diagonal(0.005),
                (100.0, 150.0),
            ),
        ];
        let pinned = [
            (
                (119.9186049310403, -29.989869852631916),
                (399.99999999999505, 599.9999999999868),
            ),
            (
                (89.77540320348389, 60.087809836039675),
                (200.000000000004, 299.9999999999946),
            ),
            (
                (-1.6419144601457878, 85.10800132287785),
                (100.00000000000097, 150.00000000000034),
            ),
        ];
        for ((projection, crval, crpix, cd, pixel), (sky, back)) in cases.into_iter().zip(pinned) {
            let wcs = WcsBuilder::new()
                .crpix(crpix, crpix)
                .crval(crval.0, crval.1)
                .cd_matrix(cd)
                .projection(projection)
                .build()
                .unwrap();
            assert_pixel_pins(&wcs, &[(pixel, sky, back)]);
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
        tan_builder().distortion(distortion).build().unwrap()
    }

    // In each of the three tests below, against 60 digits, pixel to sky is
    // within 1.53 ulp and sky to pixel within 64 ulp, at most 0.24 of the
    // summed half-ulp rounding of the steps, most of it θ's near 90°.
    #[test]
    fn test_pipeline_pinned_through_sip() {
        // SIP corrects pixels. Without a higher-order term the iterative
        // inverse would start at its answer.
        let mut sip = SipDistortion::new([512.0, 512.0], 2, 2);
        sip.set_a(2, 0, 1e-6).unwrap();
        sip.set_b(0, 2, 1e-6).unwrap();
        let wcs = tan_wcs_with(DistortionModel::Sip(sip));
        let sky = (179.98303319487104, 44.9879988880949);
        let back = (500.0000000000002, 499.99999999998664);
        assert_pixel_pins(&wcs, &[((500.0, 500.0), sky, back)]);
    }

    #[test]
    fn test_pipeline_pinned_through_tpv() {
        // TPV corrects intermediate coordinates, here ξ = x + 1e-4 x² and
        // η = y + 1e-4 y².
        let mut tpv = TpvDistortion::identity();
        tpv.set_pv1(4, 1e-4).unwrap();
        tpv.set_pv2(4, 1e-4).unwrap();
        let wcs = tan_wcs_with(DistortionModel::Tpv(Box::new(tpv)));
        let sky = (179.98303301166573, 44.98799875846779);
        let back = (499.99999999999784, 499.99999999999915);
        assert_pixel_pins(&wcs, &[((500.0, 500.0), sky, back)]);
    }

    // The projection's parameters go out as PV2_m, and one at its default
    // does not go out at all.
    #[test]
    fn test_to_keywords_writes_the_projection_parameters() {
        let wcs = tan_builder()
            .projection(Projection::cyp(1.0, 0.5))
            .build()
            .unwrap();
        let keywords = wcs.to_keywords().into_iter();
        let pv: Vec<_> = keywords.filter(|k| k.name.starts_with("PV")).collect();
        assert_eq!(pv, [WcsKeyword::real("PV2_2", 0.5)]);
    }
}
