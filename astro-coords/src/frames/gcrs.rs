use crate::{
    aberration::{apply_aberration, compute_earth_state, remove_aberration},
    frames::{CIRSPosition, ICRSPosition},
    transforms::CoordinateFrame,
    CoordError, CoordResult, Distance,
};
use astro_core::{matrix::RotationMatrix3, Angle, Vector3};
use astro_time::{transforms::NutationCalculator, TT};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GCRSPosition {
    ra: Angle,
    dec: Angle,
    epoch: TT,
    distance: Option<Distance>,
}

impl GCRSPosition {
    pub fn new(ra: Angle, dec: Angle, epoch: TT) -> CoordResult<Self> {
        let ra = ra.validate_right_ascension()?;
        let dec = dec.validate_declination(false)?;

        Ok(Self {
            ra,
            dec,
            epoch,
            distance: None,
        })
    }

    pub fn with_distance(
        ra: Angle,
        dec: Angle,
        epoch: TT,
        distance: Distance,
    ) -> CoordResult<Self> {
        let mut pos = Self::new(ra, dec, epoch)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn from_degrees(ra_deg: f64, dec_deg: f64, epoch: TT) -> CoordResult<Self> {
        Self::new(
            Angle::from_degrees(ra_deg),
            Angle::from_degrees(dec_deg),
            epoch,
        )
    }

    pub fn ra(&self) -> Angle {
        self.ra
    }

    pub fn dec(&self) -> Angle {
        self.dec
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn distance(&self) -> Option<Distance> {
        self.distance
    }

    pub fn set_distance(&mut self, distance: Distance) {
        self.distance = Some(distance);
    }

    pub fn remove_distance(&mut self) {
        self.distance = None;
    }

    pub fn unit_vector(&self) -> Vector3 {
        let (sin_dec, cos_dec) = self.dec.sin_cos();
        let (sin_ra, cos_ra) = self.ra.sin_cos();

        Vector3::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
    }

    pub fn from_unit_vector(unit: Vector3, epoch: TT) -> CoordResult<Self> {
        let r = (unit.x.powi(2) + unit.y.powi(2) + unit.z.powi(2)).sqrt();

        if r == 0.0 {
            return Err(CoordError::invalid_coordinate("Zero vector"));
        }

        let x = unit.x / r;
        let y = unit.y / r;
        let z = unit.z / r;

        let d2 = x * x + y * y;
        let ra = if d2 == 0.0 { 0.0 } else { y.atan2(x) };
        let dec = if z == 0.0 { 0.0 } else { z.atan2(d2.sqrt()) };

        Self::new(Angle::from_radians(ra), Angle::from_radians(dec), epoch)
    }

    pub fn to_cirs(&self) -> CoordResult<CIRSPosition> {
        let npb_matrix = Self::gcrs_to_cirs_matrix(&self.epoch)?;

        let gcrs_vec = self.unit_vector();
        let cirs_vec = npb_matrix * gcrs_vec;

        let mut cirs = CIRSPosition::from_unit_vector(cirs_vec, self.epoch)?;

        if let Some(distance) = self.distance {
            cirs.set_distance(distance);
        }

        Ok(cirs)
    }

    pub fn from_cirs(cirs: &CIRSPosition) -> CoordResult<Self> {
        let npb_matrix = Self::gcrs_to_cirs_matrix(&cirs.epoch())?;
        let cirs_to_gcrs = npb_matrix.transpose();

        let cirs_vec = cirs.unit_vector();
        let gcrs_vec = cirs_to_gcrs * cirs_vec;

        let mut gcrs = Self::from_unit_vector(gcrs_vec, cirs.epoch())?;

        if let Some(distance) = cirs.distance() {
            gcrs.distance = Some(distance);
        }

        Ok(gcrs)
    }

    fn gcrs_to_cirs_matrix(epoch: &TT) -> CoordResult<RotationMatrix3> {
        let jd = epoch.to_julian_date();
        let t = astro_core::utils::jd_to_centuries(jd.jd1(), jd.jd2());

        let nutation = epoch
            .nutation_iau2006a()
            .map_err(|e| CoordError::CoreError {
                message: format!("Nutation calculation failed: {}", e),
            })?;

        let precession_calc = astro_core::precession::PrecessionIAU2006::new();
        let npb_matrix = precession_calc.npb_matrix_iau2006a(
            t,
            nutation.nutation_longitude(),
            nutation.nutation_obliquity(),
        );

        let cio_solution = astro_core::CioSolution::calculate(&npb_matrix, t).map_err(|e| {
            CoordError::CoreError {
                message: format!("CIO calculation failed: {}", e),
            }
        })?;

        let c2i_matrix =
            astro_core::gcrs_to_cirs_matrix(cio_solution.cip.x, cio_solution.cip.y, cio_solution.s);

        Ok(c2i_matrix)
    }
}

impl CoordinateFrame for GCRSPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        let gcrs_vec = self.unit_vector();

        let earth_state = compute_earth_state(&self.epoch)?;
        let sun_earth_dist = earth_state.heliocentric_position.magnitude();
        let icrs_vec =
            remove_aberration(gcrs_vec, earth_state.barycentric_velocity, sun_earth_dist);

        let mut icrs = ICRSPosition::from_unit_vector(icrs_vec)?;

        if let Some(distance) = self.distance {
            icrs.set_distance(distance);
        }

        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let icrs_vec = icrs.unit_vector();

        let earth_state = compute_earth_state(epoch)?;
        let sun_earth_dist = earth_state.heliocentric_position.magnitude();
        let gcrs_vec = apply_aberration(icrs_vec, earth_state.barycentric_velocity, sun_earth_dist);

        let mut gcrs = Self::from_unit_vector(gcrs_vec, *epoch)?;

        if let Some(distance) = icrs.distance() {
            gcrs.distance = Some(distance);
        }

        Ok(gcrs)
    }
}

impl std::fmt::Display for GCRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GCRS(RA={:.6}°, Dec={:.6}°, epoch=J{:.1}",
            self.ra.degrees(),
            self.dec.degrees(),
            self.epoch.julian_year()
        )?;

        if let Some(distance) = self.distance {
            write!(f, ", d={}", distance)?;
        }

        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcrs_creation() {
        let epoch = TT::j2000();
        let pos = GCRSPosition::from_degrees(180.0, 45.0, epoch).unwrap();

        assert_eq!(pos.ra().degrees(), 180.0);
        assert_eq!(pos.dec().degrees(), 45.0);
        assert_eq!(pos.epoch(), epoch);
        assert_eq!(pos.distance(), None);
    }

    #[test]
    fn test_gcrs_with_distance() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(10.0).unwrap();
        let pos = GCRSPosition::with_distance(
            Angle::from_degrees(90.0),
            Angle::from_degrees(30.0),
            epoch,
            distance,
        )
        .unwrap();

        assert_eq!(pos.distance().unwrap(), distance);
    }

    #[test]
    fn test_unit_vector_conversion() {
        let epoch = TT::j2000();

        let vernal_equinox = GCRSPosition::from_degrees(0.0, 0.0, epoch).unwrap();
        let unit_vec = vernal_equinox.unit_vector();

        assert_eq!(unit_vec.x, 1.0);
        assert_eq!(unit_vec.y, 0.0);
        assert_eq!(unit_vec.z, 0.0);

        let recovered = GCRSPosition::from_unit_vector(unit_vec, epoch).unwrap();
        assert_eq!(recovered.ra().degrees(), vernal_equinox.ra().degrees());
        assert_eq!(recovered.dec().degrees(), vernal_equinox.dec().degrees());
    }

    #[test]
    fn test_icrs_to_gcrs_applies_aberration() {
        let epoch = TT::j2000();
        let icrs = ICRSPosition::from_degrees(90.0, 23.0).unwrap();

        let gcrs = GCRSPosition::from_icrs(&icrs, &epoch).unwrap();

        let sep_arcsec = icrs
            .angular_separation(&ICRSPosition::from_unit_vector(gcrs.unit_vector()).unwrap())
            .arcseconds();

        assert!(
            sep_arcsec > 15.0 && sep_arcsec < 25.0,
            "ICRS→GCRS aberration should be ~20 arcsec, got {:.2} arcsec",
            sep_arcsec
        );
    }

    #[test]
    fn test_gcrs_to_icrs_roundtrip() {
        let test_positions = [
            (0.0, 0.0),
            (90.0, 0.0),
            (180.0, 45.0),
            (270.0, -60.0),
            (45.0, 89.0),
        ];

        for (ra, dec) in test_positions {
            let epoch = TT::j2000();
            let icrs = ICRSPosition::from_degrees(ra, dec).unwrap();
            let gcrs = GCRSPosition::from_icrs(&icrs, &epoch).unwrap();
            let recovered = gcrs.to_icrs(&epoch).unwrap();

            // Iterative inverse aberration gives ~70 nano-arcsec precision
            let diff_arcsec = icrs.angular_separation(&recovered).arcseconds();
            assert!(
                diff_arcsec < 1e-7,
                "Roundtrip for ({}, {}) should be < 100 nano-arcsec, got {:.2e} arcsec",
                ra,
                dec,
                diff_arcsec
            );
        }
    }

    #[test]
    fn test_gcrs_to_cirs_to_gcrs_roundtrip() {
        let epoch = TT::j2000();
        let original = GCRSPosition::from_degrees(120.0, 30.0, epoch).unwrap();

        let cirs = original.to_cirs().unwrap();
        let recovered = GCRSPosition::from_cirs(&cirs).unwrap();

        // GCRS→CIRS→GCRS is just matrix multiplication (transpose is exact inverse)
        // This should be very precise - use angular separation for robustness
        let sep_arcsec = {
            let orig_vec = original.unit_vector();
            let rec_vec = recovered.unit_vector();
            let dot = orig_vec.x * rec_vec.x + orig_vec.y * rec_vec.y + orig_vec.z * rec_vec.z;
            dot.clamp(-1.0, 1.0).acos() * 206264.806247
        };
        assert!(
            sep_arcsec < 1e-10,
            "GCRS→CIRS→GCRS roundtrip should be < 0.1 nano-arcsec, got {:.2e} arcsec",
            sep_arcsec
        );
    }

    #[test]
    fn test_icrs_to_gcrs_to_cirs_chain() {
        // Note: GCRS applies only aberration, while CIRS (from ICRS) also applies
        // gravitational light deflection by the Sun. The difference between paths
        // is the light deflection effect, which can be up to ~1.75" at the solar limb.
        let epoch = TT::j2000();
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        let gcrs = GCRSPosition::from_icrs(&icrs, &epoch).unwrap();
        let cirs_via_gcrs = gcrs.to_cirs().unwrap();

        let cirs_direct = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();

        // Light deflection causes a difference between paths.
        // For typical stars not near the Sun, this is ~1-10 mas.
        let ra_diff_arcsec =
            (cirs_via_gcrs.ra().radians() - cirs_direct.ra().radians()).abs() * 206264.806247;
        let dec_diff_arcsec =
            (cirs_via_gcrs.dec().radians() - cirs_direct.dec().radians()).abs() * 206264.806247;

        // Light deflection should be < 0.1" for stars far from the Sun
        assert!(
            ra_diff_arcsec < 0.1,
            "RA difference (light deflection) should be < 0.1\", got {:.4}\"",
            ra_diff_arcsec
        );
        assert!(
            dec_diff_arcsec < 0.1,
            "Dec difference (light deflection) should be < 0.1\", got {:.4}\"",
            dec_diff_arcsec
        );
    }

    #[test]
    fn test_aberration_varies_with_epoch() {
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        let epoch_jan = TT::from_julian_date(astro_time::JulianDate::new(2451545.0, 0.0));
        let epoch_jul = TT::from_julian_date(astro_time::JulianDate::new(2451545.0, 182.5));

        let gcrs_jan = GCRSPosition::from_icrs(&icrs, &epoch_jan).unwrap();
        let gcrs_jul = GCRSPosition::from_icrs(&icrs, &epoch_jul).unwrap();

        let ra_diff_arcsec = (gcrs_jan.ra().degrees() - gcrs_jul.ra().degrees()).abs() * 3600.0;
        let dec_diff_arcsec = (gcrs_jan.dec().degrees() - gcrs_jul.dec().degrees()).abs() * 3600.0;

        assert!(
            ra_diff_arcsec > 1.0 || dec_diff_arcsec > 1.0,
            "Aberration should differ between epochs: RA={:.2}\", Dec={:.2}\"",
            ra_diff_arcsec,
            dec_diff_arcsec
        );
    }

    #[test]
    fn test_distance_preservation() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(100.0).unwrap();
        let icrs = ICRSPosition::from_degrees_with_distance(90.0, 45.0, distance).unwrap();

        let gcrs = GCRSPosition::from_icrs(&icrs, &epoch).unwrap();
        assert_eq!(gcrs.distance().unwrap(), distance);

        let cirs = gcrs.to_cirs().unwrap();
        assert_eq!(cirs.distance().unwrap(), distance);

        let recovered_gcrs = GCRSPosition::from_cirs(&cirs).unwrap();
        assert_eq!(recovered_gcrs.distance().unwrap(), distance);

        let recovered_icrs = recovered_gcrs.to_icrs(&epoch).unwrap();
        assert_eq!(recovered_icrs.distance().unwrap(), distance);
    }

    #[test]
    fn test_coordinate_validation() {
        let epoch = TT::j2000();

        assert!(GCRSPosition::from_degrees(0.0, 0.0, epoch).is_ok());
        assert!(GCRSPosition::from_degrees(359.99, 89.99, epoch).is_ok());

        assert!(GCRSPosition::from_degrees(0.0, 91.0, epoch).is_err());
        assert!(GCRSPosition::from_degrees(0.0, -91.0, epoch).is_err());
    }

    #[test]
    fn test_display_formatting() {
        let epoch = TT::j2000();
        let pos = GCRSPosition::from_degrees(123.456789, -67.123456, epoch).unwrap();
        let display = format!("{}", pos);

        assert!(display.contains("GCRS"));
        assert!(display.contains("RA=123.456789°"));
        assert!(display.contains("Dec=-67.123456°"));
        assert!(display.contains("J2000.0"));
    }

    #[test]
    fn test_from_unit_vector_north_pole() {
        let epoch = TT::j2000();
        let north_pole_vec = Vector3::new(0.0, 0.0, 1.0);
        let pos = GCRSPosition::from_unit_vector(north_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_from_unit_vector_south_pole() {
        let epoch = TT::j2000();
        let south_pole_vec = Vector3::new(0.0, 0.0, -1.0);
        let pos = GCRSPosition::from_unit_vector(south_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), -std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_pole_roundtrip() {
        let epoch = TT::j2000();

        let north_pole = GCRSPosition::from_degrees(0.0, 90.0, epoch).unwrap();
        let unit_vec = north_pole.unit_vector();
        let recovered = GCRSPosition::from_unit_vector(unit_vec, epoch).unwrap();

        assert_eq!(recovered.ra().radians(), 0.0);
        assert_eq!(recovered.dec().degrees(), 90.0);

        let south_pole = GCRSPosition::from_degrees(0.0, -90.0, epoch).unwrap();
        let unit_vec = south_pole.unit_vector();
        let recovered = GCRSPosition::from_unit_vector(unit_vec, epoch).unwrap();

        assert_eq!(recovered.ra().radians(), 0.0);
        assert_eq!(recovered.dec().degrees(), -90.0);
    }
}
