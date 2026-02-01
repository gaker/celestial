use crate::{
    aberration::{
        apply_aberration, apply_light_deflection, compute_earth_state, remove_aberration,
        remove_light_deflection,
    },
    frames::{ICRSPosition, TIRSPosition},
    transforms::CoordinateFrame,
    CoordError, CoordResult, Distance,
};
use astro_core::{matrix::RotationMatrix3, Angle, Vector3};
use astro_time::{transforms::NutationCalculator, TT};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CIRSPosition {
    ra: Angle,
    dec: Angle,
    epoch: TT,
    distance: Option<Distance>,
}

impl CIRSPosition {
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

    fn icrs_to_cirs_matrix(epoch: &TT) -> CoordResult<RotationMatrix3> {
        Self::icrs_to_cirs_matrix_with_eop(epoch, None)
    }

    fn icrs_to_cirs_matrix_with_eop(
        epoch: &TT,
        eop: Option<&crate::eop::EopParameters>,
    ) -> CoordResult<RotationMatrix3> {
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

        let (x, y) = match eop {
            Some(eop) => (
                eop.corrected_cip_x(cio_solution.cip.x),
                eop.corrected_cip_y(cio_solution.cip.y),
            ),
            None => (cio_solution.cip.x, cio_solution.cip.y),
        };

        Ok(astro_core::gcrs_to_cirs_matrix(x, y, cio_solution.s))
    }

    /// Transforms this CIRS position to Terrestrial Intermediate Reference System (TIRS).
    ///
    /// The position vector is scaled by distance (in AU) before transformation. If no distance
    /// is set, a unit vector is used. The resulting TIRS vector will have the same units (AU or
    /// dimensionless) as the input.
    pub fn to_tirs(&self, eop: &crate::eop::EopParameters) -> CoordResult<TIRSPosition> {
        let cirs_vec = if let Some(distance) = self.distance {
            self.unit_vector() * distance.au()
        } else {
            self.unit_vector()
        };

        TIRSPosition::from_cirs(cirs_vec, &self.epoch, eop)
    }

    pub fn to_hour_angle(
        &self,
        observer: &astro_core::Location,
        delta_t: f64,
    ) -> CoordResult<crate::frames::HourAnglePosition> {
        use astro_time::scales::conversions::ToUT1WithDeltaT;
        use astro_time::sidereal::GAST;

        let ut1 = self.epoch.to_ut1_with_delta_t(delta_t)?;
        let gast = GAST::from_ut1_and_tt(&ut1, &self.epoch)?;

        let last = gast.to_last(observer);

        let ha_rad = last.radians() - self.ra.radians();
        let ha = astro_core::angle::wrap_pm_pi(ha_rad);

        crate::frames::HourAnglePosition::new(
            Angle::from_radians(ha),
            self.dec,
            *observer,
            self.epoch,
        )
    }
}

impl CoordinateFrame for CIRSPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        let icrs_to_cirs = Self::icrs_to_cirs_matrix(&self.epoch)?;
        let cirs_to_icrs = icrs_to_cirs.transpose();

        // Step 1: Remove precession-nutation
        let cirs_vec = self.unit_vector();
        let apparent_vec = cirs_to_icrs * cirs_vec;

        let earth_state = compute_earth_state(&self.epoch)?;
        let sun_earth_dist = earth_state.heliocentric_position.magnitude();

        // Step 2: Remove stellar aberration
        let deflected_vec = remove_aberration(
            apparent_vec,
            earth_state.barycentric_velocity,
            sun_earth_dist,
        );

        // Sun to observer unit vector (heliocentric position normalized)
        let sun_to_earth = Vector3::new(
            earth_state.heliocentric_position.x / sun_earth_dist,
            earth_state.heliocentric_position.y / sun_earth_dist,
            earth_state.heliocentric_position.z / sun_earth_dist,
        );

        // Step 3: Remove gravitational light deflection
        let icrs_vec = remove_light_deflection(deflected_vec, sun_to_earth, sun_earth_dist);

        let mut icrs = ICRSPosition::from_unit_vector(icrs_vec)?;

        if let Some(distance) = self.distance {
            icrs.set_distance(distance);
        }

        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let icrs_to_cirs = Self::icrs_to_cirs_matrix(epoch)?;

        let icrs_vec = icrs.unit_vector();

        let earth_state = compute_earth_state(epoch)?;
        let sun_earth_dist = earth_state.heliocentric_position.magnitude();

        // Sun to observer unit vector (heliocentric position normalized)
        // heliocentric_position is Earth's position relative to Sun, so normalizing gives Sun→Earth direction
        let sun_to_earth = Vector3::new(
            earth_state.heliocentric_position.x / sun_earth_dist,
            earth_state.heliocentric_position.y / sun_earth_dist,
            earth_state.heliocentric_position.z / sun_earth_dist,
        );

        // Step 1: Apply gravitational light deflection by the Sun
        let deflected_vec = apply_light_deflection(icrs_vec, sun_to_earth, sun_earth_dist);

        // Step 2: Apply stellar aberration
        let apparent_vec = apply_aberration(
            deflected_vec,
            earth_state.barycentric_velocity,
            sun_earth_dist,
        );

        // Step 3: Apply precession-nutation (C2I matrix)
        let cirs_vec = icrs_to_cirs * apparent_vec;

        let mut cirs = Self::from_unit_vector(cirs_vec, *epoch)?;

        if let Some(distance) = icrs.distance() {
            cirs.distance = Some(distance);
        }

        Ok(cirs)
    }
}

impl std::fmt::Display for CIRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CIRS(RA={:.6}°, Dec={:.6}°, epoch=J{:.1}",
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
    fn test_cirs_creation() {
        let epoch = TT::j2000();
        let pos = CIRSPosition::from_degrees(180.0, 45.0, epoch).unwrap();

        assert_eq!(pos.ra().degrees(), 180.0);
        assert_eq!(pos.dec().degrees(), 45.0);
        assert_eq!(pos.epoch(), epoch);
        assert_eq!(pos.distance(), None);
    }

    #[test]
    fn test_cirs_with_distance() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(10.0).unwrap();
        let pos = CIRSPosition::with_distance(
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

        // Test vernal equinox direction
        let vernal_equinox = CIRSPosition::from_degrees(0.0, 0.0, epoch).unwrap();
        let unit_vec = vernal_equinox.unit_vector();

        assert_eq!(unit_vec.x, 1.0);
        assert_eq!(unit_vec.y, 0.0);
        assert_eq!(unit_vec.z, 0.0);

        // Test round-trip
        let recovered = CIRSPosition::from_unit_vector(unit_vec, epoch).unwrap();
        assert_eq!(recovered.ra().degrees(), vernal_equinox.ra().degrees());
        assert_eq!(recovered.dec().degrees(), vernal_equinox.dec().degrees());
    }

    #[test]
    fn test_icrs_to_cirs_transformation() {
        let epoch = TT::j2000();
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        // Transform to CIRS
        let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();

        // At J2000, CIRS should be very close to ICRS (small precession/nutation effects)
        // But not exactly equal due to frame bias
        assert!((cirs.ra().degrees() - icrs.ra().degrees()).abs() < 1.0);
        assert!((cirs.dec().degrees() - icrs.dec().degrees()).abs() < 1.0);
    }

    #[test]
    fn test_cirs_to_icrs_roundtrip() {
        let epoch = TT::j2000();
        let original_icrs = ICRSPosition::from_degrees(120.0, 30.0).unwrap();

        // ICRS → CIRS → ICRS
        let cirs = CIRSPosition::from_icrs(&original_icrs, &epoch).unwrap();
        let recovered_icrs = cirs.to_icrs(&epoch).unwrap();

        // Iterative inverse (aberration + light deflection) gives ~50 nano-arcsec precision
        let diff_arcsec = original_icrs
            .angular_separation(&recovered_icrs)
            .arcseconds();
        assert!(
            diff_arcsec < 1e-7,
            "CIRS roundtrip should be < 100 nano-arcsec, got {:.2e} arcsec",
            diff_arcsec
        );
    }

    #[test]
    fn test_coordinate_validation() {
        let epoch = TT::j2000();

        // Valid coordinates
        assert!(CIRSPosition::from_degrees(0.0, 0.0, epoch).is_ok());
        assert!(CIRSPosition::from_degrees(359.99, 89.99, epoch).is_ok());

        // Invalid declination
        assert!(CIRSPosition::from_degrees(0.0, 91.0, epoch).is_err());
        assert!(CIRSPosition::from_degrees(0.0, -91.0, epoch).is_err());
    }

    #[test]
    fn test_distance_preservation() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(100.0).unwrap();
        let icrs = ICRSPosition::from_degrees_with_distance(90.0, 45.0, distance).unwrap();

        // Distance should be preserved through transformation
        let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();
        assert_eq!(cirs.distance().unwrap(), distance);

        let recovered_icrs = cirs.to_icrs(&epoch).unwrap();
        assert_eq!(recovered_icrs.distance().unwrap(), distance);
    }

    #[test]
    fn test_display_formatting() {
        let epoch = TT::j2000();
        let pos = CIRSPosition::from_degrees(123.456789, -67.123456, epoch).unwrap();
        let display = format!("{}", pos);

        assert!(display.contains("CIRS"));
        assert!(display.contains("RA=123.456789°"));
        assert!(display.contains("Dec=-67.123456°"));
        assert!(display.contains("J2000.0"));
    }

    #[test]
    fn test_aberration_applied_in_transformation() {
        let epoch = TT::j2000();
        let icrs = ICRSPosition::from_degrees(90.0, 23.0).unwrap();

        let icrs_vec = icrs.unit_vector();
        let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();
        let cirs_vec = cirs.unit_vector();

        let npb_only = CIRSPosition::icrs_to_cirs_matrix(&epoch).unwrap() * icrs_vec;

        let diff_with_aberr = (cirs_vec.x - npb_only.x).powi(2)
            + (cirs_vec.y - npb_only.y).powi(2)
            + (cirs_vec.z - npb_only.z).powi(2);
        let aberr_arcsec = diff_with_aberr.sqrt() * 206264.806247;

        assert!(
            aberr_arcsec > 15.0 && aberr_arcsec < 25.0,
            "Aberration should be ~20 arcsec, got {:.2} arcsec",
            aberr_arcsec
        );
    }

    #[test]
    fn test_aberration_varies_with_epoch() {
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        let epoch_jan = TT::from_julian_date(astro_time::JulianDate::new(2451545.0, 0.0));
        let epoch_jul = TT::from_julian_date(astro_time::JulianDate::new(2451545.0, 182.5));

        let cirs_jan = CIRSPosition::from_icrs(&icrs, &epoch_jan).unwrap();
        let cirs_jul = CIRSPosition::from_icrs(&icrs, &epoch_jul).unwrap();

        let ra_diff_arcsec = (cirs_jan.ra().degrees() - cirs_jul.ra().degrees()).abs() * 3600.0;
        let dec_diff_arcsec = (cirs_jan.dec().degrees() - cirs_jul.dec().degrees()).abs() * 3600.0;

        assert!(
            ra_diff_arcsec > 1.0 || dec_diff_arcsec > 1.0,
            "Aberration should cause measurable difference between epochs: RA={:.2}\", Dec={:.2}\"",
            ra_diff_arcsec,
            dec_diff_arcsec
        );
    }

    #[test]
    fn test_aberration_roundtrip_precision() {
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
            let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();
            let recovered = cirs.to_icrs(&epoch).unwrap();

            // Iterative inverse (aberration + light deflection) gives ~70 nano-arcsec precision
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
    fn test_from_unit_vector_north_pole() {
        let epoch = TT::j2000();
        let north_pole_vec = Vector3::new(0.0, 0.0, 1.0);
        let pos = CIRSPosition::from_unit_vector(north_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_from_unit_vector_south_pole() {
        let epoch = TT::j2000();
        let south_pole_vec = Vector3::new(0.0, 0.0, -1.0);
        let pos = CIRSPosition::from_unit_vector(south_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), -std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_pole_roundtrip() {
        let epoch = TT::j2000();

        let north_pole = CIRSPosition::from_degrees(0.0, 90.0, epoch).unwrap();
        let unit_vec = north_pole.unit_vector();
        let recovered = CIRSPosition::from_unit_vector(unit_vec, epoch).unwrap();

        assert_eq!(recovered.ra().radians(), 0.0);
        assert_eq!(recovered.dec().degrees(), 90.0);

        let south_pole = CIRSPosition::from_degrees(0.0, -90.0, epoch).unwrap();
        let unit_vec = south_pole.unit_vector();
        let recovered = CIRSPosition::from_unit_vector(unit_vec, epoch).unwrap();

        assert_eq!(recovered.ra().radians(), 0.0);
        assert_eq!(recovered.dec().degrees(), -90.0);
    }
}
