use crate::astrom::earth::EarthRotation;
use crate::astrom::Astrom;
use crate::distance::Distance;
use crate::eop::record::EopParameters;
use crate::errors::CoordResult;
use crate::frames::direction::spherical_angles;
use crate::frames::icrs::ICRSPosition;
use crate::frames::tirs::TIRSPosition;
use crate::frames::topocentric::HourAnglePosition;
use crate::transforms::CoordinateFrame;
use celestial_core::location::Location;
use celestial_core::{angle::Angle, matrix::Vector3};
use celestial_time::scales::tt::TT;

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
        let (ra, dec) = spherical_angles(unit)?;
        Self::new(ra, dec, epoch)
    }

    /// Transforms this CIRS position to Terrestrial Intermediate Reference System (TIRS).
    ///
    /// The position vector is scaled by distance (in AU) before transformation. If no distance
    /// is set, a unit vector is used. The resulting TIRS vector will have the same units (AU or
    /// dimensionless) as the input.
    pub fn to_tirs(&self, eop: &EopParameters) -> CoordResult<TIRSPosition> {
        EarthRotation::new(&self.epoch, eop)?.cirs_to_tirs(self)
    }

    pub fn to_hour_angle(
        &self,
        observer: &Location,
        eop: &EopParameters,
    ) -> CoordResult<HourAnglePosition> {
        EarthRotation::new(&self.epoch, eop)?.cirs_to_hour_angle(self, observer)
    }
}

impl CoordinateFrame for CIRSPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        Astrom::new(&self.epoch)?.cirs_to_icrs(self)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        Astrom::new(epoch)?.icrs_to_cirs(icrs)
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
    use crate::frames::gcrs::GCRSPosition;
    use celestial_core::constants::{ARCSEC_PER_RAD, J2000_JD};

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

        let vernal_equinox = CIRSPosition::from_degrees(0.0, 0.0, epoch).unwrap();
        let unit_vec = vernal_equinox.unit_vector();

        assert_eq!(unit_vec.x, 1.0);
        assert_eq!(unit_vec.y, 0.0);
        assert_eq!(unit_vec.z, 0.0);

        let recovered = CIRSPosition::from_unit_vector(unit_vec, epoch).unwrap();
        assert_eq!(recovered.ra().degrees(), vernal_equinox.ra().degrees());
        assert_eq!(recovered.dec().degrees(), vernal_equinox.dec().degrees());
    }

    #[test]
    fn test_icrs_to_cirs_transformation() {
        let epoch = TT::j2000();
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();
        let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();
        // eraAtci13 with no proper motion or parallax.
        assert_eq!(
            [cirs.ra().radians(), cirs.dec().radians()],
            [3.1415883485583755, 0.7853497187153947]
        );
    }

    // The inverse iterates, so the place comes back up to 1e-13 rad off. These are eraAtic13
    // applied to eraAtci13's output, at J2000: (ra, dec) in degrees and the recovered radians.
    const ROUND_TRIPS: [(f64, f64, [f64; 2]); 6] = [
        (120.0, 30.0, [2.094395102393274, 0.5235987755982817]),
        (0.0, 0.0, [6.283185307179421, -7.180864054642573e-14]),
        (90.0, 0.0, [1.5707963267949248, -2.0438947981007283e-15]),
        (180.0, 45.0, [3.1415926535898935, 0.7853981633971286]),
        (270.0, -60.0, [4.712388980384646, -1.0471975511965956]),
        (45.0, 89.0, [0.7853981633977007, 1.5533430342749595]),
    ];

    #[test]
    fn test_cirs_to_icrs_roundtrip() {
        let epoch = TT::j2000();
        for (ra, dec, back) in ROUND_TRIPS {
            let icrs = ICRSPosition::from_degrees(ra, dec).unwrap();
            let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();
            let recovered = cirs.to_icrs(&epoch).unwrap();
            let radians = [recovered.ra().radians(), recovered.dec().radians()];
            assert_eq!(radians, back, "({ra}, {dec})");
        }
    }

    #[test]
    fn test_coordinate_validation() {
        let epoch = TT::j2000();

        assert!(CIRSPosition::from_degrees(0.0, 0.0, epoch).is_ok());
        assert!(CIRSPosition::from_degrees(359.99, 89.99, epoch).is_ok());

        assert!(CIRSPosition::from_degrees(0.0, 91.0, epoch).is_err());
        assert!(CIRSPosition::from_degrees(0.0, -91.0, epoch).is_err());
    }

    #[test]
    fn test_distance_preservation() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(100.0).unwrap();
        let icrs = ICRSPosition::from_degrees_with_distance(90.0, 45.0, distance).unwrap();

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

        let cirs = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();
        let cirs_vec = cirs.unit_vector();

        let npb_only = GCRSPosition::new(icrs.ra(), icrs.dec(), epoch)
            .unwrap()
            .to_cirs()
            .unwrap()
            .unit_vector();

        let aberr_arcsec = (cirs_vec - npb_only).magnitude() * ARCSEC_PER_RAD;

        assert!(
            aberr_arcsec > 15.0 && aberr_arcsec < 25.0,
            "Aberration should be ~20 arcsec, got {:.2} arcsec",
            aberr_arcsec
        );
    }

    #[test]
    fn test_aberration_varies_with_epoch() {
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        let epoch_jan =
            TT::from_julian_date(celestial_time::julian::JulianDate::new(J2000_JD, 0.0));
        let epoch_jul =
            TT::from_julian_date(celestial_time::julian::JulianDate::new(J2000_JD, 182.5));

        let cirs_jan = CIRSPosition::from_icrs(&icrs, &epoch_jan).unwrap();
        let cirs_jul = CIRSPosition::from_icrs(&icrs, &epoch_jul).unwrap();

        let ra_diff_arcsec = libm::fabs((cirs_jan.ra() - cirs_jul.ra()).arcseconds());
        let dec_diff_arcsec = libm::fabs((cirs_jan.dec() - cirs_jul.dec()).arcseconds());

        assert!(
            ra_diff_arcsec > 1.0 || dec_diff_arcsec > 1.0,
            "Aberration should cause measurable difference between epochs: RA={:.2}\", Dec={:.2}\"",
            ra_diff_arcsec,
            dec_diff_arcsec
        );
    }

    #[test]
    fn test_from_unit_vector_north_pole() {
        let epoch = TT::j2000();
        let north_pole_vec = Vector3::new(0.0, 0.0, 1.0);
        let pos = CIRSPosition::from_unit_vector(north_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), celestial_core::constants::HALF_PI);
    }

    #[test]
    fn test_from_unit_vector_south_pole() {
        let epoch = TT::j2000();
        let south_pole_vec = Vector3::new(0.0, 0.0, -1.0);
        let pos = CIRSPosition::from_unit_vector(south_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), -celestial_core::constants::HALF_PI);
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
