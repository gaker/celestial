use crate::astrom::Astrom;
use crate::distance::Distance;
use crate::errors::CoordResult;
use crate::frames::cirs::CIRSPosition;
use crate::frames::direction::spherical_angles;
use crate::frames::icrs::ICRSPosition;
use crate::transforms::CoordinateFrame;
use celestial_core::{angle::Angle, matrix::Vector3};
use celestial_time::scales::tt::TT;

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
        let (ra, dec) = spherical_angles(unit)?;
        Self::new(ra, dec, epoch)
    }

    pub fn to_cirs(&self) -> CoordResult<CIRSPosition> {
        Astrom::new(&self.epoch)?.gcrs_to_cirs(self)
    }

    pub fn from_cirs(cirs: &CIRSPosition) -> CoordResult<Self> {
        Astrom::new(&cirs.epoch())?.cirs_to_gcrs(cirs)
    }
}

impl CoordinateFrame for GCRSPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        Astrom::new(&self.epoch)?.gcrs_to_icrs(self)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        Astrom::new(epoch)?.icrs_to_gcrs(icrs)
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
    use celestial_core::constants::J2000_JD;
    use celestial_core::test_helpers::assert_ulp_le;

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

    // eraAtciq and eraAticq at J2000 with the bias-precession-nutation matrix set to the
    // identity, which leaves light deflection and aberration: (ra, dec) in degrees, the GCRS
    // place and the ICRS place recovered from it, in radians. The inverse iterates, so the
    // place comes back up to 1e-13 rad off.
    const ROUND_TRIPS: [(f64, f64, [f64; 2], [f64; 2]); 6] = [
        (
            90.0,
            23.0,
            [1.5709042617137516, 0.4014255858005104],
            [1.5707963267949328, 0.40142572795869574],
        ),
        (
            0.0,
            0.0,
            [6.28316855099383, -7.2646264083301395e-6],
            [6.28318530717942, -7.180864394117496e-14],
        ),
        (
            90.0,
            0.0,
            [1.5708956814203485, -7.2697790260613576e-6],
            [1.5707963267949248, -2.0438846337053613e-15],
        ),
        (
            180.0,
            45.0,
            [3.141616359684239, 0.7853227720634331],
            [3.1415926535898935, 0.7853981633971284],
        ),
        (
            270.0,
            -60.0,
            [4.712190247338286, -1.0471867059848916],
            [4.712388980384646, -1.0471975511965954],
        ),
        (
            45.0,
            89.0,
            [0.788760192471021, 1.5534249056444125],
            [0.7853981633977009, 1.5533430342749595],
        ),
    ];

    #[test]
    fn test_gcrs_to_icrs_roundtrip() {
        let epoch = TT::j2000();
        for (ra, dec, gcrs_place, back) in ROUND_TRIPS {
            let icrs = ICRSPosition::from_degrees(ra, dec).unwrap();
            let gcrs = GCRSPosition::from_icrs(&icrs, &epoch).unwrap();
            let radians = [gcrs.ra().radians(), gcrs.dec().radians()];
            assert_eq!(radians, gcrs_place, "({ra}, {dec})");
            let recovered = gcrs.to_icrs(&epoch).unwrap();
            let radians = [recovered.ra().radians(), recovered.dec().radians()];
            assert_eq!(radians, back, "({ra}, {dec})");
        }
    }

    #[test]
    fn test_gcrs_to_cirs_to_gcrs_roundtrip() {
        let epoch = TT::j2000();
        let original = GCRSPosition::from_degrees(120.0, 30.0, epoch).unwrap();

        let cirs = original.to_cirs().unwrap();
        let recovered = GCRSPosition::from_cirs(&cirs).unwrap();

        // A rotation and its transpose, with a trip through a vector on either side.
        let ra = recovered.ra().radians();
        assert_ulp_le(ra, original.ra().radians(), 2, "ra");
        assert_ulp_le(
            recovered.dec().radians(),
            original.dec().radians(),
            2,
            "dec",
        );
    }

    #[test]
    fn test_icrs_to_gcrs_to_cirs_chain() {
        // GCRS carries the light deflection and aberration, so stopping there on the way to
        // CIRS costs only the extra trip through a vector.
        let epoch = TT::j2000();
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        let gcrs = GCRSPosition::from_icrs(&icrs, &epoch).unwrap();
        let cirs_via_gcrs = gcrs.to_cirs().unwrap();

        let cirs_direct = CIRSPosition::from_icrs(&icrs, &epoch).unwrap();

        let ra = cirs_via_gcrs.ra().radians();
        assert_ulp_le(ra, cirs_direct.ra().radians(), 1, "ra");
        let dec = cirs_via_gcrs.dec().radians();
        assert_ulp_le(dec, cirs_direct.dec().radians(), 1, "dec");
    }

    #[test]
    fn test_aberration_varies_with_epoch() {
        let icrs = ICRSPosition::from_degrees(180.0, 45.0).unwrap();

        let epoch_jan =
            TT::from_julian_date(celestial_time::julian::JulianDate::new(J2000_JD, 0.0));
        let epoch_jul =
            TT::from_julian_date(celestial_time::julian::JulianDate::new(J2000_JD, 182.5));

        let gcrs_jan = GCRSPosition::from_icrs(&icrs, &epoch_jan).unwrap();
        let gcrs_jul = GCRSPosition::from_icrs(&icrs, &epoch_jul).unwrap();

        let ra_diff_arcsec = libm::fabs((gcrs_jan.ra() - gcrs_jul.ra()).arcseconds());
        let dec_diff_arcsec = libm::fabs((gcrs_jan.dec() - gcrs_jul.dec()).arcseconds());

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
        assert_eq!(pos.dec().radians(), celestial_core::constants::HALF_PI);
    }

    #[test]
    fn test_from_unit_vector_south_pole() {
        let epoch = TT::j2000();
        let south_pole_vec = Vector3::new(0.0, 0.0, -1.0);
        let pos = GCRSPosition::from_unit_vector(south_pole_vec, epoch).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), -celestial_core::constants::HALF_PI);
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
