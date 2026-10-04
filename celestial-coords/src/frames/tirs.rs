use crate::astrom::earth::EarthRotation;
use crate::eop::record::EopParameters;
use crate::errors::CoordError;
use crate::errors::CoordResult;
use crate::frames::cirs::CIRSPosition;
use crate::frames::itrs::ITRSPosition;
use celestial_core::matrix::Vector3;
use celestial_time::scales::tt::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TIRSPosition {
    x: f64,
    y: f64,
    z: f64,
    epoch: TT,
}

impl TIRSPosition {
    pub fn new(x: f64, y: f64, z: f64, epoch: TT) -> CoordResult<Self> {
        if !(x.is_finite() && y.is_finite() && z.is_finite()) {
            return Err(CoordError::invalid_coordinate(
                "TIRS position must be finite",
            ));
        }
        Ok(Self { x, y, z, epoch })
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn position_vector(&self) -> Vector3 {
        Vector3::new(self.x, self.y, self.z)
    }

    pub fn from_position_vector(pos: Vector3, epoch: TT) -> CoordResult<Self> {
        Self::new(pos.x, pos.y, pos.z, epoch)
    }

    pub fn geocentric_distance(&self) -> f64 {
        self.position_vector().magnitude()
    }

    pub fn distance_to(&self, other: &Self) -> f64 {
        (self.position_vector() - other.position_vector()).magnitude()
    }

    pub fn to_cirs(&self, eop: &EopParameters) -> CoordResult<CIRSPosition> {
        EarthRotation::new(&self.epoch, eop)?.tirs_to_cirs(self)
    }

    pub fn to_itrs(&self, eop: &EopParameters) -> CoordResult<ITRSPosition> {
        EarthRotation::new(&self.epoch, eop)?.tirs_to_itrs(self)
    }
}

impl std::fmt::Display for TIRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TIRS(X={:.3}, Y={:.3}, Z={:.3}, epoch=J{:.1})",
            self.x,
            self.y,
            self.z,
            self.epoch.julian_year()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT, WGS84_SEMI_MAJOR_AXIS};
    use celestial_core::test_helpers::assert_ulp_le;

    const MJD_J2000: f64 = J2000_JD - MJD_ZERO_POINT;

    #[test]
    fn test_tirs_creation() {
        let epoch = TT::j2000();
        let pos = TIRSPosition::new(1000000.0, 2000000.0, 3000000.0, epoch).unwrap();

        assert_eq!(pos.x(), 1000000.0);
        assert_eq!(pos.y(), 2000000.0);
        assert_eq!(pos.z(), 3000000.0);
        assert_eq!(pos.epoch(), epoch);
    }

    #[test]
    fn test_vector_operations() {
        let epoch = TT::j2000();
        let original = TIRSPosition::new(1000.0, 2000.0, 3000.0, epoch).unwrap();

        let vec = original.position_vector();
        assert_eq!(vec.x, 1000.0);
        assert_eq!(vec.y, 2000.0);
        assert_eq!(vec.z, 3000.0);

        let recovered = TIRSPosition::from_position_vector(vec, epoch).unwrap();
        assert_eq!(recovered.x(), original.x());
        assert_eq!(recovered.y(), original.y());
        assert_eq!(recovered.z(), original.z());
    }

    #[test]
    fn test_distance_calculations() {
        let epoch = TT::j2000();

        let pos1 = TIRSPosition::new(1000.0, 0.0, 0.0, epoch).unwrap();
        let pos2 = TIRSPosition::new(2000.0, 0.0, 0.0, epoch).unwrap();

        assert_eq!(pos1.distance_to(&pos2), 1000.0);
        assert_eq!(pos2.distance_to(&pos1), 1000.0);

        assert_eq!(pos1.geocentric_distance(), 1000.0);
        assert_eq!(pos2.geocentric_distance(), 2000.0);
    }

    #[test]
    fn test_itrs_transformation_roundtrip() {
        use crate::eop::record::EopRecord;

        let epoch = TT::j2000();
        let original_tirs = TIRSPosition::new(4000000.0, 3000000.0, 5000000.0, epoch).unwrap();

        // Without polar motion the two frames coincide at J2000, so give it some to undo.
        let eop = EopRecord::new(MJD_J2000, 0.1, 0.2, 0.3)
            .unwrap()
            .to_parameters();

        let itrs = original_tirs.to_itrs(&eop).unwrap();
        let recovered_tirs = itrs.to_tirs(&eop).unwrap();

        assert_ulp_le(recovered_tirs.x(), original_tirs.x(), 1, "X roundtrip");
        assert_ulp_le(recovered_tirs.y(), original_tirs.y(), 1, "Y roundtrip");
        assert_ulp_le(recovered_tirs.z(), original_tirs.z(), 1, "Z roundtrip");
    }

    #[test]
    fn test_itrs_is_tirs_without_polar_motion_at_j2000() {
        use crate::eop::record::EopRecord;

        // s' is zero at J2000, so with no polar motion the two frames coincide.
        let epoch = TT::j2000();
        let tirs = TIRSPosition::new(1000000.0, -2000000.0, 3000000.0, epoch).unwrap();
        let eop = EopRecord::new(MJD_J2000, 0.0, 0.0, 0.3)
            .unwrap()
            .to_parameters();

        let itrs = tirs.to_itrs(&eop).unwrap();

        assert_eq!(itrs.position_vector(), tirs.position_vector());
    }

    #[test]
    fn test_polar_motion_moves_the_pole_toward_x_p() {
        use crate::eop::record::EopRecord;

        // The CIP sits at (x_p, -y_p) in the ITRS.
        let epoch = TT::j2000();
        let pole = TIRSPosition::new(0.0, 0.0, 1.0, epoch).unwrap();
        let eop = EopRecord::new(MJD_J2000, 0.3, 0.0, 0.3)
            .unwrap()
            .to_parameters();
        let xp = eop.x_p * celestial_core::constants::ARCSEC_TO_RAD;

        let itrs = pole.to_itrs(&eop).unwrap();

        assert_eq!(
            itrs.position_vector(),
            Vector3::new(libm::sin(xp), 0.0, libm::cos(xp))
        );
    }

    #[test]
    fn test_display_formatting() {
        let pos = TIRSPosition::new(1234567.89, -987654.32, 555666.77, TT::j2000()).unwrap();
        assert_eq!(
            pos.to_string(),
            "TIRS(X=1234567.890, Y=-987654.320, Z=555666.770, epoch=J2000.0)"
        );
    }

    #[test]
    fn test_non_finite_position_is_err() {
        let epoch = TT::j2000();
        for (x, y, z) in [
            (f64::NAN, 0.0, 0.0),
            (0.0, f64::INFINITY, 0.0),
            (1.0, 0.0, f64::NEG_INFINITY),
        ] {
            assert!(
                TIRSPosition::new(x, y, z, epoch).is_err(),
                "({x}, {y}, {z})"
            );
        }
        let vector = Vector3::new(0.0, 0.0, f64::NAN);
        assert!(TIRSPosition::from_position_vector(vector, epoch).is_err());
    }

    #[test]
    fn test_itrs_consistency() {
        use crate::eop::record::EopRecord;

        let epoch = TT::j2000();
        let itrs_original = ITRSPosition::new(2000000.0, 1000000.0, 6000000.0, epoch).unwrap();

        let eop = EopRecord::new(MJD_J2000, 0.1, 0.2, 0.3)
            .unwrap()
            .to_parameters();

        let tirs = itrs_original.to_tirs(&eop).unwrap();
        let itrs_recovered = tirs.to_itrs(&eop).unwrap();

        assert_ulp_le(itrs_recovered.x(), itrs_original.x(), 1, "ITRS X roundtrip");
        assert_ulp_le(itrs_recovered.y(), itrs_original.y(), 1, "ITRS Y roundtrip");
        assert_ulp_le(itrs_recovered.z(), itrs_original.z(), 1, "ITRS Z roundtrip");
    }

    #[test]
    fn test_non_finite_eop_is_rejected() {
        use crate::eop::record::{EopParameters, EopRecord};
        use crate::frames::cirs::CIRSPosition;

        let epoch = TT::j2000();
        let good = EopRecord::new(MJD_J2000, 0.1, 0.2, 0.3)
            .unwrap()
            .to_parameters();
        let cirs = CIRSPosition::from_degrees(10.0, 20.0, epoch).unwrap();
        let tirs = TIRSPosition::new(1.0, 0.0, 0.0, epoch).unwrap();
        let cases = [
            (
                "mjd",
                EopParameters {
                    mjd: f64::NAN,
                    ..good.clone()
                },
            ),
            (
                "ut1_utc",
                EopParameters {
                    ut1_utc: f64::NAN,
                    ..good.clone()
                },
            ),
            (
                "x_p",
                EopParameters {
                    x_p: f64::NAN,
                    ..good.clone()
                },
            ),
            (
                "y_p",
                EopParameters {
                    y_p: f64::INFINITY,
                    ..good.clone()
                },
            ),
        ];
        for (field, eop) in cases {
            assert!(cirs.to_tirs(&eop).is_err(), "CIRS to TIRS, {field}");
            assert!(tirs.to_cirs(&eop).is_err(), "TIRS to CIRS, {field}");
            assert!(tirs.to_itrs(&eop).is_err(), "TIRS to ITRS, {field}");
        }
    }

    #[test]
    fn test_cirs_transformation() {
        use crate::eop::record::EopRecord;

        let epoch = TT::j2000();
        let eop = EopRecord::new(MJD_J2000, 0.0, 0.0, 0.3)
            .unwrap()
            .to_parameters();

        let tirs = TIRSPosition::new(WGS84_SEMI_MAJOR_AXIS, 0.0, 0.0, epoch).unwrap();

        let cirs = tirs.to_cirs(&eop).unwrap();

        assert!(cirs.ra().degrees() >= 0.0);
        assert!(cirs.ra().degrees() < 360.0);
        assert!(cirs.dec().degrees() >= -90.0);
        assert!(cirs.dec().degrees() <= 90.0);
    }
}
