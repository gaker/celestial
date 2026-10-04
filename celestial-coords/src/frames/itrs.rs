use crate::astrom::earth::EarthRotation;
use crate::eop::record::EopParameters;
use crate::errors::CoordError;
use crate::errors::CoordResult;
use crate::frames::tirs::TIRSPosition;
use celestial_core::constants::{HALF_PI, WGS84_FLATTENING, WGS84_SEMI_MAJOR_AXIS};
use celestial_core::{angle::Angle, matrix::Vector3};
use celestial_time::scales::tt::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ITRSPosition {
    x: f64,
    y: f64,
    z: f64,
    epoch: TT,
}

impl ITRSPosition {
    pub fn new(x: f64, y: f64, z: f64, epoch: TT) -> CoordResult<Self> {
        if !(x.is_finite() && y.is_finite() && z.is_finite()) {
            return Err(CoordError::invalid_coordinate(
                "ITRS position must be finite",
            ));
        }
        Ok(Self { x, y, z, epoch })
    }

    pub fn from_geodetic(
        longitude: Angle,
        latitude: Angle,
        height: f64,
        epoch: TT,
    ) -> CoordResult<Self> {
        let latitude = latitude.validate_latitude()?;
        if !(longitude.radians().is_finite() && height.is_finite()) {
            return Err(CoordError::invalid_coordinate(
                "geodetic longitude and height must be finite",
            ));
        }
        let (x, y, z) = geocentric(longitude, latitude, height);
        Self::new(x, y, z, epoch)
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

    pub fn to_geodetic(&self) -> CoordResult<(Angle, Angle, f64)> {
        let (x, y, z) = (self.x, self.y, self.z);
        if !(x.is_finite() && y.is_finite() && z.is_finite()) {
            return Err(CoordError::invalid_coordinate(
                "ITRS position must be finite",
            ));
        }
        if x == 0.0 && y == 0.0 && z == 0.0 {
            return Err(CoordError::invalid_coordinate(
                "geodetic coordinates are undefined at the geocentre",
            ));
        }
        let p2 = x * x + y * y;
        let longitude = if p2 > 0.0 { libm::atan2(y, x) } else { 0.0 };
        let (latitude, height) = geodetic_latitude_height(p2, libm::fabs(z));
        let latitude = if z < 0.0 { -latitude } else { latitude };
        Ok((
            Angle::from_radians(longitude),
            Angle::from_radians(latitude),
            height,
        ))
    }

    pub fn geocentric_distance(&self) -> f64 {
        self.position_vector().magnitude()
    }

    pub fn distance_to(&self, other: &Self) -> f64 {
        (self.position_vector() - other.position_vector()).magnitude()
    }

    pub fn to_tirs(&self, eop: &EopParameters) -> CoordResult<TIRSPosition> {
        EarthRotation::new(&self.epoch, eop)?.itrs_to_tirs(self)
    }
}

fn geocentric(longitude: Angle, latitude: Angle, height: f64) -> (f64, f64, f64) {
    let (sin_lat, cos_lat) = latitude.sin_cos();
    let (sin_lon, cos_lon) = longitude.sin_cos();
    let w = 1.0 - WGS84_FLATTENING;
    let w2 = w * w;
    let d = cos_lat * cos_lat + w2 * sin_lat * sin_lat;
    let ac = WGS84_SEMI_MAJOR_AXIS / libm::sqrt(d);
    let a_s = w2 * ac;
    let r = (ac + height) * cos_lat;
    (r * cos_lon, r * sin_lon, (a_s + height) * sin_lat)
}

// Fukushima (2006), J. Geodesy 79, 689: one Newton step with a Halley
// correction, for the latitude magnitude and height given the squared
// distance from the polar axis and |z|. On the axis the latitude is 90 deg.
fn geodetic_latitude_height(p2: f64, abs_z: f64) -> (f64, f64) {
    let a = WGS84_SEMI_MAJOR_AXIS;
    let e2 = (2.0 - WGS84_FLATTENING) * WGS84_FLATTENING;
    let ec2 = 1.0 - e2;
    let ec = libm::sqrt(ec2);
    if p2 <= a * a * 1e-32 {
        return (HALF_PI, abs_z - a * ec);
    }
    let p = libm::sqrt(p2);
    let (s1, cc) = halley_corrected(p / a, abs_z / a, e2, ec);
    let s12 = s1 * s1;
    let cc2 = cc * cc;
    let height = (p * cc + abs_z * s1 - a * libm::sqrt(ec2 * s12 + cc2)) / libm::sqrt(s12 + cc2);
    (libm::atan(s1 / cc), height)
}

// The sine and cosine numerators of the latitude, from the normalised
// distances from the axis (pn) and the equator (s0).
fn halley_corrected(pn: f64, s0: f64, e2: f64, ec: f64) -> (f64, f64) {
    let zc = ec * s0;
    let c0 = ec * pn;
    let c02 = c0 * c0;
    let c03 = c02 * c0;
    let s02 = s0 * s0;
    let s03 = s02 * s0;
    let a02 = c02 + s02;
    let a0 = libm::sqrt(a02);
    let a03 = a02 * a0;
    let d0 = zc * a03 + e2 * s03;
    let f0 = pn * a03 - e2 * c03;
    let b0 = e2 * e2 * 1.5 * s02 * c02 * pn * (a0 - ec);
    (d0 * f0 - b0 * s0, ec * (f0 * f0 - b0 * c0))
}

impl std::fmt::Display for ITRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ITRS(X={:.3}, Y={:.3}, Z={:.3}, epoch=J{:.1})",
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

    #[test]
    fn test_itrs_creation() {
        let epoch = TT::j2000();
        let pos = ITRSPosition::new(1000000.0, 2000000.0, 3000000.0, epoch).unwrap();

        assert_eq!(pos.x(), 1000000.0);
        assert_eq!(pos.y(), 2000000.0);
        assert_eq!(pos.z(), 3000000.0);
        assert_eq!(pos.epoch(), epoch);
    }

    #[test]
    fn test_vector_operations() {
        let epoch = TT::j2000();
        let original = ITRSPosition::new(1000.0, 2000.0, 3000.0, epoch).unwrap();

        let vec = original.position_vector();
        assert_eq!(vec.x, 1000.0);
        assert_eq!(vec.y, 2000.0);
        assert_eq!(vec.z, 3000.0);

        let recovered = ITRSPosition::from_position_vector(vec, epoch).unwrap();
        assert_eq!(recovered.x(), original.x());
        assert_eq!(recovered.y(), original.y());
        assert_eq!(recovered.z(), original.z());
    }

    #[test]
    fn test_geodetic_conversion_roundtrip() {
        let epoch = TT::j2000();

        // Greenwich Observatory
        let greenwich_lon = Angle::from_degrees(0.0);
        let greenwich_lat = Angle::from_degrees(51.4769);
        let greenwich_height = 47.0; // meters

        let itrs =
            ITRSPosition::from_geodetic(greenwich_lon, greenwich_lat, greenwich_height, epoch)
                .unwrap();

        let (lon, lat, height) = itrs.to_geodetic().unwrap();

        // All values equal ERFA's gd2gc and gc2gd outputs. gc2gd's single
        // Halley step recovers the latitude one ulp above the input and the
        // height 1.2 nm low.
        assert_eq!(greenwich_lat.radians(), 0.8984413937198691);
        assert_eq!(
            (itrs.x(), itrs.y(), itrs.z()),
            (3980688.823882222, 0.0, 4966798.928379789)
        );
        assert_eq!(lon.degrees(), greenwich_lon.degrees());
        assert_eq!(lat.radians(), 0.8984413937198692);
        assert_eq!(height, 46.999999998804135);
    }

    #[test]
    fn test_geodetic_conversion_equator() {
        let epoch = TT::j2000();

        let pos = ITRSPosition::from_geodetic(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            0.0,
            epoch,
        )
        .unwrap();

        assert_eq!(pos.x(), WGS84_SEMI_MAJOR_AXIS);
        assert_eq!(pos.y(), 0.0);
        assert_eq!(pos.z(), 0.0);
    }

    #[test]
    fn test_to_geodetic_at_the_poles() {
        let epoch = TT::j2000();
        let polar_radius = 6356752.314245179;
        let north = ITRSPosition::new(0.0, 0.0, polar_radius, epoch).unwrap();
        let (lon, lat, height) = north.to_geodetic().unwrap();
        assert_eq!((lon.radians(), lat.radians(), height), (0.0, HALF_PI, 0.0));

        let below_south = ITRSPosition::new(0.0, 0.0, -polar_radius + 1000.0, epoch).unwrap();
        let (_, lat, height) = below_south.to_geodetic().unwrap();
        assert_eq!((lat.radians(), height), (-HALF_PI, -1000.0));
    }

    #[test]
    fn test_to_geodetic_rejects_the_geocentre() {
        let pos = ITRSPosition::new(0.0, 0.0, 0.0, TT::j2000()).unwrap();
        assert!(pos.to_geodetic().is_err());
    }

    #[test]
    fn test_from_geodetic_rejects_invalid_input() {
        let epoch = TT::j2000();
        let geodetic = |lon: f64, lat: f64, height: f64| {
            ITRSPosition::from_geodetic(
                Angle::from_degrees(lon),
                Angle::from_degrees(lat),
                height,
                epoch,
            )
        };
        assert!(geodetic(0.0, 120.0, 0.0).is_err());
        assert!(geodetic(0.0, -90.5, 0.0).is_err());
        assert!(geodetic(f64::NAN, 45.0, 0.0).is_err());
        assert!(geodetic(0.0, 45.0, f64::INFINITY).is_err());
        assert!(geodetic(370.0, 90.0, 400_000.0).is_ok());
    }

    #[test]
    fn test_distance_calculations() {
        let epoch = TT::j2000();

        let pos1 = ITRSPosition::new(1000.0, 0.0, 0.0, epoch).unwrap();
        let pos2 = ITRSPosition::new(2000.0, 0.0, 0.0, epoch).unwrap();

        assert_eq!(pos1.distance_to(&pos2), 1000.0);
        assert_eq!(pos2.distance_to(&pos1), 1000.0);

        assert_eq!(pos1.geocentric_distance(), 1000.0);
        assert_eq!(pos2.geocentric_distance(), 2000.0);
    }

    #[test]
    fn test_display_formatting() {
        let pos = ITRSPosition::new(1234567.89, -987654.32, 555666.77, TT::j2000()).unwrap();
        assert_eq!(
            pos.to_string(),
            "ITRS(X=1234567.890, Y=-987654.320, Z=555666.770, epoch=J2000.0)"
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
                ITRSPosition::new(x, y, z, epoch).is_err(),
                "({x}, {y}, {z})"
            );
        }
        let vector = Vector3::new(0.0, 0.0, f64::NAN);
        assert!(ITRSPosition::from_position_vector(vector, epoch).is_err());
    }
}
