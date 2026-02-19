use crate::CoordResult;
use celestial_core::{Angle, Vector3};
use celestial_time::TT;

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
    pub fn new(x: f64, y: f64, z: f64, epoch: TT) -> Self {
        Self { x, y, z, epoch }
    }

    pub fn from_geodetic(
        longitude: Angle,
        latitude: Angle,
        height: f64,
        epoch: TT,
    ) -> CoordResult<Self> {
        const A: f64 = 6378137.0;
        const F: f64 = 1.0 / 298.257223563;

        let (sin_lat, cos_lat) = latitude.sin_cos();
        let (sin_lon, cos_lon) = longitude.sin_cos();

        let w = 1.0 - F;
        let w2 = w * w;
        let d = cos_lat * cos_lat + w2 * sin_lat * sin_lat;
        let ac = A / libm::sqrt(d);
        let a_s = w2 * ac;

        let r = (ac + height) * cos_lat;
        let x = r * cos_lon;
        let y = r * sin_lon;
        let z = (a_s + height) * sin_lat;

        Ok(Self::new(x, y, z, epoch))
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

    pub fn from_position_vector(pos: Vector3, epoch: TT) -> Self {
        Self::new(pos.x, pos.y, pos.z, epoch)
    }

    pub fn to_geodetic(&self) -> CoordResult<(Angle, Angle, f64)> {
        const A: f64 = 6378137.0;
        const F: f64 = 1.0 / 298.257223563;
        const B: f64 = A * (1.0 - F);
        const E2: f64 = F * (2.0 - F);
        let p = libm::sqrt(self.x * self.x + self.y * self.y);
        let longitude = libm::atan2(self.y, self.x);

        let theta = libm::atan2(self.z * A, p * B);
        let (sin_theta, cos_theta) = libm::sincos(theta);
        let ep2 = E2 / (1.0 - E2);
        let mut latitude = libm::atan2(
            self.z + ep2 * B * sin_theta.powi(3),
            p - E2 * A * cos_theta.powi(3),
        );
        let mut height = 0.0;

        for _ in 0..5 {
            let (sin_lat, cos_lat) = libm::sincos(latitude);
            let n = A / libm::sqrt(1.0 - E2 * sin_lat * sin_lat);
            height = p / cos_lat - n;
            latitude = libm::atan2(self.z, p * (1.0 - E2 * n / (n + height)));
        }

        Ok((
            Angle::from_radians(longitude),
            Angle::from_radians(latitude),
            height,
        ))
    }

    pub fn geocentric_distance(&self) -> f64 {
        libm::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }

    pub fn distance_to(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        libm::sqrt(dx * dx + dy * dy + dz * dz)
    }

    pub fn to_tirs(
        &self,
        epoch: &TT,
        eop: &crate::eop::EopParameters,
    ) -> CoordResult<crate::frames::TIRSPosition> {
        use crate::frames::TIRSPosition;
        TIRSPosition::from_itrs(self, epoch, eop)
    }
}

impl std::fmt::Display for ITRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ITRS(X={:.3}m, Y={:.3}m, Z={:.3}m, epoch=J{:.1})",
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
        let pos = ITRSPosition::new(1000000.0, 2000000.0, 3000000.0, epoch);

        assert_eq!(pos.x(), 1000000.0);
        assert_eq!(pos.y(), 2000000.0);
        assert_eq!(pos.z(), 3000000.0);
        assert_eq!(pos.epoch(), epoch);
    }

    #[test]
    fn test_vector_operations() {
        let epoch = TT::j2000();
        let original = ITRSPosition::new(1000.0, 2000.0, 3000.0, epoch);

        let vec = original.position_vector();
        assert_eq!(vec.x, 1000.0);
        assert_eq!(vec.y, 2000.0);
        assert_eq!(vec.z, 3000.0);

        let recovered = ITRSPosition::from_position_vector(vec, epoch);
        assert_eq!(recovered.x(), original.x());
        assert_eq!(recovered.y(), original.y());
        assert_eq!(recovered.z(), original.z());
    }

    #[test]
    fn test_geodetic_conversion_roundtrip() {
        let epoch = TT::j2000();

        // Test known location: Greenwich Observatory
        let greenwich_lon = Angle::from_degrees(0.0);
        let greenwich_lat = Angle::from_degrees(51.4769);
        let greenwich_height = 47.0; // meters

        let itrs =
            ITRSPosition::from_geodetic(greenwich_lon, greenwich_lat, greenwich_height, epoch)
                .unwrap();

        let (lon, lat, height) = itrs.to_geodetic().unwrap();

        // Test roundtrip accuracy (should be exact for this conversion)
        assert_eq!(lon.degrees(), greenwich_lon.degrees());
        assert_eq!(lat.degrees(), greenwich_lat.degrees());
        assert_eq!(height, greenwich_height);
    }

    #[test]
    fn test_geodetic_conversion_equator() {
        let epoch = TT::j2000();

        // Test equatorial position
        let pos = ITRSPosition::from_geodetic(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            0.0,
            epoch,
        )
        .unwrap();

        // Should be exactly at Earth's equatorial radius
        const A: f64 = 6378137.0;
        assert_eq!(pos.x(), A);
        assert_eq!(pos.y(), 0.0);
        assert_eq!(pos.z(), 0.0);
    }

    #[test]
    fn test_distance_calculations() {
        let epoch = TT::j2000();

        let pos1 = ITRSPosition::new(1000.0, 0.0, 0.0, epoch);
        let pos2 = ITRSPosition::new(2000.0, 0.0, 0.0, epoch);

        assert_eq!(pos1.distance_to(&pos2), 1000.0);
        assert_eq!(pos2.distance_to(&pos1), 1000.0);

        assert_eq!(pos1.geocentric_distance(), 1000.0);
        assert_eq!(pos2.geocentric_distance(), 2000.0);
    }

    #[test]
    fn test_display_formatting() {
        let epoch = TT::j2000();
        let pos = ITRSPosition::new(1234567.89, -987654.32, 555666.77, epoch);

        let display = format!("{}", pos);
        assert!(display.contains("ITRS"));
        assert!(display.contains("1234567.890m"));
        assert!(display.contains("-987654.320m"));
        assert!(display.contains("555666.770m"));
        assert!(display.contains("J2000.0"));
    }
}
