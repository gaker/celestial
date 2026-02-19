use crate::{lunar, transforms::CoordinateFrame, CoordResult, Distance, ICRSPosition};
use celestial_core::constants::HALF_PI;
use celestial_core::matrix::RotationMatrix3;
use celestial_core::utils::normalize_angle_to_positive;
use celestial_core::Angle;
use celestial_time::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SelenographicPosition {
    latitude: Angle,
    longitude: Angle,
    radius: Option<Distance>,
}

impl SelenographicPosition {
    pub fn new(latitude: Angle, longitude: Angle) -> CoordResult<Self> {
        let latitude = latitude.validate_latitude()?;
        let longitude = longitude.validate_longitude(true)?;

        Ok(Self {
            latitude,
            longitude,
            radius: None,
        })
    }

    pub fn with_radius(latitude: Angle, longitude: Angle, radius: Distance) -> CoordResult<Self> {
        let mut pos = Self::new(latitude, longitude)?;
        pos.radius = Some(radius);
        Ok(pos)
    }

    pub fn from_degrees(lat_deg: f64, lon_deg: f64) -> CoordResult<Self> {
        Self::new(Angle::from_degrees(lat_deg), Angle::from_degrees(lon_deg))
    }

    pub fn latitude(&self) -> Angle {
        self.latitude
    }

    pub fn longitude(&self) -> Angle {
        self.longitude
    }

    pub fn radius(&self) -> Option<Distance> {
        self.radius
    }

    pub fn set_radius(&mut self, radius: Distance) {
        self.radius = Some(radius);
    }

    pub fn sub_earth_point(epoch: &TT) -> CoordResult<Self> {
        let (lon, lat) = lunar::compute_sub_earth_point(epoch);
        Self::new(lat, lon)
    }

    pub fn nearside_center() -> Self {
        Self {
            latitude: Angle::ZERO,
            longitude: Angle::ZERO,
            radius: None,
        }
    }

    pub fn farside_center() -> Self {
        Self {
            latitude: Angle::ZERO,
            longitude: Angle::PI,
            radius: None,
        }
    }

    pub fn north_pole() -> Self {
        Self {
            latitude: Angle::HALF_PI,
            longitude: Angle::ZERO,
            radius: None,
        }
    }

    pub fn south_pole() -> Self {
        Self {
            latitude: -Angle::HALF_PI,
            longitude: Angle::ZERO,
            radius: None,
        }
    }

    pub fn angular_separation(&self, other: &Self) -> Angle {
        let (sin_lat1, cos_lat1) = self.latitude.sin_cos();
        let (sin_lat2, cos_lat2) = other.latitude.sin_cos();
        let delta_lon = (self.longitude - other.longitude).radians();

        let angle_rad = celestial_core::math::vincenty_angular_separation(
            sin_lat1, cos_lat1, sin_lat2, cos_lat2, delta_lon,
        );

        Angle::from_radians(angle_rad)
    }

    pub fn is_visible_from_earth(&self, epoch: &TT) -> bool {
        let sub_earth = Self::sub_earth_point(epoch).unwrap_or_else(|_| Self::nearside_center());
        let separation = self.angular_separation(&sub_earth);
        separation.degrees() < 90.0
    }
}

fn selenographic_to_icrs_matrix(epoch: &TT) -> CoordResult<RotationMatrix3> {
    let orientation = lunar::compute_lunar_orientation(epoch);
    let lib_lon = orientation.optical_libration.longitude.radians();
    let lib_lat = orientation.optical_libration.latitude.radians();
    let c = orientation.position_angle.radians();

    let moon_icrs = lunar::get_moon_icrs(epoch)?;
    let moon_ra = moon_icrs.ra().radians();
    let moon_dec = moon_icrs.dec().radians();

    let mut m = RotationMatrix3::identity();
    m.rotate_z(-lib_lon);
    m.rotate_y(-lib_lat);
    m.rotate_z(c);
    m.rotate_y(moon_dec - HALF_PI);
    m.rotate_z(-moon_ra);
    Ok(m)
}

impl CoordinateFrame for SelenographicPosition {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let m = selenographic_to_icrs_matrix(epoch)?;
        let (ra, dec) = m
            .transpose()
            .transform_spherical(self.longitude.radians(), self.latitude.radians());

        let mut icrs = ICRSPosition::new(
            Angle::from_radians(normalize_angle_to_positive(ra)),
            Angle::from_radians(dec),
        )?;

        if let Some(radius) = self.radius {
            icrs.set_distance(radius);
        }
        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let m = selenographic_to_icrs_matrix(epoch)?;
        let (lon, lat) = m.transform_spherical(icrs.ra().radians(), icrs.dec().radians());

        let mut pos = Self::new(
            Angle::from_radians(lat),
            Angle::from_radians(normalize_angle_to_positive(lon)),
        )?;

        if let Some(dist) = icrs.distance() {
            pos.set_radius(dist);
        }
        Ok(pos)
    }
}

impl std::fmt::Display for SelenographicPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Selenographic(lat={:.6}°, lon={:.6}°",
            self.latitude.degrees(),
            self.longitude.degrees()
        )?;

        if let Some(radius) = self.radius {
            write!(f, ", r={}", radius)?;
        }

        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_selenographic_creation() {
        let pos = SelenographicPosition::from_degrees(45.0, 30.0).unwrap();
        assert!((pos.latitude().degrees() - 45.0).abs() < 1e-12);
        assert!((pos.longitude().degrees() - 30.0).abs() < 1e-12);
        assert!(pos.radius().is_none());
    }

    #[test]
    fn test_selenographic_validation() {
        assert!(SelenographicPosition::from_degrees(0.0, 0.0).is_ok());
        assert!(SelenographicPosition::from_degrees(90.0, 180.0).is_ok());
        assert!(SelenographicPosition::from_degrees(-90.0, 359.0).is_ok());

        assert!(SelenographicPosition::from_degrees(95.0, 0.0).is_err());
        assert!(SelenographicPosition::from_degrees(-95.0, 0.0).is_err());
    }

    #[test]
    fn test_special_positions() {
        let nearside = SelenographicPosition::nearside_center();
        assert_eq!(nearside.latitude().degrees(), 0.0);
        assert_eq!(nearside.longitude().degrees(), 0.0);

        let farside = SelenographicPosition::farside_center();
        assert_eq!(farside.latitude().degrees(), 0.0);
        assert_eq!(farside.longitude().degrees(), 180.0);

        let north_pole = SelenographicPosition::north_pole();
        assert_eq!(north_pole.latitude().degrees(), 90.0);

        let south_pole = SelenographicPosition::south_pole();
        assert_eq!(south_pole.latitude().degrees(), -90.0);
    }

    #[test]
    fn test_angular_separation() {
        let nearside = SelenographicPosition::nearside_center();
        let farside = SelenographicPosition::farside_center();

        let sep = nearside.angular_separation(&farside);
        assert!((sep.degrees() - 180.0).abs() < 1e-10);

        let north = SelenographicPosition::north_pole();
        let sep_to_north = nearside.angular_separation(&north);
        assert!((sep_to_north.degrees() - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_visibility_from_earth_farside() {
        let epoch = TT::j2000();

        let farside = SelenographicPosition::farside_center();
        assert!(!farside.is_visible_from_earth(&epoch));
    }

    #[test]
    fn test_sub_earth_point() {
        let epoch = TT::j2000();
        let sub_earth = SelenographicPosition::sub_earth_point(&epoch).unwrap();

        assert!(
            sub_earth.latitude().degrees().abs() <= 7.5,
            "Sub-earth latitude = {}",
            sub_earth.latitude().degrees()
        );
        assert!(
            sub_earth.longitude().degrees() >= 0.0 && sub_earth.longitude().degrees() < 360.0,
            "Sub-earth longitude = {}",
            sub_earth.longitude().degrees()
        );
    }

    #[test]
    fn test_coordinate_frame_to_icrs() {
        let epoch = TT::j2000();
        let original = SelenographicPosition::from_degrees(0.0, 0.0).unwrap();

        let icrs = original.to_icrs(&epoch).unwrap();

        assert!(icrs.ra().degrees() >= 0.0 && icrs.ra().degrees() < 360.0);
        assert!(icrs.dec().degrees() >= -90.0 && icrs.dec().degrees() <= 90.0);
    }

    #[test]
    fn test_coordinate_frame_roundtrip() {
        let epoch = TT::j2000();
        let test_cases = [
            (0.0, 0.0),
            (5.0, 30.0),
            (-5.0, 90.0),
            (3.0, 180.0),
            (-3.0, 270.0),
        ];

        for (lat, lon) in test_cases {
            let original = SelenographicPosition::from_degrees(lat, lon).unwrap();
            let icrs = original.to_icrs(&epoch).unwrap();
            let recovered = SelenographicPosition::from_icrs(&icrs, &epoch).unwrap();

            let lat_err = (original.latitude().degrees() - recovered.latitude().degrees()).abs();
            let lon_diff = (original.longitude().radians() - recovered.longitude().radians()).abs();
            let lon_err = if lon_diff > std::f64::consts::PI {
                std::f64::consts::TAU - lon_diff
            } else {
                lon_diff
            } * celestial_core::constants::RAD_TO_DEG;

            assert!(
                lat_err < 1.0 / 3600.0,
                "({}, {}): Latitude error {:.6} arcsec",
                lat,
                lon,
                lat_err * 3600.0,
            );
            assert!(
                lon_err < 1.0 / 3600.0,
                "({}, {}): Longitude error {:.6} arcsec",
                lat,
                lon,
                lon_err * 3600.0,
            );
        }
    }

    #[test]
    fn test_with_radius() {
        let radius = Distance::from_au(0.00257).unwrap();
        let pos = SelenographicPosition::with_radius(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            radius,
        )
        .unwrap();

        assert!(pos.radius().is_some());
        assert_eq!(pos.radius().unwrap(), radius);
    }

    #[test]
    fn test_display_formatting() {
        let pos = SelenographicPosition::from_degrees(45.123456, 30.654321).unwrap();
        let display = format!("{}", pos);
        assert!(display.contains("45.123456"));
        assert!(display.contains("30.654321"));
        assert!(display.contains("Selenographic"));
    }
}
