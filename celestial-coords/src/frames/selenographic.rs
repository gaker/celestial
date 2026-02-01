use crate::{lunar, transforms::CoordinateFrame, CoordResult, Distance, ICRSPosition};
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

impl CoordinateFrame for SelenographicPosition {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let orientation = lunar::compute_lunar_orientation(epoch);
        let lib_lon = orientation.optical_libration.longitude.radians();
        let lib_lat = orientation.optical_libration.latitude.radians();
        let p = orientation.position_angle.radians();

        let lat = self.latitude.radians();
        let lon = self.longitude.radians();
        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_lon, cos_lon) = lon.sin_cos();

        let x_selen = cos_lat * cos_lon;
        let y_selen = cos_lat * sin_lon;
        let z_selen = sin_lat;

        let (sin_lib_lat, cos_lib_lat) = lib_lat.sin_cos();
        let x_lib = x_selen * cos_lib_lat + z_selen * sin_lib_lat;
        let y_lib = y_selen;
        let z_lib = -x_selen * sin_lib_lat + z_selen * cos_lib_lat;

        let (sin_lib_lon, cos_lib_lon) = lib_lon.sin_cos();
        let _x_earth = x_lib * cos_lib_lon - y_lib * sin_lib_lon;
        let y_earth = x_lib * sin_lib_lon + y_lib * cos_lib_lon;
        let z_earth = z_lib;

        let (sin_p, cos_p) = p.sin_cos();
        let y_proj = y_earth * cos_p + z_earth * sin_p;
        let z_proj = -y_earth * sin_p + z_earth * cos_p;

        let moon_icrs = lunar::get_moon_icrs(epoch)?;
        let moon_ra = moon_icrs.ra().radians();
        let moon_dec = moon_icrs.dec().radians();

        let scale = 0.0000113;
        let offset_x = -y_proj * scale;
        let offset_y = z_proj * scale;

        let ra_offset = offset_x / moon_dec.cos();
        let dec_offset = offset_y;

        let new_ra = moon_ra + ra_offset;
        let new_dec = moon_dec + dec_offset;

        let mut icrs = ICRSPosition::new(
            Angle::from_radians(normalize_angle_to_positive(new_ra)),
            Angle::from_radians(new_dec),
        )?;

        if let Some(radius) = self.radius {
            icrs.set_distance(radius);
        }

        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let orientation = lunar::compute_lunar_orientation(epoch);
        let lib_lon = orientation.optical_libration.longitude.radians();
        let lib_lat = orientation.optical_libration.latitude.radians();
        let p = orientation.position_angle.radians();

        let moon_icrs = lunar::get_moon_icrs(epoch)?;
        let moon_ra = moon_icrs.ra().radians();
        let moon_dec = moon_icrs.dec().radians();

        let ra_offset = icrs.ra().radians() - moon_ra;
        let dec_offset = icrs.dec().radians() - moon_dec;

        let scale = 0.0000113;
        let y_proj = -ra_offset * moon_dec.cos() / scale;
        let z_proj = dec_offset / scale;

        let (sin_p, cos_p) = p.sin_cos();
        let y_earth = y_proj * cos_p - z_proj * sin_p;
        let z_earth = y_proj * sin_p + z_proj * cos_p;

        let (sin_lib_lon, cos_lib_lon) = lib_lon.sin_cos();
        let x_earth = (1.0 - y_earth.powi(2) - z_earth.powi(2)).sqrt().max(0.0);

        let x_lib = x_earth * cos_lib_lon + y_earth * sin_lib_lon;
        let y_lib = -x_earth * sin_lib_lon + y_earth * cos_lib_lon;
        let z_lib = z_earth;

        let (sin_lib_lat, cos_lib_lat) = lib_lat.sin_cos();
        let x_selen = x_lib * cos_lib_lat - z_lib * sin_lib_lat;
        let y_selen = y_lib;
        let z_selen = x_lib * sin_lib_lat + z_lib * cos_lib_lat;

        let lat = z_selen.asin();
        let lon = y_selen.atan2(x_selen);

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
