mod air_mass;
mod hour_angle;
mod parallax;
pub mod refraction;

use crate::distance::Distance;
use crate::errors::CoordResult;
use celestial_core::constants::{HALF_PI, PI};
use celestial_core::{angle::Angle, location::Location};
use celestial_time::scales::tt::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TopocentricPosition {
    azimuth: Angle,
    elevation: Angle,
    observer: Location,
    epoch: TT,
    distance: Option<Distance>,
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HourAnglePosition {
    hour_angle: Angle,
    declination: Angle,
    observer: Location,
    epoch: TT,
    distance: Option<Distance>,
}

impl TopocentricPosition {
    pub fn new(
        azimuth: Angle,
        elevation: Angle,
        observer: Location,
        epoch: TT,
    ) -> CoordResult<Self> {
        let azimuth = azimuth.normalized()?;
        let elevation = elevation.validate_latitude()?;

        Ok(Self {
            azimuth,
            elevation,
            observer,
            epoch,
            distance: None,
        })
    }

    pub fn with_distance(
        azimuth: Angle,
        elevation: Angle,
        observer: Location,
        epoch: TT,
        distance: Distance,
    ) -> CoordResult<Self> {
        let mut pos = Self::new(azimuth, elevation, observer, epoch)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn from_degrees(
        az_deg: f64,
        el_deg: f64,
        observer: Location,
        epoch: TT,
    ) -> CoordResult<Self> {
        Self::new(
            Angle::from_degrees(az_deg),
            Angle::from_degrees(el_deg),
            observer,
            epoch,
        )
    }

    pub fn azimuth(&self) -> Angle {
        self.azimuth
    }

    pub fn elevation(&self) -> Angle {
        self.elevation
    }

    pub fn observer(&self) -> &Location {
        &self.observer
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

    pub fn zenith_angle(&self) -> Angle {
        Angle::HALF_PI - self.elevation
    }

    pub fn is_above_horizon(&self) -> bool {
        self.elevation.degrees() > 0.0
    }

    pub fn is_near_zenith(&self) -> bool {
        self.elevation.degrees() > 89.0
    }

    pub fn is_near_horizon(&self) -> bool {
        self.elevation.degrees() < 10.0 && self.is_above_horizon()
    }

    pub fn cardinal_direction(&self) -> &'static str {
        let az_deg = self.azimuth.degrees();
        if !(22.5..337.5).contains(&az_deg) {
            "N"
        } else if az_deg < 67.5 {
            "NE"
        } else if az_deg < 112.5 {
            "E"
        } else if az_deg < 157.5 {
            "SE"
        } else if az_deg < 202.5 {
            "S"
        } else if az_deg < 247.5 {
            "SW"
        } else if az_deg < 292.5 {
            "W"
        } else {
            "NW"
        }
    }
}

impl HourAnglePosition {
    pub fn new(
        hour_angle: Angle,
        declination: Angle,
        observer: Location,
        epoch: TT,
    ) -> CoordResult<Self> {
        let declination = declination.validate_declination(true)?; // beyond_pole for GEM pier-flips
        let (hour_angle, declination) = fold_over_pole(hour_angle, declination);
        let hour_angle = hour_angle.wrapped()?; // [-180°, +180°]

        Ok(Self {
            hour_angle,
            declination,
            observer,
            epoch,
            distance: None,
        })
    }

    pub fn with_distance(
        hour_angle: Angle,
        declination: Angle,
        observer: Location,
        epoch: TT,
        distance: Distance,
    ) -> CoordResult<Self> {
        let mut pos = Self::new(hour_angle, declination, observer, epoch)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn hour_angle(&self) -> Angle {
        self.hour_angle
    }

    pub fn declination(&self) -> Angle {
        self.declination
    }

    pub fn observer(&self) -> &Location {
        &self.observer
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
}

// Note: Topocentric coordinates cannot implement CoordinateFrame directly
// because they require both time AND observer location for transformation.
// They need specialized transformation methods.

impl std::fmt::Display for TopocentricPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Topocentric(Az={:.2}° {}, El={:.2}°",
            self.azimuth.degrees(),
            self.cardinal_direction(),
            self.elevation.degrees()
        )?;

        if let Some(distance) = self.distance {
            write!(f, ", d={}", distance)?;
        }

        write!(f, ")")
    }
}

impl std::fmt::Display for HourAnglePosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "HourAngle(HA={:.4}h, Dec={:.4}°",
            self.hour_angle.hours(),
            self.declination.degrees()
        )?;

        if let Some(distance) = self.distance {
            write!(f, ", d={}", distance)?;
        }

        write!(f, ")")
    }
}

// A declination past a pole names the same direction as the declination
// reflected back over it, twelve hours of hour angle away.
fn fold_over_pole(hour_angle: Angle, declination: Angle) -> (Angle, Angle) {
    let dec = declination.radians();
    if dec > HALF_PI {
        (hour_angle + Angle::PI, Angle::from_radians(PI - dec))
    } else if dec < -HALF_PI {
        (hour_angle + Angle::PI, Angle::from_radians(-PI - dec))
    } else {
        (hour_angle, declination)
    }
}

#[cfg(test)]
fn test_observer() -> Location {
    // Keck Observatory, Mauna Kea (4145m per keckobservatory.org)
    Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_topocentric_creation() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        assert_eq!(topo.azimuth(), Angle::from_degrees(180.0));
        assert_eq!(topo.elevation(), Angle::from_degrees(45.0));
        assert_eq!(
            topo.observer().latitude_degrees(),
            observer.latitude_degrees()
        );
        assert_eq!(topo.epoch(), epoch);
    }

    #[test]
    fn test_topocentric_validation() {
        let observer = test_observer();
        let epoch = TT::j2000();

        assert!(TopocentricPosition::from_degrees(0.0, 0.0, observer, epoch).is_ok());
        assert!(TopocentricPosition::from_degrees(359.999, 89.999, observer, epoch).is_ok());

        // eraAnp of 380 degrees: 20 degrees, less the rounding of 380 degrees in radians.
        let topo = TopocentricPosition::from_degrees(380.0, 45.0, observer, epoch).unwrap();
        assert_eq!(topo.azimuth().radians(), 0.3490658503988664);

        assert!(TopocentricPosition::from_degrees(0.0, 95.0, observer, epoch).is_err());
    }

    #[test]
    fn test_position_classification() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let above = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        assert!(above.is_above_horizon());
        assert!(!above.is_near_zenith());
        assert!(!above.is_near_horizon());

        let below = TopocentricPosition::from_degrees(0.0, -5.0, observer, epoch).unwrap();
        assert!(!below.is_above_horizon());

        let zenith = TopocentricPosition::from_degrees(0.0, 89.5, observer, epoch).unwrap();
        assert!(zenith.is_near_zenith());

        let horizon = TopocentricPosition::from_degrees(0.0, 5.0, observer, epoch).unwrap();
        assert!(horizon.is_near_horizon());
    }

    #[test]
    fn test_cardinal_directions() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let north = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        assert_eq!(north.cardinal_direction(), "N");

        let east = TopocentricPosition::from_degrees(90.0, 45.0, observer, epoch).unwrap();
        assert_eq!(east.cardinal_direction(), "E");

        let south = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        assert_eq!(south.cardinal_direction(), "S");

        let west = TopocentricPosition::from_degrees(270.0, 45.0, observer, epoch).unwrap();
        assert_eq!(west.cardinal_direction(), "W");

        let northeast = TopocentricPosition::from_degrees(45.0, 45.0, observer, epoch).unwrap();
        assert_eq!(northeast.cardinal_direction(), "NE");
    }

    #[test]
    fn test_hour_angle_creation() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let ha_pos = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        assert_eq!(ha_pos.hour_angle(), Angle::from_hours(2.0));
        assert_eq!(ha_pos.declination(), Angle::from_degrees(45.0));
    }

    #[test]
    fn test_with_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let distance = Distance::from_kilometers(384400.0).unwrap(); // Moon distance

        let topo = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        assert_eq!(topo.distance().unwrap().kilometers(), distance.kilometers());

        let ha_pos = HourAnglePosition::with_distance(
            Angle::from_hours(1.0),
            Angle::from_degrees(30.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        assert_eq!(
            ha_pos.distance().unwrap().kilometers(),
            distance.kilometers()
        );
    }

    #[test]
    fn test_topocentric_set_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let mut topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();

        assert!(topo.distance().is_none());

        let distance = Distance::from_kilometers(1000.0).unwrap();
        topo.set_distance(distance);

        assert_eq!(topo.distance(), Some(distance));
    }

    #[test]
    fn test_cardinal_directions_all() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let directions = [
            (0.0, "N"),
            (45.0, "NE"),
            (90.0, "E"),
            (135.0, "SE"),
            (180.0, "S"),
            (225.0, "SW"),
            (270.0, "W"),
            (315.0, "NW"),
        ];

        for (az, expected) in directions {
            let topo = TopocentricPosition::from_degrees(az, 45.0, observer, epoch).unwrap();
            assert_eq!(
                topo.cardinal_direction(),
                expected,
                "Failed for azimuth {}°",
                az
            );
        }
    }

    #[test]
    fn test_hour_angle_getters() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let ha_pos = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        assert_eq!(
            ha_pos.observer().latitude_degrees(),
            observer.latitude_degrees()
        );
        assert_eq!(ha_pos.epoch(), epoch);
        assert!(ha_pos.distance().is_none());
    }

    #[test]
    fn test_hour_angle_set_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let mut ha_pos = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        let distance = Distance::from_kilometers(500.0).unwrap();
        ha_pos.set_distance(distance);
        assert_eq!(ha_pos.distance(), Some(distance));
    }

    #[test]
    fn test_display_formatting() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        let display = format!("{}", topo);
        assert!(display.contains("Topocentric"));
        assert!(display.contains("180.00°"));
        assert!(display.contains("45.00°"));
        assert!(display.contains("S"));

        let distance = Distance::from_kilometers(1000.0).unwrap();
        let topo_dist = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();
        let display_dist = format!("{}", topo_dist);
        assert!(display_dist.contains("AU") || display_dist.contains("pc"));

        let ha = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();
        let ha_display = format!("{}", ha);
        assert!(ha_display.contains("HourAngle"));
        assert!(ha_display.contains("2."));
        assert!(ha_display.contains("45."));

        let ha_dist = HourAnglePosition::with_distance(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();
        let ha_display_dist = format!("{}", ha_dist);
        assert!(ha_display_dist.contains("AU") || ha_display_dist.contains("pc"));
    }
}
