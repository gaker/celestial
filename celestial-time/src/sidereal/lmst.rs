use super::angle::SiderealAngle;
use super::gmst::GMST;
use crate::scales::{TT, UT1};
use crate::TimeResult;
use celestial_core::Location;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LMST {
    angle: SiderealAngle,
    location: Location,
}

impl LMST {
    pub fn from_ut1_tt_and_location(ut1: &UT1, tt: &TT, location: &Location) -> TimeResult<Self> {
        let gmst = GMST::from_ut1_and_tt(ut1, tt)?;

        let lmst_rad = gmst.radians() + location.longitude;

        use celestial_core::angle::wrap_0_2pi;
        let lmst_normalized = wrap_0_2pi(lmst_rad);

        let angle = SiderealAngle::from_radians_exact(lmst_normalized);

        Ok(Self {
            angle,
            location: *location,
        })
    }

    pub fn from_hours(hours: f64, location: &Location) -> Self {
        Self {
            angle: SiderealAngle::from_hours(hours),
            location: *location,
        }
    }

    pub fn from_degrees(degrees: f64, location: &Location) -> Self {
        Self {
            angle: SiderealAngle::from_degrees(degrees),
            location: *location,
        }
    }

    pub fn from_radians(radians: f64, location: &Location) -> Self {
        Self {
            angle: SiderealAngle::from_radians(radians),
            location: *location,
        }
    }

    pub fn j2000(location: &Location) -> TimeResult<Self> {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        Self::from_ut1_tt_and_location(&ut1, &tt, location)
    }

    pub fn angle(&self) -> SiderealAngle {
        self.angle
    }

    pub fn location(&self) -> Location {
        self.location
    }

    pub fn hours(&self) -> f64 {
        self.angle.hours()
    }

    pub fn degrees(&self) -> f64 {
        self.angle.degrees()
    }

    pub fn radians(&self) -> f64 {
        self.angle.radians()
    }

    pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> f64 {
        self.angle.hour_angle_to_target(target_ra_hours)
    }

    pub fn to_gmst(&self) -> GMST {
        let longitude_hours = self.location.longitude * 12.0 / celestial_core::constants::PI;

        let gmst_hours = self.hours() - longitude_hours;

        GMST::from_hours(gmst_hours)
    }
}

impl std::fmt::Display for LMST {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let lat_deg = self.location.latitude * celestial_core::constants::RAD_TO_DEG;
        let lon_deg = self.location.longitude * celestial_core::constants::RAD_TO_DEG;
        write!(
            f,
            "LMST {} at ({:.4}°, {:.4}°)",
            self.angle, lat_deg, lon_deg
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mauna_kea() -> Location {
        Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap()
    }

    fn greenwich() -> Location {
        Location::greenwich()
    }

    #[test]
    fn test_lmst_at_greenwich_equals_gmst() {
        // At Greenwich (0° longitude), LMST should equal GMST
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = greenwich();

        let gmst = GMST::from_ut1_and_tt(&ut1, &tt).unwrap();
        let lmst = LMST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();

        // Should be identical (within numerical precision)
        assert!(
            (lmst.hours() - gmst.hours()).abs() < 1e-14,
            "LMST at Greenwich should equal GMST: LMST={}, GMST={}",
            lmst.hours(),
            gmst.hours()
        );
    }

    #[test]
    fn test_lmst_longitude_correction() {
        // Test longitude correction: 1 degree = 4 minutes = 1/15 hour
        let ut1 = UT1::j2000();
        let tt = TT::j2000();

        let greenwich_loc = greenwich();
        let east_15deg = Location::from_degrees(0.0, 15.0, 0.0).unwrap(); // 15°E = +1 hour
        let west_15deg = Location::from_degrees(0.0, -15.0, 0.0).unwrap(); // 15°W = -1 hour

        let lmst_greenwich = LMST::from_ut1_tt_and_location(&ut1, &tt, &greenwich_loc).unwrap();
        let lmst_east = LMST::from_ut1_tt_and_location(&ut1, &tt, &east_15deg).unwrap();
        let lmst_west = LMST::from_ut1_tt_and_location(&ut1, &tt, &west_15deg).unwrap();

        // 15°E should be +1 hour ahead of Greenwich
        let diff_east = lmst_east.hours() - lmst_greenwich.hours();
        assert!(
            (diff_east - 1.0).abs() < 1e-12,
            "15°E should be +1 hour: {}",
            diff_east
        );

        // 15°W should be -1 hour behind Greenwich
        let diff_west = lmst_west.hours() - lmst_greenwich.hours();
        assert!(
            (diff_west + 1.0).abs() < 1e-12,
            "15°W should be -1 hour: {}",
            diff_west
        );
    }

    #[test]
    fn test_lmst_mauna_kea() {
        // Mauna Kea is at -155.4783° = -10.365 hours west of Greenwich
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();

        let gmst = GMST::from_ut1_and_tt(&ut1, &tt).unwrap();
        let lmst = LMST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();

        // Expected longitude correction in hours
        let expected_offset = -155.4783 / 15.0; // degrees to hours
        let actual_offset = lmst.hours() - gmst.hours();

        assert!(
            (actual_offset - expected_offset).abs() < 1e-10,
            "Mauna Kea LMST offset incorrect: expected={}, actual={}",
            expected_offset,
            actual_offset
        );
    }

    #[test]
    fn test_lmst_j2000() {
        let location = mauna_kea();
        let lmst = LMST::j2000(&location).unwrap();

        // LMST should be in valid range
        let hours = lmst.hours();
        assert!(
            (0.0..24.0).contains(&hours),
            "LMST should be in [0, 24) hours: {}",
            hours
        );
    }

    #[test]
    fn test_lmst_hour_angle_calculation() {
        let location = mauna_kea();
        let lmst = LMST::from_hours(12.0, &location);
        let target_ra = 6.0;
        let hour_angle = lmst.hour_angle_to_target(target_ra);
        assert_eq!(hour_angle, 6.0);
    }

    #[test]
    fn test_lmst_to_gmst_roundtrip() {
        let location = mauna_kea();
        let original_gmst = GMST::from_hours(15.5);

        // Convert GMST -> LMST -> GMST
        let longitude_hours = location.longitude * 12.0 / celestial_core::constants::PI;
        let lmst_hours = original_gmst.hours() + longitude_hours;
        let lmst = LMST::from_hours(lmst_hours, &location);
        let recovered_gmst = lmst.to_gmst();

        assert!(
            (recovered_gmst.hours() - original_gmst.hours()).abs() < 1e-14,
            "GMST->LMST->GMST roundtrip failed: original={}, recovered={}",
            original_gmst.hours(),
            recovered_gmst.hours()
        );
    }

    #[test]
    fn test_lmst_from_constructors() {
        let location = mauna_kea();

        // Test all constructor methods produce equivalent results
        let hours = 14.5;
        let degrees = hours * 15.0;
        let radians = hours * celestial_core::constants::PI / 12.0;

        let lmst_hours = LMST::from_hours(hours, &location);
        let lmst_degrees = LMST::from_degrees(degrees, &location);
        let lmst_radians = LMST::from_radians(radians, &location);

        assert!((lmst_hours.hours() - lmst_degrees.hours()).abs() < 1e-14);
        assert!((lmst_hours.hours() - lmst_radians.hours()).abs() < 1e-14);
        assert_eq!(lmst_hours.location(), location);
        assert_eq!(lmst_degrees.location(), location);
        assert_eq!(lmst_radians.location(), location);
    }

    #[test]
    fn test_lmst_display() {
        let location = mauna_kea();
        let lmst = LMST::from_hours(12.5, &location);
        let display = format!("{}", lmst);

        // Should include LMST time and location coordinates
        assert!(display.contains("LMST"));
        assert!(display.contains("19.8283")); // Latitude
        assert!(display.contains("-155.4783")); // Longitude
    }

    #[test]
    fn test_lmst_location_enforcement() {
        // Test that LMST enforces location through type system
        let location = mauna_kea();
        let lmst = LMST::from_hours(12.0, &location);

        // Location is always available and cannot be None/invalid
        let stored_location = lmst.location();
        assert_eq!(stored_location.latitude, location.latitude);
        assert_eq!(stored_location.longitude, location.longitude);
        assert_eq!(stored_location.height, location.height);
    }

    #[test]
    fn test_extreme_longitudes() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();

        // Test extreme valid longitudes
        let east_extreme = Location::from_degrees(0.0, 180.0, 0.0).unwrap(); // 180°E = +12 hours
        let west_extreme = Location::from_degrees(0.0, -180.0, 0.0).unwrap(); // 180°W = -12 hours

        let gmst = GMST::from_ut1_and_tt(&ut1, &tt).unwrap();
        let lmst_east = LMST::from_ut1_tt_and_location(&ut1, &tt, &east_extreme).unwrap();
        let lmst_west = LMST::from_ut1_tt_and_location(&ut1, &tt, &west_extreme).unwrap();

        // 180°E should be +12 hours ahead
        let diff_east = lmst_east.hours() - gmst.hours();
        let expected_east = if diff_east < 0.0 {
            diff_east + 24.0
        } else {
            diff_east
        };
        assert!(
            (expected_east - 12.0).abs() < 1e-12,
            "180°E should be +12 hours"
        );

        // 180°W should be -12 hours behind (equivalent to +12 hours due to 24h wrap)
        let diff_west = lmst_west.hours() - gmst.hours();
        let expected_west = if diff_west > 12.0 {
            diff_west - 24.0
        } else {
            diff_west
        };
        assert!(
            (expected_west + 12.0).abs() < 1e-12,
            "180°W should be -12 hours"
        );
    }

    #[test]
    fn test_lmst_constructors_and_accessors() {
        let location = mauna_kea();

        // Test from_degrees constructor
        let lmst_deg = LMST::from_degrees(90.0, &location);
        assert_eq!(lmst_deg.degrees(), 90.0);
        assert_eq!(lmst_deg.hours(), 6.0);

        // Test from_radians constructor
        let lmst_rad = LMST::from_radians(celestial_core::constants::PI * 1.5, &location);
        assert!((lmst_rad.radians() - celestial_core::constants::PI * 1.5).abs() < 1e-15);
        assert_eq!(lmst_rad.hours(), 18.0);

        // Test angle() accessor
        let angle = lmst_deg.angle();
        assert_eq!(angle.degrees(), 90.0);

        // Test location() accessor
        let stored_location = lmst_deg.location();
        assert_eq!(stored_location.latitude, location.latitude);
        assert_eq!(stored_location.longitude, location.longitude);
        assert_eq!(stored_location.height, location.height);

        // Test degrees() method
        let degrees = lmst_deg.degrees();
        assert_eq!(degrees, 90.0);
    }
}
