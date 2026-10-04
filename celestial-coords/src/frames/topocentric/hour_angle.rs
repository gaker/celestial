use super::{HourAnglePosition, TopocentricPosition};
use crate::astrom::earth::EarthRotation;
use crate::eop::record::EopParameters;
use crate::errors::CoordResult;
use crate::frames::cirs::CIRSPosition;
use celestial_core::angle::Angle;
use celestial_core::constants::{HALF_PI, TWOPI};

impl TopocentricPosition {
    pub fn to_hour_angle(&self) -> CoordResult<HourAnglePosition> {
        let (sa, ca) = self.azimuth.sin_cos();
        let (se, ce) = self.elevation.sin_cos();
        let (sp, cp) = self.observer.latitude_angle().sin_cos();
        let x = -ca * ce * sp + se * cp;
        let y = -sa * ce;
        let z = ca * ce * cp + se * sp;
        let r = libm::sqrt(x * x + y * y);
        let hour_angle = if r != 0.0 { libm::atan2(y, x) } else { 0.0 };
        // atan2 already gives [-pi, +pi]; HourAnglePosition::new would turn
        // +pi into -pi.
        Ok(HourAnglePosition {
            hour_angle: Angle::from_radians(hour_angle),
            declination: Angle::from_radians(libm::atan2(z, r)),
            observer: self.observer,
            epoch: self.epoch,
            distance: self.distance,
        })
    }

    pub fn to_cirs(&self, eop: &EopParameters) -> CoordResult<CIRSPosition> {
        self.to_hour_angle()?.to_cirs(eop)
    }
}

impl HourAnglePosition {
    pub fn to_topocentric(&self) -> CoordResult<TopocentricPosition> {
        let (azimuth, elevation) = self.azimuth_elevation();
        let mut topo = TopocentricPosition::new(
            Angle::from_radians(azimuth),
            Angle::from_radians(elevation),
            self.observer,
            self.epoch,
        )?;
        if let Some(distance) = self.distance {
            topo.set_distance(distance);
        }
        Ok(topo)
    }

    // x points north along the horizon, y east and z to the zenith, so azimuth runs from north
    // through east.
    fn azimuth_elevation(&self) -> (f64, f64) {
        let (sin_ha, cos_ha) = self.hour_angle.sin_cos();
        let (sin_dec, cos_dec) = self.declination.sin_cos();
        let (sin_lat, cos_lat) = self.observer.latitude_angle().sin_cos();
        let x = -cos_ha * cos_dec * sin_lat + sin_dec * cos_lat;
        let y = -sin_ha * cos_dec;
        let z = cos_ha * cos_dec * cos_lat + sin_dec * sin_lat;
        let r = libm::sqrt(x * x + y * y);
        let azimuth = if r != 0.0 { libm::atan2(y, x) } else { 0.0 };
        let azimuth = if azimuth < 0.0 {
            azimuth + TWOPI
        } else {
            azimuth
        };
        (azimuth, libm::atan2(z, r))
    }

    pub fn parallactic_angle(&self) -> Angle {
        let (sp, cp) = self.observer.latitude_angle().sin_cos();
        let (sh, ch) = self.hour_angle.sin_cos();
        let (sd, cd) = self.declination.sin_cos();
        let sqsz = cp * sh;
        let cqsz = sp * cd - cp * sd * ch;
        if sqsz == 0.0 && cqsz == 0.0 {
            return Angle::ZERO;
        }
        Angle::from_radians(libm::atan2(sqsz, cqsz))
    }

    pub fn is_circumpolar(&self) -> bool {
        let (colatitude, hemisphere) = self.colatitude();
        self.declination.radians() * hemisphere > colatitude
    }

    pub fn never_rises(&self) -> bool {
        let (colatitude, hemisphere) = self.colatitude();
        -self.declination.radians() * hemisphere > colatitude
    }

    // The angle from the observer's visible pole to the horizon, and +1 or
    // -1 for the hemisphere that pole is in.
    fn colatitude(&self) -> (f64, f64) {
        let latitude = self.observer.latitude();
        (
            HALF_PI - libm::fabs(latitude),
            libm::copysign(1.0, latitude),
        )
    }

    pub fn to_cirs(&self, eop: &EopParameters) -> CoordResult<CIRSPosition> {
        EarthRotation::new(&self.epoch, eop)?.hour_angle_to_cirs(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::Distance;
    use crate::eop::record::EopRecord;
    use crate::frames::topocentric::test_observer;
    use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT, PI};
    use celestial_core::location::Location;
    use celestial_time::scales::tt::TT;

    #[test]
    fn test_hour_angle_to_topocentric() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let ha_pos =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(45.0), observer, epoch)
                .unwrap();

        // eraHd2ae: due north, since the declination is above the latitude.
        let topo = ha_pos.to_topocentric().unwrap();
        assert_eq!(
            [topo.azimuth().radians(), topo.elevation().radians()],
            [0.0, 1.13146728347064]
        );
    }

    #[test]
    fn test_circumpolar() {
        let observer = test_observer(); // Latitude ~20°N
        let epoch = TT::j2000();

        let high_dec =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(85.0), observer, epoch)
                .unwrap();
        assert!(high_dec.is_circumpolar());

        let low_dec =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(0.0), observer, epoch).unwrap();
        assert!(!low_dec.is_circumpolar());

        let neg_dec =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(-85.0), observer, epoch)
                .unwrap();
        assert!(neg_dec.never_rises());
    }

    #[test]
    fn test_parallactic_angle() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // On the meridian south of the zenith, the pole and the zenith lie
        // in the same direction.
        let meridian =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(0.0), observer, epoch).unwrap();
        assert_eq!(meridian.parallactic_angle(), Angle::ZERO);

        // At the zenith the angle is undefined and reported as zero.
        let zenith =
            HourAnglePosition::new(Angle::ZERO, observer.latitude_angle(), observer, epoch)
                .unwrap();
        assert_eq!(zenith.parallactic_angle(), Angle::ZERO);

        // West of the meridian (positive hour angle) the angle is positive.
        let west = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(10.0),
            observer,
            epoch,
        )
        .unwrap();
        assert!(west.parallactic_angle().radians() > 0.0);
    }

    #[test]
    fn test_circumpolar_southern_hemisphere() {
        let observer = Location::from_degrees(-30.0, -70.0, 0.0).unwrap();
        let epoch = TT::j2000();
        let at = |dec: f64| {
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(dec), observer, epoch).unwrap()
        };

        assert!(at(-70.0).is_circumpolar());
        assert!(!at(-70.0).never_rises());
        assert!(at(70.0).never_rises());
        assert!(!at(70.0).is_circumpolar());
        assert!(!at(-50.0).is_circumpolar());
        assert!(!at(50.0).never_rises());
    }

    #[test]
    fn test_beyond_pole_declination_folds_back() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let ha = Angle::from_hours(1.0);
        let flipped = (ha + Angle::PI).wrapped().unwrap();

        let north =
            HourAnglePosition::new(ha, Angle::from_degrees(120.0), observer, epoch).unwrap();
        assert_eq!(
            north.declination().radians(),
            PI - Angle::from_degrees(120.0).radians()
        );
        assert_eq!(north.hour_angle(), flipped);
        assert!(north.to_cirs(&eop_at_j2000()).is_ok());

        let south =
            HourAnglePosition::new(ha, Angle::from_degrees(-100.0), observer, epoch).unwrap();
        assert_eq!(
            south.declination().radians(),
            -PI - Angle::from_degrees(-100.0).radians()
        );
        assert_eq!(south.hour_angle(), flipped);
        assert!(south.to_cirs(&eop_at_j2000()).is_ok());
    }

    #[test]
    fn test_hour_angle_with_distance_to_topocentric() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let distance = Distance::from_kilometers(1000.0).unwrap();

        let ha_pos = HourAnglePosition::with_distance(
            Angle::ZERO,
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        let topo = ha_pos.to_topocentric().unwrap();
        assert_eq!(topo.distance(), Some(distance));
    }

    #[test]
    fn test_topocentric_to_hour_angle() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        let ha = topo.to_hour_angle().unwrap();

        // eraAe2hd: the sine of the double nearest pi puts the hour angle 1e-16 rad
        // east of the meridian.
        assert_eq!(
            [ha.hour_angle().radians(), ha.declination().radians()],
            [-9.568181415836897e-17, -0.4393290433242568]
        );
    }

    #[test]
    fn test_topocentric_hour_angle_roundtrip() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let test_cases = [
            (Angle::from_hours(0.0), Angle::from_degrees(45.0)),
            (Angle::from_hours(2.0), Angle::from_degrees(30.0)),
            (Angle::from_hours(-3.0), Angle::from_degrees(60.0)),
            (Angle::from_hours(6.0), Angle::from_degrees(0.0)),
        ];

        for (ha, dec) in test_cases {
            let original = HourAnglePosition::new(ha, dec, observer, epoch).unwrap();

            let topo = original.to_topocentric().unwrap();
            let recovered = topo.to_hour_angle().unwrap();

            // A rotation and its inverse, so a few units in the last place of a radian.
            let ha_diff = recovered.hour_angle().radians() - original.hour_angle().radians();
            let dec_diff = recovered.declination().radians() - original.declination().radians();
            let ctx = format!("HA={}h, Dec={}°", ha.hours(), dec.degrees());
            assert!(libm::fabs(ha_diff) < 1e-15, "{ctx}: {ha_diff:e} rad");
            assert!(libm::fabs(dec_diff) < 1e-15, "{ctx}: {dec_diff:e} rad");
        }
    }

    #[test]
    fn test_topocentric_to_hour_angle_distance_preservation() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let distance = Distance::from_kilometers(384400.0).unwrap();

        let topo = TopocentricPosition::with_distance(
            Angle::from_degrees(90.0),
            Angle::from_degrees(30.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        let ha = topo.to_hour_angle().unwrap();
        assert_eq!(ha.distance().unwrap().kilometers(), distance.kilometers());
    }

    #[test]
    fn test_topocentric_to_hour_angle_cardinal_points() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Due east (Az=90°): object is rising, HA should be negative (before meridian)
        let east = TopocentricPosition::from_degrees(90.0, 30.0, observer, epoch).unwrap();
        let ha_east = east.to_hour_angle().unwrap();
        assert!(
            ha_east.hour_angle().hours() < 0.0 || ha_east.hour_angle().hours() > 12.0,
            "East object should have negative or >12h hour angle, got {}h",
            ha_east.hour_angle().hours()
        );

        // Due west (Az=270°): object is setting, HA should be positive
        let west = TopocentricPosition::from_degrees(270.0, 30.0, observer, epoch).unwrap();
        let ha_west = west.to_hour_angle().unwrap();
        assert!(
            ha_west.hour_angle().hours() > 0.0 && ha_west.hour_angle().hours() < 12.0,
            "West object should have positive hour angle, got {}h",
            ha_west.hour_angle().hours()
        );
    }

    fn eop_at_j2000() -> EopParameters {
        EopRecord::new(J2000_JD - MJD_ZERO_POINT, 0.0, 0.0, 0.3)
            .unwrap()
            .to_parameters()
    }

    #[test]
    fn test_hour_angle_to_cirs() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let eop = eop_at_j2000();

        let ha = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        let cirs = ha.to_cirs(&eop).unwrap();

        assert!(cirs.ra().degrees() >= 0.0 && cirs.ra().degrees() < 360.0);
        assert_eq!(cirs.dec().degrees(), ha.declination().degrees());
    }

    #[test]
    fn test_hour_angle_cirs_roundtrip() {
        use crate::frames::cirs::CIRSPosition;

        let observer = test_observer();
        let epoch = TT::j2000();
        let eop = eop_at_j2000();

        let original_cirs = CIRSPosition::from_degrees(120.0, 35.0, epoch).unwrap();

        let ha = original_cirs.to_hour_angle(&observer, &eop).unwrap();
        let recovered_cirs = ha.to_cirs(&eop).unwrap();

        assert_eq!(recovered_cirs.ra(), original_cirs.ra());
        assert_eq!(recovered_cirs.dec(), original_cirs.dec());
    }

    #[test]
    fn test_topocentric_to_cirs() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let eop = eop_at_j2000();

        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        let cirs = topo.to_cirs(&eop).unwrap();

        assert!(cirs.ra().degrees() >= 0.0 && cirs.ra().degrees() < 360.0);
        assert!(cirs.dec().degrees() >= -90.0 && cirs.dec().degrees() <= 90.0);
    }

    #[test]
    fn test_full_reverse_chain_roundtrip() {
        use crate::frames::cirs::CIRSPosition;

        let observer = test_observer();
        let epoch = TT::j2000();
        let eop = eop_at_j2000();

        let original_cirs = CIRSPosition::from_degrees(200.0, 40.0, epoch).unwrap();

        let ha = original_cirs.to_hour_angle(&observer, &eop).unwrap();
        let topo = ha.to_topocentric().unwrap();
        let recovered_ha = topo.to_hour_angle().unwrap();
        let recovered_cirs = recovered_ha.to_cirs(&eop).unwrap();

        assert_eq!(recovered_cirs.ra(), original_cirs.ra());
        assert_eq!(recovered_cirs.dec(), original_cirs.dec());
    }
}
