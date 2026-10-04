use super::TopocentricPosition;
use crate::distance::Distance;
use crate::errors::{CoordError, CoordResult};
use celestial_core::angle::Angle;
use celestial_core::constants::WGS84_SEMI_MAJOR_AXIS;
use celestial_core::matrix::Vector3;

// Positions are handled as (north, east, up) vectors in metres, in the
// observer's geodetic horizon frame.
impl TopocentricPosition {
    // self is the geocentric place.
    pub fn diurnal_parallax(&self) -> CoordResult<Option<Angle>> {
        if self.distance.is_none() {
            return Ok(None);
        }
        let topocentric = self.with_diurnal_parallax()?;
        Ok(Some(self.elevation - topocentric.elevation))
    }

    pub fn horizontal_parallax(&self) -> CoordResult<Option<Angle>> {
        let Some(distance) = self.distance else {
            return Ok(None);
        };
        let metres = distance.kilometers() * 1000.0;
        if metres <= WGS84_SEMI_MAJOR_AXIS {
            return Err(CoordError::invalid_distance(
                "horizontal parallax needs an object outside the Earth's equatorial radius",
            ));
        }
        Ok(Some(Angle::from_radians(libm::asin(
            WGS84_SEMI_MAJOR_AXIS / metres,
        ))))
    }

    // self is the geocentric place; the result is the place seen by the observer.
    pub fn with_diurnal_parallax(&self) -> CoordResult<Self> {
        let Some(distance) = self.distance else {
            return Ok(self.clone());
        };
        let o = self.observer_vector();
        check_outside_observer(distance.kilometers() * 1000.0, o)?;
        self.at_vector(self.vector_metres(distance) - o)
    }

    // self is the place seen by the observer; the result is the geocentric place.
    pub fn without_diurnal_parallax(&self) -> CoordResult<Self> {
        let Some(distance) = self.distance else {
            return Ok(self.clone());
        };
        let o = self.observer_vector();
        let p = self.vector_metres(distance) + o;
        check_outside_observer(p.magnitude(), o)?;
        self.at_vector(p)
    }

    fn vector_metres(&self, distance: Distance) -> Vector3 {
        let direction = Vector3::from_spherical(self.azimuth.radians(), self.elevation.radians());
        direction * (distance.kilometers() * 1000.0)
    }

    // The observer's geocentric position, rotated from its meridian plane
    // (u from the axis, v from the equator) by the geodetic latitude.
    fn observer_vector(&self) -> Vector3 {
        let (u, v) = self.observer.to_geocentric_meters();
        let (sp, cp) = self.observer.latitude_angle().sin_cos();
        Vector3::new(-sp * u + cp * v, 0.0, cp * u + sp * v)
    }

    fn at_vector(&self, p: Vector3) -> CoordResult<Self> {
        let (azimuth, elevation) = p.to_spherical();
        Self::with_distance(
            Angle::from_radians(azimuth),
            Angle::from_radians(elevation),
            self.observer,
            self.epoch,
            Distance::from_kilometers(p.magnitude() / 1000.0)?,
        )
    }
}

fn check_outside_observer(geocentric_metres: f64, observer: Vector3) -> CoordResult<()> {
    if geocentric_metres <= observer.magnitude() {
        return Err(CoordError::invalid_distance(
            "diurnal parallax needs an object farther from the geocentre than the observer",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::topocentric::test_observer;
    use celestial_core::constants::HALF_PI;
    use celestial_core::location::Location;
    use celestial_time::scales::tt::TT;

    fn on_equator(azimuth: Angle, elevation: Angle, km: f64) -> TopocentricPosition {
        let observer = Location::from_degrees(0.0, 0.0, 0.0).unwrap();
        let distance = Distance::from_kilometers(km).unwrap();
        TopocentricPosition::with_distance(azimuth, elevation, observer, TT::j2000(), distance)
            .unwrap()
    }

    #[test]
    fn test_diurnal_parallax_moon() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Moon at mean distance: 384,400 km ≈ 0.00257 AU
        let moon_distance = Distance::from_kilometers(384400.0).unwrap();

        // Moon at horizon (maximum parallax)
        let moon_horizon = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(0.0),
            observer,
            epoch,
            moon_distance,
        )
        .unwrap();

        // Horizontal parallax for Moon: ~57 arcmin = 0.95°
        let h_parallax = moon_horizon.horizontal_parallax().unwrap().unwrap();
        assert!(h_parallax.degrees() > 0.9 && h_parallax.degrees() < 1.0);
        assert!(h_parallax.arcminutes() > 55.0 && h_parallax.arcminutes() < 59.0);

        // At the horizon the diurnal parallax is the horizontal parallax, apart from the
        // tilt of the geocentric zenith away from the geodetic one.
        let diurnal = moon_horizon.diurnal_parallax().unwrap().unwrap();
        assert!(libm::fabs(diurnal.degrees() - h_parallax.degrees()) < 0.001);
    }

    #[test]
    fn test_diurnal_parallax_sun() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Sun at 1 AU
        let sun_distance = Distance::from_au(1.0).unwrap();
        let sun_horizon = TopocentricPosition::with_distance(
            Angle::from_degrees(90.0),
            Angle::from_degrees(0.0),
            observer,
            epoch,
            sun_distance,
        )
        .unwrap();

        // Solar horizontal parallax: ~8.794 arcsec
        let h_parallax = sun_horizon.horizontal_parallax().unwrap().unwrap();
        assert!(h_parallax.arcseconds() > 8.7 && h_parallax.arcseconds() < 8.9);
    }

    #[test]
    fn test_diurnal_parallax_mars_opposition() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Mars at closest approach: ~0.38 AU
        let mars_distance = Distance::from_au(0.38).unwrap();
        let mars_horizon = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(0.0),
            observer,
            epoch,
            mars_distance,
        )
        .unwrap();

        // Mars horizontal parallax at opposition: ~23 arcsec
        let h_parallax = mars_horizon.horizontal_parallax().unwrap().unwrap();
        assert!(h_parallax.arcseconds() > 22.0 && h_parallax.arcseconds() < 24.0);
    }

    #[test]
    fn test_diurnal_parallax_at_various_elevations() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Moon at mean distance
        let moon_distance = Distance::from_kilometers(384400.0).unwrap();

        let elevations = vec![0.0, 30.0, 45.0, 60.0, 90.0];

        for elev in elevations {
            let pos = TopocentricPosition::with_distance(
                Angle::from_degrees(0.0),
                Angle::from_degrees(elev),
                observer,
                epoch,
                moon_distance,
            )
            .unwrap();

            let parallax = pos.diurnal_parallax().unwrap().unwrap();
            let h_par = pos.horizontal_parallax().unwrap().unwrap();

            // Off the equator the geocentric zenith is tilted from the
            // geodetic one, so the parallax stays positive even at 90°.
            if elev == 0.0 {
                assert!(libm::fabs(parallax.degrees() - h_par.degrees()) < 0.001);
            } else {
                assert!(parallax.degrees() > 0.0);
                assert!(parallax.degrees() < h_par.degrees());
            }
        }
    }

    #[test]
    fn test_diurnal_parallax_with_without() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Moon at 45° elevation
        let moon_distance = Distance::from_kilometers(384400.0).unwrap();
        let geocentric = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            moon_distance,
        )
        .unwrap();

        let topocentric = geocentric.with_diurnal_parallax().unwrap();

        assert!(topocentric.elevation().degrees() < geocentric.elevation().degrees());

        let back_to_geocentric = topocentric.without_diurnal_parallax().unwrap();
        assert_eq!(back_to_geocentric, geocentric);
    }

    #[test]
    fn test_diurnal_parallax_without_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Position without distance (star)
        let star_pos = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();

        assert_eq!(star_pos.diurnal_parallax().unwrap(), None);
        assert_eq!(star_pos.horizontal_parallax().unwrap(), None);

        let with_par = star_pos.with_diurnal_parallax().unwrap();
        assert_eq!(
            with_par.elevation().degrees(),
            star_pos.elevation().degrees()
        );

        let without_par = star_pos.without_diurnal_parallax().unwrap();
        assert_eq!(
            without_par.elevation().degrees(),
            star_pos.elevation().degrees()
        );
    }

    #[test]
    fn test_diurnal_parallax_on_equator() {
        // On the equator the observer sits a = 6378137 m straight up the
        // geocentric radius, so the parallax at the horizon is atan(a/d)
        // and at the zenith it vanishes.
        let horizon = on_equator(Angle::ZERO, Angle::ZERO, 384_400.0);
        let metres = horizon.distance().unwrap().kilometers() * 1000.0;
        let parallax = horizon.diurnal_parallax().unwrap().unwrap();
        assert_eq!(
            parallax.radians(),
            libm::atan2(WGS84_SEMI_MAJOR_AXIS, metres)
        );

        let zenith = on_equator(Angle::ZERO, Angle::from_radians(HALF_PI), 384_400.0);
        assert_eq!(zenith.diurnal_parallax().unwrap().unwrap().radians(), 0.0);
    }

    #[test]
    fn test_horizontal_parallax_uses_wgs84_radius() {
        let pos = on_equator(Angle::ZERO, Angle::from_degrees(30.0), 384_400.0);
        let metres = pos.distance().unwrap().kilometers() * 1000.0;
        let expected = libm::asin(WGS84_SEMI_MAJOR_AXIS / metres);
        assert_eq!(
            pos.horizontal_parallax().unwrap().unwrap().radians(),
            expected
        );
    }

    #[test]
    fn test_diurnal_parallax_rejects_objects_inside_the_earth() {
        let geocentric = on_equator(Angle::ZERO, Angle::from_degrees(-10.0), 5000.0);
        assert!(geocentric.with_diurnal_parallax().is_err());
        assert!(geocentric.diurnal_parallax().is_err());
        assert!(geocentric.horizontal_parallax().is_err());

        let topocentric = on_equator(Angle::ZERO, Angle::from_degrees(-60.0), 3000.0);
        assert!(topocentric.without_diurnal_parallax().is_err());
    }
}
