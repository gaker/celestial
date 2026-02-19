use celestial_coords::Vector3;
use celestial_core::AstroResult;
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

use crate::earth::Vsop2013Earth;

const DT_DAYS: f64 = 1.0 / celestial_core::constants::SECONDS_PER_DAY_F64;

pub struct Vsop2013Sun;

impl Vsop2013Sun {
    pub fn heliocentric_position(&self, _tdb: &TDB) -> AstroResult<Vector3> {
        Ok(Vector3::new(0.0, 0.0, 0.0))
    }

    pub fn geocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
        let earth = Vsop2013Earth::new();
        let earth_pos = earth.heliocentric_position(tdb)?;
        Ok(Vector3::new(-earth_pos.x, -earth_pos.y, -earth_pos.z))
    }

    pub fn geocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
        let pos = self.geocentric_position(tdb)?;
        let jd = tdb.to_julian_date();
        let t_minus = TDB::from_julian_date(JulianDate::new(jd.jd1(), jd.jd2() - DT_DAYS));
        let t_plus = TDB::from_julian_date(JulianDate::new(jd.jd1(), jd.jd2() + DT_DAYS));
        let p_minus = self.geocentric_position(&t_minus)?;
        let p_plus = self.geocentric_position(&t_plus)?;
        let inv_2dt = 1.0 / (2.0 * DT_DAYS);
        let vel = Vector3::new(
            (p_plus.x - p_minus.x) * inv_2dt,
            (p_plus.y - p_minus.y) * inv_2dt,
            (p_plus.z - p_minus.z) * inv_2dt,
        );
        Ok((pos, vel))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;
    use celestial_time::julian::JulianDate;

    #[test]
    fn sun_heliocentric_is_origin() {
        let sun = Vsop2013Sun;
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let pos = sun.heliocentric_position(&tdb).unwrap();

        assert_eq!(pos.x, 0.0);
        assert_eq!(pos.y, 0.0);
        assert_eq!(pos.z, 0.0);
    }

    #[test]
    fn sun_geocentric_is_negative_earth() {
        let sun = Vsop2013Sun;
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let sun_geo = sun.geocentric_position(&tdb).unwrap();

        let dist_au = libm::sqrt(sun_geo.x.powi(2) + sun_geo.y.powi(2) + sun_geo.z.powi(2));
        assert!(
            dist_au > 0.98 && dist_au < 1.02,
            "Sun-Earth distance {} AU should be ~1 AU",
            dist_au
        );
    }

    #[test]
    fn sun_geocentric_various_epochs() {
        let sun = Vsop2013Sun;

        for days_offset in [0, 91, 182, 273] {
            let jd = J2000_JD + days_offset as f64;
            let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
            let sun_geo = sun.geocentric_position(&tdb).unwrap();
            let dist_au = libm::sqrt(sun_geo.x.powi(2) + sun_geo.y.powi(2) + sun_geo.z.powi(2));

            assert!(
                dist_au > 0.983 && dist_au < 1.017,
                "Day {}: distance {} AU outside expected range (0.983-1.017 AU)",
                days_offset,
                dist_au
            );
        }
    }

    #[test]
    fn sun_geocentric_state_velocity() {
        let sun = Vsop2013Sun;
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let (pos, vel) = sun.geocentric_state(&tdb).unwrap();

        let dist_au = libm::sqrt(pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2));
        assert!(dist_au > 0.98 && dist_au < 1.02);

        let speed_au_day = libm::sqrt(vel.x.powi(2) + vel.y.powi(2) + vel.z.powi(2));
        assert!(
            speed_au_day > 0.016 && speed_au_day < 0.018,
            "Sun apparent speed {} AU/day should match Earth orbital speed ~0.017 AU/day",
            speed_au_day
        );
    }
}
