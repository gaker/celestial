use celestial_coords::Vector3;
use celestial_core::constants::{AU_KM, MOON_EARTH_MASS_RATIO};
use celestial_core::AstroResult;
use celestial_time::julian::JulianDate;
use celestial_time::TDB;

use crate::moon::ElpMpp02Moon;
use crate::planets::Vsop2013Emb;

const DT_DAYS: f64 = 1.0 / celestial_core::constants::SECONDS_PER_DAY_F64;

pub struct Vsop2013Earth {
    emb: Vsop2013Emb,
    moon: ElpMpp02Moon,
}

impl Default for Vsop2013Earth {
    fn default() -> Self {
        Self::new()
    }
}

impl Vsop2013Earth {
    pub fn new() -> Self {
        Self {
            emb: Vsop2013Emb,
            moon: ElpMpp02Moon::new(),
        }
    }

    pub fn heliocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
        let emb_pos = self.emb.heliocentric_position(tdb)?;
        let moon_geo_km = self.moon.geocentric_position_icrs(tdb)?;
        let moon_geo_au = [
            moon_geo_km[0] / AU_KM,
            moon_geo_km[1] / AU_KM,
            moon_geo_km[2] / AU_KM,
        ];

        Ok(Vector3::new(
            emb_pos.x - moon_geo_au[0] * MOON_EARTH_MASS_RATIO,
            emb_pos.y - moon_geo_au[1] * MOON_EARTH_MASS_RATIO,
            emb_pos.z - moon_geo_au[2] * MOON_EARTH_MASS_RATIO,
        ))
    }

    pub fn heliocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
        let pos = self.heliocentric_position(tdb)?;
        let jd = tdb.to_julian_date();
        let t_minus = TDB::from_julian_date(JulianDate::new(jd.jd1(), jd.jd2() - DT_DAYS));
        let t_plus = TDB::from_julian_date(JulianDate::new(jd.jd1(), jd.jd2() + DT_DAYS));
        let p_minus = self.heliocentric_position(&t_minus)?;
        let p_plus = self.heliocentric_position(&t_plus)?;
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
    fn earth_differs_from_emb() {
        let earth = Vsop2013Earth::new();
        let emb = Vsop2013Emb;
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let earth_pos = earth.heliocentric_position(&tdb).unwrap();
        let emb_pos = emb.heliocentric_position(&tdb).unwrap();

        let diff_au = ((earth_pos.x - emb_pos.x).powi(2)
            + (earth_pos.y - emb_pos.y).powi(2)
            + (earth_pos.z - emb_pos.z).powi(2))
        .sqrt();
        let diff_km = diff_au * AU_KM;

        println!("Earth-EMB difference at J2000: {:.1} km", diff_km);
        assert!(
            diff_km > 4000.0 && diff_km < 5000.0,
            "Earth-EMB difference {} km should be ~4670 km (Moon pulls EMB toward it)",
            diff_km
        );
    }

    #[test]
    fn earth_heliocentric_distance() {
        let earth = Vsop2013Earth::new();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let pos = earth.heliocentric_position(&tdb).unwrap();
        let dist_au = (pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2)).sqrt();

        assert!(
            dist_au > 0.98 && dist_au < 1.02,
            "Earth heliocentric distance {} AU should be ~1 AU",
            dist_au
        );
    }

    #[test]
    fn earth_heliocentric_seasonal_variation() {
        let earth = Vsop2013Earth::new();

        let perihelion_jd = 2451547.5;
        let aphelion_jd = 2451730.5;

        let tdb_peri = TDB::from_julian_date(JulianDate::new(perihelion_jd, 0.0));
        let tdb_aph = TDB::from_julian_date(JulianDate::new(aphelion_jd, 0.0));

        let pos_peri = earth.heliocentric_position(&tdb_peri).unwrap();
        let pos_aph = earth.heliocentric_position(&tdb_aph).unwrap();

        let dist_peri = (pos_peri.x.powi(2) + pos_peri.y.powi(2) + pos_peri.z.powi(2)).sqrt();
        let dist_aph = (pos_aph.x.powi(2) + pos_aph.y.powi(2) + pos_aph.z.powi(2)).sqrt();

        println!("Perihelion distance: {:.6} AU", dist_peri);
        println!("Aphelion distance: {:.6} AU", dist_aph);

        assert!(
            dist_peri < dist_aph,
            "Perihelion {} AU should be less than aphelion {} AU",
            dist_peri,
            dist_aph
        );
        assert!(
            dist_peri > 0.98 && dist_peri < 0.985,
            "Perihelion {} AU should be ~0.983 AU",
            dist_peri
        );
        assert!(
            dist_aph > 1.01 && dist_aph < 1.02,
            "Aphelion {} AU should be ~1.017 AU",
            dist_aph
        );
    }

    #[test]
    fn earth_default_impl() {
        let earth: Vsop2013Earth = Default::default();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let pos = earth.heliocentric_position(&tdb).unwrap();

        let dist_au = (pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2)).sqrt();
        assert!(dist_au > 0.98 && dist_au < 1.02);
    }

    #[test]
    fn earth_heliocentric_state_velocity() {
        let earth = Vsop2013Earth::new();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let (pos, vel) = earth.heliocentric_state(&tdb).unwrap();

        let dist_au = (pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2)).sqrt();
        assert!(dist_au > 0.98 && dist_au < 1.02);

        let speed_au_day = (vel.x.powi(2) + vel.y.powi(2) + vel.z.powi(2)).sqrt();
        assert!(
            speed_au_day > 0.016 && speed_au_day < 0.018,
            "Earth orbital speed {} AU/day should be ~0.017 AU/day (~30 km/s)",
            speed_au_day
        );
    }
}
