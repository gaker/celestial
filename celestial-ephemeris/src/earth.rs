use celestial_core::constants::MOON_EMB_MASS_RATIO;
use celestial_core::errors::AstroResult;
use celestial_core::matrix::Vector3;
use celestial_time::scales::tdb::TDB;

use crate::moon::ElpMpp02Moon;
use crate::planets::Vsop2013Emb;

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
        let emb = self.emb.heliocentric_position(tdb)?;
        let moon = self.moon.geocentric_position(tdb)?;
        Ok(emb - moon * MOON_EMB_MASS_RATIO)
    }

    pub fn heliocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
        let (emb_position, emb_velocity) = self.emb.heliocentric_state(tdb)?;
        let (moon_position, moon_velocity) = self.moon.geocentric_state(tdb)?;
        Ok((
            emb_position - moon_position * MOON_EMB_MASS_RATIO,
            emb_velocity - moon_velocity * MOON_EMB_MASS_RATIO,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;
    use celestial_core::errors::{AstroError, MathErrorKind};
    use celestial_time::julian::JulianDate;

    #[test]
    fn default_is_new() {
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        assert_eq!(
            Vsop2013Earth::default().heliocentric_state(&tdb).unwrap(),
            Vsop2013Earth::new().heliocentric_state(&tdb).unwrap()
        );
    }

    #[test]
    fn earth_is_the_emb_less_the_moon_share() {
        let tdb = TDB::from_julian_date(JulianDate::new(2445665.5, 0.45541555912031195));
        let (emb_position, emb_velocity) = Vsop2013Emb.heliocentric_state(&tdb).unwrap();
        let (moon_position, moon_velocity) = ElpMpp02Moon::new().geocentric_state(&tdb).unwrap();
        let position = emb_position - moon_position * MOON_EMB_MASS_RATIO;
        let velocity = emb_velocity - moon_velocity * MOON_EMB_MASS_RATIO;
        let earth = Vsop2013Earth::new();
        assert_eq!(
            earth.heliocentric_state(&tdb).unwrap(),
            (position, velocity)
        );
        assert_eq!(earth.heliocentric_position(&tdb).unwrap(), position);
    }

    fn error_kind<T: std::fmt::Debug>(result: AstroResult<T>) -> MathErrorKind {
        match result {
            Err(AstroError::MathError { kind, .. }) => kind,
            other => panic!("expected a MathError, got {:?}", other),
        }
    }

    #[test]
    fn non_finite_epochs_are_errors() {
        let earth = Vsop2013Earth::new();
        for (jd1, jd2) in [(f64::NAN, 0.0), (J2000_JD, f64::INFINITY)] {
            let tdb = TDB::from_julian_date(JulianDate::new(jd1, jd2));
            let kinds = [
                error_kind(earth.heliocentric_position(&tdb)),
                error_kind(earth.heliocentric_state(&tdb)),
            ];
            assert_eq!(kinds, [MathErrorKind::NotFinite; 2]);
        }
    }

    #[test]
    fn epochs_outside_the_moon_range_are_errors() {
        let earth = Vsop2013Earth::new();
        let tdb = TDB::from_julian_date(JulianDate::new(3_000_000.0, 0.0));
        assert!(Vsop2013Emb.heliocentric_position(&tdb).is_ok());
        assert_eq!(
            error_kind(earth.heliocentric_position(&tdb)),
            MathErrorKind::OutOfRange
        );
    }
}
