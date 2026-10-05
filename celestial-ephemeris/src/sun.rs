use celestial_core::errors::AstroResult;
use celestial_core::matrix::Vector3;
use celestial_time::scales::tdb::TDB;

use crate::validity::VSOP2013;

pub struct Vsop2013Sun;

impl Vsop2013Sun {
    pub fn heliocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
        VSOP2013.days_from_j2000(tdb)?;
        Ok(Vector3::zeros())
    }

    pub fn heliocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
        Ok((self.heliocentric_position(tdb)?, Vector3::zeros()))
    }

    pub fn geocentric_position(&self, tdb: &TDB, earth: &Vector3) -> AstroResult<Vector3> {
        Ok(self.heliocentric_position(tdb)? - *earth)
    }

    pub fn geocentric_state(
        &self,
        tdb: &TDB,
        (earth_position, earth_velocity): &(Vector3, Vector3),
    ) -> AstroResult<(Vector3, Vector3)> {
        let (position, velocity) = self.heliocentric_state(tdb)?;
        Ok((position - *earth_position, velocity - *earth_velocity))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;
    use celestial_core::errors::{AstroError, MathErrorKind};
    use celestial_time::julian::JulianDate;

    use crate::earth::Vsop2013Earth;

    fn tdb(jd1: f64, jd2: f64) -> TDB {
        TDB::from_julian_date(JulianDate::new(jd1, jd2))
    }

    #[test]
    fn heliocentric_state_is_the_origin_at_rest() {
        let t = tdb(J2000_JD, 0.25);
        let origin = Vector3::zeros();
        assert_eq!(Vsop2013Sun.heliocentric_position(&t).unwrap(), origin);
        assert_eq!(
            Vsop2013Sun.heliocentric_state(&t).unwrap(),
            (origin, origin)
        );
    }

    #[test]
    fn geocentric_state_is_minus_the_earth_state() {
        let earth = Vsop2013Earth::new();
        for jd2 in [0.0, 91.25, 182.5, 273.75] {
            let t = tdb(J2000_JD, jd2);
            let (position, velocity) = earth.heliocentric_state(&t).unwrap();
            let sun = Vsop2013Sun.geocentric_state(&t, &(position, velocity));
            assert_eq!(sun.unwrap(), (-position, -velocity));
            let sun = Vsop2013Sun.geocentric_position(&t, &position);
            assert_eq!(sun.unwrap(), -position);
        }
    }

    fn error_kind<T: std::fmt::Debug>(result: AstroResult<T>) -> MathErrorKind {
        match result {
            Err(AstroError::MathError { kind, .. }) => kind,
            other => panic!("expected a MathError, got {:?}", other),
        }
    }

    #[test]
    fn non_finite_epochs_are_errors() {
        let sun = Vsop2013Sun;
        let earth = (Vector3::x_axis(), Vector3::y_axis() * 0.0172);
        for (jd1, jd2) in [(f64::NAN, 0.0), (J2000_JD, f64::NEG_INFINITY)] {
            let t = tdb(jd1, jd2);
            let kinds = [
                error_kind(sun.heliocentric_position(&t)),
                error_kind(sun.heliocentric_state(&t)),
                error_kind(sun.geocentric_position(&t, &earth.0)),
                error_kind(sun.geocentric_state(&t, &earth)),
            ];
            assert_eq!(kinds, [MathErrorKind::NotFinite; 4]);
        }
    }

    #[test]
    fn epochs_outside_the_theory_range_are_errors() {
        assert_eq!(
            error_kind(Vsop2013Sun.heliocentric_position(&tdb(1e12, 0.0))),
            MathErrorKind::OutOfRange
        );
    }
}
