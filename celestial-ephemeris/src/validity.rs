use celestial_core::constants::{DAYS_PER_JULIAN_MILLENNIUM, J2000_JD};
use celestial_core::errors::{AstroError, AstroResult, MathErrorKind};
use celestial_time::scales::tdb::TDB;

pub(crate) struct Window {
    theory: &'static str,
    years: &'static str,
    first_day: f64,
    last_day: f64,
}

// Simon et al. 2013 (references/ephemeris/VSOP2013.md): -4000 to +8000, and
// 0 to 4000 for Pluto.
pub(crate) const VSOP2013: Window = Window::millennia("VSOP2013", "-4000 to +8000", -6.0, 6.0);
pub(crate) const VSOP2013_PLUTO: Window =
    Window::millennia("VSOP2013 Pluto", "0 to +4000", -2.0, 2.0);

// DE406's span: the long-span JPL integration of ELP/MPP02's era.
pub(crate) const ELPMPP02: Window = Window::millennia("ELP/MPP02", "-3000 to +3000", -5.0, 1.0);

impl Window {
    const fn millennia(theory: &'static str, years: &'static str, first: f64, last: f64) -> Self {
        Self {
            theory,
            years,
            first_day: first * DAYS_PER_JULIAN_MILLENNIUM,
            last_day: last * DAYS_PER_JULIAN_MILLENNIUM,
        }
    }

    pub(crate) fn days_from_j2000(&self, tdb: &TDB) -> AstroResult<f64> {
        let jd = tdb.to_julian_date();
        let days = (jd.jd1() - J2000_JD) + jd.jd2();
        if !days.is_finite() {
            return Err(self.error(MathErrorKind::NotFinite, "TDB epoch must be finite"));
        }
        if days >= self.first_day && days <= self.last_day {
            return Ok(days);
        }
        let reason = format!(
            "TDB JD {} is outside the theory's range, years {}",
            J2000_JD + days,
            self.years
        );
        Err(self.error(MathErrorKind::OutOfRange, &reason))
    }

    fn error(&self, kind: MathErrorKind, reason: &str) -> AstroError {
        AstroError::math_error(self.theory, kind, reason)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_time::julian::JulianDate;

    fn tdb(jd1: f64, jd2: f64) -> TDB {
        TDB::from_julian_date(JulianDate::new(jd1, jd2))
    }

    #[test]
    fn days_are_measured_from_j2000_without_summing_the_parts() {
        let days = VSOP2013.days_from_j2000(&tdb(J2000_JD, 0.25)).unwrap();
        assert_eq!(days, 0.25);
        let days = VSOP2013.days_from_j2000(&tdb(2460000.5, 1e-9)).unwrap();
        assert_eq!(days, (2460000.5 - J2000_JD) + 1e-9);
    }

    #[test]
    fn range_edges_are_inclusive() {
        let edge = |w: &Window, jd: f64| w.days_from_j2000(&tdb(jd, 0.0)).unwrap();
        assert_eq!(edge(&VSOP2013, 260045.0), -2191500.0);
        assert_eq!(edge(&VSOP2013, 4643045.0), 2191500.0);
        assert_eq!(edge(&ELPMPP02, 625295.0), -1826250.0);
        assert_eq!(edge(&ELPMPP02, 2816795.0), 365250.0);
    }

    #[test]
    fn out_of_range_error_names_the_theory_and_its_years() {
        let err = VSOP2013_PLUTO
            .days_from_j2000(&tdb(3182046.0, 0.0))
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "Math error in VSOP2013 Pluto (out of range): \
             TDB JD 3182046 is outside the theory's range, years 0 to +4000"
        );
    }

    #[test]
    fn non_finite_error_names_the_theory() {
        let err = ELPMPP02.days_from_j2000(&tdb(f64::NAN, 0.0)).unwrap_err();
        assert_eq!(
            err.to_string(),
            "Math error in ELP/MPP02 (not finite): TDB epoch must be finite"
        );
    }
}
