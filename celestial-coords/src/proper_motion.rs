use celestial_core::Angle;
use celestial_time::TT;

const MAS_PER_DEGREE: f64 = 3_600_000.0;
const DAYS_PER_YEAR: f64 = 365.25;

pub fn propagate(
    ra: Angle,
    dec: Angle,
    pm_ra_cos_dec_mas_per_year: f64,
    pm_dec_mas_per_year: f64,
    from_epoch: TT,
    to_epoch: TT,
) -> (Angle, Angle) {
    let dt_years = epoch_delta_years(from_epoch, to_epoch);
    let cos_dec = libm::cos(dec.radians());

    let delta_ra_deg = if cos_dec.abs() > f64::EPSILON {
        pm_ra_cos_dec_mas_per_year * dt_years / MAS_PER_DEGREE / cos_dec
    } else {
        0.0
    };
    let delta_dec_deg = pm_dec_mas_per_year * dt_years / MAS_PER_DEGREE;

    let ra_out = Angle::from_degrees(ra.degrees() + delta_ra_deg);
    let dec_out = Angle::from_degrees(dec.degrees() + delta_dec_deg);
    (ra_out, dec_out)
}

pub fn propagate_degrees(
    ra_deg: f64,
    dec_deg: f64,
    pm_ra_cos_dec_mas_per_year: f64,
    pm_dec_mas_per_year: f64,
    from_epoch: TT,
    to_epoch: TT,
) -> (f64, f64) {
    let (ra, dec) = propagate(
        Angle::from_degrees(ra_deg),
        Angle::from_degrees(dec_deg),
        pm_ra_cos_dec_mas_per_year,
        pm_dec_mas_per_year,
        from_epoch,
        to_epoch,
    );
    (ra.degrees(), dec.degrees())
}

fn epoch_delta_years(from_epoch: TT, to_epoch: TT) -> f64 {
    let from_jd = from_epoch.to_julian_date();
    let to_jd = to_epoch.to_julian_date();
    (to_jd - from_jd).to_f64() / DAYS_PER_YEAR
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;
    use celestial_time::JulianDate;

    fn tt_at_jd(jd: f64) -> TT {
        TT::from_julian_date(JulianDate::new(jd, 0.0))
    }

    fn j2016() -> TT {
        tt_at_jd(2_457_389.0)
    }

    #[test]
    fn zero_pm_returns_input() {
        let from = j2016();
        let to = tt_at_jd(2_457_389.0 + DAYS_PER_YEAR);
        let (ra, dec) = propagate_degrees(100.0, 45.0, 0.0, 0.0, from, to);
        assert!((ra - 100.0).abs() < 1e-12);
        assert!((dec - 45.0).abs() < 1e-12);
    }

    #[test]
    fn one_year_propagation_matches_hand_calc() {
        let from = j2016();
        let to = tt_at_jd(2_457_389.0 + DAYS_PER_YEAR);
        let (ra, dec) = propagate_degrees(100.0, 45.0, 3600.0, 3600.0, from, to);

        let expected_dec = 45.0 + 3600.0 / MAS_PER_DEGREE;
        let cos_dec = libm::cos(45.0_f64.to_radians());
        let expected_ra = 100.0 + (3600.0 / MAS_PER_DEGREE) / cos_dec;

        assert!((dec - expected_dec).abs() < 1e-12);
        assert!((ra - expected_ra).abs() < 1e-12);
    }

    #[test]
    fn negative_delta_reverses_motion() {
        let from = tt_at_jd(2_457_389.0 + DAYS_PER_YEAR);
        let to = j2016();
        let (ra, dec) = propagate_degrees(100.0, 45.0, 3600.0, 3600.0, from, to);

        let expected_dec = 45.0 - 3600.0 / MAS_PER_DEGREE;
        let cos_dec = libm::cos(45.0_f64.to_radians());
        let expected_ra = 100.0 - (3600.0 / MAS_PER_DEGREE) / cos_dec;

        assert!((dec - expected_dec).abs() < 1e-12);
        assert!((ra - expected_ra).abs() < 1e-12);
    }

    #[test]
    fn j2000_to_j2016_kochab_pm() {
        let from = tt_at_jd(J2000_JD);
        let to = j2016();
        let pm_ra_cos_dec = -32.29;
        let pm_dec = 11.91;
        let (ra, dec) = propagate_degrees(222.6764, 74.1555, pm_ra_cos_dec, pm_dec, from, to);

        let dt_years = (2_457_389.0 - J2000_JD) / DAYS_PER_YEAR;
        let cos_dec = libm::cos(74.1555_f64.to_radians());
        let expected_ra = 222.6764 + pm_ra_cos_dec * dt_years / MAS_PER_DEGREE / cos_dec;
        let expected_dec = 74.1555 + pm_dec * dt_years / MAS_PER_DEGREE;

        assert!((ra - expected_ra).abs() < 1e-10);
        assert!((dec - expected_dec).abs() < 1e-10);
    }

    #[test]
    fn pole_dec_does_not_blow_up_ra() {
        let from = j2016();
        let to = tt_at_jd(2_457_389.0 + DAYS_PER_YEAR);
        let (ra, _dec) = propagate_degrees(0.0, 90.0, 1000.0, 0.0, from, to);
        assert!(ra.is_finite());
        assert_eq!(ra, 0.0);
    }
}
