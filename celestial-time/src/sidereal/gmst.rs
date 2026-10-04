use super::lmst::LMST;
use crate::scales::tt::TT;
use crate::scales::ut1::UT1;
use crate::transforms::rotation::earth_rotation_angle;
use crate::TimeResult;
use celestial_core::angle::wrap_0_2pi;

greenwich_sidereal_time!(GMST, calculate_gmst_iau2006, LMST, to_lmst);

fn calculate_gmst_iau2006(ut1: &UT1, tt: &TT) -> TimeResult<f64> {
    let t = tt.centuries_since_j2000()?;

    let era = earth_rotation_angle(&ut1.to_julian_date())?;

    // IAU 2006 polynomial correction (Horner's method for precision)
    let polynomial_arcsec = 0.014506
        + t * (4612.156534
            + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * (-0.0000000368)))));

    let gmst = era + polynomial_arcsec * celestial_core::constants::ARCSEC_TO_RAD;

    Ok(wrap_0_2pi(gmst)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use crate::TimeError;
    use celestial_core::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD, MJD_ZERO_POINT};
    use celestial_core::location::Location;

    #[test]
    fn test_gmst_j2000_matches_erfa() {
        // eraGmst06 with UT1 = TT = J2000.0.
        let expected = 4.8949612831508285;
        let gmst = GMST::from_ut1_and_tt(&UT1::j2000(), &TT::j2000()).unwrap();
        assert_eq!(gmst.radians(), expected);
    }

    #[test]
    fn test_to_lmst_matches_erfa() {
        // anp(eraGmst06 + elong) with elong = -2.7144 rad. Each case changes in the
        // last bit if the angle goes through hours.
        let location = Location::new(0.5, -2.7144, 0.0).unwrap();
        let cases = [
            (51692.992, 51692.99279861111, 1.543180103011061),
            (51804.361, 51804.361798611106, 5.777533198070602),
            (51915.729999999996, 51915.730798611105, 3.728700986076002),
        ];
        for (ut1_mjd, tt_mjd, expected) in cases {
            let ut1 = UT1::from_julian_date(JulianDate::new(MJD_ZERO_POINT, ut1_mjd));
            let tt = TT::from_julian_date(JulianDate::new(MJD_ZERO_POINT, tt_mjd));
            let gmst = GMST::from_ut1_and_tt(&ut1, &tt).unwrap();
            let lmst = gmst.to_lmst(&location).unwrap();
            assert_eq!(lmst.radians(), expected, "UT1 MJD {ut1_mjd}");
        }
    }

    #[test]
    fn test_twenty_centuries_is_the_limit() {
        // eraGmst06 with UT1 = TT = J2000.0 ± 20 centuries.
        let span = 20.0 * DAYS_PER_JULIAN_CENTURY;
        for (jd, expected) in [
            (J2000_JD - span, 4.628864591802277),
            (J2000_JD + span, 5.166408763419673),
        ] {
            let jd = JulianDate::new(jd, 0.0);
            let gmst = gmst_at(jd, jd).unwrap();
            assert_eq!(gmst.radians(), expected, "JD {jd}");
        }

        let past = JulianDate::new(J2000_JD - 21.0 * DAYS_PER_JULIAN_CENTURY, 0.0);
        assert_eq!(
            gmst_at(past, past).unwrap_err(),
            TimeError::InvalidEpoch(
                "Epoch too far from J2000.0 for the IAU models: -21.0 centuries".into()
            )
        );
    }

    #[test]
    fn test_out_of_range_input_is_rejected() {
        let j2000 = JulianDate::new(J2000_JD, 0.0);
        let far = JulianDate::new(J2000_JD + 1e13, 0.0);
        let nan = JulianDate::new(f64::NAN, 0.0);
        assert_eq!(
            gmst_at(j2000, far).unwrap_err(),
            TimeError::InvalidEpoch(
                "Epoch too far from J2000.0 for the IAU models: 273785078.7 centuries".into()
            )
        );
        assert_eq!(
            gmst_at(j2000, nan).unwrap_err(),
            TimeError::InvalidEpoch("Julian Date (NaN, 0) is not finite".into())
        );
        assert_eq!(
            gmst_at(far, j2000).unwrap_err(),
            TimeError::CalculationError(
                "Time value out of valid range: 10000000000000 days from J2000".into()
            )
        );
    }

    fn gmst_at(ut1: JulianDate, tt: JulianDate) -> TimeResult<GMST> {
        GMST::from_ut1_and_tt(&UT1::from_julian_date(ut1), &TT::from_julian_date(tt))
    }
}
