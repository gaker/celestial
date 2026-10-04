use super::last::LAST;
use crate::scales::tt::TT;
use crate::scales::ut1::UT1;
use crate::transforms::nutation::NutationCalculator;
use crate::transforms::rotation::earth_rotation_angle;
use crate::TimeResult;
use celestial_core::angle::wrap_0_2pi;
use celestial_core::cio::CioSolution;

greenwich_sidereal_time!(GAST, calculate_gast_iau2006a, LAST, to_last);

fn calculate_gast_iau2006a(ut1: &UT1, tt: &TT) -> TimeResult<f64> {
    let era = earth_rotation_angle(&ut1.to_julian_date())?;

    let tt_centuries = tt.centuries_since_j2000()?;
    let npb_matrix = calculate_npb_matrix_iau2006a(tt, tt_centuries)?;

    let cio_solution = CioSolution::calculate(&npb_matrix, tt_centuries).map_err(|e| {
        crate::TimeError::CalculationError(format!("CIO calculation failed: {}", e))
    })?;

    let gast = era - cio_solution.equation_of_origins;

    Ok(wrap_0_2pi(gast)?)
}

fn calculate_npb_matrix_iau2006a(
    tt: &TT,
    t: f64,
) -> TimeResult<celestial_core::matrix::RotationMatrix3> {
    let nutation_result = tt.nutation_iau2006a()?;
    let dpsi = nutation_result.nutation_longitude();
    let deps = nutation_result.nutation_obliquity();

    let precession_calc = celestial_core::precession::PrecessionIAU2006::new();
    let npb_matrix = precession_calc.npb_matrix_iau2006a(t, dpsi, deps);

    Ok(npb_matrix)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use celestial_core::constants::MJD_ZERO_POINT;
    use celestial_core::location::Location;

    #[test]
    fn test_gast_j2000_matches_erfa() {
        // eraGst06a with UT1 = TT = J2000.0.
        let expected = 4.894899322716232;
        let gast = GAST::from_ut1_and_tt(&UT1::j2000(), &TT::j2000()).unwrap();
        assert_eq!(gast.radians(), expected);
    }

    #[test]
    fn test_to_last_matches_erfa() {
        // anp(eraGst06a + elong) with elong = -2.7144 rad. Each case changes in the
        // last bit if the angle goes through hours.
        let location = Location::new(0.5, -2.7144, 0.0).unwrap();
        let cases = [
            (51655.869, 51655.86979861111, 0.13165617455530532),
            (51692.992, 51692.99279861111, 1.5431071876122355),
            (51804.361, 51804.361798611106, 5.777460455752527),
        ];
        for (ut1_mjd, tt_mjd, expected) in cases {
            let ut1 = UT1::from_julian_date(JulianDate::new(MJD_ZERO_POINT, ut1_mjd));
            let tt = TT::from_julian_date(JulianDate::new(MJD_ZERO_POINT, tt_mjd));
            let gast = GAST::from_ut1_and_tt(&ut1, &tt).unwrap();
            let last = gast.to_last(&location).unwrap();
            assert_eq!(last.radians(), expected, "UT1 MJD {ut1_mjd}");
        }
    }
}
