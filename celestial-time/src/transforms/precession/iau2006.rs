use super::{PrecessionModel, PrecessionResult};
use crate::scales::tt::TT;
use crate::{TimeError, TimeResult};
use celestial_core::precession::PrecessionIAU2006 as CoreCalculator;

pub(super) fn calculate(tt: &TT) -> TimeResult<PrecessionResult> {
    // Core checks the epoch too, but reports it as a calculation error.
    tt.centuries_since_j2000()?;
    let jd = tt.to_julian_date();
    let calculator = CoreCalculator::new();
    let core_result = calculator.compute(jd.jd1(), jd.jd2()).map_err(|e| {
        TimeError::CalculationError(format!("IAU 2006 precession calculation failed: {}", e))
    })?;

    Ok(PrecessionResult {
        bias_matrix: core_result.bias_matrix,
        precession_matrix: core_result.precession_matrix,
        bias_precession_matrix: core_result.bias_precession_matrix,
        model: PrecessionModel::IAU2006,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use crate::scales::tt::TT;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_bad_epoch_is_rejected() {
        let tt = TT::from_julian_date(JulianDate::new(f64::NAN, 0.0));
        assert_eq!(
            calculate(&tt),
            Err(TimeError::InvalidEpoch(
                "Julian Date (NaN, 0) is not finite".into()
            ))
        );
    }

    #[test]
    fn test_iau2006_precession_matches_erfa() {
        // eraBp06
        let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 9131.987654321));
        let result = calculate(&tt).unwrap();
        assert_eq!(result.model, PrecessionModel::IAU2006);
        assert_eq!(
            result.bias_matrix.elements(),
            &[
                [
                    0.9999999999999941,
                    -7.078368960971556e-8,
                    8.056213977613186e-8
                ],
                [
                    7.078368694637676e-8,
                    0.9999999999999969,
                    3.3059437354321375e-8
                ],
                [
                    -8.056214211620057e-8,
                    -3.305943169218395e-8,
                    0.9999999999999962
                ],
            ]
        );
        assert_eq!(
            result.precession_matrix.elements(),
            &[
                [
                    0.9999814200444181,
                    -0.005590934815983376,
                    -0.0024292002453885774
                ],
                [
                    0.0055909348911237,
                    0.9999843705785343,
                    -6.759881166231728e-6
                ],
                [
                    0.002429200072449081,
                    -6.821744841607659e-6,
                    0.9999970494658831
                ],
            ]
        );
        assert_eq!(
            result.bias_precession_matrix.elements(),
            &[
                [
                    0.9999814198443672,
                    -0.005591005518049832,
                    -0.0024291198695787926
                ],
                [
                    0.005591005674248898,
                    0.9999843701830078,
                    -6.726371827858735e-6
                ],
                [
                    0.0024291195100617845,
                    -6.854976123460421e-6,
                    0.9999970496613553
                ],
            ]
        );
    }
}
