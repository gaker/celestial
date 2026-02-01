use super::{PrecessionModel, PrecessionResult};
use crate::{TimeError, TimeResult, TT};
use astro_core::precession::iau2006::PrecessionIAU2006 as CoreCalculator;

pub fn calculate(tt: &TT) -> TimeResult<PrecessionResult> {
    let jd = tt.to_julian_date();
    let calculator = CoreCalculator::new();
    let core_result = calculator.compute(jd.jd1(), jd.jd2()).map_err(|_| {
        TimeError::CalculationError("IAU 2006 precession calculation failed".to_string())
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
    use crate::TT;
    use astro_core::constants::J2000_JD;

    #[test]
    fn test_iau2006_precession_calculation() {
        let tt = TT::j2000();
        let result = calculate(&tt).unwrap();

        assert_eq!(result.model, PrecessionModel::IAU2006);

        let bias = result.bias_matrix.elements();
        let prec = result.precession_matrix.elements();
        let bp = result.bias_precession_matrix.elements();

        assert_eq!(bias.len(), 3);
        assert_eq!(prec.len(), 3);
        assert_eq!(bp.len(), 3);

        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                let diff = (result.bias_precession_matrix.get(i, j) - expected).abs();
                assert!(
                    diff < 1e-6,
                    "Bias-precession matrix at J2000 should be near identity"
                );
            }
        }
    }

    #[test]
    fn test_iau2006_precession_matrices_valid() {
        let tt = TT::j2000();
        let result = calculate(&tt).unwrap();

        for i in 0..3 {
            for j in 0..3 {
                assert!(result.bias_matrix.get(i, j).is_finite());
                assert!(result.precession_matrix.get(i, j).is_finite());
                assert!(result.bias_precession_matrix.get(i, j).is_finite());
            }
        }
    }

    #[test]
    fn test_iau2006_two_part_julian_date() {
        use crate::JulianDate;

        let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 0.5));
        let result = calculate(&tt).unwrap();

        assert_eq!(result.model, PrecessionModel::IAU2006);

        for i in 0..3 {
            for j in 0..3 {
                assert!(result.bias_matrix.get(i, j).is_finite());
                assert!(result.precession_matrix.get(i, j).is_finite());
                assert!(result.bias_precession_matrix.get(i, j).is_finite());
            }
        }
    }
}
