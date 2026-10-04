use super::{PrecessionModel, PrecessionResult};
use crate::scales::tt::TT;
use crate::{TimeError, TimeResult};
use celestial_core::precession::PrecessionIAU2000 as CoreCalculator;

pub(super) fn calculate(tt: &TT) -> TimeResult<PrecessionResult> {
    // Core checks the epoch too, but reports it as a calculation error.
    tt.centuries_since_j2000()?;
    let jd = tt.to_julian_date();
    let calculator = CoreCalculator::new();
    let core_result = calculator.compute(jd.jd1(), jd.jd2()).map_err(|e| {
        TimeError::CalculationError(format!("IAU 2000 precession calculation failed: {}", e))
    })?;

    Ok(PrecessionResult {
        bias_matrix: core_result.bias_matrix,
        precession_matrix: core_result.precession_matrix,
        bias_precession_matrix: core_result.bias_precession_matrix,
        model: PrecessionModel::IAU2000,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian::JulianDate;
    use crate::scales::tt::TT;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_iau2000_precession_matches_erfa() {
        // eraBp00
        let tt = TT::from_julian_date(JulianDate::new(J2000_JD, 9131.987654321));
        let result = calculate(&tt).unwrap();
        assert_eq!(result.model, PrecessionModel::IAU2000);
        assert_eq!(
            result.bias_matrix.elements(),
            &[
                [
                    0.9999999999999942,
                    -7.078279744199198e-8,
                    8.056217146976134e-8
                ],
                [
                    7.078279477857338e-8,
                    0.9999999999999969,
                    3.3060414542221364e-8
                ],
                [
                    -8.056217380986972e-8,
                    -3.306040883980552e-8,
                    0.9999999999999962
                ],
            ]
        );
        assert_eq!(
            result.precession_matrix.elements(),
            &[
                [
                    0.9999814200285311,
                    -0.005590937403741324,
                    -0.0024292008294464575
                ],
                [
                    0.005590937477350291,
                    0.9999843705640702,
                    -6.76051634278644e-6
                ],
                [
                    0.002429200660031423,
                    -6.821119224786898e-6,
                    0.9999970494644599
                ],
            ]
        );
        assert_eq!(
            result.bias_precession_matrix.elements(),
            &[
                [
                    0.9999814198284848,
                    -0.005591008104913233,
                    -0.0024291204536105297
                ],
                [
                    0.005591008259583384,
                    0.9999843701685484,
                    -6.727006026896097e-6
                ],
                [
                    0.002429120097612483,
                    -6.8543514816990364e-6,
                    0.9999970496599323
                ],
            ]
        );
    }

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
}
