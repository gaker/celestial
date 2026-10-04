use super::{NutationModel, NutationResult};
use crate::scales::tt::TT;
use crate::{TimeError, TimeResult};
use celestial_core::nutation::NutationIAU2000B as CoreCalculator;

pub(super) fn calculate(tt: &TT) -> TimeResult<NutationResult> {
    // Core checks the epoch too, but reports it as a calculation error.
    tt.centuries_since_j2000()?;
    let calculator = CoreCalculator::new();
    let jd = tt.to_julian_date();
    let core_result = calculator.compute(jd.jd1(), jd.jd2()).map_err(|e| {
        TimeError::CalculationError(format!("IAU 2000B nutation calculation failed: {}", e))
    })?;

    Ok(NutationResult::new(core_result, NutationModel::IAU2000B))
}
