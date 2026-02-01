use super::{NutationModel, NutationResult};
use crate::{TimeError, TimeResult, TT};
use astro_core::nutation::NutationIAU2000A as CoreCalculator;

pub fn calculate(tt: &TT) -> TimeResult<NutationResult> {
    let jd = tt.to_julian_date();

    let calculator = CoreCalculator::new();
    let core_result = calculator.compute(jd.jd1(), jd.jd2()).map_err(|_| {
        TimeError::CalculationError("IAU 2000A nutation calculation failed".to_string())
    })?;

    Ok(NutationResult::new(core_result, NutationModel::IAU2000A))
}
