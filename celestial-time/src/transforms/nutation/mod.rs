mod iau2000a;
mod iau2000b;
mod iau2006a;

use crate::scales::tt::TT;
use crate::TimeResult;

#[derive(Debug)]
pub struct NutationResult {
    core_result: celestial_core::nutation::NutationResult,
    model: NutationModel,
}

impl NutationResult {
    pub(crate) fn new(
        core_result: celestial_core::nutation::NutationResult,
        model: NutationModel,
    ) -> Self {
        Self { core_result, model }
    }

    pub fn nutation_longitude(&self) -> f64 {
        self.core_result.delta_psi
    }

    pub fn nutation_obliquity(&self) -> f64 {
        self.core_result.delta_eps
    }

    pub fn model(&self) -> NutationModel {
        self.model
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NutationModel {
    IAU2000A,
    IAU2000B,
    IAU2006A,
}

pub trait NutationCalculator {
    fn nutation_iau2000a(&self) -> TimeResult<NutationResult>;

    fn nutation_iau2000b(&self) -> TimeResult<NutationResult>;

    fn nutation_iau2006a(&self) -> TimeResult<NutationResult>;

    fn nutation(&self) -> TimeResult<NutationResult> {
        self.nutation_iau2006a()
    }
}

impl NutationCalculator for TT {
    fn nutation_iau2000a(&self) -> TimeResult<NutationResult> {
        iau2000a::calculate(self)
    }

    fn nutation_iau2000b(&self) -> TimeResult<NutationResult> {
        iau2000b::calculate(self)
    }

    fn nutation_iau2006a(&self) -> TimeResult<NutationResult> {
        iau2006a::calculate(self)
    }
}

#[cfg(test)]
mod integration_tests {
    use super::*;
    use crate::scales::tt::TT;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_nutation_models_match_erfa() {
        // eraNut00a, eraNut00b and eraNut06a
        let models = [
            TT::nutation_iau2000a,
            TT::nutation_iau2000b,
            TT::nutation_iau2006a,
        ];
        let cases = [
            (
                0.0,
                [
                    (-6.754422426417298e-5, -2.7970831192374137e-5),
                    (-6.754261253992235e-5, -2.7970923310985653e-5),
                    (-6.754425598969512e-5, -2.7970831192374137e-5),
                ],
            ),
            (
                9131.987654321,
                [
                    (1.2842967635988853e-6, 4.1352977363819786e-5),
                    (1.2870930153185514e-6, 4.135339367587379e-5),
                    (1.2842964750095786e-6, 4.135294864806038e-5),
                ],
            ),
        ];
        for (jd2, expected) in cases {
            let tt = TT::from_julian_date(crate::julian::JulianDate::new(J2000_JD, jd2));
            for (model, (dpsi, deps)) in models.iter().zip(expected) {
                let result = model(&tt).unwrap();
                let got = (result.nutation_longitude(), result.nutation_obliquity());
                assert_eq!(got, (dpsi, deps), "{jd2}, {:?}", result.model());
            }
        }
    }

    #[test]
    fn test_nutation_trait_methods() {
        let j2000_tt = TT::j2000();

        let default_nutation = j2000_tt.nutation().unwrap();

        assert_eq!(default_nutation.model, NutationModel::IAU2006A);
    }

    #[test]
    fn test_nutation_result_model_getter() {
        let j2000_tt = TT::j2000();
        let result_2000a = j2000_tt.nutation_iau2000a().unwrap();
        let result_2000b = j2000_tt.nutation_iau2000b().unwrap();
        let result_2006a = j2000_tt.nutation_iau2006a().unwrap();

        assert_eq!(result_2000a.model(), NutationModel::IAU2000A);
        assert_eq!(result_2000b.model(), NutationModel::IAU2000B);
        assert_eq!(result_2006a.model(), NutationModel::IAU2006A);
    }

    #[test]
    fn test_nutation_epoch_too_far_from_j2000() {
        let days = 25.0 * celestial_core::constants::DAYS_PER_JULIAN_CENTURY;
        let tt = TT::from_julian_date(crate::julian::JulianDate::new(J2000_JD, days));
        let expected = crate::TimeError::InvalidEpoch(
            "Epoch too far from J2000.0 for the IAU models: 25.0 centuries".into(),
        );
        assert_eq!(tt.nutation_iau2000a().unwrap_err(), expected);
        assert_eq!(tt.nutation_iau2000b().unwrap_err(), expected);
        assert_eq!(tt.nutation_iau2006a().unwrap_err(), expected);
    }
}
