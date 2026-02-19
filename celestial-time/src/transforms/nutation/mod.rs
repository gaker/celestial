pub mod iau2000a;
pub mod iau2000b;
pub mod iau2006a;

use crate::{TimeResult, TT};

#[derive(Debug)]
pub struct NutationResult {
    core_result: celestial_core::nutation::NutationResult,
    model: NutationModel,
}

impl NutationResult {
    pub fn new(
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
    use crate::TT;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_all_nutation_models_at_j2000() {
        let j2000_tt = TT::j2000();

        let result_2000a = j2000_tt.nutation_iau2000a().unwrap();
        let result_2000b = j2000_tt.nutation_iau2000b().unwrap();
        let result_2006a = j2000_tt.nutation_iau2006a().unwrap();

        assert!(
            result_2000a.nutation_longitude().abs() < 1e-3,
            "2000A nutation too large"
        );
        assert!(
            result_2000b.nutation_longitude().abs() < 1e-3,
            "2000B nutation too large"
        );
        assert!(
            result_2006a.nutation_longitude().abs() < 1e-3,
            "2006A nutation too large"
        );

        assert_eq!(result_2000a.model, NutationModel::IAU2000A);
        assert_eq!(result_2000b.model, NutationModel::IAU2000B);
        assert_eq!(result_2006a.model, NutationModel::IAU2006A);
    }

    #[test]
    fn test_iau2000b_is_abbreviated_version() {
        let j2000_tt = TT::j2000();

        let result_2000a = j2000_tt.nutation_iau2000a().unwrap();
        let result_2000b = j2000_tt.nutation_iau2000b().unwrap();

        let diff_psi =
            (result_2000a.nutation_longitude() - result_2000b.nutation_longitude()).abs();
        let diff_eps =
            (result_2000a.nutation_obliquity() - result_2000b.nutation_obliquity()).abs();

        assert!(
            diff_psi < 5e-9,
            "2000B differs too much from 2000A in longitude: {:.3} mas",
            diff_psi * 206264806.247
        );
        assert!(
            diff_eps < 5e-9,
            "2000B differs too much from 2000A in obliquity: {:.3} mas",
            diff_eps * 206264806.247
        );
    }

    #[test]
    fn test_iau2006a_corrections_reasonable() {
        let j2000_tt = TT::j2000();

        let result_2000a = j2000_tt.nutation_iau2000a().unwrap();
        let result_2006a = j2000_tt.nutation_iau2006a().unwrap();

        let diff_psi =
            (result_2000a.nutation_longitude() - result_2006a.nutation_longitude()).abs();
        let diff_eps =
            (result_2000a.nutation_obliquity() - result_2006a.nutation_obliquity()).abs();

        assert!(
            diff_psi < 1e-8,
            "2006A corrections too large relative to 2000A"
        );
        assert!(
            diff_eps < 1e-8,
            "2006A corrections too large relative to 2000A"
        );
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
        use crate::JulianDate;

        let far_future_jd = J2000_JD + (25.0 * celestial_core::constants::DAYS_PER_JULIAN_CENTURY);
        let far_future_tt = TT::from_julian_date(JulianDate::from_f64(far_future_jd));

        let result = far_future_tt.nutation_iau2006a();
        assert!(result.is_err());

        if let Err(crate::TimeError::InvalidEpoch(msg)) = result {
            assert!(msg.contains("Epoch too far from J2000.0"));
        } else {
            panic!("Expected InvalidEpoch error");
        }
    }
}

mod utils {
    use crate::{TimeError, TimeResult, TT};

    pub fn tt_to_centuries(tt: &TT) -> TimeResult<f64> {
        let jd = tt.to_julian_date();
        let centuries = celestial_core::utils::jd_to_centuries(jd.jd1(), jd.jd2());

        if centuries.abs() > 20.0 {
            return Err(TimeError::InvalidEpoch(format!(
                "Epoch too far from J2000.0 for nutation model: {:.1} centuries",
                centuries
            )));
        }

        Ok(centuries)
    }
}

pub(crate) use utils::tt_to_centuries;
