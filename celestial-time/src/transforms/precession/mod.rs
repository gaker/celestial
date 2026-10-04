mod iau2000;
mod iau2006;

use crate::scales::tt::TT;
use crate::TimeResult;
use celestial_core::matrix::RotationMatrix3;

#[derive(Debug, Clone, PartialEq)]
pub struct PrecessionResult {
    pub bias_matrix: RotationMatrix3,
    pub precession_matrix: RotationMatrix3,
    pub bias_precession_matrix: RotationMatrix3,
    pub model: PrecessionModel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrecessionModel {
    IAU2000,
    IAU2006,
}

pub trait PrecessionCalculator {
    fn precession_iau2000(&self) -> TimeResult<PrecessionResult>;

    fn precession_iau2006(&self) -> TimeResult<PrecessionResult>;

    fn precession(&self) -> TimeResult<PrecessionResult> {
        self.precession_iau2006()
    }

    fn bias_precession_matrix(&self) -> TimeResult<RotationMatrix3> {
        Ok(self.precession_iau2006()?.bias_precession_matrix)
    }
}

impl PrecessionCalculator for TT {
    fn precession_iau2000(&self) -> TimeResult<PrecessionResult> {
        iau2000::calculate(self)
    }

    fn precession_iau2006(&self) -> TimeResult<PrecessionResult> {
        iau2006::calculate(self)
    }
}
