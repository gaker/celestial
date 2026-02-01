pub mod polynomial;
pub mod sip;
pub mod tnx;
pub mod tpv;

pub use sip::SipDistortion;
pub use tnx::{CrossTerms, SurfaceType, TnxDistortion, TnxSurface};
pub use tpv::TpvDistortion;

use crate::error::WcsResult;

pub trait Distortion {
    fn apply(&self, x: f64, y: f64) -> (f64, f64);
    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)>;
    fn operates_on_pixels(&self) -> bool;
}

#[derive(Debug, Clone)]
pub enum DistortionModel {
    Sip(SipDistortion),
    Tpv(Box<TpvDistortion>),
    Tnx(TnxDistortion),
}

impl Distortion for DistortionModel {
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        match self {
            Self::Sip(d) => d.apply(x, y),
            Self::Tpv(d) => d.as_ref().apply(x, y),
            Self::Tnx(d) => d.apply(x, y),
        }
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        match self {
            Self::Sip(d) => d.apply_inverse(x, y),
            Self::Tpv(d) => d.as_ref().apply_inverse(x, y),
            Self::Tnx(d) => d.apply_inverse(x, y),
        }
    }

    fn operates_on_pixels(&self) -> bool {
        match self {
            Self::Sip(_) => true,
            Self::Tpv(_) => false,
            Self::Tnx(_) => false,
        }
    }
}

impl Distortion for SipDistortion {
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        SipDistortion::apply(self, x, y)
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        SipDistortion::apply_inverse(self, x, y)
    }

    fn operates_on_pixels(&self) -> bool {
        true
    }
}

impl Distortion for TpvDistortion {
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        TpvDistortion::apply(self, x, y)
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        TpvDistortion::apply_inverse(self, x, y)
    }

    fn operates_on_pixels(&self) -> bool {
        false
    }
}

impl Distortion for TnxDistortion {
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        TnxDistortion::apply(self, x, y)
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        TnxDistortion::apply_inverse(self, x, y)
    }

    fn operates_on_pixels(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distortion_model_sip() {
        let sip = SipDistortion::new([512.0, 512.0], 2, 2);
        let model = DistortionModel::Sip(sip);
        assert!(model.operates_on_pixels());
        let (x, y) = model.apply(100.0, 200.0);
        assert_eq!((x, y), (100.0, 200.0));
    }

    #[test]
    fn test_distortion_model_tpv() {
        let tpv = TpvDistortion::identity();
        let model = DistortionModel::Tpv(Box::new(tpv));
        assert!(!model.operates_on_pixels());
        let (x, y) = model.apply(0.5, 0.5);
        assert_eq!((x, y), (0.5, 0.5));
    }

    #[test]
    fn test_distortion_trait_sip() {
        let sip = SipDistortion::new([512.0, 512.0], 2, 2);
        let d: &dyn Distortion = &sip;
        assert!(d.operates_on_pixels());
    }

    #[test]
    fn test_distortion_trait_tpv() {
        let tpv = TpvDistortion::identity();
        let d: &dyn Distortion = &tpv;
        assert!(!d.operates_on_pixels());
    }
}
