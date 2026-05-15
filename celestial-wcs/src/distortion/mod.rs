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
        Self::apply(self, x, y)
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        Self::apply_inverse(self, x, y)
    }

    fn operates_on_pixels(&self) -> bool {
        true
    }
}

impl Distortion for TpvDistortion {
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        Self::apply(self, x, y)
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        Self::apply_inverse(self, x, y)
    }

    fn operates_on_pixels(&self) -> bool {
        false
    }
}

impl Distortion for TnxDistortion {
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        Self::apply(self, x, y)
    }

    fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        Self::apply_inverse(self, x, y)
    }

    fn operates_on_pixels(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distortion_model_dispatches_to_each_variant() {
        // SIP is the pixel-space distortion; the empty (degree-0 effective)
        // model leaves coordinates unchanged.
        let sip = SipDistortion::new([512.0, 512.0], 2, 2);
        let sip_model = DistortionModel::Sip(sip.clone());
        assert!(sip_model.operates_on_pixels());
        assert_eq!(sip_model.apply(100.0, 200.0), (100.0, 200.0));
        let sip_trait: &dyn Distortion = &sip;
        assert!(sip_trait.operates_on_pixels());

        // TPV is the intermediate-space distortion; identity leaves
        // intermediates unchanged.
        let tpv = TpvDistortion::identity();
        let tpv_model = DistortionModel::Tpv(Box::new(tpv.clone()));
        assert!(!tpv_model.operates_on_pixels());
        assert_eq!(tpv_model.apply(0.5, 0.5), (0.5, 0.5));
        let tpv_trait: &dyn Distortion = &tpv;
        assert!(!tpv_trait.operates_on_pixels());
    }
}
