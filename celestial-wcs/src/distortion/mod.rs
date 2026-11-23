mod polynomial;
pub mod sip;
pub mod tpv;

pub use sip::SipDistortion;
pub use tpv::TpvDistortion;

#[derive(Debug, Clone)]
pub enum DistortionModel {
    Sip(SipDistortion),
    Tpv(Box<TpvDistortion>),
}
