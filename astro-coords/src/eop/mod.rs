pub mod cache;
pub mod embedded;
pub mod interpolate;
pub mod manager;
pub mod parser;
pub mod record;

#[cfg(feature = "eop-download")]
pub mod download;

pub use manager::{EopBuilder, EopManager};
pub use record::{EopParameters, EopRecord};

pub use crate::CoordResult as EopResult;
