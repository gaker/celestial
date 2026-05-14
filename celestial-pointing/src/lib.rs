pub mod error;
pub mod observation;
pub mod commands;
pub(crate) mod diurnal;
pub mod model;
pub mod parser;
pub mod plot;
pub(crate) mod prepare;
pub mod session;
pub mod solver;
pub mod terms;
pub mod writer;

#[cfg(test)]
mod test_support;

pub use error::{Error, Result};
