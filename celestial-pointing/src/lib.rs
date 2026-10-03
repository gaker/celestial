pub mod commands;
pub(crate) mod diurnal;
pub mod error;
pub mod model;
pub mod observation;
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
