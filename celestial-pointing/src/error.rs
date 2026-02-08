use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("parse error: {0}")]
    Parse(String),

    #[error("fit error: {0}")]
    Fit(String),

    #[error("unknown term: {0}")]
    UnknownTerm(String),

    #[error("invalid harmonic specification: {0}")]
    InvalidHarmonic(String),

    #[error("no LST set - use LST command to set local sidereal time")]
    NoLst,

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, Error>;
