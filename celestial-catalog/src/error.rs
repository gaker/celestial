use celestial_coords::errors::CoordError;
use thiserror::Error;

pub type CatalogResult<T> = Result<T, CatalogError>;

#[derive(Debug, Error)]
pub enum CatalogError {
    // The I/O error is the source, so reporters such as anyhow print it once, after the
    // context.
    #[error("{context}")]
    Io {
        context: String,
        source: std::io::Error,
    },

    #[error("Invalid catalog header: {message}")]
    InvalidHeader { message: String },

    #[error("Corrupt catalog pixel table: {message}")]
    CorruptTable { message: String },

    #[error("Invalid input: {message}")]
    InvalidInput { message: String },

    #[error("Coordinate calculation failed: {0}")]
    Coord(#[from] CoordError),
}

impl CatalogError {
    pub(crate) fn invalid_header(message: impl Into<String>) -> Self {
        Self::InvalidHeader {
            message: message.into(),
        }
    }

    pub(crate) fn corrupt_table(message: impl Into<String>) -> Self {
        Self::CorruptTable {
            message: message.into(),
        }
    }

    pub(crate) fn invalid_input(message: impl Into<String>) -> Self {
        Self::InvalidInput {
            message: message.into(),
        }
    }
}
