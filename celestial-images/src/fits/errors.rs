#[derive(Debug, thiserror::Error)]
pub enum FitsError {
    #[error("Invalid FITS format: {0}")]
    InvalidFormat(String),

    #[error("Keyword {keyword} not found")]
    KeywordNotFound { keyword: String },

    #[error("Data type mismatch: expected {expected:?}, got {actual:?}")]
    TypeMismatch {
        expected: crate::core::types::BitPix,
        actual: crate::core::types::BitPix,
    },

    #[error("Invalid BITPIX value: {0}")]
    InvalidBitPix(i32),

    #[error("Invalid keyword: {0}")]
    InvalidKeyword(String),

    #[error("Invalid keyword value: {keyword} = {value}")]
    InvalidKeywordValue { keyword: String, value: String },

    #[error("Unsupported compression: {0}")]
    UnsupportedCompression(String),

    #[error("HDU not found: {0}")]
    HduNotFound(usize),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, FitsError>;
