#[derive(Debug, thiserror::Error)]
pub enum FitsError {
    #[error("Invalid FITS format: {0}")]
    InvalidFormat(String),

    #[error("Keyword {keyword} not found")]
    KeywordNotFound { keyword: String },

    #[error("Data type mismatch: expected {expected:?}, got {actual:?}")]
    TypeMismatch {
        expected: crate::core::BitPix,
        actual: crate::core::BitPix,
    },

    #[error("Invalid BITPIX value: {0}")]
    InvalidBitPix(i32),

    #[error("Header parsing error: {0}")]
    HeaderParse(String),

    #[error("Invalid keyword: {0}")]
    InvalidKeyword(String),

    #[error("Invalid keyword value: {keyword} = {value}")]
    InvalidKeywordValue { keyword: String, value: String },

    #[error("Checksum verification failed: {0}")]
    ChecksumMismatch(String),

    #[error("DATASUM verification failed: expected {expected}, computed {computed}")]
    DatasumMismatch { expected: String, computed: u32 },

    #[error("Unsupported compression: {0}")]
    UnsupportedCompression(String),

    #[error("HDU not found: {0}")]
    HduNotFound(usize),

    #[error("EOF reached unexpectedly")]
    UnexpectedEof,

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, FitsError>;
