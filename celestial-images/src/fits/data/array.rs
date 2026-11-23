use crate::core::types::{BitPix, ByteOrder};
use crate::fits::errors::{FitsError, Result};

// Callers can name `DataArray` in bounds but cannot implement it. The
// readers and writers support only the sample types listed here.
mod sealed {
    pub trait Sealed {}

    impl Sealed for u8 {}
    impl Sealed for i16 {}
    impl Sealed for i32 {}
    impl Sealed for i64 {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

pub trait DataArray: sealed::Sealed + Sized + PartialEq {
    const BITPIX: BitPix;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>>;
    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>>;
}

// Maps each `N`-byte sample of `bytes` through `sample`. A length that is
// not a whole number of samples is an error.
pub(crate) fn decode<const N: usize, T>(
    bytes: &[u8],
    sample: impl Fn([u8; N]) -> T,
) -> Result<Vec<T>> {
    let (samples, rest) = bytes.as_chunks::<N>();
    if !rest.is_empty() {
        return Err(FitsError::InvalidFormat(format!(
            "the data is {} bytes, not a whole number of {N}-byte samples",
            bytes.len()
        )));
    }
    Ok(samples.iter().map(|&s| sample(s)).collect())
}
