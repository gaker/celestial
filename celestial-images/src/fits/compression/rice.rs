use crate::fits::data::array;
use crate::fits::errors::{FitsError, Result};
use crate::ricecomp::RiceCompressible;

// FITS 4.0 Table 37 defaults, used when ZNAMEi = 'BLOCKSIZE' or 'BYTEPIX' is
// absent.
pub(super) const DEFAULT_BLOCK_SIZE: usize = 32;
const DEFAULT_BYTEPIX: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct RiceParams {
    pub(super) block_size: usize,
    pub(super) bytepix: usize,
}

impl Default for RiceParams {
    fn default() -> Self {
        Self {
            block_size: DEFAULT_BLOCK_SIZE,
            bytepix: DEFAULT_BYTEPIX,
        }
    }
}

// BYTEPIX 1 values are unsigned, the same as BITPIX 8 pixels. They decode
// through i8, whose differences wrap the same way.
pub(super) fn decode(bytes: &[u8], pixels: usize, params: RiceParams) -> Result<Vec<i64>> {
    let block = params.block_size;
    match params.bytepix {
        1 => Ok(i8::decompress(bytes, pixels, block)?
            .into_iter()
            .map(|v| i64::from(v as u8))
            .collect()),
        2 => Ok(widen(i16::decompress(bytes, pixels, block)?)),
        4 => Ok(widen(i32::decompress(bytes, pixels, block)?)),
        n => Err(FitsError::UnsupportedCompression(format!(
            "RICE_1 tiles with BYTEPIX = {n}"
        ))),
    }
}

fn widen<T: Into<i64>>(values: Vec<T>) -> Vec<i64> {
    values.into_iter().map(Into::into).collect()
}

// Every block costs at least its split-level code (3, 4 or 5 bits), so `len`
// bytes cannot encode more pixels than this.
pub(super) fn max_pixels(len: usize, params: RiceParams) -> usize {
    let code_bits = match params.bytepix {
        1 => 3,
        2 => 4,
        _ => 5,
    };
    (len.saturating_mul(8) / code_bits)
        .saturating_add(1)
        .saturating_mul(params.block_size)
}

// Compresses one tile of big-endian pixels `bytepix` bytes wide, with the
// block size that the writer records in ZVAL1.
pub(super) fn compress(data: &[u8], bytepix: usize) -> Result<Vec<u8>> {
    match bytepix {
        1 => i8::compress(&array::decode(data, i8::from_be_bytes)?, DEFAULT_BLOCK_SIZE),
        2 => i16::compress(
            &array::decode(data, i16::from_be_bytes)?,
            DEFAULT_BLOCK_SIZE,
        ),
        4 => i32::compress(
            &array::decode(data, i32::from_be_bytes)?,
            DEFAULT_BLOCK_SIZE,
        ),
        n => Err(FitsError::UnsupportedCompression(format!(
            "RICE_1 tiles with BYTEPIX = {n}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(block_size: usize, bytepix: usize) -> RiceParams {
        RiceParams {
            block_size,
            bytepix,
        }
    }

    #[test]
    fn decode_reads_each_bytepix() {
        let i32s = [100i32, 101, 102, 103, 104, 105, 106, 107, 108];
        let stream = i32::compress(&i32s, 32).unwrap();
        let expected: Vec<i64> = i32s.iter().map(|&v| v.into()).collect();
        assert_eq!(decode(&stream, 9, params(32, 4)).unwrap(), expected);

        let i16s = [1000i16, 1001, -2, 1003];
        let stream = i16::compress(&i16s, 32).unwrap();
        assert_eq!(
            decode(&stream, 4, params(32, 2)).unwrap(),
            vec![1000, 1001, -2, 1003]
        );

        let bytes = [0u8, 200, 255, 3];
        let stream = compress(&bytes, 1).unwrap();
        assert_eq!(
            decode(&stream, 4, params(32, 1)).unwrap(),
            vec![0, 200, 255, 3]
        );
    }

    #[test]
    fn decode_uses_the_block_size() {
        let values: Vec<i16> = (0..40).map(|i| 1000 + (i % 5)).collect();
        let stream = i16::compress(&values, 16).unwrap();
        let expected: Vec<i64> = values.iter().map(|&v| v.into()).collect();
        assert_eq!(
            decode(&stream, 40, params(16, 2)).unwrap(),
            expected.clone()
        );
        assert_ne!(decode(&stream, 40, params(32, 2)).ok(), Some(expected));
    }

    #[test]
    fn decode_rejects_bytepix_8() {
        assert!(matches!(
            decode(&[0; 16], 2, params(32, 8)),
            Err(FitsError::UnsupportedCompression(_))
        ));
    }

    #[test]
    fn max_pixels_bounds_a_constant_stream() {
        // A constant i32 stream: 4 first-value bytes, then 5 bits per block.
        let values = vec![7i32; 32 * 100];
        let stream = i32::compress(&values, 32).unwrap();
        assert!(max_pixels(stream.len(), params(32, 4)) >= values.len());
    }

    #[test]
    fn compress_roundtrips_each_bitpix() {
        let i32s: Vec<u8> = [5i32, 6, 7].iter().flat_map(|v| v.to_be_bytes()).collect();
        let stream = compress(&i32s, 4).unwrap();
        assert_eq!(decode(&stream, 3, params(32, 4)).unwrap(), vec![5, 6, 7]);

        let i8s: Vec<u8> = (0..25).collect();
        let stream = compress(&i8s, 1).unwrap();
        let expected: Vec<i64> = (0..25).collect();
        assert_eq!(decode(&stream, 25, params(32, 1)).unwrap(), expected);
    }
}
