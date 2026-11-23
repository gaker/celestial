mod bitbuffer;
mod compress;
mod decompress;

use crate::fits::errors::Result;

// Rice coding parameters for each pixel width. FSBITS bits hold the code of
// a block's split level, and a split level of FSMAX marks a block stored as
// raw BBITS-bit differences. Differences wrap modulo 2^BBITS, so unsigned
// pixels of the same width encode through the signed type.
pub(crate) trait RiceCompressible: Copy {
    const FSBITS: usize;
    const FSMAX: usize;
    const BBITS: usize;
    const BYTES_PER_ELEMENT: usize;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>>;
    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>>;

    fn as_signed(self) -> i32;
    fn from_signed(val: i32) -> Self;

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32 {
        Self::from_signed(next_pix.wrapping_sub(last_pix)).as_signed()
    }
}

macro_rules! rice_compressible {
    ($t:ty, $fsbits:expr, $fsmax:expr) => {
        impl RiceCompressible for $t {
            const FSBITS: usize = $fsbits;
            const FSMAX: usize = $fsmax;
            const BBITS: usize = <$t>::BITS as usize;
            const BYTES_PER_ELEMENT: usize = std::mem::size_of::<$t>();

            fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>> {
                compress::compress_rice(data, block_size)
            }

            fn decompress(
                compressed: &[u8],
                output_len: usize,
                block_size: usize,
            ) -> Result<Vec<Self>> {
                decompress::decompress_rice(compressed, output_len, block_size)
            }

            fn as_signed(self) -> i32 {
                i32::from(self)
            }

            // Truncation keeps the low bits, which wraps the value modulo
            // 2^BBITS.
            fn from_signed(val: i32) -> Self {
                val as Self
            }
        }
    };
}

rice_compressible!(i8, 3, 6);
rice_compressible!(i16, 4, 14);
rice_compressible!(i32, 5, 25);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_i32_constants() {
        assert_eq!(i32::FSBITS, 5);
        assert_eq!(i32::FSMAX, 25);
        assert_eq!(i32::BBITS, 32);
        assert_eq!(i32::BYTES_PER_ELEMENT, 4);
    }

    #[test]
    fn test_i16_constants() {
        assert_eq!(i16::FSBITS, 4);
        assert_eq!(i16::FSMAX, 14);
        assert_eq!(i16::BBITS, 16);
        assert_eq!(i16::BYTES_PER_ELEMENT, 2);
    }

    #[test]
    fn test_i8_constants() {
        assert_eq!(i8::FSBITS, 3);
        assert_eq!(i8::FSMAX, 6);
        assert_eq!(i8::BBITS, 8);
        assert_eq!(i8::BYTES_PER_ELEMENT, 1);
    }

    #[test]
    fn test_i32_as_signed() {
        let val: i32 = -42;
        assert_eq!(val.as_signed(), -42);
    }

    #[test]
    fn test_i32_from_signed() {
        assert_eq!(i32::from_signed(-42), -42);
    }

    #[test]
    fn test_i16_as_signed() {
        let val: i16 = -42;
        assert_eq!(val.as_signed(), -42);
    }

    #[test]
    fn test_i16_from_signed() {
        assert_eq!(i16::from_signed(-42), -42);
    }

    #[test]
    fn test_i8_as_signed() {
        let val: i8 = -42;
        assert_eq!(val.as_signed(), -42);
    }

    #[test]
    fn test_i8_from_signed() {
        assert_eq!(i8::from_signed(-42), -42);
    }

    #[test]
    fn test_i32_compute_difference() {
        assert_eq!(i32::compute_difference(10, 5), 5);
        assert_eq!(i32::compute_difference(5, 10), -5);
        assert_eq!(i32::compute_difference(i32::MAX, i32::MIN), -1);
    }

    #[test]
    fn test_i16_compute_difference() {
        assert_eq!(i16::compute_difference(10, 5), 5);
        assert_eq!(i16::compute_difference(5, 10), -5);
        assert_eq!(i16::compute_difference(32767, -32768), -1);
    }

    #[test]
    fn test_i8_compute_difference() {
        assert_eq!(i8::compute_difference(10, 5), 5);
        assert_eq!(i8::compute_difference(5, 10), -5);
        assert_eq!(i8::compute_difference(127, -128), -1);
    }

    #[test]
    fn test_i32_round_trip_compression() {
        let data = vec![1i32, 2, 3, 4, 5];
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_i16_round_trip_compression() {
        let data = vec![1i16, 2, 3, 4, 5];
        let compressed = i16::compress(&data, 32).unwrap();
        let decompressed = i16::decompress(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_i8_round_trip_compression() {
        let data = vec![1i8, 2, 3, 4, 5];
        let compressed = i8::compress(&data, 32).unwrap();
        let decompressed = i8::decompress(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_with_negative_values() {
        let data = vec![-10i32, -5, 0, 5, 10];
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_with_extreme_values() {
        let data = vec![i32::MIN, -1, 0, 1, i32::MAX];
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_empty() {
        let data: Vec<i32> = Vec::new();
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, 0, 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_single_element() {
        let data = vec![42i32];
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, 1, 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_repeated_values() {
        let data = vec![7i32; 100];
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, 100, 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_small_block_size() {
        let data = vec![1i32, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let compressed = i32::compress(&data, 3).unwrap();
        let decompressed = i32::decompress(&compressed, 10, 3).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_ascending_sequence() {
        let data: Vec<i32> = (0..1000).collect();
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, 1000, 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_round_trip_descending_sequence() {
        let data: Vec<i32> = (0..1000).rev().collect();
        let compressed = i32::compress(&data, 32).unwrap();
        let decompressed = i32::decompress(&compressed, 1000, 32).unwrap();
        assert_eq!(data, decompressed);
    }

    // One large difference among zeros sets a low split level, so its unary
    // prefix runs past 32 bits.
    #[test]
    fn round_trip_with_a_lone_large_difference() {
        let mut i32s = vec![0i32; 64];
        i32s[31] = 1000;
        let stream = i32::compress(&i32s, 32).unwrap();
        assert_eq!(i32::decompress(&stream, 64, 32).unwrap(), i32s);

        let mut i16s = vec![0i16; 64];
        i16s[31] = 1000;
        let stream = i16::compress(&i16s, 32).unwrap();
        assert_eq!(i16::decompress(&stream, 64, 32).unwrap(), i16s);

        let i8s: Vec<i8> = (0..64).map(|i| if i < 37 { 0 } else { -56 }).collect();
        let stream = i8::compress(&i8s, 32).unwrap();
        assert_eq!(i8::decompress(&stream, 64, 32).unwrap(), i8s);
    }
}
