mod bitbuffer;
mod compress;
mod decompress;

use crate::fits::Result;

pub trait RiceCompressible: Copy {
    const FSBITS: usize;
    const FSMAX: usize;
    const BBITS: usize;
    const BYTES_PER_ELEMENT: usize;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>>;
    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>>;

    fn as_signed(self) -> i32;
    fn from_signed(val: i32) -> Self;

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32;
}

impl RiceCompressible for i32 {
    const FSBITS: usize = 5;
    const FSMAX: usize = 25;
    const BBITS: usize = 32;
    const BYTES_PER_ELEMENT: usize = 4;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>> {
        compress::compress_rice(data, block_size)
    }

    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>> {
        decompress::decompress_rice(compressed, output_len, block_size)
    }

    fn as_signed(self) -> i32 {
        self
    }

    fn from_signed(val: i32) -> Self {
        val
    }

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32 {
        next_pix.wrapping_sub(last_pix)
    }
}

impl RiceCompressible for i16 {
    const FSBITS: usize = 4;
    const FSMAX: usize = 14;
    const BBITS: usize = 16;
    const BYTES_PER_ELEMENT: usize = 2;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>> {
        compress::compress_rice(data, block_size)
    }

    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>> {
        decompress::decompress_rice(compressed, output_len, block_size)
    }

    fn as_signed(self) -> i32 {
        self as i32
    }

    fn from_signed(val: i32) -> Self {
        val as i16
    }

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32 {
        let diff = (next_pix as i16).wrapping_sub(last_pix as i16);
        diff as i32
    }
}

impl RiceCompressible for u8 {
    const FSBITS: usize = 3;
    const FSMAX: usize = 6;
    const BBITS: usize = 8;
    const BYTES_PER_ELEMENT: usize = 1;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>> {
        let signed_data: Vec<i8> = data.iter().map(|&x| x as i8).collect();
        compress::compress_rice(&signed_data, block_size)
    }

    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>> {
        let signed_result: Vec<i8> =
            decompress::decompress_rice(compressed, output_len, block_size)?;
        Ok(signed_result.into_iter().map(|x| x as u8).collect())
    }

    fn as_signed(self) -> i32 {
        self as i8 as i32
    }

    fn from_signed(val: i32) -> Self {
        val as u8
    }

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32 {
        let diff = (next_pix as i8).wrapping_sub(last_pix as i8);
        diff as i32
    }
}

impl RiceCompressible for u16 {
    const FSBITS: usize = 4;
    const FSMAX: usize = 14;
    const BBITS: usize = 16;
    const BYTES_PER_ELEMENT: usize = 2;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>> {
        let signed_data: Vec<i16> = data.iter().map(|&x| x as i16).collect();
        compress::compress_rice(&signed_data, block_size)
    }

    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>> {
        let signed_result: Vec<i16> =
            decompress::decompress_rice(compressed, output_len, block_size)?;
        Ok(signed_result.into_iter().map(|x| x as u16).collect())
    }

    fn as_signed(self) -> i32 {
        self as i16 as i32
    }

    fn from_signed(val: i32) -> Self {
        val as u16
    }

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32 {
        let diff = (next_pix as i16).wrapping_sub(last_pix as i16);
        diff as i32
    }
}

impl RiceCompressible for i8 {
    const FSBITS: usize = 3;
    const FSMAX: usize = 6;
    const BBITS: usize = 8;
    const BYTES_PER_ELEMENT: usize = 1;

    fn compress(data: &[Self], block_size: usize) -> Result<Vec<u8>> {
        compress::compress_rice(data, block_size)
    }

    fn decompress(compressed: &[u8], output_len: usize, block_size: usize) -> Result<Vec<Self>> {
        decompress::decompress_rice(compressed, output_len, block_size)
    }

    fn as_signed(self) -> i32 {
        self as i32
    }

    fn from_signed(val: i32) -> Self {
        val as i8
    }

    fn compute_difference(next_pix: i32, last_pix: i32) -> i32 {
        let diff = (next_pix as i8).wrapping_sub(last_pix as i8);
        diff as i32
    }
}

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
    fn test_u8_constants() {
        assert_eq!(u8::FSBITS, 3);
        assert_eq!(u8::FSMAX, 6);
        assert_eq!(u8::BBITS, 8);
        assert_eq!(u8::BYTES_PER_ELEMENT, 1);
    }

    #[test]
    fn test_u16_constants() {
        assert_eq!(u16::FSBITS, 4);
        assert_eq!(u16::FSMAX, 14);
        assert_eq!(u16::BBITS, 16);
        assert_eq!(u16::BYTES_PER_ELEMENT, 2);
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
    fn test_u8_as_signed() {
        let val: u8 = 200;
        assert_eq!(val.as_signed(), -56);
    }

    #[test]
    fn test_u8_from_signed() {
        assert_eq!(u8::from_signed(-56), 200);
    }

    #[test]
    fn test_u16_as_signed() {
        let val: u16 = 40000;
        assert_eq!(val.as_signed(), -25536);
    }

    #[test]
    fn test_u16_from_signed() {
        assert_eq!(u16::from_signed(-25536), 40000);
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
    fn test_u8_compute_difference() {
        assert_eq!(u8::compute_difference(10, 5), 5);
        assert_eq!(u8::compute_difference(5, 10), -5);
        assert_eq!(u8::compute_difference(127, -128), -1);
    }

    #[test]
    fn test_u16_compute_difference() {
        assert_eq!(u16::compute_difference(10, 5), 5);
        assert_eq!(u16::compute_difference(5, 10), -5);
        assert_eq!(u16::compute_difference(32767, -32768), -1);
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
    fn test_u8_round_trip_compression() {
        let data = vec![1u8, 2, 3, 4, 5];
        let compressed = u8::compress(&data, 32).unwrap();
        let decompressed = u8::decompress(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_u16_round_trip_compression() {
        let data = vec![1u16, 2, 3, 4, 5];
        let compressed = u16::compress(&data, 32).unwrap();
        let decompressed = u16::decompress(&compressed, data.len(), 32).unwrap();
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
}
