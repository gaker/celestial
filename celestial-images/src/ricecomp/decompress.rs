use super::bitbuffer::BitReader;
use super::RiceCompressible;
use crate::fits::Result;

pub fn decompress_rice<T: RiceCompressible>(
    compressed: &[u8],
    output_len: usize,
    block_size: usize,
) -> Result<Vec<T>> {
    if output_len == 0 {
        return Ok(Vec::new());
    }

    let first_val = read_first_value_direct::<T>(compressed)?;
    let mut last_pix = first_val.as_signed();

    let rice_data = &compressed[T::BYTES_PER_ELEMENT..];
    let mut reader = BitReader::new(rice_data);

    let mut result = Vec::with_capacity(output_len);

    let mut i = 0;
    while i < output_len {
        let fs = read_split_level::<T>(&mut reader)?;

        let this_block = (output_len - i).min(block_size);

        process_block::<T>(&mut reader, &mut result, &mut last_pix, this_block, fs)?;

        i += this_block;
    }

    Ok(result)
}

fn read_first_value_direct<T: RiceCompressible>(compressed: &[u8]) -> Result<T> {
    if compressed.len() < T::BYTES_PER_ELEMENT {
        return Err(crate::fits::FitsError::InvalidFormat(
            "Compressed data too short for first pixel".to_string(),
        ));
    }

    let mut value = 0u32;
    for &byte in compressed.iter().take(T::BYTES_PER_ELEMENT) {
        value = (value << 8) | (byte as u32);
    }

    let signed_val = value as i32;
    Ok(T::from_signed(signed_val))
}

fn process_block<T: RiceCompressible>(
    reader: &mut BitReader,
    result: &mut Vec<T>,
    last_pix: &mut i32,
    block_size: usize,
    fs: i32,
) -> Result<()> {
    if fs < 0 {
        for _ in 0..block_size {
            result.push(T::from_signed(*last_pix));
        }
    } else if fs == T::FSMAX as i32 {
        high_entropy_block::<T>(reader, result, last_pix, block_size)?;
    } else {
        rice_block::<T>(reader, result, last_pix, block_size, fs as usize)?;
    }

    Ok(())
}

fn read_split_level<T: RiceCompressible>(reader: &mut BitReader) -> Result<i32> {
    let fs_bits = reader.read_bits(T::FSBITS)?;
    let fs = fs_bits as i32 - 1;
    Ok(fs)
}

fn high_entropy_block<T: RiceCompressible>(
    reader: &mut BitReader,
    result: &mut Vec<T>,
    last_pix: &mut i32,
    block_size: usize,
) -> Result<()> {
    for _ in 0..block_size {
        let diff_bits = reader.read_bits(T::BBITS)?;

        let diff = if (diff_bits & 1) == 0 {
            (diff_bits >> 1) as i32
        } else {
            !((diff_bits >> 1) as i32)
        };

        let pixel_value = diff.wrapping_add(*last_pix);
        result.push(T::from_signed(pixel_value));
        *last_pix = pixel_value;
    }
    Ok(())
}

fn rice_block<T: RiceCompressible>(
    reader: &mut BitReader,
    result: &mut Vec<T>,
    last_pix: &mut i32,
    block_size: usize,
    fs: usize,
) -> Result<()> {
    for _ in 0..block_size {
        let leading_zeros = reader.count_leading_zeros()?;
        let bottom_bits = if fs > 0 { reader.read_bits(fs)? } else { 0 };

        let diff_mapped = (leading_zeros << fs) as u32 | bottom_bits;

        let diff = if (diff_mapped & 1) == 0 {
            (diff_mapped >> 1) as i32
        } else {
            !((diff_mapped >> 1) as i32)
        };

        let pixel_value = diff.wrapping_add(*last_pix);
        result.push(T::from_signed(pixel_value));
        *last_pix = pixel_value;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ricecomp::compress;

    #[test]
    fn test_decompress_rice_empty() {
        let result = decompress_rice::<i32>(&[], 0, 32).unwrap();
        assert_eq!(result, Vec::<i32>::new());
    }

    #[test]
    fn test_read_first_value_direct_i32() {
        let data = vec![0x12, 0x34, 0x56, 0x78];
        let result = read_first_value_direct::<i32>(&data).unwrap();
        assert_eq!(result, 0x12345678);
    }

    #[test]
    fn test_read_first_value_direct_i16() {
        let data = vec![0x12, 0x34];
        let result = read_first_value_direct::<i16>(&data).unwrap();
        assert_eq!(result, 0x1234);
    }

    #[test]
    fn test_read_first_value_direct_i8() {
        let data = vec![0x12];
        let result = read_first_value_direct::<i8>(&data).unwrap();
        assert_eq!(result, 0x12);
    }

    #[test]
    fn test_read_first_value_direct_insufficient_data() {
        let data = vec![0x12];
        let result = read_first_value_direct::<i32>(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_split_level() {
        let data = vec![0b01000000];
        let mut reader = BitReader::new(&data);
        let fs = read_split_level::<i32>(&mut reader).unwrap();
        assert_eq!(fs, 7);
    }

    #[test]
    fn test_read_split_level_zero() {
        let data = vec![0b00001000];
        let mut reader = BitReader::new(&data);
        let fs = read_split_level::<i32>(&mut reader).unwrap();
        assert_eq!(fs, 0);
    }

    #[test]
    fn test_read_split_level_negative() {
        let data = vec![0b00000000];
        let mut reader = BitReader::new(&data);
        let fs = read_split_level::<i32>(&mut reader).unwrap();
        assert_eq!(fs, -1);
    }

    #[test]
    fn test_process_block_low_entropy() {
        let data = vec![0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 42;

        process_block::<i32>(&mut reader, &mut result, &mut last_pix, 3, -1).unwrap();

        assert_eq!(result.len(), 3);
        assert_eq!(result, vec![42, 42, 42]);
    }

    #[test]
    fn test_high_entropy_block() {
        let data = vec![0x24, 0x68, 0xAC, 0xF0];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        high_entropy_block::<i32>(&mut reader, &mut result, &mut last_pix, 1).unwrap();

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_rice_block() {
        let data = vec![0x80, 0x00];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        rice_block::<i32>(&mut reader, &mut result, &mut last_pix, 1, 0).unwrap();

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_rice_block_with_fs() {
        let data = vec![0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        rice_block::<i32>(&mut reader, &mut result, &mut last_pix, 1, 2).unwrap();

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_decompress_rice_single_pixel() {
        let mut data = vec![0x00, 0x00, 0x00, 0x05];
        data.extend_from_slice(&[0b00000000]);

        let result = decompress_rice::<i32>(&data, 1, 32).unwrap();
        assert_eq!(result, vec![5]);
    }

    #[test]
    fn test_process_block_high_entropy() {
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        process_block::<i32>(&mut reader, &mut result, &mut last_pix, 1, 25).unwrap();

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_process_block_normal_rice() {
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        process_block::<i32>(&mut reader, &mut result, &mut last_pix, 1, 2).unwrap();

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_decompress_rice_multiple_blocks() {
        let mut data = vec![0x00, 0x00, 0x00, 0x01];
        data.extend_from_slice(&[0b00000000, 0b00000000]);

        let result = decompress_rice::<i32>(&data, 3, 2).unwrap();
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], 1);
        assert_eq!(result[1], 1);
        assert_eq!(result[2], 1);
    }

    #[test]
    fn test_decompress_rice_i16() {
        let mut data = vec![0x12, 0x34];
        data.extend_from_slice(&[0b00000000]);

        let result = decompress_rice::<i16>(&data, 1, 32).unwrap();
        assert_eq!(result, vec![0x1234]);
    }

    #[test]
    fn test_decompress_rice_i8() {
        let mut data = vec![0x42];
        data.extend_from_slice(&[0b00010000]);

        let result = decompress_rice::<i8>(&data, 1, 32).unwrap();
        assert_eq!(result, vec![0x42]);
    }

    #[test]
    fn test_decompress_rice_u8() {
        let mut data = vec![0x42];
        data.extend_from_slice(&[0b00010000]);

        let result = u8::decompress(&data, 1, 32).unwrap();
        assert_eq!(result, vec![0x42]);
    }

    #[test]
    fn test_decompress_rice_u16() {
        let mut data = vec![0x12, 0x34];
        data.extend_from_slice(&[0b00000000]);

        let result = u16::decompress(&data, 1, 32).unwrap();
        assert_eq!(result, vec![0x1234]);
    }

    #[test]
    fn test_decompress_rice_with_real_rice_coding() {
        let original_data = vec![1i32, 3, 2, 5, 4];
        let compressed = compress::compress_rice(&original_data, 32).unwrap();
        let decompressed = decompress_rice::<i32>(&compressed, original_data.len(), 32).unwrap();
        assert_eq!(original_data, decompressed);
    }

    #[test]
    fn test_decompress_rice_error_insufficient_rice_data() {
        let mut data = vec![0x00, 0x00, 0x00, 0x01];
        data.extend_from_slice(&[0b00010000]);

        let result = decompress_rice::<i32>(&data, 2, 32);
        assert!(result.is_err());
    }

    #[test]
    fn test_high_entropy_block_edge_cases() {
        let data = vec![0x00, 0x00, 0x00, 0x01, 0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        high_entropy_block::<i32>(&mut reader, &mut result, &mut last_pix, 1).unwrap();
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_rice_block_zero_leading_zeros() {
        let data = vec![0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        rice_block::<i32>(&mut reader, &mut result, &mut last_pix, 1, 1).unwrap();
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_rice_block_max_fs() {
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 0;

        rice_block::<i32>(&mut reader, &mut result, &mut last_pix, 1, 24).unwrap();
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_process_block_boundary_fs_values() {
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        let mut result = Vec::new();
        let mut last_pix = 5;

        process_block::<i32>(&mut reader, &mut result, &mut last_pix, 2, 25).unwrap();
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_decompress_large_dataset() {
        let data = vec![1i32, 2, 3, 4, 5, 4, 3, 2, 1, 0];
        let compressed = compress::compress_rice(&data, 32).unwrap();
        let decompressed = decompress_rice::<i32>(&compressed, data.len(), 32).unwrap();
        assert_eq!(data, decompressed);
    }

    #[test]
    fn test_decompress_different_block_sizes() {
        let data = vec![10i32, 20, 30, 40, 50, 60, 70, 80];

        for &block_size in &[1, 2, 3, 4, 8, 16] {
            let compressed = compress::compress_rice(&data, block_size).unwrap();
            let decompressed = decompress_rice::<i32>(&compressed, data.len(), block_size).unwrap();
            assert_eq!(data, decompressed, "Failed for block size {}", block_size);
        }
    }

    #[test]
    fn test_read_first_value_signed_boundary() {
        let data = vec![0x80, 0x00, 0x00, 0x00];
        let result = read_first_value_direct::<i32>(&data).unwrap();
        assert_eq!(result, i32::MIN);

        let data = vec![0x7F, 0xFF, 0xFF, 0xFF];
        let result = read_first_value_direct::<i32>(&data).unwrap();
        assert_eq!(result, i32::MAX);
    }
}
