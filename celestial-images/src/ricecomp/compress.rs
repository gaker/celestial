use super::bitbuffer::BitWriter;
use super::RiceCompressible;
use crate::fits::Result;

pub fn compress_rice<T: RiceCompressible>(data: &[T], block_size: usize) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    let mut writer = BitWriter::with_capacity(data.len() + data.len() / 10);
    write_first_value::<T>(&mut writer, data[0])?;

    let mut last_pix = data[0].as_signed();
    let mut diff_buffer = vec![0u32; block_size];

    compress_blocks::<T>(
        &mut writer,
        data,
        &mut last_pix,
        &mut diff_buffer,
        block_size,
    )?;
    Ok(writer.finish())
}

fn write_first_value<T: RiceCompressible>(writer: &mut BitWriter, value: T) -> Result<()> {
    let val = value.as_signed() as u32;
    writer.output_nbits(val, T::BYTES_PER_ELEMENT * 8)
}

fn compress_blocks<T: RiceCompressible>(
    writer: &mut BitWriter,
    data: &[T],
    last_pix: &mut i32,
    diff_buffer: &mut [u32],
    block_size: usize,
) -> Result<()> {
    let mut i = 0;
    while i < data.len() {
        let this_block = (data.len() - i).min(block_size);
        let pixel_sum = compute_differences(&data[i..i + this_block], last_pix, diff_buffer)?;
        let fs = calculate_split_level(pixel_sum, this_block);

        write_block::<T>(writer, &diff_buffer[..this_block], fs, pixel_sum)?;
        i += this_block;
    }
    Ok(())
}

fn compute_differences<T: RiceCompressible>(
    data: &[T],
    last_pix: &mut i32,
    diff_buffer: &mut [u32],
) -> Result<f64> {
    let mut pixel_sum = 0.0;
    for (j, &pixel) in data.iter().enumerate() {
        let next_pix = pixel.as_signed();

        let pdiff = T::compute_difference(next_pix, *last_pix);

        let mapped_diff = map_signed_to_unsigned(pdiff);
        diff_buffer[j] = mapped_diff;
        pixel_sum += mapped_diff as f64;
        *last_pix = next_pix;
    }
    Ok(pixel_sum)
}

fn map_signed_to_unsigned(pdiff: i32) -> u32 {
    if pdiff < 0 {
        !(pdiff << 1) as u32
    } else {
        (pdiff << 1) as u32
    }
}

fn calculate_split_level(pixel_sum: f64, block_size: usize) -> usize {
    let dpsum = (pixel_sum - (block_size as f64 / 2.0) - 1.0) / block_size as f64;
    let dpsum = dpsum.max(0.0);
    let mut psum = (dpsum as u32) >> 1;
    let mut fs = 0;
    while psum > 0 {
        fs += 1;
        psum >>= 1;
    }
    fs
}

fn write_block<T: RiceCompressible>(
    writer: &mut BitWriter,
    diff_buffer: &[u32],
    fs: usize,
    pixel_sum: f64,
) -> Result<()> {
    if fs >= T::FSMAX {
        write_high_entropy_block::<T>(writer, diff_buffer)
    } else if fs == 0 && pixel_sum == 0.0 {
        write_low_entropy_block::<T>(writer)
    } else {
        write_rice_coded_block::<T>(writer, diff_buffer, fs)
    }
}

fn write_high_entropy_block<T: RiceCompressible>(
    writer: &mut BitWriter,
    diff_buffer: &[u32],
) -> Result<()> {
    writer.output_nbits((T::FSMAX + 1) as u32, T::FSBITS)?;

    for &mapped_diff in diff_buffer {
        writer.output_nbits(mapped_diff, T::BBITS)?;
    }
    Ok(())
}

fn write_low_entropy_block<T: RiceCompressible>(writer: &mut BitWriter) -> Result<()> {
    writer.output_nbits(0, T::FSBITS)
}

fn write_rice_coded_block<T: RiceCompressible>(
    writer: &mut BitWriter,
    diff_buffer: &[u32],
    fs: usize,
) -> Result<()> {
    writer.output_nbits((fs + 1) as u32, T::FSBITS)?;
    let fsmask = (1u32 << fs) - 1;

    for &v in diff_buffer.iter() {
        write_rice_value(writer, v, fs, fsmask)?;
    }
    Ok(())
}

fn write_rice_value(writer: &mut BitWriter, v: u32, fs: usize, fsmask: u32) -> Result<()> {
    let top = v >> fs;

    writer.output_nbits(1, (top + 1) as usize)?;

    if fs > 0 {
        writer.output_nbits(v & fsmask, fs)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compress_rice_empty() {
        let data: Vec<i32> = vec![];
        let result = compress_rice(&data, 32).unwrap();
        assert_eq!(result, Vec::<u8>::new());
    }

    #[test]
    fn test_compress_rice_single_element() {
        let data = vec![42i32];
        let result = compress_rice(&data, 32).unwrap();
        assert!(result.len() >= 4);
    }

    #[test]
    fn test_compress_rice_multiple_elements() {
        let data = vec![1i32, 2, 3, 4, 5];
        let result = compress_rice(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_map_signed_to_unsigned() {
        assert_eq!(map_signed_to_unsigned(0), 0);
        assert_eq!(map_signed_to_unsigned(1), 2);
        assert_eq!(map_signed_to_unsigned(-1), 1);
        assert_eq!(map_signed_to_unsigned(2), 4);
        assert_eq!(map_signed_to_unsigned(-2), 3);
    }

    #[test]
    fn test_calculate_split_level() {
        assert_eq!(calculate_split_level(0.0, 32), 0);
        assert_eq!(calculate_split_level(16.0, 32), 0);
        assert_eq!(calculate_split_level(64.0, 32), 0);
        assert_eq!(calculate_split_level(100.0, 32), 1);
    }

    #[test]
    fn test_write_first_value_i32() {
        let mut writer = BitWriter::with_capacity(10);
        write_first_value::<i32>(&mut writer, 0x12345678).unwrap();
        let result = writer.finish();
        assert_eq!(result, vec![0x12, 0x34, 0x56, 0x78]);
    }

    #[test]
    fn test_write_first_value_i16() {
        let mut writer = BitWriter::with_capacity(10);
        write_first_value::<i16>(&mut writer, 0x1234).unwrap();
        let result = writer.finish();
        assert_eq!(result, vec![0x12, 0x34]);
    }

    #[test]
    fn test_write_first_value_i8() {
        let mut writer = BitWriter::with_capacity(10);
        write_first_value::<i8>(&mut writer, 0x12).unwrap();
        let result = writer.finish();
        assert_eq!(result, vec![0x12]);
    }

    #[test]
    fn test_compute_differences() {
        let data = vec![1i32, 3, 2, 5];
        let mut last_pix = 0;
        let mut diff_buffer = vec![0u32; 4];

        let pixel_sum = compute_differences(&data, &mut last_pix, &mut diff_buffer).unwrap();

        assert_eq!(diff_buffer[0], map_signed_to_unsigned(1));
        assert_eq!(diff_buffer[1], map_signed_to_unsigned(2));
        assert_eq!(diff_buffer[2], map_signed_to_unsigned(-1));
        assert_eq!(diff_buffer[3], map_signed_to_unsigned(3));

        let expected_sum = (map_signed_to_unsigned(1)
            + map_signed_to_unsigned(2)
            + map_signed_to_unsigned(-1)
            + map_signed_to_unsigned(3)) as f64;
        assert_eq!(pixel_sum, expected_sum);
    }

    #[test]
    fn test_compress_rice_zero_differences() {
        let data = vec![5i32, 5, 5, 5];
        let result = compress_rice(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_i16() {
        let data = vec![1i16, 2, 3, 4, 5];
        let result = compress_rice(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_i8() {
        let data = vec![1i8, 2, 3, 4, 5];
        let result = compress_rice(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_u8() {
        let data = vec![1u8, 2, 3, 4, 5];
        let result = u8::compress(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_u16() {
        let data = vec![1u16, 2, 3, 4, 5];
        let result = u16::compress(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_small_block_size() {
        let data = vec![1i32, 2, 3, 4, 5];
        let result = compress_rice(&data, 2).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_large_values() {
        let data = vec![i32::MAX, i32::MIN, 0, 1000000, -1000000];
        let result = compress_rice(&data, 32).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_rice_value() {
        let mut writer = BitWriter::with_capacity(10);
        write_rice_value(&mut writer, 5, 2, 0b11).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_rice_value_zero_fs() {
        let mut writer = BitWriter::with_capacity(10);
        write_rice_value(&mut writer, 3, 0, 0).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_high_entropy_block() {
        let mut writer = BitWriter::with_capacity(100);
        let diff_buffer = vec![1u32, 2, 3, 4];
        write_high_entropy_block::<i32>(&mut writer, &diff_buffer).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_low_entropy_block() {
        let mut writer = BitWriter::with_capacity(10);
        write_low_entropy_block::<i32>(&mut writer).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_rice_coded_block() {
        let mut writer = BitWriter::with_capacity(100);
        let diff_buffer = vec![1u32, 2, 3, 4];
        write_rice_coded_block::<i32>(&mut writer, &diff_buffer, 2).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_blocks_single_block() {
        let mut writer = BitWriter::with_capacity(100);
        let data = vec![1i32, 2, 3];
        let mut last_pix = 0;
        let mut diff_buffer = vec![0u32; 32];

        compress_blocks::<i32>(&mut writer, &data, &mut last_pix, &mut diff_buffer, 32).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_blocks_multiple_blocks() {
        let mut writer = BitWriter::with_capacity(100);
        let data = vec![1i32, 2, 3, 4, 5];
        let mut last_pix = 0;
        let mut diff_buffer = vec![0u32; 2];

        compress_blocks::<i32>(&mut writer, &data, &mut last_pix, &mut diff_buffer, 2).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_block_high_entropy() {
        let mut writer = BitWriter::with_capacity(100);
        let diff_buffer = vec![1000u32, 2000, 3000, 4000];

        write_block::<i32>(&mut writer, &diff_buffer, 25, 5000.0).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_write_block_exact_low_entropy() {
        let mut writer = BitWriter::with_capacity(100);
        let diff_buffer = vec![0u32, 0, 0, 0];

        write_block::<i32>(&mut writer, &diff_buffer, 0, 0.0).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_calculate_split_level_edge_cases() {
        assert_eq!(calculate_split_level(-100.0, 32), 0);
        assert_eq!(calculate_split_level(0.5, 1), 0);
        assert_eq!(calculate_split_level(1000.0, 10), 6);
    }

    #[test]
    fn test_write_rice_value_large() {
        let mut writer = BitWriter::with_capacity(100);

        write_rice_value(&mut writer, 1023, 10, (1 << 10) - 1).unwrap();
        let result = writer.finish();
        assert!(!result.is_empty());
    }

    #[test]
    fn test_compress_rice_very_large_block() {
        let data: Vec<i32> = (0..10000).collect();
        let result = compress_rice(&data, 100).unwrap();
        assert!(!result.is_empty());
        assert!(result.len() > 4);
    }

    #[test]
    fn test_compress_with_different_data_patterns() {
        let data = vec![42i32; 1000];
        let result1 = compress_rice(&data, 32).unwrap();

        let data2: Vec<i32> = (0..1000).map(|x| x * 17 + 23).collect();
        let result2 = compress_rice(&data2, 32).unwrap();

        assert!(result1.len() < result2.len());
    }

    #[test]
    fn test_map_signed_to_unsigned_edge_cases() {
        assert_eq!(
            map_signed_to_unsigned(i32::MIN),
            (i32::MAX as u32).wrapping_mul(2).wrapping_add(1)
        );
        assert_eq!(
            map_signed_to_unsigned(i32::MAX),
            (i32::MAX as u32).wrapping_mul(2)
        );
    }
}
