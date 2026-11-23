use crate::fits::errors::{FitsError, Result};

pub(super) struct BitWriter {
    buffer: Vec<u8>,
    bit_buffer: u32,
    bits_to_go: u8,
}

impl BitWriter {
    pub(super) fn with_capacity(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            bit_buffer: 0,
            bits_to_go: 8,
        }
    }

    // `bit_buffer` holds the 8 - `bits_to_go` bits of the partial byte. Those
    // bits and `n` <= 32 new bits fit 40 bits, so a 64-bit buffer takes the
    // whole value in one shift.
    pub(super) fn output_nbits(&mut self, bits: u32, n: usize) -> Result<()> {
        if n > 32 {
            return Err(FitsError::InvalidFormat(
                "Cannot output more than 32 bits at once".to_string(),
            ));
        }
        let value = u64::from(bits) & ((1u64 << n) - 1);
        let buffer = (u64::from(self.bit_buffer) << n) | value;
        let mut bits_to_go = i32::from(self.bits_to_go) - n as i32;
        while bits_to_go <= 0 {
            self.buffer.push((buffer >> -bits_to_go) as u8);
            bits_to_go += 8;
        }
        self.bit_buffer = (buffer & 0xff) as u32;
        self.bits_to_go = bits_to_go as u8;
        Ok(())
    }

    pub(super) fn flush(&mut self) {
        if self.bits_to_go < 8 {
            self.buffer.push((self.bit_buffer << self.bits_to_go) as u8);
        }
    }

    pub(super) fn finish(mut self) -> Vec<u8> {
        self.flush();
        self.buffer
    }
}

// `bit_buffer` holds the `bits_remaining` unread bits of the bytes read so
// far, at most 8 between reads. A read of up to 32 bits fills it to at most
// 39 bits, so it is 64 bits wide.
pub(super) struct BitReader<'a> {
    data: &'a [u8],
    position: usize,
    bit_buffer: u64,
    bits_remaining: u8,
}

impl<'a> BitReader<'a> {
    pub(super) fn new(data: &'a [u8]) -> Self {
        let mut reader = Self {
            data,
            position: 0,
            bit_buffer: 0,
            bits_remaining: 0,
        };

        // Initialize with first byte if available
        if !data.is_empty() {
            reader.bit_buffer = u64::from(data[0]);
            reader.position = 1;
            reader.bits_remaining = 8;
        }

        reader
    }

    pub(super) fn read_bits(&mut self, n: usize) -> Result<u32> {
        if n > 32 {
            return Err(FitsError::InvalidFormat(
                "Cannot read more than 32 bits at once".to_string(),
            ));
        }
        while usize::from(self.bits_remaining) < n {
            let byte = self.data.get(self.position).ok_or_else(|| {
                FitsError::InvalidFormat("Unexpected end of compressed data".to_string())
            })?;
            self.bit_buffer = (self.bit_buffer << 8) | u64::from(*byte);
            self.position += 1;
            self.bits_remaining += 8;
        }
        self.bits_remaining -= n as u8;
        let result = (self.bit_buffer >> self.bits_remaining) & ((1u64 << n) - 1);
        self.bit_buffer &= (1u64 << self.bits_remaining) - 1;
        Ok(result as u32)
    }

    pub(super) fn count_leading_zeros(&mut self) -> Result<usize> {
        let mut count = 0;

        while self.read_bits(1)? == 0 {
            count += 1;
        }

        Ok(count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_writer() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(0b1010, 4).is_ok());
        let result = writer.finish();
        assert_eq!(result[0], 0b10100000);
    }

    #[test]
    fn test_bit_writer_multiple_writes() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(0b1010, 4).is_ok());
        assert!(writer.output_nbits(0b1100, 4).is_ok());
        let result = writer.finish();
        assert_eq!(result[0], 0b10101100);
    }

    #[test]
    fn test_bit_writer_cross_byte_boundary() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(0b111111, 6).is_ok());
        assert!(writer.output_nbits(0b10101010, 8).is_ok());
        let result = writer.finish();
        assert_eq!(result[0], 0b11111110);
        assert_eq!(result[1], 0b10101000); // Last 2 bits shifted to make room for next write
    }

    #[test]
    fn test_bit_writer_large_values() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(0x12345678, 32).is_ok());
        let result = writer.finish();
        assert_eq!(result, vec![0x12, 0x34, 0x56, 0x78]);
    }

    // Pending bits of a partial byte plus a wide value can exceed 32 bits.
    #[test]
    fn bit_writer_keeps_pending_bits_before_a_wide_value() {
        for pending in 1..8 {
            for n in 25..=32 {
                let mut writer = BitWriter::with_capacity(8);
                writer.output_nbits(u32::MAX, pending).unwrap();
                writer.output_nbits(1, n).unwrap();
                let total = pending + n;
                let expected = ((1u64 << pending) - 1) << n | 1;
                let expected = expected << (64 - total);
                let bytes = expected.to_be_bytes()[..total.div_ceil(8)].to_vec();
                assert_eq!(writer.finish(), bytes, "{pending} then {n}");
            }
        }
    }

    #[test]
    fn test_bit_writer_zero_bits() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(123, 0).is_ok());
        let result = writer.finish();
        assert_eq!(result, Vec::<u8>::new());
    }

    #[test]
    fn test_bit_writer_invalid_size() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(123, 33).is_err());
    }

    #[test]
    fn test_bit_writer_single_bit() {
        let mut writer = BitWriter::with_capacity(10);
        assert!(writer.output_nbits(1, 1).is_ok());
        assert!(writer.output_nbits(0, 1).is_ok());
        assert!(writer.output_nbits(1, 1).is_ok());
        let result = writer.finish();
        assert_eq!(result[0], 0b10100000);
    }

    #[test]
    fn test_bit_writer_empty() {
        let writer = BitWriter::with_capacity(10);
        let result = writer.finish();
        assert_eq!(result, Vec::<u8>::new());
    }

    #[test]
    fn test_bit_reader() {
        let data = vec![0b10101100];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.read_bits(4).unwrap(), 0b1010);
        assert_eq!(reader.read_bits(4).unwrap(), 0b1100);
    }

    #[test]
    fn test_bit_reader_cross_byte_boundary() {
        let data = vec![0b11111110, 0b10101010];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.read_bits(6).unwrap(), 0b111111);
        assert_eq!(reader.read_bits(8).unwrap(), 0b10101010);
        assert_eq!(reader.read_bits(2).unwrap(), 0b10);
    }

    #[test]
    fn test_bit_reader_32_bit_read() {
        let data = vec![0x12, 0x34, 0x56, 0x78];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.read_bits(32).unwrap(), 0x12345678);
    }

    #[test]
    fn test_bit_reader_zero_bits() {
        let data = vec![0xFF];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.read_bits(0).unwrap(), 0);
        assert_eq!(reader.read_bits(8).unwrap(), 0xFF);
    }

    #[test]
    fn test_bit_reader_empty_data() {
        let data = vec![];
        let mut reader = BitReader::new(&data);

        assert!(reader.read_bits(1).is_err());
    }

    #[test]
    fn test_bit_reader_insufficient_data() {
        let data = vec![0xFF];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.read_bits(8).unwrap(), 0xFF);
        assert!(reader.read_bits(1).is_err());
    }

    #[test]
    fn test_bit_reader_invalid_size() {
        let data = vec![0xFF];
        let mut reader = BitReader::new(&data);

        assert!(reader.read_bits(33).is_err());
    }

    #[test]
    fn test_bit_reader_count_leading_zeros() {
        let data = vec![0b00001111];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.count_leading_zeros().unwrap(), 4);
        assert_eq!(reader.read_bits(3).unwrap(), 0b111); // Only 3 bits left after consuming one bit for the terminating '1'
    }

    #[test]
    fn test_bit_reader_count_leading_zeros_cross_boundary() {
        let data = vec![0b00000000, 0b00111111];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.count_leading_zeros().unwrap(), 10);
        assert_eq!(reader.read_bits(5).unwrap(), 0b11111); // Only 5 bits left after consuming one bit for the terminating '1'
    }

    #[test]
    fn test_bit_reader_count_leading_zeros_immediate_one() {
        let data = vec![0b11111111];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.count_leading_zeros().unwrap(), 0);
        assert_eq!(reader.read_bits(7).unwrap(), 0b1111111);
    }

    #[test]
    fn test_bit_reader_incremental_large() {
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);

        assert_eq!(reader.read_bits(25).unwrap(), 0x1FFFFFF);
        assert_eq!(reader.read_bits(7).unwrap(), 0x7F);
    }

    #[test]
    fn test_round_trip_small() {
        let mut writer = BitWriter::with_capacity(10);
        writer.output_nbits(0b1010, 4).unwrap();
        writer.output_nbits(0b1100, 4).unwrap();
        let data = writer.finish();

        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(4).unwrap(), 0b1010);
        assert_eq!(reader.read_bits(4).unwrap(), 0b1100);
    }

    #[test]
    fn test_round_trip_complex() {
        let mut writer = BitWriter::with_capacity(20);
        writer.output_nbits(0b111, 3).unwrap();
        writer.output_nbits(0x1234, 16).unwrap();
        writer.output_nbits(0b10101, 5).unwrap();
        let data = writer.finish();

        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(3).unwrap(), 0b111);
        assert_eq!(reader.read_bits(16).unwrap(), 0x1234);
        assert_eq!(reader.read_bits(5).unwrap(), 0b10101);
    }

    #[test]
    fn test_bit_writer_flush() {
        let mut writer = BitWriter::with_capacity(10);
        writer.output_nbits(0b1010, 4).unwrap();
        writer.flush();
        let result = writer.finish();
        assert_eq!(result[0], 0b10100000);
    }

    #[test]
    fn test_bit_writer_partial_byte() {
        let mut writer = BitWriter::with_capacity(10);
        writer.output_nbits(0b101, 3).unwrap();
        let result = writer.finish();
        assert_eq!(result[0], 0b10100000);
    }

    #[test]
    fn test_bit_reader_read_bits_small() {
        let data = vec![0b11110000];
        let mut reader = BitReader::new(&data);

        // This should trigger the read_bits_small path for incremental reading
        assert_eq!(reader.read_bits(4).unwrap(), 0b1111);
        assert_eq!(reader.read_bits(4).unwrap(), 0b0000);
    }

    #[test]
    fn test_bit_reader_mask_edge_case() {
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF];
        let mut reader = BitReader::new(&data);

        // Test with n=32 to hit the mask edge case
        assert_eq!(reader.read_bits(32).unwrap(), 0xFFFFFFFF);
    }

    #[test]
    fn test_bit_reader_leading_zeros_end_of_data() {
        let data = vec![0b00000000];
        let mut reader = BitReader::new(&data);

        // Should fail when trying to count zeros past end of data
        assert!(reader.count_leading_zeros().is_err());
    }

    #[test]
    fn test_bit_writer_output_nbits_exact_capacity() {
        let mut writer = BitWriter::with_capacity(1);
        assert!(writer.output_nbits(0xFF, 8).is_ok());
        let result = writer.finish();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 0xFF);
    }

    #[test]
    fn test_bit_reader_multiple_refills() {
        let data = vec![0x12, 0x34, 0x56, 0x78, 0x9A];
        let mut reader = BitReader::new(&data);

        // Read in small chunks to force multiple buffer refills
        assert_eq!(reader.read_bits(4).unwrap(), 0x1);
        assert_eq!(reader.read_bits(4).unwrap(), 0x2);
        assert_eq!(reader.read_bits(8).unwrap(), 0x34);
        assert_eq!(reader.read_bits(16).unwrap(), 0x5678);
        assert_eq!(reader.read_bits(8).unwrap(), 0x9A);
    }
}
