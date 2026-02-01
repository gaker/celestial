use crate::fits::{FitsError, Result};

pub struct BitWriter {
    buffer: Vec<u8>,
    bit_buffer: u32,
    bits_to_go: u8,
}

impl BitWriter {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            bit_buffer: 0,
            bits_to_go: 8,
        }
    }

    pub fn output_nbits(&mut self, bits: u32, n: usize) -> Result<()> {
        if n > 32 {
            return Err(FitsError::InvalidFormat(
                "Cannot output more than 32 bits at once".to_string(),
            ));
        }

        const MASK: [u32; 33] = [
            0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff, 0x7ff, 0xfff, 0x1fff,
            0x3fff, 0x7fff, 0xffff, 0x1ffff, 0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff,
            0x7fffff, 0xffffff, 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff,
            0x7fffffff, 0xffffffff,
        ];

        let mut local_bit_buffer = self.bit_buffer;
        let mut local_bits_to_go = self.bits_to_go as i8;
        let mut remaining_bits = n as i8;

        if local_bits_to_go + remaining_bits > 32 {
            local_bit_buffer <<= local_bits_to_go;
            local_bit_buffer |=
                (bits >> (remaining_bits - local_bits_to_go)) & MASK[local_bits_to_go as usize];
            self.buffer.push((local_bit_buffer & 0xff) as u8);
            remaining_bits -= local_bits_to_go;
            local_bits_to_go = 8;
        }

        local_bit_buffer <<= remaining_bits;
        local_bit_buffer |= bits & MASK[remaining_bits as usize];
        local_bits_to_go -= remaining_bits;

        while local_bits_to_go <= 0 {
            self.buffer
                .push(((local_bit_buffer >> (-local_bits_to_go)) & 0xff) as u8);
            local_bits_to_go += 8;
        }

        self.bit_buffer = local_bit_buffer;
        self.bits_to_go = local_bits_to_go as u8;

        Ok(())
    }

    pub fn flush(&mut self) {
        if self.bits_to_go < 8 {
            self.buffer.push((self.bit_buffer << self.bits_to_go) as u8);
        }
    }

    pub fn finish(mut self) -> Vec<u8> {
        self.flush();
        self.buffer
    }
}

pub struct BitReader<'a> {
    data: &'a [u8],
    position: usize,
    bit_buffer: u32,
    bits_remaining: u8,
}

impl<'a> BitReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        let mut reader = Self {
            data,
            position: 0,
            bit_buffer: 0,
            bits_remaining: 0,
        };

        // Initialize with first byte if available
        if !data.is_empty() {
            reader.bit_buffer = data[0] as u32;
            reader.position = 1;
            reader.bits_remaining = 8;
        }

        reader
    }

    pub fn read_bits(&mut self, n: usize) -> Result<u32> {
        if n > 32 {
            return Err(FitsError::InvalidFormat(
                "Cannot read more than 32 bits at once".to_string(),
            ));
        }

        if n == 0 {
            return Ok(0);
        }

        // For large reads (like 32 bits), use incremental approach to avoid buffer overflow issues
        if n > 24 {
            return self.read_bits_incremental(n);
        }

        // Ensure we have enough bits
        while self.bits_remaining < n as u8 {
            if self.position >= self.data.len() {
                return Err(FitsError::InvalidFormat(
                    "Unexpected end of compressed data".to_string(),
                ));
            }

            self.bit_buffer = (self.bit_buffer << 8) | (self.data[self.position] as u32);
            self.position += 1;
            self.bits_remaining += 8;
        }

        // Extract the bits we need (high-order n bits)
        let mask = if n == 32 { 0xFFFFFFFF } else { (1u32 << n) - 1 };
        self.bits_remaining -= n as u8;
        let result = (self.bit_buffer >> self.bits_remaining) & mask;
        self.bit_buffer &= (1 << self.bits_remaining) - 1;

        Ok(result)
    }

    // CFITSIO-style incremental bit reading for large values
    fn read_bits_incremental(&mut self, n: usize) -> Result<u32> {
        let mut result = 0u32;
        let mut bits_needed = n;

        while bits_needed > 0 {
            let chunk_size = bits_needed.min(8);
            let bits = self.read_bits_small(chunk_size)?;
            result = (result << chunk_size) | bits;
            bits_needed -= chunk_size;
        }

        Ok(result)
    }

    // Read small number of bits (<=8) using the original logic
    fn read_bits_small(&mut self, n: usize) -> Result<u32> {
        if self.bits_remaining < n as u8 {
            if self.position >= self.data.len() {
                return Err(FitsError::InvalidFormat(
                    "Unexpected end of compressed data".to_string(),
                ));
            }

            self.bit_buffer = (self.bit_buffer << 8) | (self.data[self.position] as u32);
            self.position += 1;
            self.bits_remaining += 8;
        }

        self.bits_remaining -= n as u8;
        let result = self.bit_buffer >> self.bits_remaining;
        self.bit_buffer &= (1 << self.bits_remaining) - 1;

        Ok(result & ((1u32 << n) - 1))
    }

    pub fn count_leading_zeros(&mut self) -> Result<usize> {
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
