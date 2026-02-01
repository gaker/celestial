use super::BinaryTableHdu;
use crate::fits::compression::{decompress_tile, CompressionAlgorithm, DecompressionParams};
use crate::fits::{FitsError, Result};
use std::io::{Read, Seek, SeekFrom};

impl BinaryTableHdu {
    pub fn is_compressed_image(&self) -> bool {
        self.header
            .get_keyword_value("ZIMAGE")
            .and_then(|v| v.as_logical())
            .unwrap_or(false)
    }

    pub fn compression_algorithm(&self) -> Option<&str> {
        if self.is_compressed_image() {
            self.header
                .get_keyword_value("ZCMPTYPE")
                .and_then(|v| v.as_string())
        } else {
            None
        }
    }

    pub fn quantization_level(&self) -> Option<i64> {
        if self.is_compressed_image() {
            self.header
                .get_keyword_value("ZQUANTIZ")
                .and_then(|v| v.as_integer())
        } else {
            None
        }
    }

    pub fn get_compression_algorithm(&self) -> Option<CompressionAlgorithm> {
        self.compression_algorithm()
            .and_then(CompressionAlgorithm::from_fits_name)
    }

    pub fn get_tile_dimensions(&self) -> Result<(usize, usize)> {
        let znaxis1 = self
            .header
            .get_keyword_value("ZNAXIS1")
            .and_then(|v| v.as_integer())
            .ok_or_else(|| FitsError::KeywordNotFound {
                keyword: "ZNAXIS1".to_string(),
            })?;

        let znaxis2 = self
            .header
            .get_keyword_value("ZNAXIS2")
            .and_then(|v| v.as_integer())
            .ok_or_else(|| FitsError::KeywordNotFound {
                keyword: "ZNAXIS2".to_string(),
            })?;

        Ok((znaxis1 as usize, znaxis2 as usize))
    }

    pub fn get_bits_per_pixel(&self) -> Result<i32> {
        self.header
            .get_keyword_value("ZBITPIX")
            .and_then(|v| v.as_integer())
            .map(|v| v as i32)
            .ok_or_else(|| FitsError::KeywordNotFound {
                keyword: "ZBITPIX".to_string(),
            })
    }

    pub fn decompress_image_tile<R>(&self, reader: &mut R, tile_index: usize) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        if !self.is_compressed_image() {
            return Err(FitsError::InvalidFormat(
                "HDU is not a compressed image".to_string(),
            ));
        }

        let algorithm = self.get_compression_algorithm().ok_or_else(|| {
            FitsError::InvalidFormat("Unknown or unsupported compression algorithm".to_string())
        })?;

        let tile_dimensions = self.get_tile_dimensions()?;
        let bits_per_pixel = self.get_bits_per_pixel()?;
        let quantization_level = self.quantization_level();

        let params = DecompressionParams::new(
            algorithm,
            quantization_level,
            tile_dimensions,
            bits_per_pixel,
        );

        let compressed_data = self.read_tile_data(reader, tile_index)?;

        decompress_tile(&compressed_data, &params)
    }

    fn read_tile_data<R>(&self, reader: &mut R, tile_index: usize) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        let compressed_data_col = self.column_by_name("COMPRESSED_DATA").unwrap_or(0);

        let row_count = self.number_of_rows().unwrap_or(0) as usize;
        if tile_index >= row_count {
            return Err(FitsError::InvalidFormat(format!(
                "Tile index {} out of range (0..{})",
                tile_index, row_count
            )));
        }

        let col_info = self.column_info(compressed_data_col)?;
        let (data_type, _repeat) = self.parse_binary_format(&col_info.format)?;

        if data_type.starts_with('P') || data_type.starts_with('Q') {
            self.read_variable_length_tile_data(reader, compressed_data_col, tile_index)
        } else {
            self.read_fixed_length_tile_data(reader, compressed_data_col, tile_index)
        }
    }

    fn read_variable_length_tile_data<R>(
        &self,
        reader: &mut R,
        column: usize,
        tile_index: usize,
    ) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        let read_params = self.prepare_column_read(column)?;
        let descriptor_size = read_params.bytes_per_element;

        let data_start = self.info.data_start;
        let row_offset = tile_index * read_params.row_size;
        let descriptor_position = data_start + row_offset as u64 + read_params.column_offset as u64;

        reader.seek(SeekFrom::Start(descriptor_position))?;

        let mut descriptor = vec![0u8; descriptor_size];
        reader.read_exact(&mut descriptor)?;

        let (element_count, heap_offset) = if descriptor_size == 8 {
            let count =
                u32::from_be_bytes([descriptor[0], descriptor[1], descriptor[2], descriptor[3]])
                    as usize;
            let offset =
                u32::from_be_bytes([descriptor[4], descriptor[5], descriptor[6], descriptor[7]])
                    as u64;
            (count, offset)
        } else {
            let count = u64::from_be_bytes([
                descriptor[0],
                descriptor[1],
                descriptor[2],
                descriptor[3],
                descriptor[4],
                descriptor[5],
                descriptor[6],
                descriptor[7],
            ]) as usize;
            let offset = u64::from_be_bytes([
                descriptor[8],
                descriptor[9],
                descriptor[10],
                descriptor[11],
                descriptor[12],
                descriptor[13],
                descriptor[14],
                descriptor[15],
            ]);
            (count, offset)
        };

        if element_count == 0 {
            return Ok(Vec::new());
        }

        let table_rows = self.number_of_rows().unwrap_or(0) as usize;
        let heap_start = data_start + (table_rows * read_params.row_size) as u64;
        let tile_data_position = heap_start + heap_offset;

        reader.seek(SeekFrom::Start(tile_data_position))?;
        let mut tile_data = vec![0u8; element_count];
        reader.read_exact(&mut tile_data)?;

        Ok(tile_data)
    }

    fn read_fixed_length_tile_data<R>(
        &self,
        reader: &mut R,
        column: usize,
        tile_index: usize,
    ) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        let read_params = self.prepare_column_read(column)?;
        let column_bytes = read_params.width * read_params.bytes_per_element;

        let data_start = self.info.data_start;
        let row_offset = tile_index * read_params.row_size;
        let column_position = data_start + row_offset as u64 + read_params.column_offset as u64;

        reader.seek(SeekFrom::Start(column_position))?;
        let mut tile_data = vec![0u8; column_bytes];
        reader.read_exact(&mut tile_data)?;

        Ok(tile_data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Header, Keyword};
    use crate::fits::io::reader::HduInfo;
    use std::io::Cursor;

    fn create_compressed_header() -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::logical("ZIMAGE", true));
        header.add_keyword(Keyword::string("ZCMPTYPE", "GZIP_1"));
        header.add_keyword(Keyword::integer("ZQUANTIZ", 16));
        header.add_keyword(Keyword::integer("ZNAXIS1", 512));
        header.add_keyword(Keyword::integer("ZNAXIS2", 512));
        header.add_keyword(Keyword::integer("ZBITPIX", -32));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TTYPE1", "COMPRESSED_DATA"));
        header.add_keyword(Keyword::string("TFORM1", "1PB"));
        header.add_keyword(Keyword::integer("NAXIS2", 10));
        header
    }

    fn create_uncompressed_header() -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::integer("NAXIS2", 5));
        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 1000,
        }
    }

    #[test]
    fn is_compressed_image_returns_true_when_zimage_true() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert!(hdu.is_compressed_image());
    }

    #[test]
    fn is_compressed_image_returns_false_when_zimage_false() {
        let mut header = create_compressed_header();
        header.add_keyword(Keyword::logical("ZIMAGE", false));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert!(!hdu.is_compressed_image());
    }

    #[test]
    fn is_compressed_image_returns_false_when_zimage_missing() {
        let header = create_uncompressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert!(!hdu.is_compressed_image());
    }

    #[test]
    fn compression_algorithm_returns_some_when_compressed() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.compression_algorithm(), Some("GZIP_1"));
    }

    #[test]
    fn compression_algorithm_returns_none_when_not_compressed() {
        let header = create_uncompressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.compression_algorithm(), None);
    }

    #[test]
    fn quantization_level_returns_some_when_compressed() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.quantization_level(), Some(16));
    }

    #[test]
    fn quantization_level_returns_none_when_not_compressed() {
        let header = create_uncompressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.quantization_level(), None);
    }

    #[test]
    fn get_compression_algorithm_returns_gzip() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let algorithm = hdu.get_compression_algorithm();
        assert!(algorithm.is_some());
        assert_eq!(algorithm.unwrap().fits_name(), "GZIP_1");
    }

    #[test]
    fn get_compression_algorithm_returns_none_for_unknown() {
        let mut header = create_compressed_header();
        header.add_keyword(Keyword::string("ZCMPTYPE", "UNKNOWN_ALGO"));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert!(hdu.get_compression_algorithm().is_none());
    }

    #[test]
    fn get_tile_dimensions_returns_correct_values() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.get_tile_dimensions();
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), (512, 512));
    }

    #[test]
    fn get_tile_dimensions_fails_when_znaxis1_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::logical("ZIMAGE", true));
        header.add_keyword(Keyword::string("ZCMPTYPE", "GZIP_1"));
        header.add_keyword(Keyword::integer("ZQUANTIZ", 16));
        header.add_keyword(Keyword::integer("ZNAXIS2", 512));
        header.add_keyword(Keyword::integer("ZBITPIX", -32));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.get_tile_dimensions();
        assert!(result.is_err());
    }

    #[test]
    fn get_tile_dimensions_fails_when_znaxis2_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::logical("ZIMAGE", true));
        header.add_keyword(Keyword::string("ZCMPTYPE", "GZIP_1"));
        header.add_keyword(Keyword::integer("ZQUANTIZ", 16));
        header.add_keyword(Keyword::integer("ZNAXIS1", 512));
        header.add_keyword(Keyword::integer("ZBITPIX", -32));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.get_tile_dimensions();
        assert!(result.is_err());
    }

    #[test]
    fn get_bits_per_pixel_returns_correct_value() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.get_bits_per_pixel();
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), -32);
    }

    #[test]
    fn get_bits_per_pixel_fails_when_zbitpix_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::logical("ZIMAGE", true));
        header.add_keyword(Keyword::string("ZCMPTYPE", "GZIP_1"));
        header.add_keyword(Keyword::integer("ZQUANTIZ", 16));
        header.add_keyword(Keyword::integer("ZNAXIS1", 512));
        header.add_keyword(Keyword::integer("ZNAXIS2", 512));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.get_bits_per_pixel();
        assert!(result.is_err());
    }

    #[test]
    fn decompress_image_tile_fails_when_not_compressed() {
        let header = create_uncompressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 1000]);

        let result = hdu.decompress_image_tile(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn decompress_image_tile_fails_with_unknown_algorithm() {
        let mut header = create_compressed_header();
        header.add_keyword(Keyword::string("ZCMPTYPE", "UNKNOWN"));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 1000]);

        let result = hdu.decompress_image_tile(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn read_tile_data_fails_with_invalid_tile_index() {
        let header = create_compressed_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 1000]);

        let result = hdu.read_tile_data(&mut cursor, 20);
        assert!(result.is_err());
    }

    #[test]
    fn read_tile_data_with_p_format() {
        let mut header = create_compressed_header();
        header.add_keyword(Keyword::string("TFORM1", "1PB"));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 1000]);

        let result = hdu.read_tile_data(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn read_tile_data_with_q_format() {
        let mut header = create_compressed_header();
        header.add_keyword(Keyword::string("TFORM1", "1QB"));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 1000]);

        let result = hdu.read_tile_data(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn read_tile_data_with_fixed_format() {
        let mut header = create_compressed_header();
        header.add_keyword(Keyword::string("TFORM1", "10B"));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 1000]);

        let result = hdu.read_tile_data(&mut cursor, 0);
        assert!(result.is_err());
    }
}
