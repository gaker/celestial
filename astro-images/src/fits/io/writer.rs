use crate::core::ByteOrder;
use crate::fits::compression::{compress_tile, CompressionAlgorithm, CompressionParams};
use crate::fits::data::array::DataArray;
use crate::fits::header::{Header, Keyword, KeywordValue};
use crate::fits::util::checksum;
use crate::fits::Result;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const FITS_BLOCK_SIZE: usize = 2880;

struct TileExtractionParams {
    width: usize,
    col: usize,
    row: usize,
    tile_width: usize,
    tile_height: usize,
    bytes_per_pixel: usize,
}

struct CompressedHeaderParams {
    nrows: usize,
    row_size: usize,
    heap_size: usize,
}

pub struct FitsWriter {
    writer: BufWriter<File>,
}

impl FitsWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    pub fn write_primary_image<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        keywords: &[Keyword],
    ) -> Result<()>
    where
        T: DataArray + Clone,
    {
        self.write_primary_image_internal(data, dimensions, keywords, false)
    }

    pub fn write_primary_image_with_checksum<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        keywords: &[Keyword],
    ) -> Result<()>
    where
        T: DataArray + Clone,
    {
        self.write_primary_image_internal(data, dimensions, keywords, true)
    }

    fn write_primary_image_internal<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        keywords: &[Keyword],
        with_checksum: bool,
    ) -> Result<()>
    where
        T: DataArray + Clone,
    {
        let mut header = self.build_primary_header::<T>(dimensions, keywords);

        let data_bytes = if !data.is_empty() {
            T::to_bytes(data, ByteOrder::BigEndian)?
        } else {
            Vec::new()
        };

        let padded_data = checksum::pad_to_block(&data_bytes);

        if with_checksum {
            self.write_with_checksums(&mut header, &padded_data)?;
        } else {
            self.write_header_and_data(&header, &padded_data)?;
        }

        self.writer.flush()?;
        Ok(())
    }

    fn build_primary_header<T: DataArray>(
        &self,
        dimensions: &[usize],
        keywords: &[Keyword],
    ) -> Header {
        let mut header = Header::new();

        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("BITPIX", T::BITPIX.value() as i64));
        header.add_keyword(Keyword::integer("NAXIS", dimensions.len() as i64));

        for (i, &dim) in dimensions.iter().enumerate() {
            let axis_keyword = format!("NAXIS{}", i + 1);
            header.add_keyword(Keyword::integer(axis_keyword, dim as i64));
        }

        header.add_keyword(Keyword::logical("EXTEND", false));
        self.add_keywords(&mut header, keywords);

        header
    }

    fn write_with_checksums(&mut self, header: &mut Header, data_bytes: &[u8]) -> Result<()> {
        let data_sum = checksum::calculate_datasum(data_bytes);
        header.add_keyword(Keyword::string(
            "DATASUM",
            &checksum::format_datasum(data_sum),
        ));
        header.add_keyword(Keyword::string("CHECKSUM", "0000000000000000"));
        header.add_keyword(Keyword::new("END".to_string()));

        let header_bytes_with_placeholder =
            self.serialize_header_with_checksum(header, "0000000000000000")?;
        let checksum_value =
            checksum::create_checksum_card_value(&header_bytes_with_placeholder, data_sum);

        self.write_header_with_checksum(header, &checksum_value)?;
        self.writer.write_all(data_bytes)?;
        Ok(())
    }

    fn serialize_header_with_checksum(
        &self,
        header: &Header,
        checksum_value: &str,
    ) -> Result<Vec<u8>> {
        let mut header_bytes = Vec::new();
        for keyword in header.keywords() {
            let card = if keyword.name == "CHECKSUM" {
                self.format_checksum_card(checksum_value)?
            } else {
                self.format_card(keyword)?
            };
            header_bytes.extend_from_slice(&card);
        }
        self.pad_header_bytes_ref(&mut header_bytes);
        Ok(header_bytes)
    }

    fn write_header_with_checksum(&mut self, header: &Header, checksum_value: &str) -> Result<()> {
        let mut header_bytes = Vec::new();

        for keyword in header.keywords() {
            let card = if keyword.name == "CHECKSUM" {
                self.format_checksum_card(checksum_value)?
            } else {
                self.format_card(keyword)?
            };
            header_bytes.extend_from_slice(&card);
        }

        self.pad_header_bytes(&mut header_bytes);
        self.writer.write_all(&header_bytes)?;
        Ok(())
    }

    fn format_checksum_card(&self, checksum_value: &str) -> Result<[u8; 80]> {
        let mut card = [b' '; 80];
        card[0..8].copy_from_slice(b"CHECKSUM");
        card[8] = b'=';
        card[9] = b' ';
        let value_str = format!("'{}'", checksum_value);
        let value_bytes = value_str.as_bytes();
        card[10..10 + value_bytes.len()].copy_from_slice(value_bytes);
        Ok(card)
    }

    fn write_header_and_data(&mut self, header: &Header, data_bytes: &[u8]) -> Result<()> {
        let mut header_copy = header.clone();
        header_copy.add_keyword(Keyword::new("END".to_string()));
        self.write_header(&header_copy)?;
        if !data_bytes.is_empty() {
            self.writer.write_all(data_bytes)?;
        }
        Ok(())
    }

    fn pad_header_bytes(&self, header_bytes: &mut Vec<u8>) {
        let padding_needed = FITS_BLOCK_SIZE - (header_bytes.len() % FITS_BLOCK_SIZE);
        if padding_needed < FITS_BLOCK_SIZE {
            header_bytes.resize(header_bytes.len() + padding_needed, b' ');
        }
    }

    fn pad_header_bytes_ref(&self, header_bytes: &mut Vec<u8>) {
        self.pad_header_bytes(header_bytes);
    }

    fn write_header(&mut self, header: &Header) -> Result<()> {
        let mut header_bytes = Vec::new();

        for keyword in header.keywords() {
            let card = self.format_card(keyword)?;
            header_bytes.extend_from_slice(&card);
        }

        let padding_needed = FITS_BLOCK_SIZE - (header_bytes.len() % FITS_BLOCK_SIZE);
        if padding_needed < FITS_BLOCK_SIZE {
            header_bytes.resize(header_bytes.len() + padding_needed, b' ');
        }

        self.writer.write_all(&header_bytes)?;
        Ok(())
    }

    fn format_card(&self, keyword: &Keyword) -> Result<[u8; 80]> {
        let mut card = [b' '; 80];

        if keyword.name == "END" {
            card[0..3].copy_from_slice(b"END");
            return Ok(card);
        }

        let name_bytes = keyword.name.as_bytes();
        let name_len = name_bytes.len().min(8);
        card[0..name_len].copy_from_slice(&name_bytes[0..name_len]);

        // Handle HISTORY and COMMENT keywords (no value, text goes after keyword)
        if keyword.name == "HISTORY" || keyword.name == "COMMENT" {
            if let Some(comment) = &keyword.comment {
                let comment_bytes = comment.as_bytes();
                let comment_len = comment_bytes.len().min(72);
                card[8..8 + comment_len].copy_from_slice(&comment_bytes[0..comment_len]);
            }
            return Ok(card);
        }

        if let Some(value) = &keyword.value {
            card[8] = b'=';
            card[9] = b' ';

            let value_str = match value {
                KeywordValue::Logical(b) => {
                    format!("{:>20}", if *b { "T" } else { "F" })
                }
                KeywordValue::Integer(i) => {
                    format!("{:>20}", i)
                }
                KeywordValue::Real(f) => {
                    format!("{:>20}", f)
                }
                KeywordValue::String(s) => {
                    let truncated = if s.len() > 18 { &s[..18] } else { s };
                    format!("'{:<18}'", truncated)
                }
                KeywordValue::Complex(r, i) => {
                    format!("({}, {})", r, i)
                }
            };

            let value_bytes = value_str.as_bytes();
            let value_len = value_bytes.len().min(20);
            card[10..10 + value_len].copy_from_slice(&value_bytes[0..value_len]);

            if let Some(comment) = &keyword.comment {
                card[30] = b' ';
                card[31] = b'/';
                card[32] = b' ';

                let comment_bytes = comment.as_bytes();
                let comment_len = comment_bytes.len().min(47);
                card[33..33 + comment_len].copy_from_slice(&comment_bytes[0..comment_len]);
            }
        }

        Ok(card)
    }

    pub fn write_compressed_image<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        tile_size: (usize, usize),
        algorithm: CompressionAlgorithm,
        keywords: &[Keyword],
    ) -> Result<()>
    where
        T: DataArray,
    {
        self.write_minimal_primary_header()?;

        let data_bytes = T::to_bytes(data, ByteOrder::BigEndian)?;
        let compressed_tiles =
            self.compress_tiles(&data_bytes, dimensions, tile_size, T::BITPIX.value())?;

        self.write_compressed_extension(
            &compressed_tiles,
            dimensions,
            tile_size,
            T::BITPIX.value(),
            algorithm,
            keywords,
        )?;

        self.writer.flush()?;
        Ok(())
    }

    fn write_minimal_primary_header(&mut self) -> Result<()> {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::logical("EXTEND", true));
        header.add_keyword(Keyword::new("END".to_string()));
        self.write_header(&header)
    }

    fn compress_tiles(
        &self,
        data: &[u8],
        dimensions: &[usize],
        tile_size: (usize, usize),
        bitpix: i32,
    ) -> Result<Vec<Vec<u8>>> {
        let width = dimensions.first().copied().unwrap_or(1);
        let height = dimensions.get(1).copied().unwrap_or(1);
        let bytes_per_pixel = (bitpix.abs() / 8) as usize;

        let mut tiles = Vec::new();
        let mut row = 0;

        while row < height {
            let tile_height = (height - row).min(tile_size.1);
            let mut col = 0;

            while col < width {
                let tile_width = (width - col).min(tile_size.0);
                let tile_params = TileExtractionParams {
                    width,
                    col,
                    row,
                    tile_width,
                    tile_height,
                    bytes_per_pixel,
                };
                let tile_data = self.extract_tile(data, &tile_params);

                let params = CompressionParams::rice(tile_width, tile_height, bitpix);
                let compressed = compress_tile(&tile_data, &params)?;
                tiles.push(compressed);

                col += tile_size.0;
            }
            row += tile_size.1;
        }

        Ok(tiles)
    }

    fn extract_tile(&self, data: &[u8], params: &TileExtractionParams) -> Vec<u8> {
        let row_bytes = params.width * params.bytes_per_pixel;
        let mut tile_data =
            Vec::with_capacity(params.tile_width * params.tile_height * params.bytes_per_pixel);

        for tile_row in 0..params.tile_height {
            let src_row = params.row + tile_row;
            let start = src_row * row_bytes + params.col * params.bytes_per_pixel;
            let end = start + params.tile_width * params.bytes_per_pixel;
            if end <= data.len() {
                tile_data.extend_from_slice(&data[start..end]);
            }
        }

        tile_data
    }

    fn write_compressed_extension(
        &mut self,
        compressed_tiles: &[Vec<u8>],
        dimensions: &[usize],
        tile_size: (usize, usize),
        bitpix: i32,
        algorithm: CompressionAlgorithm,
        keywords: &[Keyword],
    ) -> Result<()> {
        let heap_data = self.build_heap(compressed_tiles);
        let row_size = 8;
        let nrows = compressed_tiles.len();

        let header_params = CompressedHeaderParams {
            nrows,
            row_size,
            heap_size: heap_data.len(),
        };
        let header = self.build_compressed_header(
            dimensions,
            tile_size,
            bitpix,
            algorithm,
            &header_params,
            keywords,
        );

        self.write_header(&header)?;
        self.write_binary_table_data(compressed_tiles, &heap_data, row_size)?;

        Ok(())
    }

    fn build_heap(&self, compressed_tiles: &[Vec<u8>]) -> Vec<u8> {
        let mut heap = Vec::new();
        for tile in compressed_tiles {
            heap.extend_from_slice(tile);
        }
        heap
    }

    fn build_compressed_header(
        &self,
        dimensions: &[usize],
        tile_size: (usize, usize),
        bitpix: i32,
        algorithm: CompressionAlgorithm,
        params: &CompressedHeaderParams,
        keywords: &[Keyword],
    ) -> Header {
        let mut header = Header::new();

        self.add_extension_keywords(&mut header, params.nrows, params.row_size, params.heap_size);
        self.add_column_keywords(&mut header);
        self.add_compression_keywords(&mut header, dimensions, tile_size, bitpix, algorithm);
        self.add_keywords(&mut header, keywords);

        header.add_keyword(Keyword::new("END".to_string()));
        header
    }

    fn add_extension_keywords(
        &self,
        header: &mut Header,
        nrows: usize,
        row_size: usize,
        heap_size: usize,
    ) {
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", row_size as i64));
        header.add_keyword(Keyword::integer("NAXIS2", nrows as i64));
        header.add_keyword(Keyword::integer("PCOUNT", heap_size as i64));
        header.add_keyword(Keyword::integer("GCOUNT", 1));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
    }

    fn add_column_keywords(&self, header: &mut Header) {
        header.add_keyword(Keyword::string("TTYPE1", "COMPRESSED_DATA"));
        header.add_keyword(Keyword::string("TFORM1", "1PB"));
    }

    fn add_compression_keywords(
        &self,
        header: &mut Header,
        dimensions: &[usize],
        tile_size: (usize, usize),
        bitpix: i32,
        algorithm: CompressionAlgorithm,
    ) {
        header.add_keyword(Keyword::logical("ZIMAGE", true));
        header.add_keyword(Keyword::string("ZCMPTYPE", algorithm.fits_name()));
        header.add_keyword(Keyword::integer("ZBITPIX", bitpix as i64));
        header.add_keyword(Keyword::integer("ZNAXIS", dimensions.len() as i64));

        for (i, &dim) in dimensions.iter().enumerate() {
            header.add_keyword(Keyword::integer(format!("ZNAXIS{}", i + 1), dim as i64));
        }

        header.add_keyword(Keyword::integer("ZTILE1", tile_size.0 as i64));
        header.add_keyword(Keyword::integer("ZTILE2", tile_size.1 as i64));
    }

    fn add_keywords(&self, header: &mut Header, keywords: &[Keyword]) {
        for keyword in keywords {
            // Skip mandatory keywords - they're already added by build_primary_header
            if !keyword.is_mandatory() {
                header.add_keyword(keyword.clone());
            }
        }
    }

    fn write_binary_table_data(
        &mut self,
        compressed_tiles: &[Vec<u8>],
        heap_data: &[u8],
        row_size: usize,
    ) -> Result<()> {
        let mut heap_offset: u32 = 0;

        for tile in compressed_tiles {
            let count = tile.len() as u32;
            self.writer.write_all(&count.to_be_bytes())?;
            self.writer.write_all(&heap_offset.to_be_bytes())?;
            heap_offset += count;
        }

        let table_size = compressed_tiles.len() * row_size;
        let total_size = table_size + heap_data.len();

        self.writer.write_all(heap_data)?;

        let padding_needed = (FITS_BLOCK_SIZE - (total_size % FITS_BLOCK_SIZE)) % FITS_BLOCK_SIZE;
        if padding_needed > 0 {
            let padding = vec![0u8; padding_needed];
            self.writer.write_all(&padding)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Keyword, KeywordValue};
    use tempfile::NamedTempFile;

    #[test]
    fn fits_writer_create_success() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let writer = FitsWriter::create(path);
        assert!(writer.is_ok());
    }

    #[test]
    fn write_primary_image_i16_data() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<i16> = vec![1, 2, 3, 4, 5, 6];
        let dimensions = vec![2, 3];

        let result = writer.write_primary_image(&data, &dimensions, &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn write_primary_image_with_keywords() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<u8> = vec![0, 1, 2, 3];
        let dimensions = vec![2, 2];

        let keywords = vec![
            Keyword::string("TELESCOP", "Test Telescope"),
            Keyword::real("EXPTIME", 30.5),
            Keyword::logical("ISLIGHT", true),
            Keyword::integer("COUNT", 42),
        ];

        let result = writer.write_primary_image(&data, &dimensions, &keywords);
        assert!(result.is_ok());
    }

    #[test]
    fn write_primary_image_empty_data() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<f32> = vec![];
        let dimensions = vec![0];

        let result = writer.write_primary_image(&data, &dimensions, &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn format_card_end_keyword() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let end_keyword = Keyword::new("END".to_string());

        let card = writer.format_card(&end_keyword).unwrap();
        assert_eq!(&card[0..3], b"END");
        assert_eq!(card[3], b' ');
    }

    #[test]
    fn format_card_logical_value() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let keyword = Keyword::logical("SIMPLE", true);

        let card = writer.format_card(&keyword).unwrap();
        assert_eq!(&card[0..6], b"SIMPLE");
        assert_eq!(card[8], b'=');
        assert_eq!(card[9], b' ');
        assert!(String::from_utf8_lossy(&card[10..30]).contains('T'));
    }

    #[test]
    fn format_card_integer_value() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let keyword = Keyword::integer("NAXIS", 2);

        let card = writer.format_card(&keyword).unwrap();
        assert_eq!(&card[0..5], b"NAXIS");
        assert_eq!(card[8], b'=');
        assert!(String::from_utf8_lossy(&card[10..30]).contains('2'));
    }

    #[test]
    fn format_card_real_value() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let keyword = Keyword::real("EXPTIME", 30.5);

        let card = writer.format_card(&keyword).unwrap();
        assert_eq!(&card[0..7], b"EXPTIME");
        assert_eq!(card[8], b'=');
        assert!(String::from_utf8_lossy(&card[10..30]).contains("30.5"));
    }

    #[test]
    fn format_card_string_value() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let keyword = Keyword::string("OBJECT", "M31");

        let card = writer.format_card(&keyword).unwrap();
        assert_eq!(&card[0..6], b"OBJECT");
        assert_eq!(card[8], b'=');
        assert!(String::from_utf8_lossy(&card[10..30]).contains("'M31"));
    }

    #[test]
    fn format_card_string_truncation() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let long_string = "This is a very long string that exceeds 18 characters";
        let keyword = Keyword::string("LONGNAME", long_string);

        let card = writer.format_card(&keyword).unwrap();
        let value_section = String::from_utf8_lossy(&card[10..30]);
        assert!(value_section.len() <= 20);
    }

    #[test]
    fn format_card_with_comment() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let mut keyword = Keyword::integer("BITPIX", 16);
        keyword.comment = Some("Bits per pixel".to_string());

        let card = writer.format_card(&keyword).unwrap();
        assert_eq!(card[31], b'/');
        assert!(String::from_utf8_lossy(&card[33..]).contains("Bits per pixel"));
    }

    #[test]
    fn format_card_complex_value() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let mut keyword = Keyword::new("COMPLEX".to_string());
        keyword.value = Some(KeywordValue::Complex(1.0, 2.0));

        let card = writer.format_card(&keyword).unwrap();
        let value_section = String::from_utf8_lossy(&card[10..30]);
        assert!(value_section.contains("1"));
        assert!(value_section.contains("2"));
    }

    #[test]
    fn write_header_padding() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::new("END".to_string()));

        let result = writer.write_header(&header);
        assert!(result.is_ok());
    }

    #[test]
    fn write_compressed_image_i16() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<i16> = (0..64).collect();
        let dimensions = vec![8, 8];
        let tile_size = (8, 8);

        let result = writer.write_compressed_image(
            &data,
            &dimensions,
            tile_size,
            CompressionAlgorithm::Rice,
            &[],
        );
        assert!(result.is_ok());
    }

    #[test]
    fn write_compressed_image_i32() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<i32> = (0..64).collect();
        let dimensions = vec![8, 8];
        let tile_size = (4, 4);

        let result = writer.write_compressed_image(
            &data,
            &dimensions,
            tile_size,
            CompressionAlgorithm::Rice,
            &[],
        );
        assert!(result.is_ok());
    }

    #[test]
    fn write_compressed_image_with_keywords() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<i16> = (0..100).collect();
        let dimensions = vec![10, 10];
        let tile_size = (10, 10);

        let keywords = vec![
            Keyword::string("OBJECT", "TestImage"),
            Keyword::integer("EXPTIME", 60),
        ];

        let result = writer.write_compressed_image(
            &data,
            &dimensions,
            tile_size,
            CompressionAlgorithm::Rice,
            &keywords,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn extract_tile_full_tile() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();

        let data: Vec<u8> = (0..16).collect();
        let params = TileExtractionParams {
            width: 4,
            col: 0,
            row: 0,
            tile_width: 2,
            tile_height: 2,
            bytes_per_pixel: 1,
        };
        let tile = writer.extract_tile(&data, &params);

        assert_eq!(tile, vec![0, 1, 4, 5]);
    }

    #[test]
    fn extract_tile_offset() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();

        let data: Vec<u8> = (0..16).collect();
        let params = TileExtractionParams {
            width: 4,
            col: 2,
            row: 2,
            tile_width: 2,
            tile_height: 2,
            bytes_per_pixel: 1,
        };
        let tile = writer.extract_tile(&data, &params);

        assert_eq!(tile, vec![10, 11, 14, 15]);
    }

    #[test]
    fn build_heap() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();

        let tiles = vec![vec![1u8, 2, 3], vec![4u8, 5], vec![6u8]];
        let heap = writer.build_heap(&tiles);

        assert_eq!(heap, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn add_compression_keywords_sets_zimage() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let mut header = Header::new();
        let dimensions = vec![100, 100];

        writer.add_compression_keywords(
            &mut header,
            &dimensions,
            (32, 32),
            16,
            CompressionAlgorithm::Rice,
        );

        let zimage = header
            .get_keyword_value("ZIMAGE")
            .and_then(|v| v.as_logical());
        assert_eq!(zimage, Some(true));

        let zcmptype = header
            .get_keyword_value("ZCMPTYPE")
            .and_then(|v| v.as_string());
        assert_eq!(zcmptype, Some("RICE_1"));

        let zbitpix = header
            .get_keyword_value("ZBITPIX")
            .and_then(|v| v.as_integer());
        assert_eq!(zbitpix, Some(16));
    }

    #[test]
    fn compressed_image_roundtrip_i16() {
        use crate::fits::io::reader::FitsFile;
        use crate::fits::Hdu;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let original_data: Vec<i16> = (0..64).collect();
        let dimensions = vec![8, 8];
        let tile_size = (8, 8);

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_compressed_image(
                    &original_data,
                    &dimensions,
                    tile_size,
                    CompressionAlgorithm::Rice,
                    &[],
                )
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        assert_eq!(fits.num_hdus(), 2);

        let hdu = fits.read_hdu(1).unwrap();
        if let Hdu::BinaryTable(table) = hdu {
            assert!(table.is_compressed_image());
            assert_eq!(table.compression_algorithm(), Some("RICE_1"));
            assert_eq!(table.get_bits_per_pixel().unwrap(), 16);
            assert_eq!(table.get_tile_dimensions().unwrap(), (8, 8));
        } else {
            panic!("Expected binary table HDU");
        }
    }

    #[test]
    fn compressed_image_roundtrip_i32() {
        use crate::fits::io::reader::FitsFile;
        use crate::fits::Hdu;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let original_data: Vec<i32> = (0..100).collect();
        let dimensions = vec![10, 10];
        let tile_size = (5, 5);

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_compressed_image(
                    &original_data,
                    &dimensions,
                    tile_size,
                    CompressionAlgorithm::Rice,
                    &[],
                )
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let hdu = fits.read_hdu(1).unwrap();

        if let Hdu::BinaryTable(table) = hdu {
            assert!(table.is_compressed_image());
            assert_eq!(table.get_bits_per_pixel().unwrap(), 32);
        } else {
            panic!("Expected binary table HDU");
        }
    }

    #[test]
    fn compressed_image_with_multiple_tiles() {
        use crate::fits::io::reader::FitsFile;
        use crate::fits::Hdu;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let original_data: Vec<i16> = (0..256).collect();
        let dimensions = vec![16, 16];
        let tile_size = (4, 4);

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_compressed_image(
                    &original_data,
                    &dimensions,
                    tile_size,
                    CompressionAlgorithm::Rice,
                    &[],
                )
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let hdu = fits.read_hdu(1).unwrap();

        if let Hdu::BinaryTable(table) = hdu {
            assert!(table.is_compressed_image());
            let nrows = table.number_of_rows().unwrap();
            assert_eq!(nrows, 16);
        } else {
            panic!("Expected binary table HDU");
        }
    }

    #[test]
    fn write_primary_image_with_checksum_creates_valid_file() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..100).collect();
        let dimensions = vec![10, 10];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image_with_checksum(&data, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        assert_eq!(fits.num_hdus(), 1);

        let header = fits.get_header(0).unwrap();
        assert!(header.get_keyword_value("DATASUM").is_some());
        assert!(header.get_keyword_value("CHECKSUM").is_some());
    }

    #[test]
    fn write_with_checksum_validates_on_read() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..64).collect();
        let dimensions = vec![8, 8];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image_with_checksum(&data, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let result = fits.validate_hdu_checksum(0);
        assert!(
            result.is_ok(),
            "Checksum validation call failed: {:?}",
            result.err()
        );

        let checksum_result = result.unwrap();
        assert!(checksum_result.has_checksums());
        assert!(checksum_result.datasum_valid == Some(true));
    }

    #[test]
    fn write_checksum_empty_data() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<u8> = vec![];
        let dimensions = vec![0];

        let mut writer = FitsWriter::create(path).unwrap();
        let result = writer.write_primary_image_with_checksum(&data, &dimensions, &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn format_checksum_card_correct_format() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let checksum_value = "1234567890ABCDEF";

        let card = writer.format_checksum_card(checksum_value).unwrap();
        assert_eq!(&card[0..8], b"CHECKSUM");
        assert_eq!(card[8], b'=');
        assert_eq!(card[9], b' ');
        assert!(String::from_utf8_lossy(&card).contains(checksum_value));
    }

    #[test]
    fn build_primary_header_includes_required_keywords() {
        let writer = FitsWriter::create(NamedTempFile::new().unwrap().path()).unwrap();
        let dimensions = vec![10, 20];

        let header = writer.build_primary_header::<i16>(&dimensions, &[]);

        assert!(header.get_keyword_value("SIMPLE").is_some());
        assert!(header.get_keyword_value("BITPIX").is_some());
        assert!(header.get_keyword_value("NAXIS").is_some());
        assert!(header.get_keyword_value("NAXIS1").is_some());
        assert!(header.get_keyword_value("NAXIS2").is_some());
    }

    #[test]
    fn roundtrip_2d_image_preserves_orientation() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let original: Vec<i16> = (0..100).collect();
        let dimensions = vec![10, 10];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image(&original, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let (_, read_back): (_, Vec<i16>) = fits.primary_hdu_with_data().unwrap();

        assert_eq!(original, read_back);
    }

    #[test]
    fn roundtrip_gradient_image_correct_orientation() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let mut gradient: Vec<i16> = Vec::with_capacity(25);
        for row in 0..5 {
            for col in 0..5 {
                gradient.push((row * 10 + col) as i16);
            }
        }
        let dimensions = vec![5, 5];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image(&gradient, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let (_, read_back): (_, Vec<i16>) = fits.primary_hdu_with_data().unwrap();

        assert_eq!(read_back[0], 0);
        assert_eq!(read_back[4], 4);
        assert_eq!(read_back[20], 40);
        assert_eq!(read_back[24], 44);
        assert_eq!(gradient, read_back);
    }
}
