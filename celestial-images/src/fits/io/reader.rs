use crate::fits::hdu::{Hdu, HduTrait};
use crate::fits::header::{Header, HeaderParser};
use crate::fits::util::buffer_pool::BufferPool;
use crate::fits::util::checksum;
use crate::fits::{FitsError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

const FITS_BLOCK_SIZE: usize = 2880;

#[derive(Debug, Clone)]
pub struct ChecksumResult {
    pub datasum_valid: Option<bool>,
    pub checksum_valid: Option<bool>,
}

impl ChecksumResult {
    pub fn has_checksums(&self) -> bool {
        self.datasum_valid.is_some() || self.checksum_valid.is_some()
    }

    pub fn all_valid(&self) -> bool {
        self.datasum_valid.unwrap_or(true) && self.checksum_valid.unwrap_or(true)
    }
}

fn validate_fits_block_alignment(size: u64, context: &str) -> Result<()> {
    if !size.is_multiple_of(FITS_BLOCK_SIZE as u64) {
        return Err(FitsError::InvalidFormat(format!(
            "FITS {context} not aligned to {FITS_BLOCK_SIZE}-byte blocks: {size}"
        )));
    }
    Ok(())
}

#[derive(Debug)]
pub struct FitsFile<R> {
    reader: R,
    hdus: Vec<HduInfo>,
    current_hdu: usize,
    header_cache: HashMap<usize, Header>,
    buffer_pool: BufferPool,
}

pub struct FitsReader {
    inner: BufReader<File>,
}

#[derive(Debug, Clone)]
pub struct HduInfo {
    pub index: usize,
    pub header_start: u64,
    pub header_size: usize,
    pub data_start: u64,
    pub data_size: usize,
}

impl FitsFile<FitsReader> {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let reader = FitsReader::new(file);
        let mut fits = FitsFile {
            reader,
            hdus: Vec::new(),
            current_hdu: 0,
            header_cache: HashMap::new(),
            buffer_pool: BufferPool::default(),
        };
        fits.scan_hdus()?;
        Ok(fits)
    }
}

impl<R: Read + Seek> FitsFile<R> {
    pub fn new(reader: R) -> Result<Self> {
        let mut fits = FitsFile {
            reader,
            hdus: Vec::new(),
            current_hdu: 0,
            header_cache: HashMap::new(),
            buffer_pool: BufferPool::default(),
        };
        fits.scan_hdus()?;
        Ok(fits)
    }

    pub fn num_hdus(&self) -> usize {
        self.hdus.len()
    }

    pub fn hdu_info(&self, index: usize) -> Option<&HduInfo> {
        self.hdus.get(index)
    }

    pub fn current_hdu_index(&self) -> usize {
        self.current_hdu
    }

    pub fn move_to_hdu(&mut self, index: usize) -> Result<()> {
        if index >= self.hdus.len() {
            return Err(FitsError::HduNotFound(index));
        }
        self.current_hdu = index;
        Ok(())
    }

    pub fn primary_hdu(&mut self) -> Result<Hdu> {
        self.move_to_hdu(0)?;
        self.read_current_hdu()
    }

    pub fn read_hdu(&mut self, index: usize) -> Result<Hdu> {
        self.move_to_hdu(index)?;
        self.read_current_hdu()
    }

    pub fn get_header(&mut self, index: usize) -> Result<Header> {
        let hdu = self.read_hdu(index)?;
        Ok(hdu.header().clone())
    }

    pub fn get_header_value(&mut self, index: usize, keyword: &str) -> Result<Option<String>> {
        let header = self.get_header(index)?;
        Ok(header.get_keyword_value(keyword).map(|v| v.to_string()))
    }

    pub fn clear_header_cache(&mut self) {
        self.header_cache.clear();
    }

    pub fn clear_buffer_pool(&mut self) {
        self.buffer_pool.clear();
    }

    pub fn validate_hdu_checksum(&mut self, index: usize) -> Result<ChecksumResult> {
        let hdu_info = self
            .hdus
            .get(index)
            .ok_or(FitsError::HduNotFound(index))?
            .clone();
        let (header_bytes, data_bytes) = self.read_hdu_raw_bytes(&hdu_info)?;
        let header = HeaderParser::parse_header(&header_bytes)?;

        self.validate_checksums(&header, &header_bytes, &data_bytes)
    }

    fn read_hdu_raw_bytes(&mut self, hdu_info: &HduInfo) -> Result<(Vec<u8>, Vec<u8>)> {
        self.reader.seek(SeekFrom::Start(hdu_info.header_start))?;

        let mut header_bytes = vec![0u8; hdu_info.header_size];
        self.reader.read_exact(&mut header_bytes)?;

        let mut data_bytes = vec![0u8; hdu_info.data_size];
        if !data_bytes.is_empty() {
            self.reader.read_exact(&mut data_bytes)?;
        }

        Ok((header_bytes, data_bytes))
    }

    fn validate_checksums(
        &self,
        header: &Header,
        header_bytes: &[u8],
        data_bytes: &[u8],
    ) -> Result<ChecksumResult> {
        let datasum_result = self.validate_datasum_if_present(header, data_bytes)?;
        let checksum_result =
            self.validate_checksum_if_present(header, header_bytes, data_bytes)?;

        Ok(ChecksumResult {
            datasum_valid: datasum_result,
            checksum_valid: checksum_result,
        })
    }

    fn validate_datasum_if_present(
        &self,
        header: &Header,
        data_bytes: &[u8],
    ) -> Result<Option<bool>> {
        let datasum_value = match header.get_keyword_value("DATASUM") {
            Some(v) => v.as_string(),
            None => return Ok(None),
        };

        let expected = datasum_value.unwrap_or_default().to_string();
        let computed = checksum::calculate_datasum(data_bytes);

        if !checksum::verify_datasum(data_bytes, &expected) {
            return Err(FitsError::DatasumMismatch { expected, computed });
        }

        Ok(Some(true))
    }

    fn validate_checksum_if_present(
        &self,
        header: &Header,
        header_bytes: &[u8],
        data_bytes: &[u8],
    ) -> Result<Option<bool>> {
        let checksum_value = match header.get_keyword_value("CHECKSUM") {
            Some(v) => v.as_string(),
            None => return Ok(None),
        };

        if checksum_value.is_none() {
            return Ok(None);
        }

        let valid = checksum::verify_hdu_checksum(header_bytes, data_bytes);
        Ok(Some(valid))
    }

    pub fn get_image_info(&mut self, index: usize) -> Result<(Vec<usize>, crate::core::BitPix)> {
        let hdu = self.read_hdu(index)?;
        let (dimensions, bitpix) = match hdu {
            Hdu::Primary(primary) => (primary.data_dimensions(), primary.bitpix()),
            Hdu::Image(image) => (image.data_dimensions(), image.bitpix()),
            _ => {
                return Err(FitsError::InvalidFormat(
                    "HDU type does not support image data".to_string(),
                ))
            }
        };

        let bitpix = bitpix.ok_or_else(|| FitsError::KeywordNotFound {
            keyword: "BITPIX".to_string(),
        })?;

        Ok((dimensions, bitpix))
    }

    pub fn primary_hdu_with_data<T>(&mut self) -> Result<(crate::fits::hdu::PrimaryHdu, Vec<T>)>
    where
        T: crate::fits::data::array::DataArray,
    {
        let primary = self.primary_hdu()?;
        match primary {
            Hdu::Primary(hdu) => {
                let data = hdu.read_data(&mut self.reader)?;
                Ok((hdu, data))
            }
            _ => unreachable!("primary_hdu always returns Primary"),
        }
    }

    pub fn read_hdu_with_data<T>(&mut self, index: usize) -> Result<(Hdu, Vec<T>)>
    where
        T: crate::fits::data::array::DataArray,
    {
        let hdu = self.read_hdu(index)?;
        let data = match &hdu {
            Hdu::Primary(primary) => primary.read_data(&mut self.reader)?,
            Hdu::Image(image) => image.read_data(&mut self.reader)?,
            _ => {
                return Err(FitsError::InvalidFormat(
                    "HDU type does not support image data reading".to_string(),
                ))
            }
        };
        Ok((hdu, data))
    }

    fn read_current_hdu(&mut self) -> Result<Hdu> {
        let hdu_info = self.hdus[self.current_hdu].clone();

        let header = if let Some(cached_header) = self.header_cache.get(&self.current_hdu) {
            cached_header.clone()
        } else {
            self.reader.seek(SeekFrom::Start(hdu_info.header_start))?;

            let mut header_data = self.buffer_pool.get(hdu_info.header_size);
            self.reader.read_exact(&mut header_data)?;

            let parsed_header = HeaderParser::parse_header(&header_data)?;
            self.header_cache
                .insert(self.current_hdu, parsed_header.clone());

            self.buffer_pool.return_buffer(header_data);

            parsed_header
        };

        if self.current_hdu == 0 {
            if !header.is_primary() {
                return Err(FitsError::InvalidFormat(
                    "First HDU must be a primary HDU".to_string(),
                ));
            }
            Ok(Hdu::Primary(crate::fits::hdu::PrimaryHdu::new(
                header, hdu_info,
            )))
        } else {
            if !header.is_extension() {
                return Err(FitsError::InvalidFormat(
                    "Non-primary HDUs must be extensions".to_string(),
                ));
            }
            self.create_extension_hdu(header, hdu_info)
        }
    }

    fn scan_hdus(&mut self) -> Result<()> {
        self.reader.seek(SeekFrom::Start(0))?;
        let mut position = 0u64;
        let mut hdu_index = 0;

        loop {
            let hdu_info = match self.scan_single_hdu(position, hdu_index) {
                Ok(info) => info,
                Err(_) => {
                    if !self.hdus.is_empty() {
                        break;
                    }
                    return Err(FitsError::InvalidFormat(
                        "Failed to scan primary HDU".to_string(),
                    ));
                }
            };

            let should_stop = self.should_stop_scanning(&hdu_info, position)?;
            position = Self::calculate_next_position(&hdu_info);

            self.hdus.push(hdu_info);

            if should_stop {
                break;
            }
            hdu_index += 1;
        }

        Ok(())
    }

    fn scan_single_hdu(&mut self, position: u64, hdu_index: usize) -> Result<HduInfo> {
        self.reader.seek(SeekFrom::Start(position))?;

        let header_size = self.determine_header_size(position)?;
        validate_fits_block_alignment(header_size as u64, "header")?;
        let header_start = position;

        let mut header_data = self.buffer_pool.get(header_size);
        self.reader.read_exact(&mut header_data)?;

        let header = HeaderParser::parse_header(&header_data)?;
        self.buffer_pool.return_buffer(header_data);

        let data_start = position + header_size as u64;
        validate_fits_block_alignment(data_start, "data start position")?;
        let data_size = self.calculate_data_size(&header)?;
        if data_size > 0 {
            validate_fits_block_alignment(data_size as u64, "data size")?;
        }

        Ok(HduInfo {
            index: hdu_index,
            header_start,
            header_size,
            data_start,
            data_size,
        })
    }

    fn calculate_next_position(hdu_info: &HduInfo) -> u64 {
        let position = hdu_info.data_start + hdu_info.data_size as u64;
        Self::align_to_block(position)
    }

    fn should_stop_scanning(&mut self, hdu_info: &HduInfo, position: u64) -> Result<bool> {
        if hdu_info.data_size == 0 && hdu_info.index > 0 {
            return Ok(true);
        }

        if hdu_info.index == 0 {
            let mut header_data = self.buffer_pool.get(hdu_info.header_size);
            let current_pos = self.reader.stream_position()?;
            self.reader.seek(SeekFrom::Start(hdu_info.header_start))?;
            self.reader.read_exact(&mut header_data)?;
            self.reader.seek(SeekFrom::Start(current_pos))?;

            let header = HeaderParser::parse_header(&header_data)?;
            self.buffer_pool.return_buffer(header_data);

            if let Some(extend_value) = header.get_keyword_value("EXTEND") {
                if let Some(logical_val) = extend_value.as_logical() {
                    if !logical_val {
                        return Ok(true);
                    }
                }
            }
        }

        self.check_end_of_file(position)
    }

    fn check_end_of_file(&mut self, position: u64) -> Result<bool> {
        if self.reader.seek(SeekFrom::Start(position)).is_err() {
            return Ok(true);
        }

        let mut test_buf = [0u8; 8];
        if self.reader.read_exact(&mut test_buf).is_err() {
            return Ok(true);
        }

        Ok(test_buf.iter().all(|&b| b == 0))
    }

    fn determine_header_size(&mut self, start_position: u64) -> Result<usize> {
        const MAX_HEADER_BLOCKS: usize = 1000;
        const CARD_SIZE: usize = 80;

        let current_pos = self.reader.stream_position()?;
        self.reader.seek(SeekFrom::Start(start_position))?;

        let mut block_buffer = vec![0u8; FITS_BLOCK_SIZE];
        let mut blocks_read = 0;

        loop {
            if blocks_read >= MAX_HEADER_BLOCKS {
                self.reader.seek(SeekFrom::Start(current_pos))?;
                return Err(FitsError::InvalidFormat(format!(
                    "Header exceeds maximum size of {} blocks ({} bytes)",
                    MAX_HEADER_BLOCKS,
                    MAX_HEADER_BLOCKS * FITS_BLOCK_SIZE
                )));
            }

            self.reader.read_exact(&mut block_buffer).map_err(|e| {
                let _ = self.reader.seek(SeekFrom::Start(current_pos));
                FitsError::InvalidFormat(format!(
                    "Unexpected end of file while scanning header at block {}: {}",
                    blocks_read, e
                ))
            })?;

            blocks_read += 1;

            for chunk in block_buffer.chunks_exact(CARD_SIZE) {
                let keyword_part = std::str::from_utf8(&chunk[0..8])
                    .map_err(|_| FitsError::InvalidFormat("Invalid UTF-8 in header".to_string()))?;

                if keyword_part.trim() == "END" {
                    self.reader.seek(SeekFrom::Start(current_pos))?;
                    return Ok(blocks_read * FITS_BLOCK_SIZE);
                }
            }
        }
    }

    fn calculate_data_size(&self, header: &Header) -> Result<usize> {
        let naxis = header
            .get_keyword_value("NAXIS")
            .and_then(|v| v.as_integer())
            .unwrap_or(0) as usize;

        if naxis == 0 {
            return Ok(0);
        }

        let bitpix = header
            .get_keyword_value("BITPIX")
            .and_then(|v| v.as_integer())
            .ok_or_else(|| FitsError::KeywordNotFound {
                keyword: "BITPIX".to_string(),
            })? as i32;

        let bytes_per_pixel = crate::core::BitPix::from_value(bitpix)
            .ok_or(FitsError::InvalidBitPix(bitpix))?
            .bytes_per_pixel();

        let mut total_pixels = 1usize;
        for i in 1..=naxis {
            let axis_name = format!("NAXIS{}", i);
            let axis_size = header
                .get_keyword_value(&axis_name)
                .and_then(|v| v.as_integer())
                .unwrap_or(1) as usize;
            total_pixels = total_pixels
                .checked_mul(axis_size)
                .ok_or_else(|| FitsError::InvalidFormat("Data dimensions too large".to_string()))?;
        }

        let data_size = total_pixels * bytes_per_pixel;
        Ok(Self::align_to_block(data_size as u64) as usize)
    }

    fn align_to_block(size: u64) -> u64 {
        size.div_ceil(FITS_BLOCK_SIZE as u64) * FITS_BLOCK_SIZE as u64
    }

    fn create_extension_hdu(&self, header: Header, hdu_info: HduInfo) -> Result<Hdu> {
        let xtension = header
            .get_keyword_value("XTENSION")
            .and_then(|v| v.as_string())
            .ok_or_else(|| FitsError::KeywordNotFound {
                keyword: "XTENSION".to_string(),
            })?;

        match xtension {
            "IMAGE" => Ok(Hdu::Image(Box::new(crate::fits::hdu::ImageHdu::new(
                header, hdu_info,
            )))),
            "TABLE" => Ok(Hdu::AsciiTable(crate::fits::hdu::AsciiTableHdu::new(
                header, hdu_info,
            ))),
            "BINTABLE" => Ok(Hdu::BinaryTable(crate::fits::hdu::BinaryTableHdu::new(
                header, hdu_info,
            ))),
            "A3DTABLE" => Ok(Hdu::RandomGroups(crate::fits::hdu::RandomGroupsHdu::new(
                header, hdu_info,
            ))),
            _ => Err(FitsError::InvalidFormat(format!(
                "Unsupported XTENSION type: {}",
                xtension
            ))),
        }
    }
}

impl FitsReader {
    pub fn new(file: File) -> Self {
        Self {
            inner: BufReader::new(file),
        }
    }
}

impl Read for FitsReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.inner.read(buf)
    }
}

impl Seek for FitsReader {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        self.inner.seek(pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Keyword, KeywordValue};
    use crate::test_utils::*;
    use std::io::Cursor;

    #[test]
    fn fits_file_new_with_valid_data() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.num_hdus(), 1);
    }

    #[test]
    fn fits_file_new_with_invalid_data() {
        let invalid_data = vec![0x00, 0x01, 0x02];
        let cursor = Cursor::new(invalid_data);
        let result = FitsFile::new(cursor);

        assert!(result.is_err());
    }

    #[test]
    fn scan_single_hdu_valid() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.hdus.len(), 1);

        let hdu_info = &fits_file.hdus[0];
        assert_eq!(hdu_info.index, 0);
        assert_eq!(hdu_info.header_start, 0);
        assert!(hdu_info.header_size > 0);
        assert!(hdu_info.data_start >= hdu_info.header_size as u64);
    }

    #[test]
    fn calculate_next_position() {
        let hdu_info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 100,
        };

        let next_pos = FitsFile::<Cursor<Vec<u8>>>::calculate_next_position(&hdu_info);

        assert_eq!(next_pos % FITS_BLOCK_SIZE as u64, 0);
        assert!(next_pos >= hdu_info.data_start + hdu_info.data_size as u64);
    }

    #[test]
    fn calculate_next_position_already_aligned() {
        let hdu_info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 2880,
        };

        let next_pos = FitsFile::<Cursor<Vec<u8>>>::calculate_next_position(&hdu_info);
        assert_eq!(next_pos, 2880 + 2880);
    }

    #[test]
    fn should_stop_scanning_conditions() {
        let fits_data = MockFitsBuilder::new()
            .card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", "8", "Bits per pixel")
            .card("NAXIS", "0", "Number of axes")
            .card("EXTEND", "T", "Has extensions")
            .build_memory();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let hdu_zero_data = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 0,
        };

        let should_stop = fits_file
            .should_stop_scanning(&hdu_zero_data, 5760)
            .unwrap();
        assert!(should_stop);

        let primary_zero_data = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let should_stop = fits_file
            .should_stop_scanning(&primary_zero_data, 0)
            .unwrap();
        assert!(!should_stop);
    }

    #[test]
    fn align_to_block_various_sizes() {
        let test_cases = [
            (0, 0),
            (1, FITS_BLOCK_SIZE as u64),
            (2880, 2880),
            (2881, 5760),
            (5760, 5760),
            (5761, 8640),
        ];

        for (input, expected) in test_cases {
            let aligned = FitsFile::<std::io::Cursor<Vec<u8>>>::align_to_block(input);
            assert_eq!(aligned, expected);
            assert_eq!(aligned % FITS_BLOCK_SIZE as u64, 0);
        }
    }

    #[test]
    fn primary_hdu_access() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let primary = fits_file.primary_hdu().unwrap();
        assert!(matches!(primary, Hdu::Primary(_)));
    }

    #[test]
    fn read_hdu_by_index() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let hdu = fits_file.read_hdu(0).unwrap();
        assert!(matches!(hdu, Hdu::Primary(_)));

        let result = fits_file.read_hdu(999);
        assert!(result.is_err());
    }

    #[test]
    fn check_end_of_file_conditions() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let is_eof = fits_file.check_end_of_file(0).unwrap();
        assert!(!is_eof);

        let is_eof = fits_file.check_end_of_file(u64::MAX).unwrap();
        assert!(is_eof);
    }

    #[test]
    fn fits_file_with_data() {
        let test_data: Vec<i16> = (0..100).collect();
        let fits_data = create_image_fits(16, &[10, 10], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.num_hdus(), 1);

        let (_primary, data): (_, Vec<i16>) = fits_file.primary_hdu_with_data().unwrap();
        assert_eq!(data.len(), 100);
        assert_eq!(data[0], 90);
        assert_eq!(data[9], 99);
        assert_eq!(data[90], 0);
        assert_eq!(data[99], 9);
    }

    #[test]
    fn error_propagation() {
        let malformed = create_malformed_fits("truncated");
        let cursor = Cursor::new(malformed);
        let result = FitsFile::new(cursor);

        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn multiple_hdu_handling() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.num_hdus(), 1);
    }

    #[test]
    fn fits_file_open_success() {
        use std::fs::File;
        use std::io::Write;

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_fits_open.fits");

        let fits_data = create_minimal_fits();
        let mut file = File::create(&temp_path).unwrap();
        file.write_all(&fits_data).unwrap();

        let result = FitsFile::open(&temp_path);
        assert!(result.is_ok());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn hdu_info_accessor() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let hdu_info = fits_file.hdu_info(0);
        assert!(hdu_info.is_some());
        assert_eq!(hdu_info.unwrap().index, 0);

        let hdu_info_invalid = fits_file.hdu_info(999);
        assert!(hdu_info_invalid.is_none());
    }

    #[test]
    fn current_hdu_index_accessor() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.current_hdu_index(), 0);
    }

    #[test]
    fn move_to_hdu_success_and_error() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.move_to_hdu(0);
        assert!(result.is_ok());
        assert_eq!(fits_file.current_hdu_index(), 0);

        let result = fits_file.move_to_hdu(999);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::HduNotFound(999)));
    }

    #[test]
    fn get_header_success() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let header = fits_file.get_header(0).unwrap();
        assert!(header.is_primary());
    }

    #[test]
    fn get_header_value_success() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let simple_value = fits_file.get_header_value(0, "SIMPLE").unwrap();
        assert!(simple_value.is_some());
        assert_eq!(simple_value.unwrap(), "T");

        let nonexistent_value = fits_file.get_header_value(0, "NONEXIST").unwrap();
        assert!(nonexistent_value.is_none());
    }

    #[test]
    fn get_image_info_success() {
        let test_data: Vec<i16> = (0..100).collect();
        let fits_data = create_image_fits(16, &[10, 10], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let (dimensions, bitpix) = fits_file.get_image_info(0).unwrap();
        assert_eq!(dimensions, vec![10, 10]);
        assert_eq!(bitpix, crate::core::BitPix::I16);
    }

    #[test]
    fn primary_hdu_with_data_success() {
        let test_data: Vec<i16> = (0..100).collect();
        let fits_data = create_image_fits(16, &[10, 10], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let (hdu, data): (_, Vec<i16>) = fits_file.primary_hdu_with_data().unwrap();
        assert!(matches!(hdu, crate::fits::hdu::PrimaryHdu { .. }));
        assert_eq!(data.len(), 100);
    }

    #[test]
    fn read_hdu_with_data_success() {
        let test_data: Vec<i16> = (0..100).collect();
        let fits_data = create_image_fits(16, &[10, 10], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let (hdu, data): (_, Vec<i16>) = fits_file.read_hdu_with_data(0).unwrap();
        match hdu {
            Hdu::Primary(_) => {
                assert_eq!(data.len(), 100);
            }
            _ => panic!("Expected primary HDU"),
        }
    }

    #[test]
    fn read_current_hdu_primary() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let hdu = fits_file.read_current_hdu().unwrap();
        assert!(matches!(hdu, Hdu::Primary(_)));
    }

    #[test]
    fn read_current_hdu_invalid_primary() {
        let fits_data = MockFitsBuilder::new()
            .card("XTENSION", "IMAGE", "Image extension")
            .card("BITPIX", "8", "Bits per pixel")
            .card("NAXIS", "0", "Number of axes")
            .build_memory();

        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();
        let result = fits_file.read_current_hdu();
        assert!(result.is_err());
    }

    #[test]
    fn read_current_hdu_extension() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let extension_hdu_info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 0,
        };

        fits_file.hdus.push(extension_hdu_info);
        fits_file.current_hdu = 1;

        let result = fits_file.read_current_hdu();
        assert!(result.is_err());
    }

    #[test]
    fn scan_hdus_error_handling() {
        let invalid_data = vec![0u8; 100];
        let cursor = Cursor::new(invalid_data);
        let result = FitsFile::new(cursor);
        assert!(result.is_err());
    }

    #[test]
    fn scan_single_hdu_components() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.hdus.len(), 1);
        let hdu_info = &fits_file.hdus[0];

        assert_eq!(hdu_info.index, 0);
        assert_eq!(hdu_info.header_start, 0);
        assert!(hdu_info.header_size > 0);
        assert!(hdu_info.data_start >= hdu_info.header_size as u64);
    }

    #[test]
    fn should_stop_scanning_extend_false() {
        let fits_data = MockFitsBuilder::new()
            .card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", "8", "Bits per pixel")
            .card("NAXIS", "0", "Number of axes")
            .card("EXTEND", "F", "No extensions")
            .build_memory();

        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let primary_hdu = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let should_stop = fits_file.should_stop_scanning(&primary_hdu, 0).unwrap();
        assert!(should_stop);
    }

    #[test]
    fn check_end_of_file_seek_error() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let is_eof = fits_file.check_end_of_file(u64::MAX).unwrap();
        assert!(is_eof);
    }

    #[test]
    fn check_end_of_file_read_error() {
        let fits_data = create_minimal_fits();
        let file_size = fits_data.len() as u64;
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();
        let is_eof = fits_file.check_end_of_file(file_size).unwrap();
        assert!(is_eof);
    }

    #[test]
    fn determine_header_size_success() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let header_size = fits_file.determine_header_size(0).unwrap();
        assert_eq!(header_size, 2880);
    }

    #[test]
    fn calculate_data_size_zero_naxis() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        fits_file.reader.seek(SeekFrom::Start(0)).unwrap();
        let mut header_data = vec![0u8; 2880];
        fits_file.reader.read_exact(&mut header_data).unwrap();
        let header = HeaderParser::parse_header(&header_data).unwrap();

        let data_size = fits_file.calculate_data_size(&header).unwrap();
        assert_eq!(data_size, 0);
    }

    #[test]
    fn calculate_data_size_missing_bitpix() {
        use crate::fits::header::{Keyword, KeywordValue};

        let mut header = Header::new();
        header.add_keyword(Keyword::new("NAXIS".to_string()).with_value(KeywordValue::Integer(2)));
        header
            .add_keyword(Keyword::new("NAXIS1".to_string()).with_value(KeywordValue::Integer(10)));
        header
            .add_keyword(Keyword::new("NAXIS2".to_string()).with_value(KeywordValue::Integer(10)));

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.calculate_data_size(&header);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            FitsError::KeywordNotFound { .. }
        ));
    }

    #[test]
    fn calculate_data_size_invalid_bitpix() {
        use crate::fits::header::{Keyword, KeywordValue};

        let mut header = Header::new();
        header.add_keyword(Keyword::new("NAXIS".to_string()).with_value(KeywordValue::Integer(2)));
        header
            .add_keyword(Keyword::new("NAXIS1".to_string()).with_value(KeywordValue::Integer(10)));
        header
            .add_keyword(Keyword::new("NAXIS2".to_string()).with_value(KeywordValue::Integer(10)));
        header
            .add_keyword(Keyword::new("BITPIX".to_string()).with_value(KeywordValue::Integer(999)));

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.calculate_data_size(&header);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidBitPix(999)));
    }

    #[test]
    fn calculate_data_size_overflow() {
        use crate::fits::header::{Keyword, KeywordValue};

        let mut header = Header::new();
        header.add_keyword(Keyword::new("NAXIS".to_string()).with_value(KeywordValue::Integer(2)));
        header.add_keyword(
            Keyword::new("NAXIS1".to_string()).with_value(KeywordValue::Integer(i64::MAX)),
        );
        header.add_keyword(
            Keyword::new("NAXIS2".to_string()).with_value(KeywordValue::Integer(i64::MAX)),
        );
        header.add_keyword(Keyword::new("BITPIX".to_string()).with_value(KeywordValue::Integer(8)));

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.calculate_data_size(&header);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn calculate_data_size_alignment() {
        use crate::fits::header::{Keyword, KeywordValue};

        let mut header = Header::new();
        header.add_keyword(Keyword::new("NAXIS".to_string()).with_value(KeywordValue::Integer(2)));
        header
            .add_keyword(Keyword::new("NAXIS1".to_string()).with_value(KeywordValue::Integer(10)));
        header
            .add_keyword(Keyword::new("NAXIS2".to_string()).with_value(KeywordValue::Integer(10)));
        header.add_keyword(Keyword::new("BITPIX".to_string()).with_value(KeywordValue::Integer(8)));

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let data_size = fits_file.calculate_data_size(&header).unwrap();
        assert_eq!(data_size % FITS_BLOCK_SIZE, 0);
        assert!(data_size >= 100);
    }

    #[test]
    fn fits_reader_trait_implementations() {
        use std::fs::File;
        use std::io::Write;

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_fits_reader_traits.fits");

        let fits_data = create_minimal_fits();
        let mut file = File::create(&temp_path).unwrap();
        file.write_all(&fits_data).unwrap();

        let file = File::open(&temp_path).unwrap();
        let mut fits_reader = FitsReader::new(file);

        let mut buffer = [0u8; 8];
        let bytes_read = fits_reader.read(&mut buffer).unwrap();
        assert_eq!(bytes_read, 8);

        let pos = fits_reader.seek(SeekFrom::Start(0)).unwrap();
        assert_eq!(pos, 0);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn comprehensive_workflow_test() {
        let test_data: Vec<i16> = (0..400).collect();
        let fits_data = create_image_fits(16, &[20, 20], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.num_hdus(), 1);
        assert!(fits_file.hdu_info(0).is_some());
        assert_eq!(fits_file.current_hdu_index(), 0);

        let header = fits_file.get_header(0).unwrap();
        assert!(header.is_primary());

        let simple_value = fits_file.get_header_value(0, "SIMPLE").unwrap();
        assert!(simple_value.is_some());

        let (dimensions, bitpix) = fits_file.get_image_info(0).unwrap();
        assert_eq!(dimensions, vec![20, 20]);
        assert_eq!(bitpix, crate::core::BitPix::I16);

        let (hdu, data): (_, Vec<i16>) = fits_file.read_hdu_with_data(0).unwrap();
        assert_eq!(data.len(), 400);
        assert!(matches!(hdu, Hdu::Primary(_)));
    }

    #[test]
    fn get_image_info_valid_header() {
        let fits_data = MockFitsBuilder::new()
            .card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", "16", "Bits per pixel")
            .card("NAXIS", "2", "Number of axes")
            .card("NAXIS1", "10", "Axis 1 size")
            .card("NAXIS2", "10", "Axis 2 size")
            .build_memory();

        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.get_image_info(0);
        assert!(result.is_ok());
        let (dimensions, bitpix) = result.unwrap();
        assert_eq!(dimensions, vec![10, 10]);
        assert_eq!(bitpix, crate::core::BitPix::I16);
    }

    #[test]
    fn scan_hdus_graceful_stop() {
        let fits_data = MockFitsBuilder::new()
            .card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", "8", "Bits per pixel")
            .card("NAXIS", "0", "Number of axes")
            .card("EXTEND", "F", "No extensions")
            .build_memory();

        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        assert_eq!(fits_file.num_hdus(), 1);
    }

    #[test]
    fn calculate_next_position_aligned_data() {
        let hdu_info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 2880,
        };

        let next_pos = FitsFile::<std::io::Cursor<Vec<u8>>>::calculate_next_position(&hdu_info);

        assert_eq!(next_pos, 5760);
        assert_eq!(next_pos % FITS_BLOCK_SIZE as u64, 0);
    }

    #[test]
    fn validate_fits_block_alignment_success() {
        assert!(validate_fits_block_alignment(0, "test").is_ok());
        assert!(validate_fits_block_alignment(2880, "test").is_ok());
        assert!(validate_fits_block_alignment(5760, "test").is_ok());
    }

    #[test]
    fn validate_fits_block_alignment_failure() {
        let result = validate_fits_block_alignment(1, "test");
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));

        let result = validate_fits_block_alignment(2879, "test");
        assert!(result.is_err());

        let result = validate_fits_block_alignment(2881, "test");
        assert!(result.is_err());
    }

    #[test]
    fn block_validation_prevents_corruption() {
        let mut bad_header = Header::new();
        bad_header
            .add_keyword(Keyword::new("NAXIS".to_string()).with_value(KeywordValue::Integer(1)));
        bad_header
            .add_keyword(Keyword::new("NAXIS1".to_string()).with_value(KeywordValue::Integer(13)));
        bad_header
            .add_keyword(Keyword::new("BITPIX".to_string()).with_value(KeywordValue::Integer(8)));

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.calculate_data_size(&bad_header);
        assert!(result.is_ok());
        let data_size = result.unwrap();
        assert_eq!(data_size % FITS_BLOCK_SIZE, 0);
    }

    #[test]
    fn create_extension_hdu_image_type() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "IMAGE"));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 10));
        header.add_keyword(Keyword::integer("NAXIS2", 10));

        let hdu_info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 200,
        };

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.create_extension_hdu(header, hdu_info);
        assert!(result.is_ok());
        assert!(matches!(result.unwrap(), Hdu::Image(_)));
    }

    #[test]
    fn create_extension_hdu_binary_table_type() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 80));
        header.add_keyword(Keyword::integer("NAXIS2", 100));
        header.add_keyword(Keyword::integer("TFIELDS", 3));

        let hdu_info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 8000,
        };

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.create_extension_hdu(header, hdu_info);
        assert!(result.is_ok());
        assert!(matches!(result.unwrap(), Hdu::BinaryTable(_)));
    }

    #[test]
    fn create_extension_hdu_ascii_table_type() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 80));
        header.add_keyword(Keyword::integer("NAXIS2", 100));
        header.add_keyword(Keyword::integer("TFIELDS", 4));

        let hdu_info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 8000,
        };

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.create_extension_hdu(header, hdu_info);
        assert!(result.is_ok());
        assert!(matches!(result.unwrap(), Hdu::AsciiTable(_)));
    }

    #[test]
    fn create_extension_hdu_unsupported_type() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "UNKNOWN"));

        let hdu_info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 0,
        };

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.create_extension_hdu(header, hdu_info);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn create_extension_hdu_missing_xtension() {
        let mut header = Header::new();
        header.add_keyword(Keyword::integer("BITPIX", 8));

        let hdu_info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 0,
        };

        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.create_extension_hdu(header, hdu_info);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            FitsError::KeywordNotFound { .. }
        ));
    }

    #[test]
    fn validate_hdu_checksum_no_checksums() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.validate_hdu_checksum(0);
        assert!(result.is_ok());

        let checksum_result = result.unwrap();
        assert!(!checksum_result.has_checksums());
        assert!(checksum_result.all_valid());
    }

    #[test]
    fn validate_hdu_checksum_invalid_index() {
        let fits_data = create_minimal_fits();
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.validate_hdu_checksum(999);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::HduNotFound(999)));
    }

    #[test]
    fn checksum_result_methods() {
        let result_with_checksums = ChecksumResult {
            datasum_valid: Some(true),
            checksum_valid: Some(true),
        };
        assert!(result_with_checksums.has_checksums());
        assert!(result_with_checksums.all_valid());

        let result_no_checksums = ChecksumResult {
            datasum_valid: None,
            checksum_valid: None,
        };
        assert!(!result_no_checksums.has_checksums());
        assert!(result_no_checksums.all_valid());

        let result_partial = ChecksumResult {
            datasum_valid: Some(true),
            checksum_valid: None,
        };
        assert!(result_partial.has_checksums());
        assert!(result_partial.all_valid());
    }

    #[test]
    fn pixinsight_large_header_fits_file() {
        let test_file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_data")
            .join("cygnus.fit");

        if !test_file.exists() {
            eprintln!("Skipping test: cygnus.fit not found at {:?}", test_file);
            return;
        }

        let result = FitsFile::open(&test_file);
        assert!(
            result.is_ok(),
            "Failed to open PixInsight FITS file with large header: {:?}",
            result.err()
        );

        let fits = result.unwrap();
        assert!(fits.num_hdus() >= 1, "Expected at least one HDU");

        let hdu_info = fits.hdu_info(0).expect("Primary HDU info should exist");
        assert!(
            hdu_info.header_size > FITS_BLOCK_SIZE * 4,
            "Expected header larger than 4 blocks (11520 bytes), got {} bytes",
            hdu_info.header_size
        );
    }
}
