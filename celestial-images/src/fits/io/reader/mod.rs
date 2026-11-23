use crate::fits::block::FITS_BLOCK_SIZE;
use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::ascii_table::AsciiTableHdu;
use crate::fits::hdu::binary_table::BinaryTableHdu;
use crate::fits::hdu::image::ImageHdu;
use crate::fits::hdu::primary::PrimaryHdu;
use crate::fits::hdu::random_groups::RandomGroupsHdu;
use crate::fits::hdu::Hdu;
use crate::fits::header::parser::{Header, HeaderParser};
use crate::fits::util::buffer_pool::BufferPool;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

mod checksum;
mod data;
pub(crate) mod data_size;
mod scan;
#[cfg(test)]
mod test_support;

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
    scan_error: Option<FitsError>,
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
        let mut fits = Self {
            reader,
            hdus: Vec::new(),
            current_hdu: 0,
            header_cache: HashMap::new(),
            buffer_pool: BufferPool::default(),
            scan_error: None,
        };
        fits.scan_hdus()?;
        Ok(fits)
    }
}

impl<R: Read + Seek> FitsFile<R> {
    pub fn new(reader: R) -> Result<Self> {
        let mut fits = Self {
            reader,
            hdus: Vec::new(),
            current_hdu: 0,
            header_cache: HashMap::new(),
            buffer_pool: BufferPool::default(),
            scan_error: None,
        };
        fits.scan_hdus()?;
        Ok(fits)
    }

    // Why the scan stopped before the end of the file, when it did. The
    // HDUs before the failure stay readable.
    pub fn scan_error(&self) -> Option<&FitsError> {
        self.scan_error.as_ref()
    }

    pub(crate) fn take_scan_error(&mut self) -> Option<FitsError> {
        self.scan_error.take()
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

    // Returns the header without classifying the HDU, so an extension of a
    // type this crate cannot read still has a readable header.
    pub fn get_header(&mut self, index: usize) -> Result<Header> {
        self.move_to_hdu(index)?;
        self.load_header(index)
    }

    pub fn get_header_value(&mut self, index: usize, keyword: &str) -> Result<Option<String>> {
        let header = self.get_header(index)?;
        Ok(header.get(keyword).value().map(|v| v.to_string()))
    }

    pub fn clear_header_cache(&mut self) {
        self.header_cache.clear();
    }

    pub fn clear_buffer_pool(&mut self) {
        self.buffer_pool.clear();
    }

    fn load_header(&mut self, index: usize) -> Result<Header> {
        if let Some(cached_header) = self.header_cache.get(&index) {
            return Ok(cached_header.clone());
        }
        let hdu_info = self.hdus.get(index).ok_or(FitsError::HduNotFound(index))?;
        let header = self.read_header_at(hdu_info.header_start, hdu_info.header_size)?;
        self.header_cache.insert(index, header.clone());
        Ok(header)
    }

    fn read_header_at(&mut self, position: u64, header_size: usize) -> Result<Header> {
        self.reader.seek(SeekFrom::Start(position))?;
        let mut header_data = self.buffer_pool.get(header_size);
        self.reader.read_exact(&mut header_data)?;
        let header = HeaderParser::parse_header(&header_data);
        self.buffer_pool.return_buffer(header_data);
        header
    }

    fn read_current_hdu(&mut self) -> Result<Hdu> {
        let hdu_info = self
            .hdus
            .get(self.current_hdu)
            .ok_or(FitsError::HduNotFound(self.current_hdu))?
            .clone();
        let header = self.load_header(self.current_hdu)?;

        if self.current_hdu == 0 {
            if !header.is_primary() {
                return Err(FitsError::InvalidFormat(
                    "First HDU must be a primary HDU".to_string(),
                ));
            }
            if is_random_groups(&header) {
                return Ok(Hdu::RandomGroups(RandomGroupsHdu::new(header, hdu_info)));
            }
            Ok(Hdu::Primary(PrimaryHdu::new(header, hdu_info)))
        } else {
            if !header.is_extension() {
                return Err(FitsError::InvalidFormat(
                    "Non-primary HDUs must be extensions".to_string(),
                ));
            }
            self.create_extension_hdu(header, hdu_info)
        }
    }

    fn create_extension_hdu(&self, header: Header, hdu_info: HduInfo) -> Result<Hdu> {
        let xtension =
            header
                .get("XTENSION")
                .as_str()
                .ok_or_else(|| FitsError::KeywordNotFound {
                    keyword: "XTENSION".to_string(),
                })?;

        match xtension {
            "IMAGE" => Ok(Hdu::Image(Box::new(ImageHdu::new(header, hdu_info)))),
            "TABLE" => Ok(Hdu::AsciiTable(AsciiTableHdu::new(header, hdu_info)?)),
            // A3DTABLE is the name binary tables had before the standard
            // named them BINTABLE.
            "BINTABLE" | "A3DTABLE" => Ok(Hdu::BinaryTable(BinaryTableHdu::new(header, hdu_info)?)),
            _ => Err(FitsError::InvalidFormat(format!(
                "Unsupported XTENSION type: {}",
                xtension
            ))),
        }
    }
}

// FITS 4.0 §6: GROUPS = T with NAXIS1 = 0 marks random groups.
fn is_random_groups(header: &Header) -> bool {
    header.get("GROUPS").as_bool() == Some(true) && header.get("NAXIS1").as_i64() == Some(0)
}

impl FitsReader {
    pub(crate) fn new(file: File) -> Self {
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
    use super::test_support::{
        create_image_fits, create_malformed_fits, create_minimal_fits, MockFitsBuilder,
    };
    use super::*;
    use crate::fits::header::keywords::{Keyword, KeywordValue};
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

        assert!(fits_file.check_end_of_file(5760).unwrap());
        assert!(!fits_file.check_end_of_file(0).unwrap());
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
        assert_eq!(data, test_data);
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
    fn get_header_value_returns_fits_value_text() {
        let mut fits = open_bytes(hdu_bytes(
            &[
                "SIMPLE  = T",
                "BITPIX  = 8",
                "NAXIS   = 0",
                "CRPIX1  = 512.0",
                "OBSERVER= 'O''Brien'",
            ],
            &[],
        ))
        .unwrap();
        assert_eq!(
            fits.get_header_value(0, "CRPIX1").unwrap().as_deref(),
            Some("512.0")
        );
        assert_eq!(
            fits.get_header_value(0, "OBSERVER").unwrap().as_deref(),
            Some("'O''Brien'")
        );
    }

    #[test]
    fn get_image_info_success() {
        let test_data: Vec<i16> = (0..100).collect();
        let fits_data = create_image_fits(16, &[10, 10], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let (dimensions, bitpix) = fits_file.get_image_info(0).unwrap();
        assert_eq!(dimensions, vec![10, 10]);
        assert_eq!(bitpix, crate::core::types::BitPix::I16);
    }

    #[test]
    fn primary_hdu_with_data_success() {
        let test_data: Vec<i16> = (0..100).collect();
        let fits_data = create_image_fits(16, &[10, 10], &test_data);
        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let (hdu, data): (_, Vec<i16>) = fits_file.primary_hdu_with_data().unwrap();
        assert!(matches!(hdu, crate::fits::hdu::primary::PrimaryHdu { .. }));
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
        assert!(matches!(
            FitsFile::new(cursor),
            Err(FitsError::InvalidFormat(_))
        ));
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
    fn read_current_hdu_past_the_last_hdu_is_not_found() {
        let mut fits_file = FitsFile::new(Cursor::new(create_minimal_fits())).unwrap();
        fits_file.current_hdu = 5;
        let result = fits_file.read_current_hdu();
        assert!(
            matches!(result, Err(FitsError::HduNotFound(5))),
            "{result:?}"
        );
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
        use crate::fits::header::keywords::{Keyword, KeywordValue};

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
        use crate::fits::header::keywords::{Keyword, KeywordValue};

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
        use crate::fits::header::keywords::{Keyword, KeywordValue};

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
        use crate::fits::header::keywords::{Keyword, KeywordValue};

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
    fn workflow_test() {
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
        assert_eq!(bitpix, crate::core::types::BitPix::I16);

        let (hdu, data): (_, Vec<i16>) = fits_file.read_hdu_with_data(0).unwrap();
        assert_eq!(data.len(), 400);
        assert!(matches!(hdu, Hdu::Primary(_)));
    }

    #[test]
    fn get_image_info_valid_header() {
        let mut fits_data = MockFitsBuilder::new()
            .card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", "16", "Bits per pixel")
            .card("NAXIS", "2", "Number of axes")
            .card("NAXIS1", "10", "Axis 1 size")
            .card("NAXIS2", "10", "Axis 2 size")
            .build_memory();
        fits_data.resize(fits_data.len() + FITS_BLOCK_SIZE, 0);

        let cursor = Cursor::new(fits_data);
        let mut fits_file = FitsFile::new(cursor).unwrap();

        let result = fits_file.get_image_info(0);
        assert!(result.is_ok());
        let (dimensions, bitpix) = result.unwrap();
        assert_eq!(dimensions, vec![10, 10]);
        assert_eq!(bitpix, crate::core::types::BitPix::I16);
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
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "80A"));

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
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "A80"));
        header.add_keyword(Keyword::integer("TBCOL1", 1));

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

    fn header_blocks(cards: &[&str]) -> Vec<u8> {
        let mut bytes: Vec<u8> = cards
            .iter()
            .chain(std::iter::once(&"END"))
            .flat_map(|c| format!("{:<80}", c).into_bytes())
            .collect();
        bytes.resize(
            bytes.len().div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE,
            b' ',
        );
        bytes
    }

    fn hdu_bytes(cards: &[&str], data: &[u8]) -> Vec<u8> {
        let mut bytes = header_blocks(cards);
        bytes.extend_from_slice(data);
        bytes.resize(bytes.len().div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE, 0);
        bytes
    }

    const EMPTY_PRIMARY: &[&str] = &["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 0", "EXTEND  = T"];

    const I16_PAIR_EXTENSION: &[&str] = &[
        "XTENSION= 'IMAGE   '",
        "BITPIX  = 16",
        "NAXIS   = 1",
        "NAXIS1  = 2",
        "PCOUNT  = 0",
        "GCOUNT  = 1",
    ];

    fn open_bytes(bytes: Vec<u8>) -> Result<FitsFile<Cursor<Vec<u8>>>> {
        FitsFile::new(Cursor::new(bytes))
    }

    fn i16_data(fits: &mut FitsFile<Cursor<Vec<u8>>>, index: usize) -> Vec<i16> {
        let (_, data): (_, Vec<i16>) = fits.read_hdu_with_data(index).unwrap();
        data
    }

    #[test]
    fn scan_continues_past_empty_extension() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &[
                "XTENSION= 'IMAGE   '",
                "BITPIX  = 16",
                "NAXIS   = 0",
                "PCOUNT  = 0",
                "GCOUNT  = 1",
            ],
            &[],
        ));
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 5, 0, 6]));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 3);
        assert_eq!(i16_data(&mut fits, 2), vec![5, 6]);
    }

    #[test]
    fn scan_skips_heap_counted_by_pcount() {
        let mut data = vec![0, 0, 0, 1];
        data.extend(std::iter::repeat_n(0, FITS_BLOCK_SIZE));
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &[
                "XTENSION= 'BINTABLE'",
                "BITPIX  = 8",
                "NAXIS   = 2",
                "NAXIS1  = 4",
                "NAXIS2  = 1",
                "PCOUNT  = 2880",
                "GCOUNT  = 1",
                "TFIELDS = 1",
                "TFORM1  = '1J      '",
            ],
            &data,
        ));
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 3);
        assert_eq!(fits.hdu_info(1).unwrap().data_size, 2 * FITS_BLOCK_SIZE);
        assert_eq!(i16_data(&mut fits, 2), vec![7, 8]);
    }

    #[test]
    fn scan_multiplies_data_size_by_gcount() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &[
                "XTENSION= 'IMAGE   '",
                "BITPIX  = 8",
                "NAXIS   = 1",
                "NAXIS1  = 2000",
                "PCOUNT  = 0",
                "GCOUNT  = 2",
            ],
            &[1; 4000],
        ));
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 3);
        assert_eq!(fits.hdu_info(1).unwrap().data_size, 2 * FITS_BLOCK_SIZE);
        assert_eq!(i16_data(&mut fits, 2), vec![7, 8]);
    }

    #[test]
    fn scan_sizes_random_groups_without_naxis1() {
        // 3 groups of 2 parameters and a 1x2 array, 4 bytes per value.
        let mut bytes = hdu_bytes(
            &[
                "SIMPLE  = T",
                "BITPIX  = -32",
                "NAXIS   = 3",
                "NAXIS1  = 0",
                "NAXIS2  = 1",
                "NAXIS3  = 2",
                "GROUPS  = T",
                "PCOUNT  = 2",
                "GCOUNT  = 3",
            ],
            &[1; 48],
        );
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 2);
        assert_eq!(fits.hdu_info(0).unwrap().data_size, FITS_BLOCK_SIZE);
        assert_eq!(
            fits.hdu_info(1).unwrap().header_start,
            2 * FITS_BLOCK_SIZE as u64
        );
        assert_eq!(i16_data(&mut fits, 1), vec![7, 8]);
    }

    #[test]
    fn scan_reads_extensions_when_extend_is_false() {
        let mut bytes = hdu_bytes(
            &["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 0", "EXTEND  = F"],
            &[],
        );
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 5, 0, 6]));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 2);
        assert_eq!(i16_data(&mut fits, 1), vec![5, 6]);
    }

    fn scan_error_after_primary(bytes: Vec<u8>) -> FitsError {
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 1);
        fits.take_scan_error().unwrap()
    }

    #[test]
    fn scan_records_bad_block_after_last_hdu() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(b"garbage!".repeat(FITS_BLOCK_SIZE / 8));
        assert!(matches!(
            scan_error_after_primary(bytes),
            FitsError::InvalidFormat(_)
        ));
    }

    #[test]
    fn scan_records_partial_block_after_last_hdu() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(b"0123456789");
        assert!(matches!(
            scan_error_after_primary(bytes),
            FitsError::InvalidFormat(_)
        ));
    }

    #[test]
    fn scan_records_error_in_later_extension_header() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &["XTENSION= 'IMAGE   '", "BITPIX  = 16", "NAXIS   = 1"],
            &[],
        ));
        assert!(matches!(
            scan_error_after_primary(bytes),
            FitsError::KeywordNotFound { .. }
        ));
    }

    #[test]
    fn scan_keeps_primary_before_truncated_extension() {
        let mut bytes = hdu_bytes(
            &["SIMPLE  = T", "BITPIX  = 16", "NAXIS   = 1", "NAXIS1  = 2"],
            &[0, 5, 0, 6],
        );
        bytes.extend(header_blocks(&[
            "XTENSION= 'IMAGE   '",
            "BITPIX  = 16",
            "NAXIS   = 1",
            "NAXIS1  = 2000",
        ]));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 1);
        assert_eq!(i16_data(&mut fits, 0), vec![5, 6]);
        assert!(matches!(
            fits.take_scan_error(),
            Some(FitsError::InvalidFormat(_))
        ));
        assert!(fits.take_scan_error().is_none());
    }

    #[test]
    fn scan_has_no_error_for_intact_file() {
        let mut fits = open_bytes(hdu_bytes(EMPTY_PRIMARY, &[])).unwrap();
        assert!(fits.take_scan_error().is_none());
    }

    #[test]
    fn scan_rejects_primary_whose_first_card_is_not_simple() {
        let bytes = hdu_bytes(
            &["BITPIX  = 16", "NAXIS   = 1", "NAXIS1  = 2", "SIMPLE  = T"],
            &[0, 5, 0, 6],
        );
        assert!(matches!(
            open_bytes(bytes),
            Err(FitsError::InvalidFormat(_))
        ));
    }

    #[test]
    fn scan_stops_at_extension_whose_first_card_is_not_xtension() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &["SIMPLE  = T", "BITPIX  = 16", "NAXIS   = 1", "NAXIS1  = 2"],
            &[0, 5, 0, 6],
        ));
        assert!(matches!(
            scan_error_after_primary(bytes),
            FitsError::InvalidFormat(_)
        ));
    }

    // A header block whose END card was damaged, so the header runs on into
    // the next block.
    fn block_without_end(cards: &[&str]) -> Vec<u8> {
        let mut bytes: Vec<u8> = cards
            .iter()
            .flat_map(|c| format!("{:<80}", c).into_bytes())
            .collect();
        bytes.resize(FITS_BLOCK_SIZE, b' ');
        bytes
    }

    // Shape of fuzz case 00562: an empty primary header with no END card,
    // followed by the header and data of an image extension.
    #[test]
    fn scan_rejects_primary_merged_with_extension_header() {
        let mut bytes = block_without_end(EMPTY_PRIMARY);
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 5, 0, 6]));
        assert!(matches!(
            open_bytes(bytes),
            Err(FitsError::InvalidFormat(ref m)) if m.contains("BITPIX")
        ));
    }

    // Shape of fuzz case 01381: two copies of the primary header, then a
    // tile-compressed table header, all read as one header.
    #[test]
    fn scan_rejects_primary_merged_with_two_more_headers() {
        let primary = ["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 0", "EXTEND  = T"];
        let mut bytes = block_without_end(&primary);
        bytes.extend(block_without_end(&primary));
        bytes.extend(hdu_bytes(
            &[
                "XTENSION= 'BINTABLE'",
                "BITPIX  = 8",
                "NAXIS   = 2",
                "NAXIS1  = 8",
                "NAXIS2  = 1",
                "PCOUNT  = 0",
                "GCOUNT  = 1",
                "TFIELDS = 1",
                "TFORM1  = '1PB(1)  '",
                "ZIMAGE  = T",
            ],
            &[0; 8],
        ));
        assert!(matches!(
            open_bytes(bytes),
            Err(FitsError::InvalidFormat(ref m)) if m.contains("SIMPLE")
        ));
    }

    #[test]
    fn scan_rejects_primary_with_bitpix_twice() {
        let bytes = hdu_bytes(
            &[
                "SIMPLE  = T",
                "BITPIX  = 16",
                "BITPIX  = 16",
                "NAXIS   = 1",
                "NAXIS1  = 2",
            ],
            &[0, 5, 0, 6],
        );
        assert!(matches!(
            open_bytes(bytes),
            Err(FitsError::InvalidFormat(ref m)) if m.contains("BITPIX")
        ));
    }

    #[test]
    fn scan_stops_at_extension_with_repeated_structural_keyword() {
        for repeated in ["NAXIS1  = 2", "TFORM1  = '1J      '"] {
            let mut cards = I16_PAIR_EXTENSION.to_vec();
            cards.extend(["BZERO   = 0", "TFORM1  = '1J      '", repeated]);
            let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
            bytes.extend(hdu_bytes(&cards, &[0, 5, 0, 6]));
            assert!(
                matches!(scan_error_after_primary(bytes), FitsError::InvalidFormat(_)),
                "{}",
                repeated
            );
        }
    }

    #[test]
    fn scan_keeps_first_value_of_repeated_user_keyword() {
        let mut fits = open_bytes(hdu_bytes(
            &[
                "SIMPLE  = T",
                "BITPIX  = 8",
                "NAXIS   = 0",
                "OBJECT  = 'first'",
                "OBJECT  = 'second'",
                "NAXIS10 = 1",
                "TFORMAT = 1",
                "TFORMAT = 2",
            ],
            &[],
        ))
        .unwrap();
        assert!(fits.take_scan_error().is_none());
        let header = fits.get_header(0).unwrap();
        assert_eq!(header.get("OBJECT").as_str(), Some("first"));
        assert_eq!(header.get("TFORMAT").as_i64(), Some(1));
    }

    #[test]
    fn get_header_reads_extension_of_unknown_type() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &[
                "XTENSION= 'FOREIGN '",
                "BITPIX  = 16",
                "NAXIS   = 1",
                "NAXIS1  = 2",
                "PCOUNT  = 0",
                "GCOUNT  = 1",
            ],
            &[0, 1, 0, 2],
        ));
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 2);
        let header = fits.get_header(1).unwrap();
        assert_eq!(header.get("XTENSION").as_str(), Some("FOREIGN"));
    }

    #[test]
    fn scan_ignores_nul_blocks_after_last_hdu() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend([0; 2 * FITS_BLOCK_SIZE]);
        assert_eq!(open_bytes(bytes).unwrap().num_hdus(), 1);
    }

    #[test]
    fn scan_accepts_data_without_final_padding() {
        let mut bytes =
            header_blocks(&["SIMPLE  = T", "BITPIX  = 16", "NAXIS   = 1", "NAXIS1  = 2"]);
        bytes.extend([0, 5, 0, 6]);
        let mut fits = open_bytes(bytes).unwrap();
        assert_eq!(fits.num_hdus(), 1);
        assert_eq!(i16_data(&mut fits, 0), vec![5, 6]);
    }

    #[test]
    fn scan_rejects_data_unit_longer_than_file() {
        let mut bytes = header_blocks(&[
            "SIMPLE  = T",
            "BITPIX  = 16",
            "NAXIS   = 2",
            "NAXIS1  = 100",
            "NAXIS2  = 100",
        ]);
        bytes.extend([0; FITS_BLOCK_SIZE]);
        assert!(open_bytes(bytes).is_err());
    }

    #[test]
    fn scan_rejects_huge_axes_before_allocating() {
        let mut bytes = header_blocks(&[
            "SIMPLE  = T",
            "BITPIX  = 16",
            "NAXIS   = 2",
            "NAXIS1  = 1000000000",
            "NAXIS2  = 1000000000",
        ]);
        bytes.extend([0; FITS_BLOCK_SIZE]);
        assert!(open_bytes(bytes).is_err());
    }

    // The old scan looped NAXIS = -1 as usize times, so the call runs on a
    // thread and a hang fails the test after the timeout.
    fn scan_within_timeout(cards: &'static [&'static str]) -> Result<usize> {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let result = open_bytes(header_blocks(cards)).map(|f| f.num_hdus());
            tx.send(result).ok();
        });
        rx.recv_timeout(std::time::Duration::from_secs(5))
            .expect("the scan did not return within 5 s")
    }

    #[test]
    fn scan_rejects_negative_naxis() {
        let result = scan_within_timeout(&["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = -1"]);
        assert!(matches!(
            result,
            Err(FitsError::InvalidKeywordValue { ref keyword, .. }) if keyword == "NAXIS"
        ));
    }

    #[test]
    fn scan_rejects_naxis_above_999() {
        let result = scan_within_timeout(&["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 1000"]);
        assert!(result.is_err());
    }

    #[test]
    fn scan_rejects_negative_axis_length() {
        let result = open_bytes(header_blocks(&[
            "SIMPLE  = T",
            "BITPIX  = 8",
            "NAXIS   = 2",
            "NAXIS1  = -1",
            "NAXIS2  = 1",
        ]));
        assert!(matches!(
            result,
            Err(FitsError::InvalidKeywordValue { ref keyword, .. }) if keyword == "NAXIS1"
        ));
    }

    #[test]
    fn scan_rejects_bitpix_outside_i32() {
        let result = open_bytes(hdu_bytes(
            &[
                "SIMPLE  = T",
                "BITPIX  = 4294967304",
                "NAXIS   = 1",
                "NAXIS1  = 4",
            ],
            &[1, 2, 3, 4],
        ));
        assert!(matches!(
            result,
            Err(FitsError::InvalidKeywordValue { ref keyword, .. }) if keyword == "BITPIX"
        ));
    }

    #[test]
    fn scan_rejects_missing_axis_length() {
        let result = open_bytes(hdu_bytes(
            &["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 2", "NAXIS1  = 4"],
            &[1, 2, 3, 4],
        ));
        assert!(matches!(
            result,
            Err(FitsError::KeywordNotFound { ref keyword }) if keyword == "NAXIS2"
        ));
        let result = open_bytes(header_blocks(&[
            "SIMPLE  = T",
            "BITPIX  = 8",
            "NAXIS   = 999",
        ]));
        assert!(matches!(
            result,
            Err(FitsError::KeywordNotFound { ref keyword }) if keyword == "NAXIS1"
        ));
    }

    #[test]
    fn scan_rejects_real_axis_length() {
        let result = open_bytes(hdu_bytes(
            &["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 1", "NAXIS1  = 4.0"],
            &[1, 2, 3, 4],
        ));
        assert!(result.is_err());
    }

    #[test]
    fn scan_rejects_negative_pcount() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(hdu_bytes(
            &[
                "XTENSION= 'IMAGE   '",
                "BITPIX  = 16",
                "NAXIS   = 1",
                "NAXIS1  = 2",
                "PCOUNT  = -4",
                "GCOUNT  = 1",
            ],
            &[0, 5, 0, 6],
        ));
        assert!(matches!(
            scan_error_after_primary(bytes),
            FitsError::InvalidKeywordValue { ref keyword, .. } if keyword == "PCOUNT"
        ));
    }

    #[test]
    fn scan_without_end_card_names_end() {
        let bytes = format!("{:<2880}", "SIMPLE  = T").into_bytes();
        let message = open_bytes(bytes).unwrap_err().to_string();
        assert!(message.contains("END"), "{}", message);
    }

    #[test]
    fn scan_rejects_text_after_end_keyword() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        let end = EMPTY_PRIMARY.len() * 80;
        bytes[end + 20] = b'0';
        bytes.extend(hdu_bytes(I16_PAIR_EXTENSION, &[0, 1, 0, 2]));
        assert!(matches!(
            open_bytes(bytes),
            Err(FitsError::InvalidFormat(_))
        ));
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

    #[test]
    fn scan_keeps_first_value_of_repeated_optional_structural_keyword() {
        let repeats = [
            ("BZERO   = 0", "BZERO   = 100", KeywordValue::Integer(0)),
            ("BSCALE  = 1", "BSCALE  = 2", KeywordValue::Integer(1)),
            ("BLANK   = -1", "BLANK   = -2", KeywordValue::Integer(-1)),
            ("EXTEND  = F", "EXTEND  = T", KeywordValue::Logical(false)),
            ("ZIMAGE  = F", "ZIMAGE  = T", KeywordValue::Logical(false)),
            ("TBCOL1  = 1", "TBCOL1  = 2", KeywordValue::Integer(1)),
            (
                "TTYPE1  = 'a'",
                "TTYPE1  = 'b'",
                KeywordValue::String("a".into()),
            ),
        ];
        for (first, last, expected) in repeats {
            let mut cards = vec!["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 0"];
            cards.extend([first, last]);
            let mut fits = open_bytes(hdu_bytes(&cards, &[])).unwrap();
            assert!(fits.take_scan_error().is_none(), "{}", last);
            let header = fits.get_header(0).unwrap();
            assert_eq!(
                header.get(last[..8].trim_end()).value(),
                Some(&expected),
                "{}",
                last
            );
        }
    }

    #[test]
    fn scan_rejects_primary_whose_naxis_holds_non_ascii() {
        let mut bytes = hdu_bytes(&["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 0"], &[]);
        bytes[2 * 80 + 10] = 0xE9;
        let result = open_bytes(bytes);
        assert!(
            matches!(
                &result,
                Err(FitsError::InvalidKeywordValue { keyword, value })
                    if keyword == "NAXIS" && value == "'?'"
            ),
            "{:?}",
            result.err()
        );
    }

    #[test]
    fn scan_error_reports_why_the_scan_stopped_without_clearing_it() {
        let mut bytes = hdu_bytes(EMPTY_PRIMARY, &[]);
        bytes.extend(b"garbage!".repeat(FITS_BLOCK_SIZE / 8));
        let fits = open_bytes(bytes).unwrap();
        for _ in 0..2 {
            assert_eq!(
                fits.scan_error().map(ToString::to_string).as_deref(),
                Some(
                    "Invalid FITS format: Header at byte 2880 has no END card in its 1 \
                     complete blocks: failed to fill whole buffer"
                )
            );
        }
        let intact = open_bytes(hdu_bytes(EMPTY_PRIMARY, &[])).unwrap();
        assert!(intact.scan_error().is_none());
    }
}
