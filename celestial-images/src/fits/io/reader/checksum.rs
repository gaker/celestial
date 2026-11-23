use super::{ChecksumResult, FitsFile, HduInfo};
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::{Header, HeaderParser};
use crate::fits::util::checksum;
use std::io::{Read, Seek, SeekFrom};

impl<R: Read + Seek> FitsFile<R> {
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

    // A mismatch is reported as Some(false) for either card. Changed data
    // fails both, since CHECKSUM covers the data unit too.
    fn validate_checksums(
        &self,
        header: &Header,
        header_bytes: &[u8],
        data_bytes: &[u8],
    ) -> Result<ChecksumResult> {
        Ok(ChecksumResult {
            datasum_valid: validate_datasum_if_present(header, data_bytes),
            checksum_valid: validate_checksum_if_present(header, header_bytes, data_bytes),
        })
    }
}

// A DATASUM card that is not a string of digits cannot match.
fn validate_datasum_if_present(header: &Header, data_bytes: &[u8]) -> Option<bool> {
    let value = header.get("DATASUM").value()?;
    let expected = value.as_string().unwrap_or_default();
    Some(checksum::verify_datasum(data_bytes, expected))
}

fn validate_checksum_if_present(
    header: &Header,
    header_bytes: &[u8],
    data_bytes: &[u8],
) -> Option<bool> {
    header.get("CHECKSUM").as_str()?;
    Some(checksum::verify_hdu_checksum(header_bytes, data_bytes))
}
