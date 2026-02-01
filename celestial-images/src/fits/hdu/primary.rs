use super::{HduTrait, HduType};
use crate::fits::header::Header;
use crate::fits::io::reader::HduInfo;

#[derive(Debug)]
pub struct PrimaryHdu {
    header: Header,
    info: HduInfo,
}

impl PrimaryHdu {
    pub fn new(header: Header, info: HduInfo) -> Self {
        Self { header, info }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn info(&self) -> &HduInfo {
        &self.info
    }

    pub fn has_data(&self) -> bool {
        self.header
            .get_keyword_value("NAXIS")
            .and_then(|v| v.as_integer())
            .unwrap_or(0)
            > 0
    }

    pub fn data_dimensions(&self) -> Vec<usize> {
        let naxis = self
            .header
            .get_keyword_value("NAXIS")
            .and_then(|v| v.as_integer())
            .unwrap_or(0) as usize;

        let mut dims = Vec::with_capacity(naxis);
        for i in 1..=naxis {
            let axis_name = format!("NAXIS{}", i);
            let axis_size = self
                .header
                .get_keyword_value(&axis_name)
                .and_then(|v| v.as_integer())
                .unwrap_or(0) as usize;
            dims.push(axis_size);
        }
        dims
    }

    pub fn bitpix(&self) -> Option<crate::core::BitPix> {
        self.header
            .get_keyword_value("BITPIX")
            .and_then(|v| v.as_integer())
            .and_then(|i| crate::core::BitPix::from_value(i as i32))
    }
}

impl HduTrait for PrimaryHdu {
    fn header(&self) -> &Header {
        &self.header
    }

    fn info(&self) -> &HduInfo {
        &self.info
    }

    fn hdu_type(&self) -> HduType {
        HduType::Primary
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::BitPix;
    use crate::fits::header::{Header, Keyword};
    use crate::fits::io::reader::HduInfo;
    use crate::fits::{FitsError, Result};
    use std::io::Cursor;

    fn create_test_header(naxis: i64, bitpix: i32) -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", naxis));
        header.add_keyword(Keyword::integer("BITPIX", bitpix as i64));
        header.add_keyword(Keyword::logical("EXTEND", false));

        if naxis > 0 {
            for i in 1..=naxis {
                let key = format!("NAXIS{}", i);
                header.add_keyword(Keyword::integer(key, 10));
            }
        }

        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 200,
        }
    }

    #[test]
    fn new_creates_primary_hdu() {
        let header = create_test_header(2, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.info.index, 0);
        assert_eq!(
            hdu.header
                .get_keyword_value("NAXIS")
                .unwrap()
                .as_integer()
                .unwrap(),
            2
        );
        assert!(hdu
            .header
            .get_keyword_value("SIMPLE")
            .unwrap()
            .as_logical()
            .unwrap());
    }

    #[test]
    fn header_returns_header_reference() {
        let header = create_test_header(2, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let header_ref = hdu.header();
        assert_eq!(
            header_ref
                .get_keyword_value("NAXIS")
                .unwrap()
                .as_integer()
                .unwrap(),
            2
        );
        assert!(header_ref
            .get_keyword_value("SIMPLE")
            .unwrap()
            .as_logical()
            .unwrap());
    }

    #[test]
    fn info_returns_info_reference() {
        let header = create_test_header(2, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let info_ref = hdu.info();
        assert_eq!(info_ref.index, 0);
        assert_eq!(info_ref.data_start, 2880);
    }

    #[test]
    fn has_data_true_when_naxis_greater_than_zero() {
        let header = create_test_header(2, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert!(hdu.has_data());
    }

    #[test]
    fn has_data_false_when_naxis_zero() {
        let header = create_test_header(0, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert!(!hdu.has_data());
    }

    #[test]
    fn has_data_false_when_naxis_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert!(!hdu.has_data());
    }

    #[test]
    fn data_dimensions_returns_correct_dimensions() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 1024));
        header.add_keyword(Keyword::integer("NAXIS2", 512));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.data_dimensions(), vec![1024, 512]);
    }

    #[test]
    fn data_dimensions_returns_empty_for_zero_naxis() {
        let header = create_test_header(0, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.data_dimensions(), vec![]);
    }

    #[test]
    fn data_dimensions_handles_missing_axis_sizes() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 3));
        header.add_keyword(Keyword::integer("NAXIS1", 100));
        header.add_keyword(Keyword::integer("NAXIS3", 50));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.data_dimensions(), vec![100, 0, 50]);
    }

    #[test]
    fn bitpix_returns_correct_value() {
        let header = create_test_header(2, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.bitpix(), Some(BitPix::I16));
    }

    #[test]
    fn bitpix_returns_correct_values_for_all_types() {
        for (bitpix_val, expected) in [
            (8, BitPix::U8),
            (16, BitPix::I16),
            (32, BitPix::I32),
            (64, BitPix::I64),
            (-32, BitPix::F32),
            (-64, BitPix::F64),
        ] {
            let header = create_test_header(1, bitpix_val);
            let info = create_test_hdu_info();
            let hdu = PrimaryHdu::new(header, info);

            assert_eq!(hdu.bitpix(), Some(expected));
        }
    }

    #[test]
    fn bitpix_returns_none_for_invalid_value() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("BITPIX", 99));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.bitpix(), None);
    }

    #[test]
    fn bitpix_returns_none_when_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.bitpix(), None);
    }

    #[test]
    fn read_data_returns_empty_when_no_data() {
        let header = create_test_header(0, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let mut cursor = Cursor::new(vec![]);
        let result: Result<Vec<i16>> = hdu.read_data(&mut cursor);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), Vec::<i16>::new());
    }

    #[test]
    fn read_data_fails_when_bitpix_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 2));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let mut cursor = Cursor::new(vec![]);
        let result: Result<Vec<i16>> = hdu.read_data(&mut cursor);
        assert!(matches!(result, Err(FitsError::KeywordNotFound { .. })));
    }

    #[test]
    fn read_data_fails_on_type_mismatch() {
        let header = create_test_header(1, 32);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let mut cursor = Cursor::new(vec![]);
        let result: Result<Vec<i16>> = hdu.read_data(&mut cursor);
        assert!(matches!(result, Err(FitsError::TypeMismatch { .. })));
    }

    #[test]
    fn read_data_success_with_matching_types() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 2));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let mut info = create_test_hdu_info();
        info.data_start = 0;
        let hdu = PrimaryHdu::new(header, info);

        let data = vec![0x01, 0x23, 0x45, 0x67];
        let mut cursor = Cursor::new(data);
        let result: Result<Vec<i16>> = hdu.read_data(&mut cursor);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), vec![0x0123, 0x4567]);
    }

    #[test]
    fn calculate_data_size_returns_zero_for_no_dimensions() {
        let header = create_test_header(0, 16);
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.calculate_data_size().unwrap(), 0);
    }

    #[test]
    fn calculate_data_size_calculates_correctly() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 10));
        header.add_keyword(Keyword::integer("NAXIS2", 5));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.calculate_data_size().unwrap(), 10 * 5 * 2);
    }

    #[test]
    fn calculate_data_size_for_different_bitpix_values() {
        for (bitpix_val, bytes_per_pixel) in [(8, 1), (16, 2), (32, 4), (64, 8), (-32, 4), (-64, 8)]
        {
            let mut header = Header::new();
            header.add_keyword(Keyword::logical("SIMPLE", true));
            header.add_keyword(Keyword::integer("NAXIS", 1));
            header.add_keyword(Keyword::integer("NAXIS1", 100));
            header.add_keyword(Keyword::integer("BITPIX", bitpix_val));
            let info = create_test_hdu_info();
            let hdu = PrimaryHdu::new(header, info);

            assert_eq!(hdu.calculate_data_size().unwrap(), 100 * bytes_per_pixel);
        }
    }

    #[test]
    fn calculate_data_size_fails_on_missing_bitpix() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 10));
        header.add_keyword(Keyword::integer("NAXIS2", 5));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let result = hdu.calculate_data_size();
        assert!(matches!(result, Err(FitsError::KeywordNotFound { .. })));
    }

    #[test]
    fn calculate_data_size_handles_overflow() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", usize::MAX as i64));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        let result = hdu.calculate_data_size();
        assert!(matches!(result, Err(FitsError::InvalidFormat(_))));
    }

    #[test]
    fn calculate_data_size_three_dimensional() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 3));
        header.add_keyword(Keyword::integer("NAXIS1", 10));
        header.add_keyword(Keyword::integer("NAXIS2", 20));
        header.add_keyword(Keyword::integer("NAXIS3", 5));
        header.add_keyword(Keyword::integer("BITPIX", -32));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert_eq!(hdu.calculate_data_size().unwrap(), 10 * 20 * 5 * 4);
    }

    #[test]
    fn all_methods_work_together_primary_hdu() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 100));
        header.add_keyword(Keyword::integer("NAXIS2", 50));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::logical("EXTEND", true));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert!(hdu.has_data());
        assert_eq!(hdu.data_dimensions(), vec![100, 50]);
        assert_eq!(hdu.bitpix(), Some(BitPix::U8));
        assert_eq!(hdu.calculate_data_size().unwrap(), 100 * 50);
        assert!(hdu
            .header()
            .get_keyword_value("SIMPLE")
            .unwrap()
            .as_logical()
            .unwrap());
        assert!(hdu
            .header()
            .get_keyword_value("EXTEND")
            .unwrap()
            .as_logical()
            .unwrap());
    }

    #[test]
    fn primary_hdu_minimal_valid() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::logical("EXTEND", false));
        let info = create_test_hdu_info();
        let hdu = PrimaryHdu::new(header, info);

        assert!(!hdu.has_data());
        assert_eq!(hdu.data_dimensions(), vec![]);
        assert_eq!(hdu.bitpix(), Some(BitPix::U8));
        assert_eq!(hdu.calculate_data_size().unwrap(), 0);
    }
}
