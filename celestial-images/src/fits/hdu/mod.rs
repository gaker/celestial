pub mod ascii_table;
pub mod binary_table;
pub mod image;
pub mod primary;
pub mod random_groups;

pub use ascii_table::{AsciiTableHdu, AsciiTableRowIterator};
pub use binary_table::{BinaryTableHdu, BinaryTableRowIterator};
pub use image::ImageHdu;
pub use primary::PrimaryHdu;
pub use random_groups::RandomGroupsHdu;

use crate::core::{BitPix, ByteOrder};
use crate::fits::data::array::DataArray;
use crate::fits::header::Header;
use crate::fits::io::reader::HduInfo;
use crate::fits::{FitsError, Result};
use std::io::{Read, Seek, SeekFrom};

#[derive(Debug, Clone, PartialEq)]
pub struct ColumnInfo {
    pub index: usize,
    pub name: Option<String>,
    pub format: String,
    pub unit: Option<String>,
    pub null_value: Option<String>,
    pub scale: Option<f64>,
    pub zero_offset: Option<f64>,
    pub display_format: Option<String>,
    pub coordinate_type: Option<String>,
    pub coordinate_reference_pixel: Option<f64>,
    pub coordinate_reference_value: Option<f64>,
    pub coordinate_increment: Option<f64>,
}

impl ColumnInfo {
    pub fn new(index: usize, format: String) -> Self {
        Self {
            index,
            name: None,
            format,
            unit: None,
            null_value: None,
            scale: None,
            zero_offset: None,
            display_format: None,
            coordinate_type: None,
            coordinate_reference_pixel: None,
            coordinate_reference_value: None,
            coordinate_increment: None,
        }
    }

    pub fn with_name(mut self, name: String) -> Self {
        self.name = Some(name);
        self
    }

    pub fn with_unit(mut self, unit: String) -> Self {
        self.unit = Some(unit);
        self
    }

    pub fn with_null_value(mut self, null_value: String) -> Self {
        self.null_value = Some(null_value);
        self
    }

    pub fn with_scale(mut self, scale: f64) -> Self {
        self.scale = Some(scale);
        self
    }

    pub fn with_zero_offset(mut self, zero: f64) -> Self {
        self.zero_offset = Some(zero);
        self
    }
}

pub trait HduTrait: std::fmt::Debug {
    fn header(&self) -> &Header;
    fn info(&self) -> &HduInfo;
    fn hdu_type(&self) -> HduType;

    fn bzero(&self) -> f64 {
        self.header()
            .get_keyword_value("BZERO")
            .and_then(|v| v.as_real())
            .unwrap_or(0.0)
    }

    fn bscale(&self) -> f64 {
        self.header()
            .get_keyword_value("BSCALE")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0)
    }

    fn needs_scaling(&self) -> bool {
        let bzero = self.bzero();
        let bscale = self.bscale();
        bzero != 0.0 || bscale != 1.0
    }

    fn has_data(&self) -> bool {
        self.header()
            .get_keyword_value("NAXIS")
            .and_then(|v| v.as_integer())
            .unwrap_or(0)
            > 0
    }

    fn data_dimensions(&self) -> Vec<usize> {
        let naxis = self
            .header()
            .get_keyword_value("NAXIS")
            .and_then(|v| v.as_integer())
            .unwrap_or(0) as usize;

        let mut dims = Vec::with_capacity(naxis);
        for i in 1..=naxis {
            let axis_name = format!("NAXIS{}", i);
            let axis_size = self
                .header()
                .get_keyword_value(&axis_name)
                .and_then(|v| v.as_integer())
                .unwrap_or(0) as usize;
            dims.push(axis_size);
        }
        dims
    }

    fn bitpix(&self) -> Option<BitPix> {
        self.header()
            .get_keyword_value("BITPIX")
            .and_then(|v| v.as_integer())
            .and_then(|i| BitPix::from_value(i as i32))
    }

    fn calculate_data_size(&self) -> Result<usize> {
        let dimensions = self.data_dimensions();
        if dimensions.is_empty() {
            return Ok(0);
        }

        let total_pixels = dimensions
            .iter()
            .try_fold(1usize, |acc, &dim| acc.checked_mul(dim))
            .ok_or_else(|| FitsError::InvalidFormat("Data dimensions too large".to_string()))?;
        let bitpix = self.bitpix().ok_or_else(|| FitsError::KeywordNotFound {
            keyword: "BITPIX".to_string(),
        })?;

        Ok(total_pixels * bitpix.bytes_per_pixel())
    }

    fn read_raw_data<R>(&self, reader: &mut R) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        if !self.has_data() {
            return Ok(Vec::new());
        }

        let data_size = self.calculate_data_size()?;

        reader.seek(SeekFrom::Start(self.info().data_start))?;

        let mut buffer = vec![0u8; data_size];
        reader.read_exact(&mut buffer)?;
        Ok(buffer)
    }

    fn write_raw_data<W>(&self, writer: &mut W, data: &[u8]) -> Result<()>
    where
        W: std::io::Write,
    {
        writer.write_all(data)?;
        Ok(())
    }

    fn read_data<T, R>(&self, reader: &mut R) -> Result<Vec<T>>
    where
        T: DataArray,
        R: Read + Seek,
    {
        if !self.has_data() {
            return Ok(Vec::new());
        }

        let expected_bitpix = T::BITPIX;
        let actual_bitpix = self.bitpix().ok_or_else(|| FitsError::KeywordNotFound {
            keyword: "BITPIX".to_string(),
        })?;

        if expected_bitpix != actual_bitpix {
            return Err(FitsError::TypeMismatch {
                expected: expected_bitpix,
                actual: actual_bitpix,
            });
        }

        let data_size = self.calculate_data_size()?;

        reader.seek(SeekFrom::Start(self.info().data_start))?;

        let mut buffer = vec![0u8; data_size];
        reader.read_exact(&mut buffer)?;

        let mut data = T::from_bytes(&buffer, ByteOrder::BigEndian)?;

        if self.needs_scaling() && (T::BITPIX == BitPix::F32 || T::BITPIX == BitPix::F64) {
            T::apply_scaling(&mut data, self.bscale(), self.bzero());
        }

        Ok(data)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HduType {
    Primary,
    Image,
    AsciiTable,
    BinaryTable,
    RandomGroups,
    Unknown(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LogicalHduType {
    Image,
    AsciiTable,
    BinaryTable,
    CompressedImage,
    Unknown(String),
}

#[derive(Debug)]
pub enum Hdu {
    Primary(PrimaryHdu),
    Image(Box<ImageHdu>),
    AsciiTable(AsciiTableHdu),
    BinaryTable(BinaryTableHdu),
    RandomGroups(RandomGroupsHdu),
}

impl Hdu {
    pub fn hdu_type(&self) -> HduType {
        match self {
            Hdu::Primary(_) => HduType::Primary,
            Hdu::Image(_) => HduType::Image,
            Hdu::AsciiTable(_) => HduType::AsciiTable,
            Hdu::BinaryTable(_) => HduType::BinaryTable,
            Hdu::RandomGroups(_) => HduType::RandomGroups,
        }
    }

    pub fn header(&self) -> &Header {
        match self {
            Hdu::Primary(hdu) => hdu.header(),
            Hdu::Image(hdu) => hdu.header(),
            Hdu::AsciiTable(hdu) => hdu.header(),
            Hdu::BinaryTable(hdu) => hdu.header(),
            Hdu::RandomGroups(hdu) => hdu.header(),
        }
    }

    pub fn info(&self) -> &HduInfo {
        match self {
            Hdu::Primary(hdu) => hdu.info(),
            Hdu::Image(hdu) => hdu.info(),
            Hdu::AsciiTable(hdu) => hdu.info(),
            Hdu::BinaryTable(hdu) => hdu.info(),
            Hdu::RandomGroups(hdu) => hdu.info(),
        }
    }

    pub fn logical_type(&self) -> LogicalHduType {
        match self {
            Hdu::Primary(_) => LogicalHduType::Image,
            Hdu::Image(_) => LogicalHduType::Image,
            Hdu::AsciiTable(_) => LogicalHduType::AsciiTable,
            Hdu::BinaryTable(hdu) => {
                if hdu.is_compressed_image() {
                    LogicalHduType::CompressedImage
                } else {
                    LogicalHduType::BinaryTable
                }
            }
            Hdu::RandomGroups(_) => LogicalHduType::Image,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Header, Keyword};
    use crate::fits::io::reader::HduInfo;
    use std::io::Cursor;

    fn create_primary_hdu() -> PrimaryHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 10));
        header.add_keyword(Keyword::integer("NAXIS2", 10));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 200,
        };

        PrimaryHdu::new(header, info)
    }

    fn create_image_hdu() -> ImageHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "IMAGE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 5));
        header.add_keyword(Keyword::integer("NAXIS2", 5));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::string("EXTNAME", "SCI"));

        let info = HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 25,
        };

        ImageHdu::new(header, info)
    }

    fn create_binary_table_hdu(compressed: bool) -> BinaryTableHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 20));
        header.add_keyword(Keyword::integer("NAXIS2", 100));
        header.add_keyword(Keyword::integer("TFIELDS", 2));

        if compressed {
            header.add_keyword(Keyword::logical("ZIMAGE", true));
            header.add_keyword(Keyword::string("ZCMPTYPE", "RICE_1"));
        }

        let info = HduInfo {
            index: 2,
            header_start: 5760,
            header_size: 2880,
            data_start: 8640,
            data_size: 2000,
        };

        BinaryTableHdu::new(header, info)
    }

    fn create_ascii_table_hdu() -> AsciiTableHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 80));
        header.add_keyword(Keyword::integer("NAXIS2", 50));
        header.add_keyword(Keyword::integer("TFIELDS", 3));

        let info = HduInfo {
            index: 3,
            header_start: 8640,
            header_size: 2880,
            data_start: 11520,
            data_size: 4000,
        };

        AsciiTableHdu::new(header, info)
    }

    fn create_random_groups_hdu() -> RandomGroupsHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 3));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("NAXIS2", 10));
        header.add_keyword(Keyword::integer("NAXIS3", 10));
        header.add_keyword(Keyword::integer("BITPIX", -32));
        header.add_keyword(Keyword::logical("GROUPS", true));
        header.add_keyword(Keyword::integer("GCOUNT", 5));
        header.add_keyword(Keyword::integer("PCOUNT", 2));

        let info = HduInfo {
            index: 4,
            header_start: 11520,
            header_size: 2880,
            data_start: 14400,
            data_size: 2000,
        };

        RandomGroupsHdu::new(header, info)
    }

    #[test]
    fn hdu_type_enum_values() {
        assert_eq!(HduType::Primary, HduType::Primary);
        assert_ne!(HduType::Primary, HduType::Image);
        assert_eq!(
            HduType::Unknown("CUSTOM".to_string()),
            HduType::Unknown("CUSTOM".to_string())
        );
        assert_ne!(
            HduType::Unknown("CUSTOM1".to_string()),
            HduType::Unknown("CUSTOM2".to_string())
        );
    }

    #[test]
    fn logical_hdu_type_enum_values() {
        assert_eq!(LogicalHduType::Image, LogicalHduType::Image);
        assert_ne!(LogicalHduType::Image, LogicalHduType::BinaryTable);
        assert_eq!(
            LogicalHduType::CompressedImage,
            LogicalHduType::CompressedImage
        );
        assert_eq!(
            LogicalHduType::Unknown("TEST".to_string()),
            LogicalHduType::Unknown("TEST".to_string())
        );
    }

    #[test]
    fn hdu_enum_hdu_type_primary() {
        let primary = create_primary_hdu();
        let hdu = Hdu::Primary(primary);

        assert_eq!(hdu.hdu_type(), HduType::Primary);
    }

    #[test]
    fn hdu_enum_hdu_type_image() {
        let image = create_image_hdu();
        let hdu = Hdu::Image(Box::new(image));

        assert_eq!(hdu.hdu_type(), HduType::Image);
    }

    #[test]
    fn hdu_enum_hdu_type_binary_table() {
        let table = create_binary_table_hdu(false);
        let hdu = Hdu::BinaryTable(table);

        assert_eq!(hdu.hdu_type(), HduType::BinaryTable);
    }

    #[test]
    fn hdu_enum_hdu_type_ascii_table() {
        let table = create_ascii_table_hdu();
        let hdu = Hdu::AsciiTable(table);

        assert_eq!(hdu.hdu_type(), HduType::AsciiTable);
    }

    #[test]
    fn hdu_enum_hdu_type_random_groups() {
        let groups = create_random_groups_hdu();
        let hdu = Hdu::RandomGroups(groups);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
    }

    #[test]
    fn hdu_enum_logical_type_primary_is_image() {
        let primary = create_primary_hdu();
        let hdu = Hdu::Primary(primary);

        assert_eq!(hdu.logical_type(), LogicalHduType::Image);
    }

    #[test]
    fn hdu_enum_logical_type_image_is_image() {
        let image = create_image_hdu();
        let hdu = Hdu::Image(Box::new(image));

        assert_eq!(hdu.logical_type(), LogicalHduType::Image);
    }

    #[test]
    fn hdu_enum_logical_type_binary_table_is_binary_table() {
        let table = create_binary_table_hdu(false);
        let hdu = Hdu::BinaryTable(table);

        assert_eq!(hdu.logical_type(), LogicalHduType::BinaryTable);
    }

    #[test]
    fn hdu_enum_logical_type_compressed_binary_table_is_compressed_image() {
        let table = create_binary_table_hdu(true);
        let hdu = Hdu::BinaryTable(table);

        assert_eq!(hdu.logical_type(), LogicalHduType::CompressedImage);
    }

    #[test]
    fn hdu_enum_logical_type_ascii_table_is_ascii_table() {
        let table = create_ascii_table_hdu();
        let hdu = Hdu::AsciiTable(table);

        assert_eq!(hdu.logical_type(), LogicalHduType::AsciiTable);
    }

    #[test]
    fn hdu_enum_logical_type_random_groups_is_image() {
        let groups = create_random_groups_hdu();
        let hdu = Hdu::RandomGroups(groups);

        assert_eq!(hdu.logical_type(), LogicalHduType::Image);
    }

    #[test]
    fn hdu_enum_header_access_all_variants() {
        let primary = create_primary_hdu();
        let image = create_image_hdu();
        let binary_table = create_binary_table_hdu(false);
        let ascii_table = create_ascii_table_hdu();
        let random_groups = create_random_groups_hdu();

        let hdu_primary = Hdu::Primary(primary);
        let hdu_image = Hdu::Image(Box::new(image));
        let hdu_binary = Hdu::BinaryTable(binary_table);
        let hdu_ascii = Hdu::AsciiTable(ascii_table);
        let hdu_random = Hdu::RandomGroups(random_groups);

        assert!(hdu_primary.header().get_keyword_value("SIMPLE").is_some());
        assert!(hdu_image.header().get_keyword_value("XTENSION").is_some());
        assert!(hdu_binary.header().get_keyword_value("TFIELDS").is_some());
        assert!(hdu_ascii.header().get_keyword_value("TFIELDS").is_some());
        assert!(hdu_random.header().get_keyword_value("GROUPS").is_some());
    }

    #[test]
    fn hdu_enum_info_access_all_variants() {
        let primary = create_primary_hdu();
        let image = create_image_hdu();
        let binary_table = create_binary_table_hdu(false);
        let ascii_table = create_ascii_table_hdu();
        let random_groups = create_random_groups_hdu();

        let hdu_primary = Hdu::Primary(primary);
        let hdu_image = Hdu::Image(Box::new(image));
        let hdu_binary = Hdu::BinaryTable(binary_table);
        let hdu_ascii = Hdu::AsciiTable(ascii_table);
        let hdu_random = Hdu::RandomGroups(random_groups);

        assert_eq!(hdu_primary.info().index, 0);
        assert_eq!(hdu_image.info().index, 1);
        assert_eq!(hdu_binary.info().index, 2);
        assert_eq!(hdu_ascii.info().index, 3);
        assert_eq!(hdu_random.info().index, 4);
    }

    #[test]
    fn hdu_trait_has_data_true_when_naxis_greater_than_zero() {
        let primary = create_primary_hdu();
        assert!(primary.has_data());
    }

    #[test]
    fn hdu_trait_has_data_false_when_naxis_zero() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert!(!primary.has_data());
    }

    #[test]
    fn hdu_trait_data_dimensions_returns_correct_values() {
        let primary = create_primary_hdu();
        assert_eq!(primary.data_dimensions(), vec![10, 10]);
    }

    #[test]
    fn hdu_trait_data_dimensions_returns_empty_for_zero_naxis() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert_eq!(primary.data_dimensions(), vec![]);
    }

    #[test]
    fn hdu_trait_bitpix_returns_correct_value() {
        let primary = create_primary_hdu();
        assert_eq!(primary.bitpix(), Some(crate::core::BitPix::I16));
    }

    #[test]
    fn hdu_trait_bitpix_returns_none_when_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert_eq!(primary.bitpix(), None);
    }

    #[test]
    fn hdu_trait_calculate_data_size_returns_correct_value() {
        let primary = create_primary_hdu();
        assert_eq!(primary.calculate_data_size().unwrap(), 10 * 10 * 2);
    }

    #[test]
    fn hdu_trait_calculate_data_size_returns_zero_for_no_data() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert_eq!(primary.calculate_data_size().unwrap(), 0);
    }

    #[test]
    fn hdu_trait_read_raw_data_returns_empty_when_no_data() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        let mut cursor = Cursor::new(vec![]);
        let result = primary.read_raw_data(&mut cursor).unwrap();
        assert_eq!(result, Vec::<u8>::new());
    }

    #[test]
    fn hdu_trait_read_raw_data_success() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 2));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 0,
            data_size: 4,
        };

        let primary = PrimaryHdu::new(header, info);
        let data = vec![0x01, 0x23, 0x45, 0x67];
        let mut cursor = Cursor::new(data.clone());
        let result = primary.read_raw_data(&mut cursor).unwrap();
        assert_eq!(result, data);
    }

    #[test]
    fn hdu_trait_write_raw_data_success() {
        let primary = create_primary_hdu();
        let mut buffer = Vec::new();
        let data = vec![0x01, 0x02, 0x03, 0x04];

        primary.write_raw_data(&mut buffer, &data).unwrap();
        assert_eq!(buffer, data);
    }

    #[test]
    fn hdu_trait_bzero_default_value() {
        let primary = create_primary_hdu();
        assert_eq!(primary.bzero(), 0.0);
    }

    #[test]
    fn hdu_trait_bzero_from_header() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::real("BZERO", 32768.0));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert_eq!(primary.bzero(), 32768.0);
    }

    #[test]
    fn hdu_trait_bscale_default_value() {
        let primary = create_primary_hdu();
        assert_eq!(primary.bscale(), 1.0);
    }

    #[test]
    fn hdu_trait_bscale_from_header() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::real("BSCALE", 2.5));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert_eq!(primary.bscale(), 2.5);
    }

    #[test]
    fn hdu_trait_needs_scaling_false_default() {
        let primary = create_primary_hdu();
        assert!(!primary.needs_scaling());
    }

    #[test]
    fn hdu_trait_needs_scaling_true_with_bzero() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::real("BZERO", 32768.0));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert!(primary.needs_scaling());
    }

    #[test]
    fn hdu_trait_needs_scaling_true_with_bscale() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::real("BSCALE", 2.0));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert!(primary.needs_scaling());
    }

    #[test]
    fn hdu_trait_read_data_with_scaling_integer_returns_raw() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 2));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::real("BSCALE", 2.0));
        header.add_keyword(Keyword::real("BZERO", 100.0));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 0,
            data_size: 4,
        };

        let primary = PrimaryHdu::new(header, info);
        let data = vec![0x00, 0x0A, 0x00, 0x14];
        let mut cursor = Cursor::new(data);
        let result: Vec<i16> = primary.read_data(&mut cursor).unwrap();
        assert_eq!(result, vec![10, 20]);
    }

    #[test]
    fn hdu_trait_read_data_without_scaling() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 2));
        header.add_keyword(Keyword::integer("BITPIX", 16));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 0,
            data_size: 4,
        };

        let primary = PrimaryHdu::new(header, info);
        let data = vec![0x00, 0x0A, 0x00, 0x14];
        let mut cursor = Cursor::new(data);
        let result: Vec<i16> = primary.read_data(&mut cursor).unwrap();
        assert_eq!(result, vec![10, 20]);
    }

    #[test]
    fn hdu_trait_bzero_integer_value_converted_to_real() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::integer("BZERO", 32768));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };

        let primary = PrimaryHdu::new(header, info);
        assert_eq!(primary.bzero(), 32768.0);
    }
}
