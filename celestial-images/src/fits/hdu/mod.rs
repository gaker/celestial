pub mod ascii_table;
pub mod binary_table;
pub mod image;
pub mod primary;
pub mod random_groups;
pub mod table;
#[cfg(test)]
mod table_tests;

use ascii_table::AsciiTableHdu;
use binary_table::BinaryTableHdu;
use image::ImageHdu;
use primary::PrimaryHdu;
use random_groups::RandomGroupsHdu;

use crate::core::types::{BitPix, ByteOrder};
use crate::fits::data::array::DataArray;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;
use crate::fits::io::reader::data_size::data_unit_len;
use crate::fits::io::reader::HduInfo;
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
    pub(crate) fn new(index: usize, format: String) -> Self {
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
}

pub(crate) trait HduTrait: std::fmt::Debug {
    fn header(&self) -> &Header;
    fn info(&self) -> &HduInfo;
    fn hdu_type(&self) -> HduType;

    fn has_data(&self) -> bool {
        self.header().get("NAXIS").as_i64().unwrap_or(0) > 0
    }

    // A missing NAXIS or NAXISn reads as 0. NAXIS outside 0..=999 and
    // negative axis lengths are errors (FITS 4.0 §4.4.1.1).
    fn data_dimensions(&self) -> Result<Vec<usize>> {
        let naxis = self.header().get("NAXIS").as_i64().unwrap_or(0);
        if !(0..=999).contains(&naxis) {
            return Err(invalid_value("NAXIS", naxis));
        }
        (1..=naxis)
            .map(|i| {
                let name = format!("NAXIS{i}");
                let length = self.header().get(&name).as_i64().unwrap_or(0);
                usize::try_from(length).map_err(|_| invalid_value(&name, length))
            })
            .collect()
    }

    fn bitpix(&self) -> Option<BitPix> {
        self.header()
            .get("BITPIX")
            .as_i64()
            .and_then(|i| i32::try_from(i).ok())
            .and_then(BitPix::from_value)
    }

    // The whole data unit before padding, with PCOUNT and GCOUNT (FITS 4.0
    // eq. 2): a binary table's heap, and every random group.
    fn calculate_data_size(&self) -> Result<usize> {
        data_unit_len(self.header())
    }

    fn checked_data_size(&self) -> Result<usize> {
        within_data_unit(self, self.calculate_data_size()?)
    }

    fn read_raw_data<R>(&self, reader: &mut R) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        if !self.has_data() {
            return Ok(Vec::new());
        }
        read_bytes(self, reader, self.checked_data_size()?)
    }

    // Stored values for every BITPIX: the image array of a primary or IMAGE
    // HDU, and the whole data unit of other HDUs. BSCALE and BZERO are
    // applied by the caller, which can pick a type that holds the physical
    // values.
    fn read_data<T, R>(&self, reader: &mut R) -> Result<Vec<T>>
    where
        T: DataArray,
        R: Read + Seek,
    {
        if !self.has_data() {
            return Ok(Vec::new());
        }
        stored_as(self, T::BITPIX)?;
        let bytes = match self.hdu_type() {
            HduType::Primary | HduType::Image => read_array(self, reader)?,
            _ => self.read_raw_data(reader)?,
        };
        T::from_bytes(&bytes, ByteOrder::BigEndian)
    }
}

// The image array at the start of a primary or IMAGE data unit: NAXIS1 x
// ... x NAXISn values of |BITPIX|/8 bytes (FITS 4.0 eq. 1). PCOUNT and
// GCOUNT do not change it.
pub(crate) fn read_array<H, R>(hdu: &H, reader: &mut R) -> Result<Vec<u8>>
where
    H: HduTrait + ?Sized,
    R: Read + Seek,
{
    if !hdu.has_data() {
        return Ok(Vec::new());
    }
    let bitpix = hdu.bitpix().ok_or_else(|| FitsError::KeywordNotFound {
        keyword: "BITPIX".to_string(),
    })?;
    let size = hdu
        .data_dimensions()?
        .iter()
        .try_fold(bitpix.bytes_per_pixel(), |n, &d| n.checked_mul(d))
        .ok_or_else(|| FitsError::InvalidFormat("Data dimensions too large".to_string()))?;
    read_bytes(hdu, reader, within_data_unit(hdu, size)?)
}

// The scan checked `info().data_size` against the file length, so a larger
// size here comes from a header that does not match the scan.
fn within_data_unit<H: HduTrait + ?Sized>(hdu: &H, size: usize) -> Result<usize> {
    let info = hdu.info();
    if size > info.data_size {
        return Err(FitsError::InvalidFormat(format!(
            "HDU {} needs {} data bytes, but its data unit holds {}",
            info.index, size, info.data_size
        )));
    }
    Ok(size)
}

fn read_bytes<H, R>(hdu: &H, reader: &mut R, size: usize) -> Result<Vec<u8>>
where
    H: HduTrait + ?Sized,
    R: Read + Seek,
{
    reader.seek(SeekFrom::Start(hdu.info().data_start))?;
    let mut buffer = vec![0u8; size];
    reader.read_exact(&mut buffer)?;
    Ok(buffer)
}

fn stored_as<H: HduTrait + ?Sized>(hdu: &H, expected: BitPix) -> Result<()> {
    let actual = hdu.bitpix().ok_or_else(|| FitsError::KeywordNotFound {
        keyword: "BITPIX".to_string(),
    })?;
    if actual != expected {
        return Err(FitsError::TypeMismatch { expected, actual });
    }
    Ok(())
}

fn invalid_value(keyword: &str, value: i64) -> FitsError {
    FitsError::InvalidKeywordValue {
        keyword: keyword.to_string(),
        value: value.to_string(),
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
            Self::Primary(_) => HduType::Primary,
            Self::Image(_) => HduType::Image,
            Self::AsciiTable(_) => HduType::AsciiTable,
            Self::BinaryTable(_) => HduType::BinaryTable,
            Self::RandomGroups(_) => HduType::RandomGroups,
        }
    }

    pub fn header(&self) -> &Header {
        match self {
            Self::Primary(hdu) => hdu.header(),
            Self::Image(hdu) => hdu.header(),
            Self::AsciiTable(hdu) => hdu.header(),
            Self::BinaryTable(hdu) => hdu.header(),
            Self::RandomGroups(hdu) => hdu.header(),
        }
    }

    pub fn info(&self) -> &HduInfo {
        match self {
            Self::Primary(hdu) => hdu.info(),
            Self::Image(hdu) => hdu.info(),
            Self::AsciiTable(hdu) => hdu.info(),
            Self::BinaryTable(hdu) => hdu.info(),
            Self::RandomGroups(hdu) => hdu.info(),
        }
    }

    pub fn logical_type(&self) -> LogicalHduType {
        match self {
            Self::Primary(_) => LogicalHduType::Image,
            Self::Image(_) => LogicalHduType::Image,
            Self::AsciiTable(_) => LogicalHduType::AsciiTable,
            Self::BinaryTable(hdu) => {
                if hdu.is_compressed_image() {
                    LogicalHduType::CompressedImage
                } else {
                    LogicalHduType::BinaryTable
                }
            }
            Self::RandomGroups(_) => LogicalHduType::Image,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::Keyword;
    use crate::fits::header::parser::Header;
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
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::string("TFORM2", "4E"));

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

        BinaryTableHdu::new(header, info).unwrap()
    }

    fn create_ascii_table_hdu() -> AsciiTableHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 80));
        header.add_keyword(Keyword::integer("NAXIS2", 50));
        header.add_keyword(Keyword::integer("TFIELDS", 3));
        for (n, col) in [(1, 1), (2, 21), (3, 41)] {
            header.add_keyword(Keyword::string(format!("TFORM{n}"), "A20".to_string()));
            header.add_keyword(Keyword::integer(format!("TBCOL{n}"), col));
        }

        let info = HduInfo {
            index: 3,
            header_start: 8640,
            header_size: 2880,
            data_start: 11520,
            data_size: 4000,
        };

        AsciiTableHdu::new(header, info).unwrap()
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

        assert!(hdu_primary.header().get("SIMPLE").exists());
        assert!(hdu_image.header().get("XTENSION").exists());
        assert!(hdu_binary.header().get("TFIELDS").exists());
        assert!(hdu_ascii.header().get("TFIELDS").exists());
        assert!(hdu_random.header().get("GROUPS").exists());
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
        assert_eq!(primary.data_dimensions().unwrap(), vec![10, 10]);
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
        assert_eq!(primary.data_dimensions().unwrap(), vec![]);
    }

    #[test]
    fn hdu_trait_bitpix_returns_correct_value() {
        let primary = create_primary_hdu();
        assert_eq!(primary.bitpix(), Some(crate::core::types::BitPix::I16));
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
    fn hdu_trait_read_data_with_scaling_float_returns_stored_values() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 2));
        header.add_keyword(Keyword::integer("BITPIX", -32));
        header.add_keyword(Keyword::real("BSCALE", 2.0));
        header.add_keyword(Keyword::real("BZERO", 1.0));

        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 0,
            data_size: 8,
        };

        let primary = PrimaryHdu::new(header, info);
        let data: Vec<u8> = [2.5f32, -1.0]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let mut cursor = Cursor::new(data);
        let result: Vec<f32> = primary.read_data(&mut cursor).unwrap();
        assert_eq!(result, vec![2.5, -1.0]);
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

    fn huge_primary_hdu() -> PrimaryHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 1_000_000_000));
        header.add_keyword(Keyword::integer("NAXIS2", 1_000_000_000));
        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 0,
            data_size: 2880,
        };
        PrimaryHdu::new(header, info)
    }

    #[test]
    fn hdu_trait_read_data_rejects_size_beyond_data_unit() {
        let primary = huge_primary_hdu();
        let result: Result<Vec<i16>> = primary.read_data(&mut Cursor::new(vec![0u8; 2880]));
        assert!(matches!(result, Err(FitsError::InvalidFormat(_))));
        let result = primary.read_raw_data(&mut Cursor::new(vec![0u8; 2880]));
        assert!(matches!(result, Err(FitsError::InvalidFormat(_))));
    }

    fn primary_with(cards: &[(&str, i64)]) -> PrimaryHdu {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        for &(name, value) in cards {
            header.add_keyword(Keyword::integer(name, value));
        }
        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };
        PrimaryHdu::new(header, info)
    }

    #[test]
    fn hdu_trait_data_dimensions_rejects_a_negative_axis_length() {
        let cards = [("BITPIX", 16), ("NAXIS", 2), ("NAXIS1", 3), ("NAXIS2", -1)];
        assert!(matches!(
            primary_with(&cards).data_dimensions(),
            Err(FitsError::InvalidKeywordValue { ref keyword, .. }) if keyword == "NAXIS2"
        ));
    }

    #[test]
    fn hdu_trait_data_dimensions_rejects_naxis_outside_0_to_999() {
        for naxis in [-1, 1000] {
            assert!(matches!(
                primary_with(&[("BITPIX", 16), ("NAXIS", naxis)]).data_dimensions(),
                Err(FitsError::InvalidKeywordValue { ref keyword, .. }) if keyword == "NAXIS"
            ));
        }
    }

    #[test]
    fn hdu_trait_bitpix_rejects_value_outside_i32() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("BITPIX", (1 << 32) + 8));
        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 0,
        };
        let primary = PrimaryHdu::new(header, info);
        assert_eq!(HduTrait::bitpix(&primary), None);
        assert_eq!(primary.bitpix(), None);
    }
}
