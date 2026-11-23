use super::{HduTrait, HduType};
use crate::core::types::BitPix;
use crate::fits::data::scaling::{integer_value, Linear};
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;
use crate::fits::io::reader::HduInfo;
use std::io::{Read, Seek, SeekFrom};

#[derive(Debug)]
pub struct RandomGroupsHdu {
    header: Header,
    info: HduInfo,
}

// One group of FITS 4.0 §6: the parameters scaled by PSCALn and PZEROn,
// and the array scaled by BSCALE and BZERO. Integer array values equal to
// BLANK are NaN.
#[derive(Debug, Clone, PartialEq)]
pub struct RandomGroup {
    pub parameters: Vec<f64>,
    pub array: Vec<f64>,
}

// The stored sizes of a random-groups data unit, in values.
struct Shape {
    bitpix: BitPix,
    parameters: usize,
    array: usize,
    groups: usize,
}

impl RandomGroupsHdu {
    pub(crate) fn new(header: Header, info: HduInfo) -> Self {
        Self { header, info }
    }

    pub fn group_count(&self) -> Option<i64> {
        self.header.get("GCOUNT").as_i64()
    }

    pub fn parameter_count(&self) -> Option<i64> {
        self.header.get("PCOUNT").as_i64()
    }

    pub fn extension_name(&self) -> Option<&str> {
        self.header.get("EXTNAME").as_str()
    }

    pub fn extension_version(&self) -> Option<i64> {
        self.header.get("EXTVER").as_i64()
    }

    pub fn read_groups<R: Read + Seek>(&self, reader: &mut R) -> Result<Vec<RandomGroup>> {
        let shape = self.shape()?;
        let size = self.checked_data_size()?;
        if shape.groups == 0 {
            return Ok(Vec::new());
        }
        // With at least one group, the size check bounds PCOUNT by the data
        // unit, so the scalings below cost no more than the data does.
        let group_len = match shape.parameters.checked_add(shape.array) {
            Some(0) | None => {
                return Err(invalid(
                    "random groups hold no values, or too many to count",
                ))
            }
            Some(n) => n,
        };
        let scaling = GroupScaling::from_header(&self.header, shape.parameters)?;
        reader.seek(SeekFrom::Start(self.info.data_start))?;
        let mut bytes = vec![0u8; size];
        reader.read_exact(&mut bytes)?;
        Ok(stored(&bytes, shape.bitpix)
            .chunks_exact(group_len)
            .map(|group| scaling.group(group))
            .collect())
    }

    // NAXIS1 = 0, and NAXIS2 ... NAXISn give the shape of each array.
    fn shape(&self) -> Result<Shape> {
        let count = |name: &str| super::table::layout::count(&self.header, name);
        let axes = count("NAXIS")?;
        let array = (2..=axes)
            .try_fold(1usize, |n, i| {
                n.checked_mul(count(&format!("NAXIS{i}")).ok()?)
            })
            .ok_or_else(|| invalid("random-groups axes are missing or too large"))?;
        let bitpix = self.bitpix().ok_or_else(|| FitsError::KeywordNotFound {
            keyword: "BITPIX".to_string(),
        })?;
        Ok(Shape {
            bitpix,
            parameters: count("PCOUNT")?,
            array,
            groups: count("GCOUNT")?,
        })
    }
}

// PSCALn and PZEROn for each parameter, and BSCALE, BZERO and BLANK for
// the array.
struct GroupScaling {
    parameters: Vec<Linear>,
    array: Linear,
    blank: Option<i64>,
}

impl GroupScaling {
    fn from_header(header: &Header, parameters: usize) -> Result<Self> {
        Ok(Self {
            parameters: (1..=parameters)
                .map(|n| Linear::from_header(header, &format!("PSCAL{n}"), &format!("PZERO{n}")))
                .collect::<Result<_>>()?,
            array: Linear::from_header(header, "BSCALE", "BZERO")?,
            blank: header.get_last("BLANK").value().and_then(integer_value),
        })
    }

    // `values` holds every parameter, then the array.
    fn group(&self, values: &[Stored]) -> RandomGroup {
        let (parameters, array) = values.split_at(self.parameters.len().min(values.len()));
        RandomGroup {
            parameters: parameters
                .iter()
                .zip(&self.parameters)
                .map(|(v, s)| v.physical(s, None))
                .collect(),
            array: array
                .iter()
                .map(|v| v.physical(&self.array, self.blank))
                .collect(),
        }
    }
}

#[derive(Clone, Copy)]
enum Stored {
    Integer(i64),
    Real(f64),
}

impl Stored {
    fn physical(self, scale: &Linear, blank: Option<i64>) -> f64 {
        match self {
            Self::Integer(v) if Some(v) == blank => f64::NAN,
            Self::Integer(v) => scale.apply(v as f64),
            Self::Real(v) => scale.apply(v),
        }
    }
}

fn stored(bytes: &[u8], bitpix: BitPix) -> Vec<Stored> {
    let size = bitpix.bytes_per_pixel();
    bytes
        .chunks_exact(size)
        .map(|b| {
            let mut be = [0u8; 8];
            be[8 - size..].copy_from_slice(b);
            let wide = u64::from_be_bytes(be);
            match bitpix {
                BitPix::U8 => Stored::Integer(wide as i64),
                BitPix::I16 => Stored::Integer(i64::from(wide as u16 as i16)),
                BitPix::I32 => Stored::Integer(i64::from(wide as u32 as i32)),
                BitPix::I64 => Stored::Integer(wide as i64),
                BitPix::F32 => Stored::Real(f64::from(f32::from_bits(wide as u32))),
                BitPix::F64 => Stored::Real(f64::from_bits(wide)),
            }
        })
        .collect()
}

fn invalid(message: &str) -> FitsError {
    FitsError::InvalidFormat(message.to_string())
}

impl HduTrait for RandomGroupsHdu {
    fn header(&self) -> &Header {
        &self.header
    }

    fn info(&self) -> &HduInfo {
        &self.info
    }

    fn hdu_type(&self) -> HduType {
        HduType::RandomGroups
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::Keyword;
    use crate::fits::header::parser::Header;
    use crate::fits::io::reader::HduInfo;

    fn create_test_header(extname: Option<&str>, groups: bool) -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 3));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("NAXIS2", 10));
        header.add_keyword(Keyword::integer("NAXIS3", 5));
        header.add_keyword(Keyword::integer("BITPIX", -32));

        if groups {
            header.add_keyword(Keyword::logical("GROUPS", true));
            header.add_keyword(Keyword::integer("GCOUNT", 10));
            header.add_keyword(Keyword::integer("PCOUNT", 3));
        }

        if let Some(name) = extname {
            header.add_keyword(Keyword::string("EXTNAME", name));
            header.add_keyword(Keyword::integer("EXTVER", 1));
        }

        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 1200,
        }
    }

    #[test]
    fn new_creates_random_groups_hdu() {
        let header = create_test_header(Some("UVDATA"), true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.info.index, 0);
        assert_eq!(hdu.group_count(), Some(10));
    }

    #[test]
    fn header_returns_header_reference() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        let header_ref = hdu.header();
        assert!(header_ref.get("SIMPLE").as_bool().unwrap());
    }

    #[test]
    fn info_returns_info_reference() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        let info_ref = hdu.info();
        assert_eq!(info_ref.index, 0);
        assert_eq!(info_ref.data_start, 2880);
    }

    #[test]
    fn hdu_type_returns_random_groups() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
    }

    #[test]
    fn group_count_returns_gcount_value() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.group_count(), Some(10));
    }

    #[test]
    fn group_count_returns_none_when_missing() {
        let header = create_test_header(None, false);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.group_count(), None);
    }

    #[test]
    fn parameter_count_returns_pcount_value() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.parameter_count(), Some(3));
    }

    #[test]
    fn parameter_count_returns_none_when_missing() {
        let header = create_test_header(None, false);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.parameter_count(), None);
    }

    #[test]
    fn extension_name_returns_extname_value() {
        let header = create_test_header(Some("VISIBILITY"), true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_name(), Some("VISIBILITY"));
    }

    #[test]
    fn extension_name_returns_none_when_missing() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_name(), None);
    }

    #[test]
    fn extension_version_returns_extver_value() {
        let header = create_test_header(Some("TEST"), true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_version(), Some(1));
    }

    #[test]
    fn extension_version_returns_none_when_missing() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_version(), None);
    }

    #[test]
    fn all_methods_work_together() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 4));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        header.add_keyword(Keyword::integer("NAXIS3", 1024));
        header.add_keyword(Keyword::integer("NAXIS4", 1024));
        header.add_keyword(Keyword::integer("BITPIX", -64));
        header.add_keyword(Keyword::logical("GROUPS", true));
        header.add_keyword(Keyword::integer("GCOUNT", 100));
        header.add_keyword(Keyword::integer("PCOUNT", 5));
        header.add_keyword(Keyword::string("EXTNAME", "RADIODATA"));
        header.add_keyword(Keyword::integer("EXTVER", 2));
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
        assert_eq!(hdu.group_count(), Some(100));
        assert_eq!(hdu.parameter_count(), Some(5));
        assert_eq!(hdu.extension_name(), Some("RADIODATA"));
        assert_eq!(hdu.extension_version(), Some(2));
    }

    #[test]
    fn minimal_valid_random_groups() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::logical("GROUPS", true));
        header.add_keyword(Keyword::integer("GCOUNT", 1));
        header.add_keyword(Keyword::integer("PCOUNT", 0));
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
        assert_eq!(hdu.group_count(), Some(1));
        assert_eq!(hdu.parameter_count(), Some(0));
        assert_eq!(hdu.extension_name(), None);
        assert_eq!(hdu.extension_version(), None);
    }
}
