mod column_ops;
mod compression;
mod format_parsing;

pub use column_ops::BinaryTableRowIterator;

use super::{HduTrait, HduType};
use crate::fits::header::Header;
use crate::fits::io::reader::HduInfo;
use std::collections::HashMap;
use std::sync::OnceLock;

#[derive(Debug)]
pub struct BinaryTableHdu {
    header: Header,
    info: HduInfo,
    column_name_index: OnceLock<HashMap<String, usize>>,
}

impl BinaryTableHdu {
    pub fn new(header: Header, info: HduInfo) -> Self {
        Self {
            header,
            info,
            column_name_index: OnceLock::new(),
        }
    }

    pub fn number_of_fields(&self) -> Option<i64> {
        self.header
            .get_keyword_value("TFIELDS")
            .and_then(|v| v.as_integer())
    }

    pub fn number_of_rows(&self) -> Option<i64> {
        self.header
            .get_keyword_value("NAXIS2")
            .and_then(|v| v.as_integer())
    }

    pub fn extension_name(&self) -> Option<&str> {
        self.header
            .get_keyword_value("EXTNAME")
            .and_then(|v| v.as_string())
    }

    pub fn extension_version(&self) -> Option<i64> {
        self.header
            .get_keyword_value("EXTVER")
            .and_then(|v| v.as_integer())
    }
}

impl HduTrait for BinaryTableHdu {
    fn header(&self) -> &Header {
        &self.header
    }

    fn info(&self) -> &HduInfo {
        &self.info
    }

    fn hdu_type(&self) -> HduType {
        HduType::BinaryTable
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::Keyword;

    fn create_test_header(extname: Option<&str>, compressed: bool) -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 100));
        header.add_keyword(Keyword::integer("NAXIS2", 50));
        header.add_keyword(Keyword::integer("TFIELDS", 3));

        if let Some(name) = extname {
            header.add_keyword(Keyword::string("EXTNAME", name));
            header.add_keyword(Keyword::integer("EXTVER", 1));
        }

        if compressed {
            header.add_keyword(Keyword::logical("ZIMAGE", true));
            header.add_keyword(Keyword::string("ZCMPTYPE", "RICE_1"));
            header.add_keyword(Keyword::integer("ZQUANTIZ", 16));
        }

        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 5000,
        }
    }

    #[test]
    fn new_creates_binary_table_hdu() {
        let header = create_test_header(Some("TEST"), false);
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.info.index, 1);
        assert_eq!(hdu.number_of_fields(), Some(3));
    }

    #[test]
    fn header_returns_header_reference() {
        let header = create_test_header(None, false);
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let header_ref = hdu.header();
        assert_eq!(
            header_ref
                .get_keyword_value("XTENSION")
                .unwrap()
                .as_string()
                .unwrap(),
            "BINTABLE"
        );
    }

    #[test]
    fn hdu_type_returns_binary_table() {
        let header = create_test_header(None, false);
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::BinaryTable);
    }
}
