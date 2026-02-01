mod column_ops;
mod compression;
mod format_parsing;

pub use column_ops::BinaryTableRowIterator;

#[cfg(test)]
mod tests;

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
