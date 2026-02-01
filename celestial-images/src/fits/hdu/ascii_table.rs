use super::{ColumnInfo, HduTrait, HduType};
use crate::core::ByteOrder;
use crate::fits::data::array::{DataArray, DataValue, TableValue};
use crate::fits::header::Header;
use crate::fits::io::reader::HduInfo;
use crate::fits::{FitsError, Result};
use std::io::{Read, Seek, SeekFrom};
use std::marker::PhantomData;

#[derive(Debug)]
pub struct AsciiTableHdu {
    header: Header,
    info: HduInfo,
}

impl AsciiTableHdu {
    pub fn new(header: Header, info: HduInfo) -> Self {
        Self { header, info }
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

    pub fn column_count(&self) -> Result<usize> {
        self.header
            .get_keyword_value("TFIELDS")
            .and_then(|v| v.as_integer())
            .map(|n| n as usize)
            .ok_or_else(|| FitsError::KeywordNotFound {
                keyword: "TFIELDS".to_string(),
            })
    }

    pub fn column_info(&self, column: usize) -> Result<ColumnInfo> {
        let column_count = self.column_count()?;
        if column >= column_count {
            return Err(FitsError::InvalidFormat(format!(
                "Column index {} out of range (0..{})",
                column, column_count
            )));
        }

        let column_index = column + 1;

        let format_key = format!("TFORM{}", column_index);
        let format = self
            .header
            .get_keyword_value(&format_key)
            .and_then(|v| v.as_string())
            .ok_or(FitsError::KeywordNotFound {
                keyword: format_key,
            })?;

        let mut info = ColumnInfo::new(column, format.to_string());

        if let Some(name) = self
            .header
            .get_keyword_value(&format!("TTYPE{}", column_index))
            .and_then(|v| v.as_string())
        {
            info = info.with_name(name.to_string());
        }

        if let Some(unit) = self
            .header
            .get_keyword_value(&format!("TUNIT{}", column_index))
            .and_then(|v| v.as_string())
        {
            info = info.with_unit(unit.to_string());
        }

        if let Some(null_val) = self
            .header
            .get_keyword_value(&format!("TNULL{}", column_index))
            .and_then(|v| v.as_string())
        {
            info = info.with_null_value(null_val.to_string());
        }

        if let Some(scale) = self
            .header
            .get_keyword_value(&format!("TSCAL{}", column_index))
            .and_then(|v| v.as_real())
        {
            info = info.with_scale(scale);
        }

        if let Some(zero) = self
            .header
            .get_keyword_value(&format!("TZERO{}", column_index))
            .and_then(|v| v.as_real())
        {
            info = info.with_zero_offset(zero);
        }

        if let Some(disp) = self
            .header
            .get_keyword_value(&format!("TDISP{}", column_index))
            .and_then(|v| v.as_string())
        {
            info.display_format = Some(disp.to_string());
        }

        Ok(info)
    }

    pub fn column_by_name(&self, name: &str) -> Result<usize> {
        let column_count = self.column_count()?;

        for i in 0..column_count {
            if let Ok(info) = self.column_info(i) {
                if let Some(col_name) = &info.name {
                    if col_name == name {
                        return Ok(i);
                    }
                }
            }
        }

        Err(FitsError::InvalidFormat(format!(
            "Column '{}' not found",
            name
        )))
    }

    pub fn all_column_info(&self) -> Result<Vec<ColumnInfo>> {
        let column_count = self.column_count()?;
        let mut columns = Vec::with_capacity(column_count);

        for i in 0..column_count {
            columns.push(self.column_info(i)?);
        }

        Ok(columns)
    }

    pub fn read_column_raw<R>(&self, reader: &mut R, column: usize) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        let info = self.column_info(column)?;
        let row_count = self.number_of_rows().unwrap_or(0) as usize;

        if row_count == 0 {
            return Ok(Vec::new());
        }

        let (_data_type, width) = self.parse_ascii_format(&info.format)?;
        let total_bytes = row_count * width;

        let column_offset = self.calculate_column_offset(column)?;
        let row_size = self.get_row_size()?;

        let data_start = self.info.data_start;
        let mut result = Vec::with_capacity(total_bytes);

        for row in 0..row_count {
            let row_start = data_start + (row * row_size) as u64;
            let column_position = row_start + column_offset as u64;

            reader.seek(SeekFrom::Start(column_position))?;

            let mut column_data = vec![0u8; width];
            reader.read_exact(&mut column_data)?;

            result.extend_from_slice(&column_data);
        }

        Ok(result)
    }

    fn calculate_column_offset(&self, column: usize) -> Result<usize> {
        let column_count = self.column_count()?;
        if column >= column_count {
            return Err(FitsError::InvalidFormat(format!(
                "Column index {} out of range (0..{})",
                column, column_count
            )));
        }

        if column == 0 {
            return Ok(0);
        }

        let mut offset = 0;
        for i in 0..column {
            let info = self.column_info(i)?;
            let (_, width) = self.parse_ascii_format(&info.format)?;
            offset += width + 1;
        }

        Ok(offset)
    }

    fn get_row_size(&self) -> Result<usize> {
        let column_count = self.column_count()?;
        let mut row_size = 0;

        for i in 0..column_count {
            let info = self.column_info(i)?;
            let (_, width) = self.parse_ascii_format(&info.format)?;
            row_size += width;
            if i < column_count - 1 {
                row_size += 1;
            }
        }

        row_size += 1;
        Ok(row_size)
    }

    fn parse_ascii_format(&self, format: &str) -> Result<(String, usize)> {
        if format.is_empty() {
            return Err(FitsError::InvalidFormat("Empty column format".to_string()));
        }

        let mut data_type = String::new();
        let mut width_str = String::new();
        let mut in_width = false;

        for ch in format.chars() {
            if ch.is_ascii_alphabetic() && !in_width {
                data_type.push(ch);
                in_width = true;
            } else if ch.is_ascii_digit() && in_width {
                width_str.push(ch);
            } else if ch == '.' {
                break;
            }
        }

        let width = if width_str.is_empty() {
            1
        } else {
            width_str.parse().unwrap_or(1)
        };

        Ok((data_type, width))
    }

    pub fn read_column_with_nulls<T, R>(
        &self,
        reader: &mut R,
        column: usize,
    ) -> Result<Vec<DataValue<T>>>
    where
        T: DataArray,
        R: Read + Seek,
    {
        let column_info = self.column_info(column)?;
        let raw_data = self.read_column_raw(reader, column)?;

        let null_value = match &column_info.null_value {
            Some(null_str) => Some(T::parse_null_value(null_str)?),
            None => None,
        };

        T::from_bytes_with_null(&raw_data, ByteOrder::BigEndian, null_value)
    }

    pub fn get_row<R>(&self, reader: &mut R, row_index: usize) -> Result<Vec<TableValue>>
    where
        R: Read + Seek,
    {
        let row_count = self.number_of_rows().unwrap_or(0) as usize;
        if row_index >= row_count {
            return Err(FitsError::InvalidFormat(format!(
                "Row index {} out of range (0..{})",
                row_index, row_count
            )));
        }

        let column_count = self.column_count()?;
        let row_bytes = self.read_row_raw(reader, row_index)?;
        self.parse_row_values(&row_bytes, column_count)
    }

    fn read_row_raw<R>(&self, reader: &mut R, row_index: usize) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        let row_size = self.get_row_size()?;
        let data_start = self.info.data_start;
        let row_offset = row_index * row_size;

        reader.seek(SeekFrom::Start(data_start + row_offset as u64))?;

        let mut row_buffer = vec![0u8; row_size];
        reader.read_exact(&mut row_buffer)?;
        Ok(row_buffer)
    }

    fn parse_row_values(&self, row_bytes: &[u8], column_count: usize) -> Result<Vec<TableValue>> {
        let mut values = Vec::with_capacity(column_count);
        let mut offset = 0;

        for col_idx in 0..column_count {
            let info = self.column_info(col_idx)?;
            let (data_type, width) = self.parse_ascii_format(&info.format)?;
            let col_data = &row_bytes[offset..offset + width];
            let value = self.parse_ascii_value(&data_type, col_data)?;
            values.push(value);
            offset += width;
            if col_idx < column_count - 1 {
                offset += 1;
            }
        }

        Ok(values)
    }

    fn parse_ascii_value(&self, data_type: &str, bytes: &[u8]) -> Result<TableValue> {
        let text = String::from_utf8_lossy(bytes).trim().to_string();
        if text.is_empty() {
            return Ok(TableValue::Null);
        }

        match data_type {
            "A" => Ok(TableValue::String(text)),
            "I" => self.parse_ascii_integer(&text),
            "F" | "E" | "D" => self.parse_ascii_float(&text),
            _ => Ok(TableValue::String(text)),
        }
    }

    fn parse_ascii_integer(&self, text: &str) -> Result<TableValue> {
        match text.parse::<i64>() {
            Ok(v) => Ok(TableValue::I64(v)),
            Err(_) => Ok(TableValue::Null),
        }
    }

    fn parse_ascii_float(&self, text: &str) -> Result<TableValue> {
        let normalized = text.replace('D', "E").replace('d', "e");
        match normalized.parse::<f64>() {
            Ok(v) => Ok(TableValue::F64(v)),
            Err(_) => Ok(TableValue::Null),
        }
    }

    pub fn get_column_by_name<R>(&self, reader: &mut R, name: &str) -> Result<Vec<TableValue>>
    where
        R: Read + Seek,
    {
        let col_idx = self.column_by_name(name)?;
        self.get_column_values(reader, col_idx)
    }

    pub fn get_column_values<R>(&self, reader: &mut R, column: usize) -> Result<Vec<TableValue>>
    where
        R: Read + Seek,
    {
        let row_count = self.number_of_rows().unwrap_or(0) as usize;
        if row_count == 0 {
            return Ok(Vec::new());
        }

        let info = self.column_info(column)?;
        let (data_type, width) = self.parse_ascii_format(&info.format)?;
        let raw_data = self.read_column_raw(reader, column)?;

        self.convert_column_to_values(&data_type, &raw_data, width, row_count)
    }

    fn convert_column_to_values(
        &self,
        data_type: &str,
        raw_data: &[u8],
        width: usize,
        row_count: usize,
    ) -> Result<Vec<TableValue>> {
        let mut values = Vec::with_capacity(row_count);

        for row in 0..row_count {
            let start = row * width;
            let end = start + width;
            let row_bytes = &raw_data[start..end];
            let value = self.parse_ascii_value(data_type, row_bytes)?;
            values.push(value);
        }

        Ok(values)
    }

    pub fn row_count(&self) -> usize {
        self.number_of_rows().unwrap_or(0) as usize
    }

    pub fn iter_rows<'a, R>(&'a self, reader: &'a mut R) -> AsciiTableRowIterator<'a, R>
    where
        R: Read + Seek,
    {
        AsciiTableRowIterator::new(self, reader)
    }
}

pub struct AsciiTableRowIterator<'a, R>
where
    R: Read + Seek,
{
    hdu: &'a AsciiTableHdu,
    reader: &'a mut R,
    current_row: usize,
    row_count: usize,
    _phantom: PhantomData<R>,
}

impl<'a, R> AsciiTableRowIterator<'a, R>
where
    R: Read + Seek,
{
    fn new(hdu: &'a AsciiTableHdu, reader: &'a mut R) -> Self {
        let row_count = hdu.row_count();
        Self {
            hdu,
            reader,
            current_row: 0,
            row_count,
            _phantom: PhantomData,
        }
    }
}

impl<'a, R> Iterator for AsciiTableRowIterator<'a, R>
where
    R: Read + Seek,
{
    type Item = Result<Vec<TableValue>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_row >= self.row_count {
            return None;
        }

        let result = self.hdu.get_row(self.reader, self.current_row);
        self.current_row += 1;
        Some(result)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.row_count - self.current_row;
        (remaining, Some(remaining))
    }
}

impl<'a, R> ExactSizeIterator for AsciiTableRowIterator<'a, R> where R: Read + Seek {}

impl HduTrait for AsciiTableHdu {
    fn header(&self) -> &Header {
        &self.header
    }

    fn info(&self) -> &HduInfo {
        &self.info
    }

    fn hdu_type(&self) -> HduType {
        HduType::AsciiTable
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Header, Keyword};
    use crate::fits::io::reader::HduInfo;

    fn create_test_header(extname: Option<&str>) -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 80));
        header.add_keyword(Keyword::integer("NAXIS2", 100));
        header.add_keyword(Keyword::integer("TFIELDS", 4));

        if let Some(name) = extname {
            header.add_keyword(Keyword::string("EXTNAME", name));
            header.add_keyword(Keyword::integer("EXTVER", 1));
        }

        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 2,
            header_start: 5760,
            header_size: 2880,
            data_start: 8640,
            data_size: 8000,
        }
    }

    #[test]
    fn new_creates_ascii_table_hdu() {
        let header = create_test_header(Some("CATALOG"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.info.index, 2);
        assert_eq!(hdu.number_of_fields(), Some(4));
    }

    #[test]
    fn header_returns_header_reference() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let header_ref = hdu.header();
        assert_eq!(
            header_ref
                .get_keyword_value("XTENSION")
                .unwrap()
                .as_string()
                .unwrap(),
            "TABLE"
        );
    }

    #[test]
    fn info_returns_info_reference() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let info_ref = hdu.info();
        assert_eq!(info_ref.index, 2);
        assert_eq!(info_ref.data_start, 8640);
    }

    #[test]
    fn hdu_type_returns_ascii_table() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::AsciiTable);
    }

    #[test]
    fn number_of_fields_returns_tfields_value() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.number_of_fields(), Some(4));
    }

    #[test]
    fn number_of_fields_returns_none_when_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.number_of_fields(), None);
    }

    #[test]
    fn number_of_rows_returns_naxis2_value() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.number_of_rows(), Some(100));
    }

    #[test]
    fn number_of_rows_returns_none_when_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.number_of_rows(), None);
    }

    #[test]
    fn extension_name_returns_extname_value() {
        let header = create_test_header(Some("PHOTOMETRY"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.extension_name(), Some("PHOTOMETRY"));
    }

    #[test]
    fn extension_name_returns_none_when_missing() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.extension_name(), None);
    }

    #[test]
    fn extension_version_returns_extver_value() {
        let header = create_test_header(Some("TEST"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.extension_version(), Some(1));
    }

    #[test]
    fn extension_version_returns_none_when_missing() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.extension_version(), None);
    }

    #[test]
    fn all_methods_work_together() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 120));
        header.add_keyword(Keyword::integer("NAXIS2", 500));
        header.add_keyword(Keyword::integer("TFIELDS", 6));
        header.add_keyword(Keyword::string("EXTNAME", "STARCAT"));
        header.add_keyword(Keyword::integer("EXTVER", 3));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::AsciiTable);
        assert_eq!(hdu.number_of_fields(), Some(6));
        assert_eq!(hdu.number_of_rows(), Some(500));
        assert_eq!(hdu.extension_name(), Some("STARCAT"));
        assert_eq!(hdu.extension_version(), Some(3));
    }

    #[test]
    fn minimal_valid_ascii_table() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 40));
        header.add_keyword(Keyword::integer("NAXIS2", 0));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::AsciiTable);
        assert_eq!(hdu.number_of_fields(), Some(1));
        assert_eq!(hdu.number_of_rows(), Some(0));
        assert_eq!(hdu.extension_name(), None);
        assert_eq!(hdu.extension_version(), None);
    }

    #[test]
    fn parse_ascii_format_handles_types() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(
            hdu.parse_ascii_format("A10").unwrap(),
            ("A".to_string(), 10)
        );
        assert_eq!(hdu.parse_ascii_format("I6").unwrap(), ("I".to_string(), 6));
        assert_eq!(
            hdu.parse_ascii_format("F12.3").unwrap(),
            ("F".to_string(), 12)
        );
        assert_eq!(
            hdu.parse_ascii_format("E15.7").unwrap(),
            ("E".to_string(), 15)
        );
    }

    #[test]
    fn parse_ascii_format_fails_for_empty_string() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert!(matches!(
            hdu.parse_ascii_format(""),
            Err(FitsError::InvalidFormat(_))
        ));
    }

    #[test]
    fn read_column_raw_returns_ascii_data() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let file_data = vec![0u8; 10000];
        let mut cursor = std::io::Cursor::new(file_data);
        let result = hdu.read_column_raw(&mut cursor, 0).unwrap();

        assert_eq!(result.len(), 100 * 10);
        assert!(result.iter().all(|&b| b == 0));
    }

    #[test]
    fn read_column_raw_returns_empty_for_zero_rows() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        header.add_keyword(Keyword::integer("NAXIS2", 0));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let mut cursor = std::io::Cursor::new(vec![]);
        let result = hdu.read_column_raw(&mut cursor, 0).unwrap();

        assert_eq!(result.len(), 0);
    }

    #[test]
    fn column_info_with_all_metadata() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        header.add_keyword(Keyword::string("TTYPE1", "NAME"));
        header.add_keyword(Keyword::string("TUNIT1", "meter"));
        header.add_keyword(Keyword::string("TNULL1", "NULL"));
        header.add_keyword(Keyword::real("TSCAL1", 2.5));
        header.add_keyword(Keyword::real("TZERO1", 100.0));
        header.add_keyword(Keyword::string("TDISP1", "A8"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let col_info = hdu.column_info(0).unwrap();
        assert_eq!(col_info.index, 0);
        assert_eq!(col_info.name, Some("NAME".to_string()));
        assert_eq!(col_info.format, "A10");
        assert_eq!(col_info.unit, Some("meter".to_string()));
        assert_eq!(col_info.null_value, Some("NULL".to_string()));
        assert_eq!(col_info.scale, Some(2.5));
        assert_eq!(col_info.zero_offset, Some(100.0));
        assert_eq!(col_info.display_format, Some("A8".to_string()));
    }

    #[test]
    fn column_info_minimal_metadata() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let col_info = hdu.column_info(0).unwrap();
        assert_eq!(col_info.index, 0);
        assert_eq!(col_info.name, None);
        assert_eq!(col_info.format, "I6");
        assert_eq!(col_info.unit, None);
        assert_eq!(col_info.null_value, None);
        assert_eq!(col_info.scale, None);
        assert_eq!(col_info.zero_offset, None);
        assert_eq!(col_info.display_format, None);
    }

    #[test]
    fn column_by_name_success() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 3));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        header.add_keyword(Keyword::string("TTYPE1", "ID"));
        header.add_keyword(Keyword::string("TFORM2", "F12.3"));
        header.add_keyword(Keyword::string("TTYPE2", "FLUX"));
        header.add_keyword(Keyword::string("TFORM3", "A20"));
        header.add_keyword(Keyword::string("TTYPE3", "NAME"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.column_by_name("ID").unwrap(), 0);
        assert_eq!(hdu.column_by_name("FLUX").unwrap(), 1);
        assert_eq!(hdu.column_by_name("NAME").unwrap(), 2);
    }

    #[test]
    fn column_by_name_not_found() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        header.add_keyword(Keyword::string("TTYPE1", "ID"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let result = hdu.column_by_name("NONEXISTENT");
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn column_by_name_no_column_names() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 2));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        header.add_keyword(Keyword::string("TFORM2", "F12.3"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let result = hdu.column_by_name("ANYTHING");
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn all_column_info_success() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 2));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        header.add_keyword(Keyword::string("TTYPE1", "ID"));
        header.add_keyword(Keyword::string("TFORM2", "F12.3"));
        header.add_keyword(Keyword::string("TTYPE2", "FLUX"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let columns = hdu.all_column_info().unwrap();
        assert_eq!(columns.len(), 2);
        assert_eq!(columns[0].index, 0);
        assert_eq!(columns[0].name, Some("ID".to_string()));
        assert_eq!(columns[0].format, "I6");
        assert_eq!(columns[1].index, 1);
        assert_eq!(columns[1].name, Some("FLUX".to_string()));
        assert_eq!(columns[1].format, "F12.3");
    }

    #[test]
    fn parse_ascii_format_comprehensive() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.parse_ascii_format("A").unwrap(), ("A".to_string(), 1));
        assert_eq!(hdu.parse_ascii_format("A1").unwrap(), ("A".to_string(), 1));
        assert_eq!(
            hdu.parse_ascii_format("A10").unwrap(),
            ("A".to_string(), 10)
        );
        assert_eq!(hdu.parse_ascii_format("I6").unwrap(), ("I".to_string(), 6));
        assert_eq!(
            hdu.parse_ascii_format("F12.3").unwrap(),
            ("F".to_string(), 12)
        );
        assert_eq!(
            hdu.parse_ascii_format("E15.7").unwrap(),
            ("E".to_string(), 15)
        );
        assert_eq!(
            hdu.parse_ascii_format("D25.17").unwrap(),
            ("D".to_string(), 25)
        );
    }

    #[test]
    fn calculate_column_offset_first_column() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 3));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        header.add_keyword(Keyword::string("TFORM2", "I6"));
        header.add_keyword(Keyword::string("TFORM3", "F12.3"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.calculate_column_offset(0).unwrap(), 0);
    }

    #[test]
    fn calculate_column_offset_subsequent_columns() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 3));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        header.add_keyword(Keyword::string("TFORM2", "I6"));
        header.add_keyword(Keyword::string("TFORM3", "F12.3"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.calculate_column_offset(1).unwrap(), 11);
        assert_eq!(hdu.calculate_column_offset(2).unwrap(), 18);
    }

    #[test]
    fn calculate_column_offset_out_of_range() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 2));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        header.add_keyword(Keyword::string("TFORM2", "I6"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let result = hdu.calculate_column_offset(5);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn get_row_size_calculation() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 3));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        header.add_keyword(Keyword::string("TFORM2", "I6"));
        header.add_keyword(Keyword::string("TFORM3", "F12.3"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.get_row_size().unwrap(), 31);
    }

    #[test]
    fn read_column_with_nulls_functionality() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::integer("NAXIS2", 0));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        header.add_keyword(Keyword::string("TNULL1", "-999"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let mut cursor = std::io::Cursor::new(vec![]);
        let result = hdu
            .read_column_with_nulls::<i32, _>(&mut cursor, 0)
            .unwrap();
        assert_eq!(result.len(), 0);

        let col_info = hdu.column_info(0).unwrap();
        assert_eq!(col_info.null_value, Some("-999".to_string()));
    }

    #[test]
    fn get_row_fails_for_invalid_index() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);
        let mut cursor = std::io::Cursor::new(vec![0u8; 100]);

        let result = hdu.get_row(&mut cursor, 200);
        assert!(result.is_err());
    }

    #[test]
    fn get_row_parses_integer() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 7));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 7,
        };
        let hdu = AsciiTableHdu::new(header, info);

        let data = b"   42 \n";
        let mut cursor = std::io::Cursor::new(data.to_vec());

        let result = hdu.get_row(&mut cursor, 0);
        assert!(result.is_ok());
        let row = result.unwrap();
        assert_eq!(row.len(), 1);
        assert_eq!(row[0], TableValue::I64(42));
    }

    #[test]
    fn get_row_parses_float() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 10));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "F9.3"));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 10,
        };
        let hdu = AsciiTableHdu::new(header, info);

        let data = b"   1.2345\n";
        let mut cursor = std::io::Cursor::new(data.to_vec());

        let result = hdu.get_row(&mut cursor, 0);
        assert!(result.is_ok());
        let row = result.unwrap();
        assert_eq!(row.len(), 1);
        if let TableValue::F64(v) = row[0] {
            assert!((v - 1.2345).abs() < 1e-6);
        } else {
            panic!("Expected F64");
        }
    }

    #[test]
    fn get_row_parses_string() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 11));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "A10"));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 11,
        };
        let hdu = AsciiTableHdu::new(header, info);

        let data = b"Hello     \n";
        let mut cursor = std::io::Cursor::new(data.to_vec());

        let result = hdu.get_row(&mut cursor, 0);
        assert!(result.is_ok());
        let row = result.unwrap();
        assert_eq!(row.len(), 1);
        assert_eq!(row[0], TableValue::String("Hello".to_string()));
    }

    #[test]
    fn get_column_by_name_success() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 7));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        header.add_keyword(Keyword::string("TTYPE1", "ID"));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 14,
        };
        let hdu = AsciiTableHdu::new(header, info);

        let data = b"   10 \n   20 \n";
        let mut cursor = std::io::Cursor::new(data.to_vec());

        let result = hdu.get_column_by_name(&mut cursor, "ID");
        assert!(result.is_ok());
        let column = result.unwrap();
        assert_eq!(column.len(), 2);
        assert_eq!(column[0], TableValue::I64(10));
        assert_eq!(column[1], TableValue::I64(20));
    }

    #[test]
    fn get_column_by_name_not_found() {
        let mut header = create_test_header(None);
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);
        let mut cursor = std::io::Cursor::new(vec![0u8; 100]);

        let result = hdu.get_column_by_name(&mut cursor, "NONEXISTENT");
        assert!(result.is_err());
    }

    #[test]
    fn row_count_returns_correct_value() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.row_count(), 100);
    }

    #[test]
    fn row_count_returns_zero_when_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        assert_eq!(hdu.row_count(), 0);
    }

    #[test]
    fn iter_rows_returns_correct_count() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 7));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 14,
        };
        let hdu = AsciiTableHdu::new(header, info);

        let data = b"   10 \n   20 \n";
        let mut cursor = std::io::Cursor::new(data.to_vec());

        let iter = hdu.iter_rows(&mut cursor);
        assert_eq!(iter.len(), 2);
    }

    #[test]
    fn parse_ascii_value_empty_is_null() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let result = hdu.parse_ascii_value("I", b"      ").unwrap();
        assert_eq!(result, TableValue::Null);
    }

    #[test]
    fn parse_ascii_float_with_d_exponent() {
        let header = create_test_header(None);
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);

        let result = hdu.parse_ascii_float("1.5D+02");
        assert!(result.is_ok());
        if let TableValue::F64(v) = result.unwrap() {
            assert!((v - 150.0).abs() < 1e-6);
        } else {
            panic!("Expected F64");
        }
    }

    #[test]
    fn get_column_values_empty_table() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "TABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 7));
        header.add_keyword(Keyword::integer("NAXIS2", 0));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "I6"));
        let info = create_test_hdu_info();
        let hdu = AsciiTableHdu::new(header, info);
        let mut cursor = std::io::Cursor::new(vec![]);

        let result = hdu.get_column_values(&mut cursor, 0);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }
}
