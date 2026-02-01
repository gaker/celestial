use super::BinaryTableHdu;
use crate::core::ByteOrder;
use crate::fits::data::array::{DataArray, DataValue, TableValue};
use crate::fits::hdu::ColumnInfo;
use crate::fits::{FitsError, Result};
use std::io::{Read, Seek, SeekFrom};
use std::marker::PhantomData;

#[derive(Debug)]
pub(super) struct ColumnReadParams {
    pub width: usize,
    pub bytes_per_element: usize,
    pub column_offset: usize,
    pub row_size: usize,
}

impl BinaryTableHdu {
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
        let format = self.get_column_format(column_index)?;
        let mut info = ColumnInfo::new(column, format);

        self.apply_column_metadata(&mut info, column_index);
        Ok(info)
    }

    fn get_column_format(&self, column_index: usize) -> Result<String> {
        let format_key = format!("TFORM{}", column_index);
        self.header
            .get_keyword_value(&format_key)
            .and_then(|v| v.as_string())
            .map(|s| s.to_string())
            .ok_or(FitsError::KeywordNotFound {
                keyword: format_key,
            })
    }

    fn apply_column_metadata(&self, info: &mut ColumnInfo, column_index: usize) {
        self.set_column_name(info, column_index);
        self.set_column_unit(info, column_index);
        self.set_column_null_value(info, column_index);
        self.set_column_scale(info, column_index);
        self.set_column_zero_offset(info, column_index);
        self.set_column_display_format(info, column_index);
    }

    fn set_column_name(&self, info: &mut ColumnInfo, column_index: usize) {
        if let Some(name) = self
            .header
            .get_keyword_value(&format!("TTYPE{}", column_index))
            .and_then(|v| v.as_string())
        {
            *info = info.clone().with_name(name.to_string());
        }
    }

    fn set_column_unit(&self, info: &mut ColumnInfo, column_index: usize) {
        if let Some(unit) = self
            .header
            .get_keyword_value(&format!("TUNIT{}", column_index))
            .and_then(|v| v.as_string())
        {
            *info = info.clone().with_unit(unit.to_string());
        }
    }

    fn set_column_null_value(&self, info: &mut ColumnInfo, column_index: usize) {
        if let Some(null_val) = self
            .header
            .get_keyword_value(&format!("TNULL{}", column_index))
            .and_then(|v| v.as_string())
        {
            *info = info.clone().with_null_value(null_val.to_string());
        }
    }

    fn set_column_scale(&self, info: &mut ColumnInfo, column_index: usize) {
        if let Some(scale) = self
            .header
            .get_keyword_value(&format!("TSCAL{}", column_index))
            .and_then(|v| v.as_real())
        {
            *info = info.clone().with_scale(scale);
        }
    }

    fn set_column_zero_offset(&self, info: &mut ColumnInfo, column_index: usize) {
        if let Some(zero) = self
            .header
            .get_keyword_value(&format!("TZERO{}", column_index))
            .and_then(|v| v.as_real())
        {
            *info = info.clone().with_zero_offset(zero);
        }
    }

    fn set_column_display_format(&self, info: &mut ColumnInfo, column_index: usize) {
        if let Some(disp) = self
            .header
            .get_keyword_value(&format!("TDISP{}", column_index))
            .and_then(|v| v.as_string())
        {
            info.display_format = Some(disp.to_string());
        }
    }

    fn build_column_index(&self) -> Result<std::collections::HashMap<String, usize>> {
        let column_count = self.column_count()?;
        let mut index = std::collections::HashMap::new();

        for i in 0..column_count {
            if let Ok(info) = self.column_info(i) {
                if let Some(col_name) = &info.name {
                    index.insert(col_name.clone(), i);
                }
            }
        }

        Ok(index)
    }

    pub fn column_by_name(&self, name: &str) -> Result<usize> {
        let index = self
            .column_name_index
            .get_or_init(|| self.build_column_index().unwrap_or_default());

        index
            .get(name)
            .copied()
            .ok_or_else(|| FitsError::InvalidFormat(format!("Column '{}' not found", name)))
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
        let row_count = self.number_of_rows().unwrap_or(0) as usize;
        if row_count == 0 {
            return Ok(Vec::new());
        }

        let read_params = self.prepare_column_read(column)?;
        self.read_column_data(reader, &read_params, row_count)
    }

    pub(super) fn prepare_column_read(&self, column: usize) -> Result<ColumnReadParams> {
        let info = self.column_info(column)?;
        let (data_type, width) = self.parse_binary_format(&info.format)?;
        let bytes_per_element = self.get_element_size(&data_type)?;
        let column_offset = self.calculate_column_offset(column)?;
        let row_size = self.get_row_size()?;

        Ok(ColumnReadParams {
            width,
            bytes_per_element,
            column_offset,
            row_size,
        })
    }

    fn read_column_data<R>(
        &self,
        reader: &mut R,
        params: &ColumnReadParams,
        row_count: usize,
    ) -> Result<Vec<u8>>
    where
        R: Read + Seek,
    {
        let row_size = params.row_size;
        let column_bytes_per_row = params.width * params.bytes_per_element;
        let total_column_bytes = row_count * column_bytes_per_row;
        let mut result = Vec::with_capacity(total_column_bytes);

        let data_start = self.info.data_start;
        reader.seek(SeekFrom::Start(data_start))?;

        let mut row_buffer = vec![0u8; row_size];
        for _ in 0..row_count {
            reader.read_exact(&mut row_buffer)?;

            let column_start = params.column_offset;
            let column_end = column_start + column_bytes_per_row;
            result.extend_from_slice(&row_buffer[column_start..column_end]);
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
            let (data_type, width) = self.parse_binary_format(&info.format)?;
            let element_size = self.get_element_size(&data_type)?;
            offset += width * element_size;
        }

        Ok(offset)
    }

    pub(super) fn get_row_size(&self) -> Result<usize> {
        let column_count = self.column_count()?;
        let mut row_size = 0;

        for i in 0..column_count {
            let info = self.column_info(i)?;
            let (data_type, width) = self.parse_binary_format(&info.format)?;
            let element_size = self.get_element_size(&data_type)?;
            row_size += width * element_size;
        }

        let padding = (8 - (row_size % 8)) % 8;
        Ok(row_size + padding)
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

    pub fn read_column_i16<R>(&self, reader: &mut R, column: usize) -> Result<Vec<DataValue<i16>>>
    where
        R: Read + Seek,
    {
        self.read_column_with_nulls(reader, column)
    }

    pub fn read_column_i32<R>(&self, reader: &mut R, column: usize) -> Result<Vec<DataValue<i32>>>
    where
        R: Read + Seek,
    {
        self.read_column_with_nulls(reader, column)
    }

    pub fn read_column_f32<R>(&self, reader: &mut R, column: usize) -> Result<Vec<DataValue<f32>>>
    where
        R: Read + Seek,
    {
        self.read_column_with_nulls(reader, column)
    }

    pub fn read_column_f64<R>(&self, reader: &mut R, column: usize) -> Result<Vec<DataValue<f64>>>
    where
        R: Read + Seek,
    {
        self.read_column_with_nulls(reader, column)
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
            let (data_type, repeat) = self.parse_binary_format(&info.format)?;
            let elem_size = self.get_element_size(&data_type)?;
            let col_bytes = repeat * elem_size;

            let col_data = &row_bytes[offset..offset + col_bytes];
            let value = self.parse_column_value(&data_type, col_data, repeat)?;
            values.push(value);
            offset += col_bytes;
        }

        Ok(values)
    }

    fn parse_column_value(
        &self,
        data_type: &str,
        bytes: &[u8],
        repeat: usize,
    ) -> Result<TableValue> {
        let type_char = data_type.chars().next().unwrap_or('X');

        match type_char {
            'L' => self.parse_logical(bytes),
            'X' => Ok(TableValue::Byte(bytes.first().copied().unwrap_or(0))),
            'B' => self.parse_byte(bytes, repeat),
            'I' => self.parse_i16(bytes, repeat),
            'J' => self.parse_i32(bytes, repeat),
            'K' => self.parse_i64(bytes, repeat),
            'A' => self.parse_string(bytes),
            'E' => self.parse_f32(bytes, repeat),
            'D' => self.parse_f64(bytes, repeat),
            'C' => self.parse_complex32(bytes),
            'M' => self.parse_complex64(bytes),
            _ => Ok(TableValue::Null),
        }
    }

    fn parse_logical(&self, bytes: &[u8]) -> Result<TableValue> {
        match bytes.first() {
            Some(b'T') | Some(1) => Ok(TableValue::Logical(true)),
            Some(b'F') | Some(0) => Ok(TableValue::Logical(false)),
            _ => Ok(TableValue::Null),
        }
    }

    fn parse_byte(&self, bytes: &[u8], repeat: usize) -> Result<TableValue> {
        if repeat == 1 {
            Ok(TableValue::Byte(bytes.first().copied().unwrap_or(0)))
        } else {
            Ok(TableValue::String(
                String::from_utf8_lossy(bytes).into_owned(),
            ))
        }
    }

    fn parse_i16(&self, bytes: &[u8], repeat: usize) -> Result<TableValue> {
        if repeat == 1 && bytes.len() >= 2 {
            Ok(TableValue::I16(i16::from_be_bytes([bytes[0], bytes[1]])))
        } else {
            self.parse_first_i16(bytes)
        }
    }

    fn parse_first_i16(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 2 {
            Ok(TableValue::I16(i16::from_be_bytes([bytes[0], bytes[1]])))
        } else {
            Ok(TableValue::Null)
        }
    }

    fn parse_i32(&self, bytes: &[u8], repeat: usize) -> Result<TableValue> {
        if repeat == 1 && bytes.len() >= 4 {
            Ok(TableValue::I32(i32::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3],
            ])))
        } else {
            self.parse_first_i32(bytes)
        }
    }

    fn parse_first_i32(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 4 {
            Ok(TableValue::I32(i32::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3],
            ])))
        } else {
            Ok(TableValue::Null)
        }
    }

    fn parse_i64(&self, bytes: &[u8], repeat: usize) -> Result<TableValue> {
        if repeat == 1 && bytes.len() >= 8 {
            let val = i64::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
            ]);
            Ok(TableValue::I64(val))
        } else {
            self.parse_first_i64(bytes)
        }
    }

    fn parse_first_i64(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 8 {
            let val = i64::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
            ]);
            Ok(TableValue::I64(val))
        } else {
            Ok(TableValue::Null)
        }
    }

    fn parse_string(&self, bytes: &[u8]) -> Result<TableValue> {
        let s = String::from_utf8_lossy(bytes).trim_end().to_string();
        Ok(TableValue::String(s))
    }

    fn parse_f32(&self, bytes: &[u8], repeat: usize) -> Result<TableValue> {
        if repeat == 1 && bytes.len() >= 4 {
            let val = f32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
            Ok(TableValue::F32(val))
        } else {
            self.parse_first_f32(bytes)
        }
    }

    fn parse_first_f32(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 4 {
            let val = f32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
            Ok(TableValue::F32(val))
        } else {
            Ok(TableValue::Null)
        }
    }

    fn parse_f64(&self, bytes: &[u8], repeat: usize) -> Result<TableValue> {
        if repeat == 1 && bytes.len() >= 8 {
            let val = f64::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
            ]);
            Ok(TableValue::F64(val))
        } else {
            self.parse_first_f64(bytes)
        }
    }

    fn parse_first_f64(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 8 {
            let val = f64::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
            ]);
            Ok(TableValue::F64(val))
        } else {
            Ok(TableValue::Null)
        }
    }

    fn parse_complex32(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 8 {
            let real = f32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
            let imag = f32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
            Ok(TableValue::Complex32(real, imag))
        } else {
            Ok(TableValue::Null)
        }
    }

    fn parse_complex64(&self, bytes: &[u8]) -> Result<TableValue> {
        if bytes.len() >= 16 {
            let real = f64::from_be_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
            ]);
            let imag = f64::from_be_bytes([
                bytes[8], bytes[9], bytes[10], bytes[11], bytes[12], bytes[13], bytes[14],
                bytes[15],
            ]);
            Ok(TableValue::Complex64(real, imag))
        } else {
            Ok(TableValue::Null)
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
        let (data_type, repeat) = self.parse_binary_format(&info.format)?;
        let raw_data = self.read_column_raw(reader, column)?;

        self.convert_column_to_values(&data_type, &raw_data, repeat, row_count)
    }

    fn convert_column_to_values(
        &self,
        data_type: &str,
        raw_data: &[u8],
        repeat: usize,
        row_count: usize,
    ) -> Result<Vec<TableValue>> {
        let elem_size = self.get_element_size(data_type)?;
        let bytes_per_row = repeat * elem_size;
        let mut values = Vec::with_capacity(row_count);

        for row in 0..row_count {
            let start = row * bytes_per_row;
            let end = start + bytes_per_row;
            let row_bytes = &raw_data[start..end];
            let value = self.parse_column_value(data_type, row_bytes, repeat)?;
            values.push(value);
        }

        Ok(values)
    }

    pub fn row_count(&self) -> usize {
        self.number_of_rows().unwrap_or(0) as usize
    }

    pub fn iter_rows<'a, R>(&'a self, reader: &'a mut R) -> BinaryTableRowIterator<'a, R>
    where
        R: Read + Seek,
    {
        BinaryTableRowIterator::new(self, reader)
    }
}

pub struct BinaryTableRowIterator<'a, R>
where
    R: Read + Seek,
{
    hdu: &'a BinaryTableHdu,
    reader: &'a mut R,
    current_row: usize,
    row_count: usize,
    _phantom: PhantomData<R>,
}

impl<'a, R> BinaryTableRowIterator<'a, R>
where
    R: Read + Seek,
{
    fn new(hdu: &'a BinaryTableHdu, reader: &'a mut R) -> Self {
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

impl<'a, R> Iterator for BinaryTableRowIterator<'a, R>
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

impl<'a, R> ExactSizeIterator for BinaryTableRowIterator<'a, R> where R: Read + Seek {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Header, Keyword};
    use crate::fits::io::reader::HduInfo;
    use celestial_core::constants::PI;
    use std::io::Cursor;

    fn create_test_header() -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", 20));
        header.add_keyword(Keyword::integer("NAXIS2", 3));
        header.add_keyword(Keyword::integer("TFIELDS", 3));

        header.add_keyword(Keyword::string("TTYPE1", "COL1"));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::string("TUNIT1", "meters"));
        header.add_keyword(Keyword::string("TNULL1", "-999"));
        header.add_keyword(Keyword::real("TSCAL1", 2.0));
        header.add_keyword(Keyword::real("TZERO1", 100.0));
        header.add_keyword(Keyword::string("TDISP1", "I8"));

        header.add_keyword(Keyword::string("TTYPE2", "COL2"));
        header.add_keyword(Keyword::string("TFORM2", "2I"));

        header.add_keyword(Keyword::string("TTYPE3", "COL3"));
        header.add_keyword(Keyword::string("TFORM3", "4E"));

        header
    }

    fn create_minimal_header() -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 1,
            header_start: 2880,
            header_size: 2880,
            data_start: 5760,
            data_size: 60,
        }
    }

    #[test]
    fn column_count_returns_tfields_value() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_count();
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 3);
    }

    #[test]
    fn column_count_fails_when_tfields_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_count();
        assert!(result.is_err());
    }

    #[test]
    fn column_info_returns_complete_metadata() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_info(0);
        assert!(result.is_ok());

        let col_info = result.unwrap();
        assert_eq!(col_info.index, 0);
        assert_eq!(col_info.format, "1J");
        assert_eq!(col_info.name, Some("COL1".to_string()));
        assert_eq!(col_info.unit, Some("meters".to_string()));
        assert_eq!(col_info.null_value, Some("-999".to_string()));
        assert_eq!(col_info.scale, Some(2.0));
        assert_eq!(col_info.zero_offset, Some(100.0));
        assert_eq!(col_info.display_format, Some("I8".to_string()));
    }

    #[test]
    fn column_info_returns_minimal_metadata() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_info(1);
        assert!(result.is_ok());

        let col_info = result.unwrap();
        assert_eq!(col_info.index, 1);
        assert_eq!(col_info.format, "2I");
        assert_eq!(col_info.name, Some("COL2".to_string()));
        assert!(col_info.unit.is_none());
        assert!(col_info.null_value.is_none());
        assert!(col_info.scale.is_none());
        assert!(col_info.zero_offset.is_none());
        assert!(col_info.display_format.is_none());
    }

    #[test]
    fn column_info_fails_for_invalid_index() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_info(10);
        assert!(result.is_err());
    }

    #[test]
    fn column_info_fails_when_format_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_info(0);
        assert!(result.is_err());
    }

    #[test]
    fn column_by_name_finds_existing_column() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_by_name("COL2");
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 1);
    }

    #[test]
    fn column_by_name_fails_for_nonexistent_column() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_by_name("NONEXISTENT");
        assert!(result.is_err());
    }

    #[test]
    fn column_by_name_uses_caching() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result1 = hdu.column_by_name("COL1");
        assert!(result1.is_ok());
        assert_eq!(result1.unwrap(), 0);

        let result2 = hdu.column_by_name("COL1");
        assert!(result2.is_ok());
        assert_eq!(result2.unwrap(), 0);
    }

    #[test]
    fn all_column_info_returns_all_columns() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.all_column_info();
        assert!(result.is_ok());

        let columns = result.unwrap();
        assert_eq!(columns.len(), 3);
        assert_eq!(columns[0].name, Some("COL1".to_string()));
        assert_eq!(columns[1].name, Some("COL2".to_string()));
        assert_eq!(columns[2].name, Some("COL3".to_string()));
    }

    #[test]
    fn read_column_raw_with_zero_rows() {
        let mut header = create_minimal_header();
        header.add_keyword(Keyword::integer("NAXIS2", 0));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.read_column_raw(&mut cursor, 0);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }

    #[test]
    fn read_column_raw_fails_for_invalid_column() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.read_column_raw(&mut cursor, 10);
        assert!(result.is_err());
    }

    #[test]
    fn prepare_column_read_calculates_parameters() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.prepare_column_read(0);
        assert!(result.is_ok());

        let params = result.unwrap();
        assert_eq!(params.width, 1);
        assert_eq!(params.bytes_per_element, 4);
        assert_eq!(params.column_offset, 0);
        assert_eq!(params.row_size, 24);
    }

    #[test]
    fn read_column_with_nulls_i16() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.read_column_i16(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn read_column_with_nulls_i32() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.read_column_i32(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn read_column_with_nulls_f32() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.read_column_f32(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn read_column_with_nulls_f64() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.read_column_f64(&mut cursor, 0);
        assert!(result.is_err());
    }

    #[test]
    fn get_row_size_with_padding() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.get_row_size();
        assert!(result.is_ok());

        assert_eq!(result.unwrap(), 24);
    }

    #[test]
    fn column_metadata_methods_coverage() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_info(2);
        assert!(result.is_ok());

        let col_info = result.unwrap();
        assert_eq!(col_info.name, Some("COL3".to_string()));
        assert!(col_info.unit.is_none());
        assert!(col_info.null_value.is_none());
        assert!(col_info.scale.is_none());
        assert!(col_info.zero_offset.is_none());
        assert!(col_info.display_format.is_none());
    }

    #[test]
    fn build_column_index_handles_missing_names() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 2));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::string("TFORM2", "1I"));
        header.add_keyword(Keyword::string("TTYPE1", "NAMED_COL"));

        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.column_by_name("NAMED_COL");
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 0);

        let result = hdu.column_by_name("UNNAMED");
        assert!(result.is_err());
    }

    #[test]
    fn get_row_fails_for_invalid_index() {
        let mut header = create_minimal_header();
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.get_row(&mut cursor, 10);
        assert!(result.is_err());
    }

    #[test]
    fn get_row_reads_data_correctly() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::integer("NAXIS2", 1));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 8,
        };
        let hdu = BinaryTableHdu::new(header, info);

        let data = vec![0x00, 0x00, 0x00, 0x2A, 0x00, 0x00, 0x00, 0x00];
        let mut cursor = Cursor::new(data);

        let result = hdu.get_row(&mut cursor, 0);
        assert!(result.is_ok());
        let row = result.unwrap();
        assert_eq!(row.len(), 1);
        assert_eq!(row[0], TableValue::I32(42));
    }

    #[test]
    fn get_column_by_name_success() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::string("TTYPE1", "VALUES"));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 16,
        };
        let hdu = BinaryTableHdu::new(header, info);

        let data = vec![
            0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00,
            0x00, 0x00,
        ];
        let mut cursor = Cursor::new(data);

        let result = hdu.get_column_by_name(&mut cursor, "VALUES");
        assert!(result.is_ok());
        let column = result.unwrap();
        assert_eq!(column.len(), 2);
        assert_eq!(column[0], TableValue::I32(1));
        assert_eq!(column[1], TableValue::I32(2));
    }

    #[test]
    fn get_column_by_name_not_found() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![0u8; 100]);

        let result = hdu.get_column_by_name(&mut cursor, "NONEXISTENT");
        assert!(result.is_err());
    }

    #[test]
    fn row_count_returns_correct_value() {
        let header = create_test_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.row_count(), 3);
    }

    #[test]
    fn row_count_returns_zero_when_missing() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(hdu.row_count(), 0);
    }

    #[test]
    fn iter_rows_returns_correct_count() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 16,
        };
        let hdu = BinaryTableHdu::new(header, info);

        let data = vec![
            0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00,
            0x00, 0x00,
        ];
        let mut cursor = Cursor::new(data);

        let iter = hdu.iter_rows(&mut cursor);
        assert_eq!(iter.len(), 2);
    }

    #[test]
    fn parse_column_value_logical() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        assert_eq!(
            hdu.parse_column_value("L", b"T", 1).unwrap(),
            TableValue::Logical(true)
        );
        assert_eq!(
            hdu.parse_column_value("L", b"F", 1).unwrap(),
            TableValue::Logical(false)
        );
        assert_eq!(
            hdu.parse_column_value("L", &[1], 1).unwrap(),
            TableValue::Logical(true)
        );
        assert_eq!(
            hdu.parse_column_value("L", &[0], 1).unwrap(),
            TableValue::Logical(false)
        );
    }

    #[test]
    fn parse_column_value_string() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let result = hdu.parse_column_value("A", b"Hello   ", 8).unwrap();
        assert_eq!(result, TableValue::String("Hello".to_string()));
    }

    #[test]
    fn parse_column_value_f64() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let pi_bytes = PI.to_be_bytes();
        let result = hdu.parse_column_value("D", &pi_bytes, 1).unwrap();
        if let TableValue::F64(v) = result {
            assert!((v - PI).abs() < f64::EPSILON);
        } else {
            panic!("Expected F64");
        }
    }

    #[test]
    fn parse_column_value_complex32() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&1.0f32.to_be_bytes());
        bytes.extend_from_slice(&2.0f32.to_be_bytes());

        let result = hdu.parse_column_value("C", &bytes, 1).unwrap();
        assert_eq!(result, TableValue::Complex32(1.0, 2.0));
    }

    #[test]
    fn parse_column_value_complex64() {
        let header = create_minimal_header();
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&3.0f64.to_be_bytes());
        bytes.extend_from_slice(&4.0f64.to_be_bytes());

        let result = hdu.parse_column_value("M", &bytes, 1).unwrap();
        assert_eq!(result, TableValue::Complex64(3.0, 4.0));
    }

    #[test]
    fn get_column_values_empty_table() {
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TFORM1", "1J"));
        header.add_keyword(Keyword::integer("NAXIS2", 0));
        let info = create_test_hdu_info();
        let hdu = BinaryTableHdu::new(header, info);
        let mut cursor = Cursor::new(vec![]);

        let result = hdu.get_column_values(&mut cursor, 0);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }
}
