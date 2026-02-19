use crate::core::{BitPix, ByteOrder};
use crate::fits::{FitsError, Result};
use byteorder::{BigEndian, LittleEndian, ReadBytesExt, WriteBytesExt};
use std::io::Cursor;

#[derive(Debug, Clone, PartialEq)]
pub enum DataValue<T> {
    Value(T),
    Null,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TableValue {
    Logical(bool),
    Byte(u8),
    I16(i16),
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
    String(String),
    Complex32(f32, f32),
    Complex64(f64, f64),
    Null,
}

impl TableValue {
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            TableValue::Byte(v) => Some(*v as i64),
            TableValue::I16(v) => Some(*v as i64),
            TableValue::I32(v) => Some(*v as i64),
            TableValue::I64(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            TableValue::Byte(v) => Some(*v as f64),
            TableValue::I16(v) => Some(*v as f64),
            TableValue::I32(v) => Some(*v as f64),
            TableValue::I64(v) => Some(*v as f64),
            TableValue::F32(v) => Some(*v as f64),
            TableValue::F64(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_string(&self) -> Option<&str> {
        match self {
            TableValue::String(s) => Some(s),
            _ => None,
        }
    }

    pub fn is_null(&self) -> bool {
        matches!(self, TableValue::Null)
    }
}

impl<T> DataValue<T> {
    pub fn is_null(&self) -> bool {
        matches!(self, DataValue::Null)
    }

    pub fn value(&self) -> Option<&T> {
        match self {
            DataValue::Value(v) => Some(v),
            DataValue::Null => None,
        }
    }

    pub fn unwrap_or(self, default: T) -> T {
        match self {
            DataValue::Value(v) => v,
            DataValue::Null => default,
        }
    }
}

pub trait DataArray: Sized + PartialEq {
    const BITPIX: BitPix;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>>;
    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>>;

    fn from_bytes_with_null(
        bytes: &[u8],
        byte_order: ByteOrder,
        null_value: Option<Self>,
    ) -> Result<Vec<DataValue<Self>>> {
        let raw_data = Self::from_bytes(bytes, byte_order)?;
        Ok(raw_data
            .into_iter()
            .map(|value| match &null_value {
                Some(null_val) if value == *null_val => DataValue::Null,
                _ => DataValue::Value(value),
            })
            .collect())
    }

    fn parse_null_value(null_str: &str) -> Result<Self>;

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64);
}

impl DataArray for u8 {
    const BITPIX: BitPix = BitPix::U8;

    fn from_bytes(bytes: &[u8], _byte_order: ByteOrder) -> Result<Vec<Self>> {
        Ok(bytes.to_vec())
    }

    fn to_bytes(data: &[Self], _byte_order: ByteOrder) -> Result<Vec<u8>> {
        Ok(data.to_vec())
    }

    fn parse_null_value(null_str: &str) -> Result<Self> {
        null_str.parse::<u8>().map_err(|_| {
            FitsError::InvalidFormat(format!("Invalid NULL value for u8: {}", null_str))
        })
    }

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64) {
        for val in data.iter_mut() {
            let scaled = (*val as f64 * bscale) + bzero;
            *val = libm::round(scaled).clamp(0.0, 255.0) as u8;
        }
    }
}

impl DataArray for i16 {
    const BITPIX: BitPix = BitPix::I16;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>> {
        let mut cursor = Cursor::new(bytes);
        let mut result = Vec::with_capacity(bytes.len() / 2);

        while cursor.position() < bytes.len() as u64 {
            let value = match byte_order {
                ByteOrder::BigEndian => cursor.read_i16::<BigEndian>(),
                ByteOrder::LittleEndian => cursor.read_i16::<LittleEndian>(),
            };

            match value {
                Ok(v) => result.push(v),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(FitsError::Io(e)),
            }
        }

        Ok(result)
    }

    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>> {
        let mut result = Vec::with_capacity(data.len() * 2);

        for &value in data {
            match byte_order {
                ByteOrder::BigEndian => result.write_i16::<BigEndian>(value)?,
                ByteOrder::LittleEndian => result.write_i16::<LittleEndian>(value)?,
            }
        }

        Ok(result)
    }

    fn parse_null_value(null_str: &str) -> Result<Self> {
        null_str.parse::<i16>().map_err(|_| {
            FitsError::InvalidFormat(format!("Invalid NULL value for i16: {}", null_str))
        })
    }

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64) {
        for val in data.iter_mut() {
            let scaled = (*val as f64 * bscale) + bzero;
            *val = libm::round(scaled).clamp(i16::MIN as f64, i16::MAX as f64) as i16;
        }
    }
}

impl DataArray for i32 {
    const BITPIX: BitPix = BitPix::I32;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>> {
        let mut cursor = Cursor::new(bytes);
        let mut result = Vec::with_capacity(bytes.len() / 4);

        while cursor.position() < bytes.len() as u64 {
            let value = match byte_order {
                ByteOrder::BigEndian => cursor.read_i32::<BigEndian>(),
                ByteOrder::LittleEndian => cursor.read_i32::<LittleEndian>(),
            };

            match value {
                Ok(v) => result.push(v),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(FitsError::Io(e)),
            }
        }

        Ok(result)
    }

    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>> {
        let mut result = Vec::with_capacity(data.len() * 4);

        for &value in data {
            match byte_order {
                ByteOrder::BigEndian => result.write_i32::<BigEndian>(value)?,
                ByteOrder::LittleEndian => result.write_i32::<LittleEndian>(value)?,
            }
        }

        Ok(result)
    }

    fn parse_null_value(null_str: &str) -> Result<Self> {
        null_str.parse::<i32>().map_err(|_| {
            FitsError::InvalidFormat(format!("Invalid NULL value for i32: {}", null_str))
        })
    }

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64) {
        for val in data.iter_mut() {
            let scaled = (*val as f64 * bscale) + bzero;
            *val = libm::round(scaled).clamp(i32::MIN as f64, i32::MAX as f64) as i32;
        }
    }
}

impl DataArray for i64 {
    const BITPIX: BitPix = BitPix::I64;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>> {
        let mut cursor = Cursor::new(bytes);
        let mut result = Vec::with_capacity(bytes.len() / 8);

        while cursor.position() < bytes.len() as u64 {
            let value = match byte_order {
                ByteOrder::BigEndian => cursor.read_i64::<BigEndian>(),
                ByteOrder::LittleEndian => cursor.read_i64::<LittleEndian>(),
            };

            match value {
                Ok(v) => result.push(v),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(FitsError::Io(e)),
            }
        }

        Ok(result)
    }

    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>> {
        let mut result = Vec::with_capacity(data.len() * 8);

        for &value in data {
            match byte_order {
                ByteOrder::BigEndian => result.write_i64::<BigEndian>(value)?,
                ByteOrder::LittleEndian => result.write_i64::<LittleEndian>(value)?,
            }
        }

        Ok(result)
    }

    fn parse_null_value(null_str: &str) -> Result<Self> {
        null_str.parse::<i64>().map_err(|_| {
            FitsError::InvalidFormat(format!("Invalid NULL value for i64: {}", null_str))
        })
    }

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64) {
        for val in data.iter_mut() {
            let scaled = (*val as f64 * bscale) + bzero;
            *val = libm::round(scaled).clamp(i64::MIN as f64, i64::MAX as f64) as i64;
        }
    }
}

impl DataArray for f32 {
    const BITPIX: BitPix = BitPix::F32;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>> {
        let mut cursor = Cursor::new(bytes);
        let mut result = Vec::with_capacity(bytes.len() / 4);

        while cursor.position() < bytes.len() as u64 {
            let value = match byte_order {
                ByteOrder::BigEndian => cursor.read_f32::<BigEndian>(),
                ByteOrder::LittleEndian => cursor.read_f32::<LittleEndian>(),
            };

            match value {
                Ok(v) => result.push(v),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(FitsError::Io(e)),
            }
        }

        Ok(result)
    }

    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>> {
        let mut result = Vec::with_capacity(data.len() * 4);

        for &value in data {
            match byte_order {
                ByteOrder::BigEndian => result.write_f32::<BigEndian>(value)?,
                ByteOrder::LittleEndian => result.write_f32::<LittleEndian>(value)?,
            }
        }

        Ok(result)
    }

    fn parse_null_value(null_str: &str) -> Result<Self> {
        null_str.parse::<f32>().map_err(|_| {
            FitsError::InvalidFormat(format!("Invalid NULL value for f32: {}", null_str))
        })
    }

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64) {
        let bscale_f32 = bscale as f32;
        let bzero_f32 = bzero as f32;
        for val in data.iter_mut() {
            *val = (*val * bscale_f32) + bzero_f32;
        }
    }
}

impl DataArray for f64 {
    const BITPIX: BitPix = BitPix::F64;

    fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>> {
        let mut cursor = Cursor::new(bytes);
        let mut result = Vec::with_capacity(bytes.len() / 8);

        while cursor.position() < bytes.len() as u64 {
            let value = match byte_order {
                ByteOrder::BigEndian => cursor.read_f64::<BigEndian>(),
                ByteOrder::LittleEndian => cursor.read_f64::<LittleEndian>(),
            };

            match value {
                Ok(v) => result.push(v),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(FitsError::Io(e)),
            }
        }

        Ok(result)
    }

    fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>> {
        let mut result = Vec::with_capacity(data.len() * 8);

        for &value in data {
            match byte_order {
                ByteOrder::BigEndian => result.write_f64::<BigEndian>(value)?,
                ByteOrder::LittleEndian => result.write_f64::<LittleEndian>(value)?,
            }
        }

        Ok(result)
    }

    fn parse_null_value(null_str: &str) -> Result<Self> {
        null_str.parse::<f64>().map_err(|_| {
            FitsError::InvalidFormat(format!("Invalid NULL value for f64: {}", null_str))
        })
    }

    fn apply_scaling(data: &mut [Self], bscale: f64, bzero: f64) {
        for val in data.iter_mut() {
            *val = (*val * bscale) + bzero;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::ByteOrder;

    #[test]
    fn u8_bitpix_constant() {
        assert_eq!(u8::BITPIX, BitPix::U8);
    }

    #[test]
    fn u8_from_bytes_identity() {
        let bytes = vec![0, 1, 255, 128];
        let result = u8::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, bytes);
    }

    #[test]
    fn u8_from_bytes_ignores_byte_order() {
        let bytes = vec![1, 2, 3];
        let big_endian = u8::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        let little_endian = u8::from_bytes(&bytes, ByteOrder::LittleEndian).unwrap();
        assert_eq!(big_endian, little_endian);
    }

    #[test]
    fn u8_to_bytes_identity() {
        let data = vec![0, 42, 255];
        let result = u8::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, data);
    }

    #[test]
    fn u8_roundtrip() {
        let original = vec![10, 20, 30, 255, 0];
        let bytes = u8::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = u8::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn i16_bitpix_constant() {
        assert_eq!(i16::BITPIX, BitPix::I16);
    }

    #[test]
    fn i16_from_bytes_big_endian() {
        let bytes = vec![0x01, 0x23, 0xFF, 0xFF];
        let result = i16::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![0x0123, -1]);
    }

    #[test]
    fn i16_from_bytes_little_endian() {
        let bytes = vec![0x23, 0x01, 0xFF, 0xFF];
        let result = i16::from_bytes(&bytes, ByteOrder::LittleEndian).unwrap();
        assert_eq!(result, vec![0x0123, -1]);
    }

    #[test]
    fn i16_from_bytes_partial() {
        let bytes = vec![0x01];
        let result = i16::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![]);
    }

    #[test]
    fn i16_to_bytes_big_endian() {
        let data = vec![0x0123, -1];
        let result = i16::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![0x01, 0x23, 0xFF, 0xFF]);
    }

    #[test]
    fn i16_to_bytes_little_endian() {
        let data = vec![0x0123, -1];
        let result = i16::to_bytes(&data, ByteOrder::LittleEndian).unwrap();
        assert_eq!(result, vec![0x23, 0x01, 0xFF, 0xFF]);
    }

    #[test]
    fn i16_roundtrip_big_endian() {
        let original = vec![0, 1, -1, 32767, -32768];
        let bytes = i16::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = i16::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn i16_roundtrip_little_endian() {
        let original = vec![0, 1, -1, 32767, -32768];
        let bytes = i16::to_bytes(&original, ByteOrder::LittleEndian).unwrap();
        let recovered = i16::from_bytes(&bytes, ByteOrder::LittleEndian).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn i32_bitpix_constant() {
        assert_eq!(i32::BITPIX, BitPix::I32);
    }

    #[test]
    fn i32_from_bytes_big_endian() {
        let bytes = vec![0x01, 0x23, 0x45, 0x67, 0xFF, 0xFF, 0xFF, 0xFF];
        let result = i32::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![0x01234567, -1]);
    }

    #[test]
    fn i32_from_bytes_little_endian() {
        let bytes = vec![0x67, 0x45, 0x23, 0x01, 0xFF, 0xFF, 0xFF, 0xFF];
        let result = i32::from_bytes(&bytes, ByteOrder::LittleEndian).unwrap();
        assert_eq!(result, vec![0x01234567, -1]);
    }

    #[test]
    fn i32_to_bytes_big_endian() {
        let data = vec![0x01234567, -1];
        let result = i32::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![0x01, 0x23, 0x45, 0x67, 0xFF, 0xFF, 0xFF, 0xFF]);
    }

    #[test]
    fn i32_roundtrip() {
        let original = vec![0, 1, -1, 2147483647, -2147483648];
        let bytes = i32::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = i32::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn i64_bitpix_constant() {
        assert_eq!(i64::BITPIX, BitPix::I64);
    }

    #[test]
    fn i64_from_bytes_big_endian() {
        let bytes = vec![
            0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
            0xFF, 0xFF,
        ];
        let result = i64::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![0x0123456789ABCDEF, -1]);
    }

    #[test]
    fn i64_to_bytes_big_endian() {
        let data = vec![0x0123456789ABCDEF, -1];
        let result = i64::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        assert_eq!(
            result,
            vec![
                0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
                0xFF, 0xFF
            ]
        );
    }

    #[test]
    fn i64_roundtrip() {
        let original = vec![0, 1, -1, 9223372036854775807, -9223372036854775808];
        let bytes = i64::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = i64::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn f32_bitpix_constant() {
        assert_eq!(f32::BITPIX, BitPix::F32);
    }

    #[test]
    fn f32_from_bytes_big_endian() {
        let bytes = vec![0x3F, 0x80, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00];
        let result = f32::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert!((result[0] - 1.0).abs() < f32::EPSILON);
        assert!((result[1] - 2.0).abs() < f32::EPSILON);
    }

    #[test]
    fn f32_to_bytes_big_endian() {
        let data = vec![1.0, 2.0];
        let result = f32::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, vec![0x3F, 0x80, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00]);
    }

    #[test]
    fn f32_roundtrip() {
        let original = vec![
            0.0,
            1.0,
            -1.0,
            std::f32::consts::PI,
            f32::INFINITY,
            f32::NEG_INFINITY,
        ];
        let bytes = f32::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = f32::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();

        for (orig, recov) in original.iter().zip(recovered.iter()) {
            if orig.is_nan() {
                assert!(recov.is_nan());
            } else {
                assert_eq!(orig, recov);
            }
        }
    }

    #[test]
    fn f32_nan_handling() {
        let original = vec![f32::NAN];
        let bytes = f32::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = f32::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert!(recovered[0].is_nan());
    }

    #[test]
    fn f64_bitpix_constant() {
        assert_eq!(f64::BITPIX, BitPix::F64);
    }

    #[test]
    fn f64_from_bytes_big_endian() {
        let bytes = vec![
            0x3F, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00,
        ];
        let result = f64::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert!((result[0] - 1.0).abs() < f64::EPSILON);
        assert!((result[1] - 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn f64_to_bytes_big_endian() {
        let data = vec![1.0, 2.0];
        let result = f64::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        assert_eq!(
            result,
            vec![
                0x3F, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00
            ]
        );
    }

    #[test]
    fn f64_roundtrip() {
        let original = vec![
            0.0,
            1.0,
            -1.0,
            celestial_core::constants::PI,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ];
        let bytes = f64::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = f64::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();

        for (orig, recov) in original.iter().zip(recovered.iter()) {
            if orig.is_nan() {
                assert!(recov.is_nan());
            } else {
                assert_eq!(orig, recov);
            }
        }
    }

    #[test]
    fn f64_nan_handling() {
        let original = vec![f64::NAN];
        let bytes = f64::to_bytes(&original, ByteOrder::BigEndian).unwrap();
        let recovered = f64::from_bytes(&bytes, ByteOrder::BigEndian).unwrap();
        assert!(recovered[0].is_nan());
    }

    #[test]
    fn empty_input_handling() {
        let empty_bytes: Vec<u8> = vec![];

        assert_eq!(
            u8::from_bytes(&empty_bytes, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );
        assert_eq!(
            i16::from_bytes(&empty_bytes, ByteOrder::BigEndian).unwrap(),
            Vec::<i16>::new()
        );
        assert_eq!(
            i32::from_bytes(&empty_bytes, ByteOrder::BigEndian).unwrap(),
            Vec::<i32>::new()
        );
        assert_eq!(
            i64::from_bytes(&empty_bytes, ByteOrder::BigEndian).unwrap(),
            Vec::<i64>::new()
        );
        assert_eq!(
            f32::from_bytes(&empty_bytes, ByteOrder::BigEndian).unwrap(),
            Vec::<f32>::new()
        );
        assert_eq!(
            f64::from_bytes(&empty_bytes, ByteOrder::BigEndian).unwrap(),
            Vec::<f64>::new()
        );
    }

    #[test]
    fn empty_data_to_bytes() {
        let empty_data: Vec<u8> = vec![];
        assert_eq!(
            u8::to_bytes(&empty_data, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );

        let empty_i16: Vec<i16> = vec![];
        assert_eq!(
            i16::to_bytes(&empty_i16, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );

        let empty_i32: Vec<i32> = vec![];
        assert_eq!(
            i32::to_bytes(&empty_i32, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );

        let empty_i64: Vec<i64> = vec![];
        assert_eq!(
            i64::to_bytes(&empty_i64, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );

        let empty_f32: Vec<f32> = vec![];
        assert_eq!(
            f32::to_bytes(&empty_f32, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );

        let empty_f64: Vec<f64> = vec![];
        assert_eq!(
            f64::to_bytes(&empty_f64, ByteOrder::BigEndian).unwrap(),
            Vec::<u8>::new()
        );
    }

    #[test]
    fn partial_bytes_handling() {
        let partial_i16 = vec![0x01];
        let result = i16::from_bytes(&partial_i16, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, Vec::<i16>::new());

        let partial_i32 = vec![0x01, 0x02, 0x03];
        let result = i32::from_bytes(&partial_i32, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, Vec::<i32>::new());

        let partial_i64 = vec![0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07];
        let result = i64::from_bytes(&partial_i64, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, Vec::<i64>::new());

        let partial_f32 = vec![0x01, 0x02, 0x03];
        let result = f32::from_bytes(&partial_f32, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, Vec::<f32>::new());

        let partial_f64 = vec![0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07];
        let result = f64::from_bytes(&partial_f64, ByteOrder::BigEndian).unwrap();
        assert_eq!(result, Vec::<f64>::new());
    }

    #[test]
    fn data_value_creation() {
        let val = DataValue::Value(42i16);
        let null = DataValue::<i16>::Null;

        assert!(!val.is_null());
        assert!(null.is_null());
        assert_eq!(val.value(), Some(&42));
        assert_eq!(null.value(), None);
    }

    #[test]
    fn data_value_unwrap_or() {
        let val = DataValue::Value(42i16);
        let null = DataValue::<i16>::Null;

        assert_eq!(val.unwrap_or(0), 42);
        assert_eq!(null.unwrap_or(0), 0);
    }

    #[test]
    fn parse_null_value_i16() {
        assert_eq!(i16::parse_null_value("42").unwrap(), 42);
        assert_eq!(i16::parse_null_value("-32768").unwrap(), -32768);
        assert!(i16::parse_null_value("invalid").is_err());
        assert!(i16::parse_null_value("100000").is_err());
    }

    #[test]
    fn parse_null_value_f32() {
        assert_eq!(f32::parse_null_value("2.75").unwrap(), 2.75);
        assert_eq!(f32::parse_null_value("-1.0").unwrap(), -1.0);
        assert!(f32::parse_null_value("NaN").unwrap().is_nan());
        assert_eq!(f32::parse_null_value("inf").unwrap(), f32::INFINITY);
        assert!(f32::parse_null_value("invalid").is_err());
    }

    #[test]
    fn from_bytes_with_null_i16() {
        let bytes = vec![0x00, 0x01, 0x00, 0x02, 0xFF, 0xFF];
        let null_value = Some(-1i16);

        let result = i16::from_bytes_with_null(&bytes, ByteOrder::BigEndian, null_value).unwrap();

        assert_eq!(result.len(), 3);
        assert_eq!(result[0], DataValue::Value(1));
        assert_eq!(result[1], DataValue::Value(2));
        assert_eq!(result[2], DataValue::Null);
    }

    #[test]
    fn from_bytes_with_null_no_nulls() {
        let bytes = vec![0x00, 0x01, 0x00, 0x02];
        let null_value = Some(-1i16);

        let result = i16::from_bytes_with_null(&bytes, ByteOrder::BigEndian, null_value).unwrap();

        assert_eq!(result.len(), 2);
        assert_eq!(result[0], DataValue::Value(1));
        assert_eq!(result[1], DataValue::Value(2));
    }

    #[test]
    fn from_bytes_with_null_disabled() {
        let bytes = vec![0x00, 0x01, 0xFF, 0xFF];
        let null_value = None;

        let result = i16::from_bytes_with_null(&bytes, ByteOrder::BigEndian, null_value).unwrap();

        assert_eq!(result.len(), 2);
        assert_eq!(result[0], DataValue::Value(1));
        assert_eq!(result[1], DataValue::Value(-1));
    }

    #[test]
    fn from_bytes_with_null_f64() {
        let data = vec![1.0, f64::NAN, 3.0];
        let bytes = f64::to_bytes(&data, ByteOrder::BigEndian).unwrap();
        let null_value = Some(f64::NAN);

        let result = f64::from_bytes_with_null(&bytes, ByteOrder::BigEndian, null_value).unwrap();

        assert_eq!(result.len(), 3);
        assert_eq!(result[0], DataValue::Value(1.0));
        assert!(matches!(result[1], DataValue::Value(v) if v.is_nan()));
        assert_eq!(result[2], DataValue::Value(3.0));
    }

    #[test]
    fn parse_null_value_all_integer_types() {
        assert_eq!(u8::parse_null_value("0").unwrap(), 0);
        assert_eq!(u8::parse_null_value("255").unwrap(), 255);
        assert!(u8::parse_null_value("256").is_err());
        assert!(u8::parse_null_value("-1").is_err());

        assert_eq!(i32::parse_null_value("0").unwrap(), 0);
        assert_eq!(i32::parse_null_value("-2147483648").unwrap(), -2147483648);
        assert_eq!(i32::parse_null_value("2147483647").unwrap(), 2147483647);
        assert!(i32::parse_null_value("9999999999").is_err());

        assert_eq!(i64::parse_null_value("0").unwrap(), 0);
        assert_eq!(
            i64::parse_null_value("-9223372036854775808").unwrap(),
            -9223372036854775808
        );
        assert_eq!(
            i64::parse_null_value("9223372036854775807").unwrap(),
            9223372036854775807
        );
        assert!(i64::parse_null_value("invalid").is_err());
    }

    #[test]
    fn parse_null_value_all_float_types() {
        assert_eq!(f64::parse_null_value("0.0").unwrap(), 0.0);
        assert_eq!(f64::parse_null_value("-1.5").unwrap(), -1.5);
        assert_eq!(f64::parse_null_value("1.23456").unwrap(), 1.23456);
        assert_eq!(f64::parse_null_value("inf").unwrap(), f64::INFINITY);
        assert_eq!(f64::parse_null_value("-inf").unwrap(), f64::NEG_INFINITY);
        assert!(f64::parse_null_value("NaN").unwrap().is_nan());
        assert!(f64::parse_null_value("not_a_number").is_err());
    }

    #[test]
    fn from_bytes_with_null_all_types() {
        let u8_bytes = vec![1, 255, 0, 128];
        let u8_result =
            u8::from_bytes_with_null(&u8_bytes, ByteOrder::BigEndian, Some(255)).unwrap();
        assert_eq!(u8_result[0], DataValue::Value(1));
        assert_eq!(u8_result[1], DataValue::Null);
        assert_eq!(u8_result[2], DataValue::Value(0));
        assert_eq!(u8_result[3], DataValue::Value(128));

        let i32_bytes = vec![0x00, 0x00, 0x00, 0x01, 0x80, 0x00, 0x00, 0x00];
        let i32_result =
            i32::from_bytes_with_null(&i32_bytes, ByteOrder::BigEndian, Some(-2147483648)).unwrap();
        assert_eq!(i32_result[0], DataValue::Value(1));
        assert_eq!(i32_result[1], DataValue::Null);

        let i64_bytes = vec![
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00,
        ];
        let i64_result =
            i64::from_bytes_with_null(&i64_bytes, ByteOrder::BigEndian, Some(-9223372036854775808))
                .unwrap();
        assert_eq!(i64_result[0], DataValue::Value(1));
        assert_eq!(i64_result[1], DataValue::Null);

        let f32_bytes = vec![0x3F, 0x80, 0x00, 0x00, 0xFF, 0x80, 0x00, 0x01];
        let null_nan = f32::from_bits(0xFF800001);
        let f32_result =
            f32::from_bytes_with_null(&f32_bytes, ByteOrder::BigEndian, Some(null_nan)).unwrap();
        assert_eq!(f32_result[0], DataValue::Value(1.0));
        assert!(matches!(f32_result[1], DataValue::Value(_)));
    }

    #[test]
    fn data_value_operations() {
        let value = DataValue::Value(42);
        let null = DataValue::<i32>::Null;

        assert!(!value.is_null());
        assert!(null.is_null());

        assert_eq!(value.value(), Some(&42));
        assert_eq!(null.value(), None);

        let cloned_value = value.clone();
        assert_eq!(cloned_value, DataValue::Value(42));

        let cloned_null = null.clone();
        assert_eq!(cloned_null, DataValue::<i32>::Null);

        assert_eq!(value.unwrap_or(-1), 42);
        assert_eq!(null.clone().unwrap_or(-1), -1);
        assert_eq!(null.unwrap_or(999), 999);
    }

    #[test]
    fn from_bytes_with_null_edge_cases() {
        let empty_result = i16::from_bytes_with_null(&[], ByteOrder::BigEndian, Some(-1)).unwrap();
        assert_eq!(empty_result.len(), 0);

        let all_null_bytes = vec![0xFF, 0xFF, 0xFF, 0xFF];
        let all_null_result =
            i16::from_bytes_with_null(&all_null_bytes, ByteOrder::BigEndian, Some(-1)).unwrap();
        assert_eq!(all_null_result.len(), 2);
        assert_eq!(all_null_result[0], DataValue::Null);
        assert_eq!(all_null_result[1], DataValue::Null);

        let no_null_result =
            i16::from_bytes_with_null(&all_null_bytes, ByteOrder::BigEndian, None).unwrap();
        assert_eq!(no_null_result.len(), 2);
        assert_eq!(no_null_result[0], DataValue::Value(-1));
        assert_eq!(no_null_result[1], DataValue::Value(-1));
    }

    #[test]
    fn parse_null_value_error_cases() {
        assert!(u8::parse_null_value("abc").is_err());
        assert!(u8::parse_null_value("3.14").is_err());

        assert!(i16::parse_null_value("not_a_number").is_err());
        assert!(i16::parse_null_value("100000").is_err());

        assert!(i32::parse_null_value("10000000000").is_err());
        assert!(i32::parse_null_value("3.14159").is_err());

        assert!(i64::parse_null_value("18446744073709551616").is_err());

        assert!(f32::parse_null_value("not_float").is_err());
        assert!(f64::parse_null_value("also_not_float").is_err());
    }

    #[test]
    fn apply_scaling_u8_identity() {
        let mut data = vec![0u8, 100, 200, 255];
        u8::apply_scaling(&mut data, 1.0, 0.0);
        assert_eq!(data, vec![0, 100, 200, 255]);
    }

    #[test]
    fn apply_scaling_u8_with_bzero() {
        let mut data = vec![0u8, 50, 100, 150];
        u8::apply_scaling(&mut data, 1.0, 10.0);
        assert_eq!(data, vec![10, 60, 110, 160]);
    }

    #[test]
    fn apply_scaling_u8_with_bscale() {
        let mut data = vec![0u8, 50, 100];
        u8::apply_scaling(&mut data, 2.0, 0.0);
        assert_eq!(data, vec![0, 100, 200]);
    }

    #[test]
    fn apply_scaling_u8_clamping() {
        let mut data = vec![200u8, 250];
        u8::apply_scaling(&mut data, 2.0, 0.0);
        assert_eq!(data, vec![255, 255]);
    }

    #[test]
    fn apply_scaling_i16_identity() {
        let mut data = vec![-100i16, 0, 100, 1000];
        i16::apply_scaling(&mut data, 1.0, 0.0);
        assert_eq!(data, vec![-100, 0, 100, 1000]);
    }

    #[test]
    fn apply_scaling_i16_unsigned_to_signed() {
        let mut data = vec![0i16, 32767, -32768, -1];
        i16::apply_scaling(&mut data, 1.0, -32768.0);
        assert_eq!(data, vec![-32768, -1, -32768, -32768]);
    }

    #[test]
    fn apply_scaling_i16_with_bscale() {
        let mut data = vec![10i16, 20, 30];
        i16::apply_scaling(&mut data, 100.0, 0.0);
        assert_eq!(data, vec![1000, 2000, 3000]);
    }

    #[test]
    fn apply_scaling_i32_bzero_bscale() {
        let mut data = vec![0i32, 1, 2, 3];
        i32::apply_scaling(&mut data, 2.0, 100.0);
        assert_eq!(data, vec![100, 102, 104, 106]);
    }

    #[test]
    fn apply_scaling_i64_identity() {
        let mut data = vec![i64::MIN, 0, i64::MAX];
        i64::apply_scaling(&mut data, 1.0, 0.0);
        assert_eq!(data, vec![i64::MIN, 0, i64::MAX]);
    }

    #[test]
    fn apply_scaling_f32_identity() {
        let mut data = vec![1.0f32, 2.0, 3.0];
        f32::apply_scaling(&mut data, 1.0, 0.0);
        assert_eq!(data, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn apply_scaling_f32_bzero_bscale() {
        let mut data = vec![0.0f32, 0.5, 1.0];
        f32::apply_scaling(&mut data, 2.0, 10.0);
        assert_eq!(data, vec![10.0, 11.0, 12.0]);
    }

    #[test]
    fn apply_scaling_f64_identity() {
        let mut data = vec![1.0f64, 2.0, 3.0];
        f64::apply_scaling(&mut data, 1.0, 0.0);
        assert_eq!(data, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn apply_scaling_f64_bzero_bscale() {
        let mut data = vec![0.0f64, 0.5, 1.0];
        f64::apply_scaling(&mut data, 2.0, 10.0);
        assert_eq!(data, vec![10.0, 11.0, 12.0]);
    }

    #[test]
    fn apply_scaling_f64_astronomical_values() {
        let mut data = vec![32768.0f64];
        f64::apply_scaling(&mut data, 1.0, -32768.0);
        assert_eq!(data, vec![0.0]);
    }

    #[test]
    fn table_value_as_i64_integer_types() {
        assert_eq!(TableValue::Byte(42).as_i64(), Some(42));
        assert_eq!(TableValue::I16(1000).as_i64(), Some(1000));
        assert_eq!(TableValue::I32(-500).as_i64(), Some(-500));
        assert_eq!(TableValue::I64(123456789).as_i64(), Some(123456789));
    }

    #[test]
    fn table_value_as_i64_non_integer_types() {
        assert_eq!(TableValue::F32(1.5).as_i64(), None);
        assert_eq!(TableValue::F64(2.75).as_i64(), None);
        assert_eq!(TableValue::String("test".to_string()).as_i64(), None);
        assert_eq!(TableValue::Null.as_i64(), None);
        assert_eq!(TableValue::Logical(true).as_i64(), None);
    }

    #[test]
    fn table_value_as_f64_numeric_types() {
        assert_eq!(TableValue::Byte(42).as_f64(), Some(42.0));
        assert_eq!(TableValue::I16(1000).as_f64(), Some(1000.0));
        assert_eq!(TableValue::I32(-500).as_f64(), Some(-500.0));
        assert_eq!(TableValue::I64(123456789).as_f64(), Some(123456789.0));
        assert_eq!(TableValue::F32(1.5).as_f64(), Some(1.5));
        assert_eq!(TableValue::F64(2.75).as_f64(), Some(2.75));
    }

    #[test]
    fn table_value_as_f64_non_numeric_types() {
        assert_eq!(TableValue::String("test".to_string()).as_f64(), None);
        assert_eq!(TableValue::Null.as_f64(), None);
        assert_eq!(TableValue::Logical(true).as_f64(), None);
        assert_eq!(TableValue::Complex32(1.0, 2.0).as_f64(), None);
    }

    #[test]
    fn table_value_as_string() {
        assert_eq!(
            TableValue::String("hello".to_string()).as_string(),
            Some("hello")
        );
        assert_eq!(TableValue::I32(42).as_string(), None);
        assert_eq!(TableValue::Null.as_string(), None);
    }

    #[test]
    fn table_value_is_null() {
        assert!(TableValue::Null.is_null());
        assert!(!TableValue::I32(0).is_null());
        assert!(!TableValue::String("".to_string()).is_null());
    }

    #[test]
    fn table_value_equality() {
        assert_eq!(TableValue::I32(42), TableValue::I32(42));
        assert_ne!(TableValue::I32(42), TableValue::I32(0));
        assert_ne!(TableValue::I32(42), TableValue::I64(42));
        assert_eq!(TableValue::Null, TableValue::Null);
    }

    #[test]
    fn table_value_clone() {
        let original = TableValue::String("test".to_string());
        let cloned = original.clone();
        assert_eq!(original, cloned);
    }

    #[test]
    fn table_value_complex_types() {
        let c32 = TableValue::Complex32(1.0, 2.0);
        let c64 = TableValue::Complex64(3.0, 4.0);

        assert_eq!(c32, TableValue::Complex32(1.0, 2.0));
        assert_eq!(c64, TableValue::Complex64(3.0, 4.0));
        assert_ne!(c32, c64);
    }
}
