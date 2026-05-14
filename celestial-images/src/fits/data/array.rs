use super::value::DataValue;
use crate::core::{BitPix, ByteOrder};
use crate::fits::Result;

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
