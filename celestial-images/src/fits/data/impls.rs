use super::array::{decode, DataArray};
use crate::core::types::{BitPix, ByteOrder};
use crate::fits::errors::Result;

macro_rules! data_array {
    ($($sample:ty => $bitpix:expr),* $(,)?) => {$(
        impl DataArray for $sample {
            const BITPIX: BitPix = $bitpix;

            fn from_bytes(bytes: &[u8], byte_order: ByteOrder) -> Result<Vec<Self>> {
                match byte_order {
                    ByteOrder::BigEndian => decode(bytes, <$sample>::from_be_bytes),
                    ByteOrder::LittleEndian => decode(bytes, <$sample>::from_le_bytes),
                }
            }

            fn to_bytes(data: &[Self], byte_order: ByteOrder) -> Result<Vec<u8>> {
                Ok(match byte_order {
                    ByteOrder::BigEndian => encode(data, |v| v.to_be_bytes()),
                    ByteOrder::LittleEndian => encode(data, |v| v.to_le_bytes()),
                })
            }
        }
    )*};
}

data_array!(
    u8 => BitPix::U8,
    i16 => BitPix::I16,
    i32 => BitPix::I32,
    i64 => BitPix::I64,
    f32 => BitPix::F32,
    f64 => BitPix::F64,
);

fn encode<const N: usize, T: Copy>(data: &[T], bytes: impl Fn(T) -> [u8; N]) -> Vec<u8> {
    let mut out = vec![0u8; data.len() * N];
    for (sample, &value) in out.as_chunks_mut::<N>().0.iter_mut().zip(data) {
        *sample = bytes(value);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::ByteOrder;

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
        assert_eq!(result, vec![1.0f32, 2.0]);
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
        assert_eq!(result, vec![1.0f64, 2.0]);
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
    fn from_bytes_rejects_a_partial_sample() {
        fn shown<T>(result: Result<Vec<T>>) -> Option<String> {
            result.err().map(|e| e.to_string())
        }
        let message = |len: usize, size: usize| {
            Some(format!(
                "Invalid FITS format: the data is {len} bytes, not a whole number of {size}-byte samples"
            ))
        };
        let bytes = [1u8, 2, 3, 4, 5, 6, 7];
        for order in [ByteOrder::BigEndian, ByteOrder::LittleEndian] {
            assert_eq!(shown(i16::from_bytes(&bytes[..1], order)), message(1, 2));
            assert_eq!(shown(i16::from_bytes(&bytes, order)), message(7, 2));
            assert_eq!(shown(i32::from_bytes(&bytes, order)), message(7, 4));
            assert_eq!(shown(i64::from_bytes(&bytes, order)), message(7, 8));
            assert_eq!(shown(f32::from_bytes(&bytes[..3], order)), message(3, 4));
            assert_eq!(shown(f64::from_bytes(&bytes, order)), message(7, 8));
        }
    }
}
