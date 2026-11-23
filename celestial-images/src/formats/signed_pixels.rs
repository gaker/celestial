use super::pixel_data::PixelData;
use crate::core::errors::{ImageError, Result};

impl PixelData {
    // Adds the FITS 4.0 Table 11 offset, 2^(n-1) for n-bit samples, which maps
    // the signed range onto the unsigned one. Flipping the sign bit is that
    // addition modulo 2^n.
    pub fn signed_to_offset_unsigned(&self) -> Result<PixelData> {
        match self {
            Self::I8(v) => Ok(Self::U8(v.iter().map(|&x| x as u8 ^ 0x80).collect())),
            Self::I16(v) => Ok(Self::U16(v.iter().map(|&x| x as u16 ^ 0x8000).collect())),
            Self::I32(v) => Ok(Self::U32(v.iter().map(|&x| x as u32 ^ (1 << 31)).collect())),
            Self::I64(v) => Ok(Self::U64(v.iter().map(|&x| x as u64 ^ (1 << 63)).collect())),
            _ => Err(not_signed(self)),
        }
    }

    // f32 holds every i8 and i16 exactly, and f64 every i32. No float type
    // holds every i64, so I64 is refused.
    pub fn signed_to_float(&self) -> Result<PixelData> {
        match self {
            Self::I8(v) => Ok(Self::F32(v.iter().map(|&x| f32::from(x)).collect())),
            Self::I16(v) => Ok(Self::F32(v.iter().map(|&x| f32::from(x)).collect())),
            Self::I32(v) => Ok(Self::F64(v.iter().map(|&x| f64::from(x)).collect())),
            Self::I64(_) => Err(ImageError::PixelConversion(
                "I64 samples have no exact float form; use \
                 PixelData::signed_to_offset_unsigned"
                    .to_string(),
            )),
            _ => Err(not_signed(self)),
        }
    }
}

fn not_signed(pixels: &PixelData) -> ImageError {
    let held = match pixels {
        PixelData::U8(_) => "U8",
        PixelData::U16(_) => "U16",
        PixelData::U32(_) => "U32",
        PixelData::U64(_) => "U64",
        PixelData::F32(_) => "F32",
        PixelData::F64(_) => "F64",
        PixelData::I8(_) | PixelData::I16(_) | PixelData::I32(_) | PixelData::I64(_) => "signed",
    };
    ImageError::PixelConversion(format!(
        "the conversion takes I8, I16, I32 or I64 samples, not {held}"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn offset_maps_each_signed_type_onto_its_unsigned_range() {
        let i8s = PixelData::I8(vec![i8::MIN, -1, 0, 1, i8::MAX]);
        let i16s = PixelData::I16(vec![i16::MIN, -1, 0, 1, i16::MAX]);
        let i32s = PixelData::I32(vec![i32::MIN, -1, 0, 1, i32::MAX]);
        let i64s = PixelData::I64(vec![i64::MIN, -1, 0, 1, i64::MAX]);
        assert_eq!(
            i8s.signed_to_offset_unsigned().unwrap().as_u8(),
            Some(&vec![0, 127, 128, 129, 255])
        );
        assert_eq!(
            i16s.signed_to_offset_unsigned().unwrap().as_u16(),
            Some(&vec![0, 32767, 32768, 32769, 65535])
        );
        assert_eq!(
            i32s.signed_to_offset_unsigned().unwrap().as_u32(),
            Some(&vec![0, 2147483647, 2147483648, 2147483649, 4294967295])
        );
        assert_eq!(
            i64s.signed_to_offset_unsigned().unwrap().as_u64(),
            Some(&vec![
                0,
                9223372036854775807,
                9223372036854775808,
                9223372036854775809,
                18446744073709551615
            ])
        );
    }

    #[test]
    fn float_conversion_is_exact_at_the_type_extremes() {
        let i8s = PixelData::I8(vec![i8::MIN, -1, 0, i8::MAX]);
        let i16s = PixelData::I16(vec![i16::MIN, -1, 0, i16::MAX]);
        let i32s = PixelData::I32(vec![i32::MIN, -1, 0, i32::MAX]);
        assert_eq!(
            i8s.signed_to_float().unwrap().as_f32(),
            Some(&vec![-128.0, -1.0, 0.0, 127.0])
        );
        assert_eq!(
            i16s.signed_to_float().unwrap().as_f32(),
            Some(&vec![-32768.0, -1.0, 0.0, 32767.0])
        );
        assert_eq!(
            i32s.signed_to_float().unwrap().as_f64(),
            Some(&vec![-2147483648.0, -1.0, 0.0, 2147483647.0])
        );
    }

    #[test]
    fn float_conversion_refuses_i64() {
        let err = PixelData::I64(vec![i64::MAX])
            .signed_to_float()
            .unwrap_err();
        assert!(matches!(err, ImageError::PixelConversion(_)));
        assert!(err
            .to_string()
            .contains("PixelData::signed_to_offset_unsigned"));
    }

    #[test]
    fn conversions_refuse_samples_that_are_not_signed() {
        for pixels in [
            PixelData::U8(vec![1]),
            PixelData::U16(vec![1]),
            PixelData::U32(vec![1]),
            PixelData::U64(vec![1]),
            PixelData::F32(vec![1.0]),
            PixelData::F64(vec![1.0]),
        ] {
            assert!(matches!(
                pixels.signed_to_offset_unsigned(),
                Err(ImageError::PixelConversion(_))
            ));
            assert!(matches!(
                pixels.signed_to_float(),
                Err(ImageError::PixelConversion(_))
            ));
        }
    }
}
