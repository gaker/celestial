use crate::xisf::header::sample_format::SampleFormat;

// Callers can name `XisfDataType` in bounds but cannot implement it. The
// reader and writer support only the sample types listed here.
mod sealed {
    pub trait Sealed {}

    impl Sealed for u8 {}
    impl Sealed for u16 {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

pub trait XisfDataType: sealed::Sealed + Copy {
    const SAMPLE_FORMAT: SampleFormat;
    fn to_le_bytes(data: &[Self]) -> Vec<u8>;
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self>;
    fn default_bounds() -> (f64, f64);
    fn calculate_bounds(data: &[Self]) -> (f64, f64);
}

impl XisfDataType for u8 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt8;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        data.to_vec()
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes.to_vec()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 255.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for u16 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt16;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 2);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .as_chunks::<2>()
            .0
            .iter()
            .map(|&c| Self::from_le_bytes(c))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 65535.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for u32 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt32;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 4);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|&c| Self::from_le_bytes(c))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 4294967295.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for u64 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt64;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 8);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|&c| Self::from_le_bytes(c))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 18446744073709551615.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for f32 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::Float32;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 4);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|&c| Self::from_le_bytes(c))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 1.0)
    }
    fn calculate_bounds(data: &[Self]) -> (f64, f64) {
        float_bounds(data.iter().map(|&v| f64::from(v)))
    }
}

impl XisfDataType for f64 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::Float64;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 8);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|&c| Self::from_le_bytes(c))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 1.0)
    }
    fn calculate_bounds(data: &[Self]) -> (f64, f64) {
        float_bounds(data.iter().copied())
    }
}

// The range of the finite samples. XISF 1.0 section 11.5.1 makes the bounds
// a black point and a white point, so a range of one value, or of none, is
// widened to include [0, 1].
fn float_bounds(samples: impl Iterator<Item = f64>) -> (f64, f64) {
    let (min, max) = samples
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), v| {
            (min.min(v), max.max(v))
        });
    if min < max {
        (min, max)
    } else {
        (min.min(0.0), max.max(1.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f32::consts::PI as F32_PI;
    use core::f64::consts::{E as F64_E, PI as F64_PI};

    #[test]
    fn sample_format_constants_match_each_impl() {
        assert!(matches!(
            <u8 as XisfDataType>::SAMPLE_FORMAT,
            SampleFormat::UInt8
        ));
        assert!(matches!(
            <u16 as XisfDataType>::SAMPLE_FORMAT,
            SampleFormat::UInt16
        ));
        assert!(matches!(
            <u32 as XisfDataType>::SAMPLE_FORMAT,
            SampleFormat::UInt32
        ));
        assert!(matches!(
            <u64 as XisfDataType>::SAMPLE_FORMAT,
            SampleFormat::UInt64
        ));
        assert!(matches!(
            <f32 as XisfDataType>::SAMPLE_FORMAT,
            SampleFormat::Float32
        ));
        assert!(matches!(
            <f64 as XisfDataType>::SAMPLE_FORMAT,
            SampleFormat::Float64
        ));
    }

    #[test]
    fn u64_roundtrip_preserves_values_and_encodes_little_endian() {
        let original: Vec<u64> = vec![0, (1 << 53) + 1, 0x0123_4567_89AB_CDEF, u64::MAX];
        let bytes = <u64 as XisfDataType>::to_le_bytes(&original);
        assert_eq!(&bytes[16..24], &0x0123_4567_89AB_CDEFu64.to_le_bytes());
        assert_eq!(<u64 as XisfDataType>::from_le_bytes(&bytes), original);
        assert_eq!(
            <u64 as XisfDataType>::default_bounds(),
            (0.0, 18446744073709551615.0)
        );
    }

    #[test]
    fn u8_roundtrip_preserves_every_byte() {
        let original: Vec<u8> = (0u8..=255).collect();
        let bytes = <u8 as XisfDataType>::to_le_bytes(&original);
        assert_eq!(bytes, original);
        let restored = <u8 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(restored, original);
    }

    #[test]
    fn u16_roundtrip_preserves_values_and_encodes_little_endian() {
        let original: Vec<u16> = vec![0, 1, 0x00FF, 0xFF00, 0x1234, u16::MAX];
        let bytes = <u16 as XisfDataType>::to_le_bytes(&original);
        assert_eq!(bytes.len(), original.len() * 2);
        assert_eq!(&bytes[0..2], &0u16.to_le_bytes());
        assert_eq!(&bytes[8..10], &0x1234u16.to_le_bytes());
        let restored = <u16 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(restored, original);
    }

    #[test]
    fn u32_roundtrip_preserves_values_and_encodes_little_endian() {
        let original: Vec<u32> = vec![0, 1, 0xDEADBEEF, 0x12345678, u32::MAX];
        let bytes = <u32 as XisfDataType>::to_le_bytes(&original);
        assert_eq!(bytes.len(), original.len() * 4);
        assert_eq!(&bytes[8..12], &0xDEADBEEFu32.to_le_bytes());
        let restored = <u32 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(restored, original);
    }

    #[test]
    fn f32_roundtrip_preserves_finite_values() {
        let original: Vec<f32> = vec![0.0, 1.0, -1.5, f32::MIN, f32::MAX, F32_PI];
        let bytes = <f32 as XisfDataType>::to_le_bytes(&original);
        assert_eq!(bytes.len(), original.len() * 4);
        let restored = <f32 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(restored, original);
    }

    #[test]
    fn f64_roundtrip_preserves_finite_values() {
        let original: Vec<f64> = vec![0.0, 1.0, -1.5, f64::MIN, f64::MAX, F64_E];
        let bytes = <f64 as XisfDataType>::to_le_bytes(&original);
        assert_eq!(bytes.len(), original.len() * 8);
        let restored = <f64 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(restored, original);
    }

    #[test]
    fn empty_roundtrips_produce_empty_output() {
        assert!(<u8 as XisfDataType>::to_le_bytes(&[]).is_empty());
        assert!(<u16 as XisfDataType>::to_le_bytes(&[]).is_empty());
        assert!(<u32 as XisfDataType>::to_le_bytes(&[]).is_empty());
        assert!(<f32 as XisfDataType>::to_le_bytes(&[]).is_empty());
        assert!(<f64 as XisfDataType>::to_le_bytes(&[]).is_empty());
    }

    #[test]
    fn from_le_bytes_drops_trailing_partial_chunks() {
        let v = <u16 as XisfDataType>::from_le_bytes(&[0x01, 0x02, 0x03]);
        assert_eq!(v, vec![u16::from_le_bytes([0x01, 0x02])]);

        let v = <u32 as XisfDataType>::from_le_bytes(&[0, 0, 0, 0, 1]);
        assert_eq!(v, vec![0u32]);

        let v = <f32 as XisfDataType>::from_le_bytes(&[0; 5]);
        assert_eq!(v, vec![0.0f32]);

        let v = <f64 as XisfDataType>::from_le_bytes(&[0; 9]);
        assert_eq!(v, vec![0.0f64]);
    }

    #[test]
    fn default_bounds_per_type() {
        assert_eq!(<u8 as XisfDataType>::default_bounds(), (0.0, 255.0));
        assert_eq!(<u16 as XisfDataType>::default_bounds(), (0.0, 65535.0));
        assert_eq!(<u32 as XisfDataType>::default_bounds(), (0.0, 4294967295.0));
        assert_eq!(<f32 as XisfDataType>::default_bounds(), (0.0, 1.0));
        assert_eq!(<f64 as XisfDataType>::default_bounds(), (0.0, 1.0));
    }

    #[test]
    fn integer_calculate_bounds_always_returns_default() {
        assert_eq!(<u8 as XisfDataType>::calculate_bounds(&[]), (0.0, 255.0));
        assert_eq!(<u8 as XisfDataType>::calculate_bounds(&[42]), (0.0, 255.0));

        assert_eq!(<u16 as XisfDataType>::calculate_bounds(&[]), (0.0, 65535.0));
        assert_eq!(
            <u16 as XisfDataType>::calculate_bounds(&[7]),
            (0.0, 65535.0)
        );

        assert_eq!(
            <u32 as XisfDataType>::calculate_bounds(&[]),
            (0.0, 4294967295.0)
        );
        assert_eq!(
            <u32 as XisfDataType>::calculate_bounds(&[1, 2]),
            (0.0, 4294967295.0)
        );
    }

    #[test]
    fn f32_calculate_bounds_empty_returns_default() {
        assert_eq!(<f32 as XisfDataType>::calculate_bounds(&[]), (0.0, 1.0));
    }

    #[test]
    fn f32_calculate_bounds_finite_range() {
        let data = [0.0f32, -5.0, 12.5, 3.25];
        let (min, max) = <f32 as XisfDataType>::calculate_bounds(&data);
        assert_eq!(min, -5.0);
        assert_eq!(max, 12.5);
    }

    #[test]
    fn f32_calculate_bounds_all_nan_returns_default() {
        let data = [f32::NAN, f32::NAN, f32::NAN];
        assert_eq!(<f32 as XisfDataType>::calculate_bounds(&data), (0.0, 1.0));
    }

    #[test]
    fn f32_calculate_bounds_mixed_finite_and_nan_ignores_nan() {
        let data = [
            1.0f32,
            f32::NAN,
            -2.0,
            f32::INFINITY,
            f32::NEG_INFINITY,
            3.0,
        ];
        let (min, max) = <f32 as XisfDataType>::calculate_bounds(&data);
        assert_eq!(min, -2.0);
        assert_eq!(max, 3.0);
    }

    #[test]
    fn f64_calculate_bounds_empty_returns_default() {
        assert_eq!(<f64 as XisfDataType>::calculate_bounds(&[]), (0.0, 1.0));
    }

    #[test]
    fn f64_calculate_bounds_finite_range() {
        let data = [0.0, -100.0, 50.0, F64_PI, -0.5];
        let (min, max) = <f64 as XisfDataType>::calculate_bounds(&data);
        assert_eq!(min, -100.0);
        assert_eq!(max, 50.0);
    }

    #[test]
    fn f64_calculate_bounds_all_nan_returns_default() {
        let data = [f64::NAN, f64::NAN];
        assert_eq!(<f64 as XisfDataType>::calculate_bounds(&data), (0.0, 1.0));
    }

    #[test]
    fn f64_calculate_bounds_mixed_finite_and_nonfinite() {
        let data = [f64::NAN, 1.0, f64::INFINITY, f64::NEG_INFINITY, -4.5, 8.25];
        let (min, max) = <f64 as XisfDataType>::calculate_bounds(&data);
        assert_eq!(min, -4.5);
        assert_eq!(max, 8.25);
    }
}
