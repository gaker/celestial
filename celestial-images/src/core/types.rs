#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BitPix {
    U8 = 8,
    I16 = 16,
    I32 = 32,
    I64 = 64,
    F32 = -32,
    F64 = -64,
}

impl BitPix {
    pub fn from_value(value: i32) -> Option<Self> {
        match value {
            8 => Some(Self::U8),
            16 => Some(Self::I16),
            32 => Some(Self::I32),
            64 => Some(Self::I64),
            -32 => Some(Self::F32),
            -64 => Some(Self::F64),
            _ => None,
        }
    }

    pub fn value(self) -> i32 {
        self as i32
    }

    pub fn bytes_per_pixel(self) -> usize {
        match self {
            Self::U8 => 1,
            Self::I16 => 2,
            Self::I32 | Self::F32 => 4,
            Self::I64 | Self::F64 => 8,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ByteOrder {
    #[default]
    BigEndian,
    LittleEndian,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bitpix_valid_values() {
        let valid_cases = [
            (8, BitPix::U8),
            (16, BitPix::I16),
            (32, BitPix::I32),
            (64, BitPix::I64),
            (-32, BitPix::F32),
            (-64, BitPix::F64),
        ];

        for (input, expected) in valid_cases {
            assert_eq!(BitPix::from_value(input), Some(expected));
            assert_eq!(expected.value(), input);
        }
    }

    #[test]
    fn bitpix_invalid_values() {
        let invalid_values = [
            0,
            1,
            -1,
            7,
            9,
            15,
            17,
            31,
            33,
            63,
            65,
            -8,
            -16,
            -31,
            -33,
            -63,
            -65,
            i32::MIN,
            i32::MAX,
            24,
            48,
            128,
            256,
        ];

        for &invalid in &invalid_values {
            assert_eq!(BitPix::from_value(invalid), None);
        }
    }

    #[test]
    fn bitpix_bytes_per_pixel() {
        assert_eq!(BitPix::U8.bytes_per_pixel(), 1);
        assert_eq!(BitPix::I16.bytes_per_pixel(), 2);
        assert_eq!(BitPix::I32.bytes_per_pixel(), 4);
        assert_eq!(BitPix::I64.bytes_per_pixel(), 8);
        assert_eq!(BitPix::F32.bytes_per_pixel(), 4);
        assert_eq!(BitPix::F64.bytes_per_pixel(), 8);
    }

    #[test]
    fn bitpix_repr_c_layout() {
        assert_eq!(BitPix::U8 as i32, 8);
        assert_eq!(BitPix::I16 as i32, 16);
        assert_eq!(BitPix::I32 as i32, 32);
        assert_eq!(BitPix::I64 as i32, 64);
        assert_eq!(BitPix::F32 as i32, -32);
        assert_eq!(BitPix::F64 as i32, -64);
    }

    #[test]
    fn bitpix_equality_and() {
        assert_eq!(BitPix::U8, BitPix::U8);
        assert_ne!(BitPix::U8, BitPix::I16);

        assert!(format!("{:?}", BitPix::F32).contains("F32"));
        assert!(format!("{:?}", BitPix::I64).contains("I64"));
    }

    #[test]
    fn bitpix_clone_and_copy() {
        let original = BitPix::F64;
        let copied = original;

        assert_eq!(BitPix::F64, copied);
    }

    #[test]
    fn byteorder_default() {
        assert_eq!(ByteOrder::default(), ByteOrder::BigEndian);
    }

    #[test]
    fn byteorder_equality_and() {
        assert_eq!(ByteOrder::BigEndian, ByteOrder::BigEndian);
        assert_ne!(ByteOrder::BigEndian, ByteOrder::LittleEndian);

        assert!(format!("{:?}", ByteOrder::BigEndian).contains("BigEndian"));
        assert!(format!("{:?}", ByteOrder::LittleEndian).contains("LittleEndian"));
    }

    #[test]
    fn byteorder_clone_and_copy() {
        let original = ByteOrder::LittleEndian;
        let cloned = original;
        let copied = original;

        assert_eq!(original, cloned);
        assert_eq!(original, copied);
    }

    #[test]
    fn bitpix_roundtrip_conversion() {
        let all_variants = [
            BitPix::U8,
            BitPix::I16,
            BitPix::I32,
            BitPix::I64,
            BitPix::F32,
            BitPix::F64,
        ];

        for &variant in &all_variants {
            let value = variant.value();
            let reconstructed = BitPix::from_value(value).unwrap();
            assert_eq!(variant, reconstructed);
        }
    }

    #[test]
    fn bitpix_data_size_calculations() {
        let test_cases = [
            (BitPix::U8, 1000, 1000),
            (BitPix::I16, 1000, 2000),
            (BitPix::I32, 1000, 4000),
            (BitPix::I64, 1000, 8000),
            (BitPix::F32, 1000, 4000),
            (BitPix::F64, 1000, 8000),
        ];

        for (bitpix, num_pixels, expected_bytes) in test_cases {
            let actual_bytes = num_pixels * bitpix.bytes_per_pixel();
            assert_eq!(actual_bytes, expected_bytes);
        }
    }

    #[test]
    fn extreme_pixel_counts() {
        let large_pixel_count = usize::MAX / 8;

        for bitpix in [
            BitPix::U8,
            BitPix::I16,
            BitPix::I32,
            BitPix::I64,
            BitPix::F32,
            BitPix::F64,
        ] {
            let bytes_per_pixel = bitpix.bytes_per_pixel();

            if large_pixel_count <= usize::MAX / bytes_per_pixel {
                let _total_bytes = large_pixel_count * bytes_per_pixel;
            }
        }
    }
}
