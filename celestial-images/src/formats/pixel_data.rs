use crate::core::types::BitPix;

mod normalize;

#[derive(Debug, Clone)]
pub enum PixelData {
    U8(Vec<u8>),
    I8(Vec<i8>),
    U16(Vec<u16>),
    I16(Vec<i16>),
    U32(Vec<u32>),
    I32(Vec<i32>),
    U64(Vec<u64>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
}

// Runs `$body` with `$data` bound to the sample vector of whichever variant
// `$pixels` holds. Every per-variant function uses this, so adding a variant
// means adding one arm here.
macro_rules! with_samples {
    ($pixels:expr, $data:ident => $body:expr) => {
        match $pixels {
            $crate::formats::pixel_data::PixelData::U8($data) => $body,
            $crate::formats::pixel_data::PixelData::I8($data) => $body,
            $crate::formats::pixel_data::PixelData::U16($data) => $body,
            $crate::formats::pixel_data::PixelData::I16($data) => $body,
            $crate::formats::pixel_data::PixelData::U32($data) => $body,
            $crate::formats::pixel_data::PixelData::I32($data) => $body,
            $crate::formats::pixel_data::PixelData::U64($data) => $body,
            $crate::formats::pixel_data::PixelData::I64($data) => $body,
            $crate::formats::pixel_data::PixelData::F32($data) => $body,
            $crate::formats::pixel_data::PixelData::F64($data) => $body,
        }
    };
}
pub(crate) use with_samples;

// Like `with_samples!`, but wraps the result back in the same variant.
macro_rules! map_samples {
    ($pixels:expr, $data:ident => $body:expr) => {
        match $pixels {
            $crate::formats::pixel_data::PixelData::U8($data) => {
                $crate::formats::pixel_data::PixelData::U8($body)
            }
            $crate::formats::pixel_data::PixelData::I8($data) => {
                $crate::formats::pixel_data::PixelData::I8($body)
            }
            $crate::formats::pixel_data::PixelData::U16($data) => {
                $crate::formats::pixel_data::PixelData::U16($body)
            }
            $crate::formats::pixel_data::PixelData::I16($data) => {
                $crate::formats::pixel_data::PixelData::I16($body)
            }
            $crate::formats::pixel_data::PixelData::U32($data) => {
                $crate::formats::pixel_data::PixelData::U32($body)
            }
            $crate::formats::pixel_data::PixelData::I32($data) => {
                $crate::formats::pixel_data::PixelData::I32($body)
            }
            $crate::formats::pixel_data::PixelData::U64($data) => {
                $crate::formats::pixel_data::PixelData::U64($body)
            }
            $crate::formats::pixel_data::PixelData::I64($data) => {
                $crate::formats::pixel_data::PixelData::I64($body)
            }
            $crate::formats::pixel_data::PixelData::F32($data) => {
                $crate::formats::pixel_data::PixelData::F32($body)
            }
            $crate::formats::pixel_data::PixelData::F64($data) => {
                $crate::formats::pixel_data::PixelData::F64($body)
            }
        }
    };
}

pub(crate) trait Sample: Copy + Default {
    fn to_f64(self) -> f64;
    fn to_f32(self) -> f32;
    // None for a float sample.
    fn to_i128(self) -> Option<i128>;
}

macro_rules! impl_sample {
    ($to_i128:expr; $($t:ty),*) => {$(
        impl Sample for $t {
            fn to_f64(self) -> f64 {
                self as f64
            }
            fn to_f32(self) -> f32 {
                self as f32
            }
            fn to_i128(self) -> Option<i128> {
                $to_i128(self)
            }
        }
    )*};
}
impl_sample!(|v| Some(i128::from(v)); u8, i8, u16, i16, u32, i32, u64, i64);
impl_sample!(|_| None; f32, f64);

macro_rules! variant_access {
    ($($variant:ident, $t:ty, $as_ref:ident, $as_mut:ident;)*) => {$(
        impl PixelData {
            pub fn $as_ref(&self) -> Option<&Vec<$t>> {
                match self {
                    Self::$variant(v) => Some(v),
                    _ => None,
                }
            }

            pub fn $as_mut(&mut self) -> Option<&mut Vec<$t>> {
                match self {
                    Self::$variant(v) => Some(v),
                    _ => None,
                }
            }
        }

        impl From<&[$t]> for PixelData {
            fn from(data: &[$t]) -> Self {
                Self::$variant(data.to_vec())
            }
        }

        impl From<Vec<$t>> for PixelData {
            fn from(data: Vec<$t>) -> Self {
                Self::$variant(data)
            }
        }
    )*};
}
variant_access! {
    U8, u8, as_u8, as_u8_mut;
    I8, i8, as_i8, as_i8_mut;
    U16, u16, as_u16, as_u16_mut;
    I16, i16, as_i16, as_i16_mut;
    U32, u32, as_u32, as_u32_mut;
    I32, i32, as_i32, as_i32_mut;
    U64, u64, as_u64, as_u64_mut;
    I64, i64, as_i64, as_i64_mut;
    F32, f32, as_f32, as_f32_mut;
    F64, f64, as_f64, as_f64_mut;
}

impl PixelData {
    pub fn len(&self) -> usize {
        with_samples!(self, v => v.len())
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn bitpix(&self) -> BitPix {
        match self {
            Self::U8(_) | Self::I8(_) => BitPix::U8,
            Self::U16(_) | Self::I16(_) => BitPix::I16,
            Self::U32(_) | Self::I32(_) => BitPix::I32,
            Self::U64(_) | Self::I64(_) => BitPix::I64,
            Self::F32(_) => BitPix::F32,
            Self::F64(_) => BitPix::F64,
        }
    }

    pub fn to_f32(&self) -> Vec<f32> {
        with_samples!(self, v => v.iter().map(|&x| x.to_f32()).collect())
    }

    pub fn convert_to_f32(&mut self) {
        let converted = self.to_f32_normalized();
        *self = Self::F32(converted);
    }

    /// Keep the first `keep` of every `stored` interleaved samples and drop
    /// the rest. Loaders use this to strip an alpha (or other extra) channel
    /// so in-memory buffers are always a clean 1- or 3-channel image.
    ///
    /// Returns `self` unchanged when there is nothing to drop.
    pub fn keep_leading_samples(self, stored: usize, keep: usize) -> Self {
        // Compacting in place keeps one buffer alive. Collecting into a new
        // vector holds both, and the new one grows by reallocation.
        fn trim<T: Copy>(mut data: Vec<T>, stored: usize, keep: usize) -> Vec<T> {
            let mut kept = 0;
            for start in (0..data.len()).step_by(stored) {
                let end = data.len().min(start + keep);
                data.copy_within(start..end, kept);
                kept += end - start;
            }
            data.truncate(kept);
            data.shrink_to_fit();
            data
        }
        if keep == 0 || keep >= stored {
            return self;
        }
        map_samples!(self, d => trim(d, stored, keep))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn i8_and_u32_variants_match_the_other_integer_variants() {
        let i8s = PixelData::from(vec![-128i8, 0, 127]);
        assert_eq!(i8s.len(), 3);
        assert_eq!(i8s.bitpix(), BitPix::U8);
        assert_eq!(i8s.as_i8(), Some(&vec![-128, 0, 127]));
        assert!(i8s.as_u8().is_none());
        assert_eq!(i8s.to_f32(), vec![-128.0, 0.0, 127.0]);
        assert_eq!(i8s.to_f64_normalized(), vec![0.0, 128.0 / 255.0, 1.0]);
        assert_eq!(i8s.to_f32_normalized(), vec![0.0, 128.0 / 255.0, 1.0]);

        let u32s = PixelData::from(&[0u32, 7, u32::MAX][..]);
        assert_eq!(u32s.len(), 3);
        assert_eq!(u32s.bitpix(), BitPix::I32);
        assert_eq!(u32s.as_u32(), Some(&vec![0, 7, u32::MAX]));
        assert!(u32s.as_i32().is_none());
        assert_eq!(u32s.to_f64_normalized(), vec![0.0, 7.0 / 4294967295.0, 1.0]);
        assert_eq!(u32s.to_f32_normalized()[2], 1.0);

        let mut pd = PixelData::U32(vec![]);
        pd.as_u32_mut().unwrap().push(3000000000);
        assert_eq!(
            pd.keep_leading_samples(1, 1).as_u32(),
            Some(&vec![3000000000])
        );
    }

    #[test]
    fn keep_leading_samples_drops_trailing_channel() {
        // RGBA → RGB: keep the first 3 of every 4. One pixel: [R,G,B,A].
        let rgba = PixelData::U16(vec![10, 20, 30, 999, 40, 50, 60, 888]);
        let rgb = rgba.keep_leading_samples(4, 3);
        assert_eq!(rgb.as_u16().unwrap(), &vec![10, 20, 30, 40, 50, 60]);
    }

    #[test]
    fn keep_leading_samples_is_a_no_op_when_nothing_to_drop() {
        let gray = PixelData::U8(vec![1, 2, 3, 4]);
        assert_eq!(
            gray.clone().keep_leading_samples(1, 1).as_u8().unwrap(),
            &vec![1, 2, 3, 4]
        );
        assert_eq!(
            gray.clone().keep_leading_samples(3, 4).as_u8().unwrap(),
            &vec![1, 2, 3, 4]
        );
        assert_eq!(
            gray.keep_leading_samples(4, 0).as_u8().unwrap(),
            &vec![1, 2, 3, 4]
        );
    }

    #[test]
    fn len_and_is_empty_across_variants() {
        assert_eq!(PixelData::U8(vec![]).len(), 0);
        assert!(PixelData::U8(vec![]).is_empty());

        assert_eq!(PixelData::U16(vec![1, 2, 3]).len(), 3);
        assert!(!PixelData::U16(vec![1, 2, 3]).is_empty());

        assert_eq!(PixelData::I16(vec![1; 10]).len(), 10);
        assert_eq!(PixelData::I32(vec![1; 5]).len(), 5);
        assert_eq!(PixelData::F32(vec![1.0; 7]).len(), 7);
        assert_eq!(PixelData::F64(vec![1.0; 2]).len(), 2);
    }

    #[test]
    fn bitpix_maps_variants_correctly() {
        assert_eq!(PixelData::U8(vec![]).bitpix(), BitPix::U8);
        assert_eq!(PixelData::U16(vec![]).bitpix(), BitPix::I16);
        assert_eq!(PixelData::I16(vec![]).bitpix(), BitPix::I16);
        assert_eq!(PixelData::I32(vec![]).bitpix(), BitPix::I32);
        assert_eq!(PixelData::F32(vec![]).bitpix(), BitPix::F32);
        assert_eq!(PixelData::F64(vec![]).bitpix(), BitPix::F64);
    }

    #[test]
    fn as_ref_accessors_return_some_only_for_matching_variant() {
        let pd = PixelData::U8(vec![1, 2, 3]);
        assert_eq!(pd.as_u8(), Some(&vec![1, 2, 3]));
        assert!(pd.as_u16().is_none());
        assert!(pd.as_i16().is_none());
        assert!(pd.as_i32().is_none());
        assert!(pd.as_f32().is_none());
        assert!(pd.as_f64().is_none());

        let pd = PixelData::U16(vec![100, 200]);
        assert_eq!(pd.as_u16(), Some(&vec![100u16, 200]));
        assert!(pd.as_u8().is_none());

        let pd = PixelData::I16(vec![-5, 5]);
        assert_eq!(pd.as_i16(), Some(&vec![-5i16, 5]));

        let pd = PixelData::I32(vec![1000]);
        assert_eq!(pd.as_i32(), Some(&vec![1000i32]));

        let pd = PixelData::F32(vec![1.5]);
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let pd = PixelData::F64(vec![2.5]);
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }

    #[test]
    fn as_mut_accessors_allow_modification() {
        let mut pd = PixelData::U8(vec![1, 2, 3]);
        pd.as_u8_mut().unwrap().push(4);
        assert_eq!(pd.as_u8(), Some(&vec![1u8, 2, 3, 4]));
        assert!(pd.as_u16_mut().is_none());

        let mut pd = PixelData::U16(vec![1]);
        pd.as_u16_mut().unwrap().push(2);
        assert_eq!(pd.as_u16(), Some(&vec![1u16, 2]));
        assert!(pd.as_u8_mut().is_none());

        let mut pd = PixelData::I16(vec![]);
        pd.as_i16_mut().unwrap().push(-1);
        assert_eq!(pd.as_i16(), Some(&vec![-1i16]));

        let mut pd = PixelData::I32(vec![]);
        pd.as_i32_mut().unwrap().push(100);
        assert_eq!(pd.as_i32(), Some(&vec![100i32]));

        let mut pd = PixelData::F32(vec![]);
        pd.as_f32_mut().unwrap().push(1.5);
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let mut pd = PixelData::F64(vec![]);
        pd.as_f64_mut().unwrap().push(2.5);
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }

    #[test]
    fn to_f32_from_every_variant() {
        assert_eq!(
            PixelData::U8(vec![0, 128, 255]).to_f32(),
            vec![0.0, 128.0, 255.0]
        );
        assert_eq!(PixelData::U16(vec![0, 65535]).to_f32(), vec![0.0, 65535.0]);
        assert_eq!(PixelData::I16(vec![-1, 1]).to_f32(), vec![-1.0, 1.0]);
        assert_eq!(PixelData::I32(vec![-42, 42]).to_f32(), vec![-42.0, 42.0]);
        assert_eq!(PixelData::F32(vec![1.5, 2.5]).to_f32(), vec![1.5, 2.5]);
        assert_eq!(PixelData::F64(vec![1.5, 2.5]).to_f32(), vec![1.5f32, 2.5]);
    }

    #[test]
    fn convert_to_f32_replaces_variant() {
        let mut pd = PixelData::U8(vec![0, 255]);
        pd.convert_to_f32();
        assert!(matches!(pd, PixelData::F32(_)));
        assert_eq!(pd.as_f32().unwrap(), &vec![0.0f32, 1.0]);

        let mut pd = PixelData::F64(vec![0.25]);
        pd.convert_to_f32();
        assert!(matches!(pd, PixelData::F32(_)));
    }

    #[test]
    fn from_slice_conversions() {
        let pd: PixelData = (&[1u8, 2, 3][..]).into();
        assert_eq!(pd.as_u8(), Some(&vec![1u8, 2, 3]));

        let pd: PixelData = (&[100u16, 200][..]).into();
        assert_eq!(pd.as_u16(), Some(&vec![100u16, 200]));

        let pd: PixelData = (&[-1i16, 1][..]).into();
        assert_eq!(pd.as_i16(), Some(&vec![-1i16, 1]));

        let pd: PixelData = (&[42i32][..]).into();
        assert_eq!(pd.as_i32(), Some(&vec![42i32]));

        let pd: PixelData = (&[1.5f32][..]).into();
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let pd: PixelData = (&[2.5f64][..]).into();
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }

    #[test]
    fn from_vec_conversions() {
        let pd: PixelData = vec![1u8, 2].into();
        assert_eq!(pd.as_u8(), Some(&vec![1u8, 2]));

        let pd: PixelData = vec![100u16].into();
        assert_eq!(pd.as_u16(), Some(&vec![100u16]));

        let pd: PixelData = vec![-1i16].into();
        assert_eq!(pd.as_i16(), Some(&vec![-1i16]));

        let pd: PixelData = vec![42i32].into();
        assert_eq!(pd.as_i32(), Some(&vec![42i32]));

        let pd: PixelData = vec![1.5f32].into();
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let pd: PixelData = vec![2.5f64].into();
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }
}
