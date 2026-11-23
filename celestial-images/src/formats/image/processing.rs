use super::Image;
use crate::core::errors::Result;
use crate::core::sample_order::{to_interleaved, to_planar};
use crate::debayer::{debayer_bilinear_u16, debayer_bilinear_u8, BayerPattern};
use crate::formats::pixel_data::{with_samples, PixelData, Sample};

impl Image {
    // Only a three-axis image with three planes is RGB. Images with more axes
    // keep their planes in stored order.
    pub(super) fn has_rgb_planes(&self) -> bool {
        self.dimensions.len() == 3 && self.channels() == 3
    }

    pub fn interleaved_to_planar(&mut self) {
        if !self.has_rgb_planes() {
            return;
        }
        with_samples!(&mut self.pixels, data => *data = to_planar(data, 3));
    }

    pub fn planar_to_interleaved(&mut self) {
        if !self.has_rgb_planes() {
            return;
        }
        with_samples!(&mut self.pixels, data => *data = to_interleaved(data, 3));
    }

    // F32 and F64 keep their type. Integer variants become F64, because the
    // scaled values lie in [0, 1] and an integer type would round them to 0 or 1.
    pub fn normalize(&mut self) {
        self.blank_to_nan();
        let (min, range) = self.min_and_range();
        let scale = |v: f64| (v - min) / range;
        match &mut self.pixels {
            PixelData::F32(data) => {
                for v in data.iter_mut() {
                    *v = scale(*v as f64) as f32;
                }
            }
            PixelData::F64(data) => {
                for v in data.iter_mut() {
                    *v = scale(*v);
                }
            }
            other => {
                let scaled =
                    with_samples!(other, data => data.iter().map(|v| scale(v.to_f64())).collect());
                self.pixels = PixelData::F64(scaled);
            }
        }
    }

    pub fn normalize_to_f32(&mut self) {
        self.blank_to_nan();
        let (min, range) = self.min_and_range();
        let normalized: Vec<f32> = with_samples!(&self.pixels, data => data
            .iter()
            .map(|v| ((v.to_f64() - min) / range) as f32)
            .collect());
        self.pixels = PixelData::F32(normalized);
    }

    // A BLANK pixel has no measured value, and a float array marks that with
    // NaN. The card is dropped because a float array has no use for BLANK.
    fn blank_to_nan(&mut self) {
        let Some(flags) = self.blank_flags() else {
            return;
        };
        let masked = with_samples!(&self.pixels, data => data
            .iter()
            .zip(&flags)
            .map(|(v, &blank)| if blank { f64::NAN } else { v.to_f64() })
            .collect());
        self.pixels = PixelData::F64(masked);
        self.remove_keyword("BLANK");
    }

    // One flag per sample, set where an integer sample equals BLANK. None
    // when the image has no BLANK card or holds float samples.
    pub(crate) fn blank_flags(&self) -> Option<Vec<bool>> {
        let blank = i128::from(self.get_last("BLANK").as_i64()?);
        with_samples!(&self.pixels, data => equals_blank(data, blank))
    }

    // NaN samples are skipped. A constant image gets the smallest positive
    // range, so it scales to 0 instead of dividing by zero.
    fn min_and_range(&self) -> (f64, f64) {
        let (min, max) = with_samples!(&self.pixels, data => data.iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(lo, hi), v| (libm::fmin(lo, v.to_f64()), libm::fmax(hi, v.to_f64())),
        ));
        (min, libm::fmax(max - min, f64::MIN_POSITIVE))
    }

    // Images with more than one channel, and sample types other than U8 and
    // U16, are left unchanged.
    pub fn debayer(&mut self, pattern: BayerPattern) -> Result<()> {
        if self.channels() != 1 {
            return Ok(());
        }
        let width = self.width();
        let height = self.height();
        self.pixels = match &self.pixels {
            PixelData::U8(data) => {
                PixelData::U8(debayer_bilinear_u8(data, width, height, pattern)?)
            }
            PixelData::U16(data) => {
                PixelData::U16(debayer_bilinear_u16(data, width, height, pattern)?)
            }
            _ => return Ok(()),
        };
        self.dimensions = vec![width, height, 3];
        Ok(())
    }
}

// The comparison runs in i128 so that 64-bit values above 2^53 match exactly.
// Float types have no integer form, so they return None before any sample
// is read, even when `data` is empty.
fn equals_blank<T: Sample>(data: &[T], blank: i128) -> Option<Vec<bool>> {
    T::default().to_i128()?;
    Some(data.iter().map(|v| v.to_i128() == Some(blank)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rgb_interleaved_u8() -> Image {
        // Two pixels, RGB values (1,2,3) and (4,5,6) — 2x1 RGB.
        Image::new(PixelData::U8(vec![1, 2, 3, 4, 5, 6]), vec![2usize, 1, 3])
    }

    fn rgb_planar_u8() -> Image {
        // Two pixels in planar: R=[1,4], G=[2,5], B=[3,6]
        Image::new(PixelData::U8(vec![1, 4, 2, 5, 3, 6]), vec![2usize, 1, 3])
    }

    #[test]
    fn interleaved_to_planar_rearranges_channels() {
        let mut img = rgb_interleaved_u8();
        img.interleaved_to_planar();
        assert_eq!(img.pixels.as_u8().unwrap(), &vec![1u8, 4, 2, 5, 3, 6]);
    }

    #[test]
    fn planar_to_interleaved_rearranges_channels() {
        let mut img = rgb_planar_u8();
        img.planar_to_interleaved();
        assert_eq!(img.pixels.as_u8().unwrap(), &vec![1u8, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn interleaved_planar_roundtrip() {
        let original = rgb_interleaved_u8();
        let mut img = original.clone();
        img.interleaved_to_planar();
        img.planar_to_interleaved();
        assert_eq!(
            img.pixels.as_u8().unwrap(),
            original.pixels.as_u8().unwrap()
        );
    }

    #[test]
    fn interleaved_to_planar_skips_non_rgb() {
        let mut mono = Image::new(PixelData::U8(vec![1, 2, 3, 4]), vec![2usize, 2]);
        let before = mono.pixels.as_u8().unwrap().clone();
        mono.interleaved_to_planar();
        assert_eq!(mono.pixels.as_u8().unwrap(), &before);
    }

    #[test]
    fn planar_to_interleaved_skips_non_rgb() {
        let mut mono = Image::new(PixelData::U8(vec![1, 2, 3, 4]), vec![2usize, 2]);
        let before = mono.pixels.as_u8().unwrap().clone();
        mono.planar_to_interleaved();
        assert_eq!(mono.pixels.as_u8().unwrap(), &before);
    }

    #[test]
    fn interleaved_to_planar_handles_all_variants() {
        let shapes = vec![2usize, 1, 3];
        for pd in [
            PixelData::U8(vec![1, 2, 3, 4, 5, 6]),
            PixelData::I8(vec![1, 2, 3, 4, 5, 6]),
            PixelData::U16(vec![1, 2, 3, 4, 5, 6]),
            PixelData::I16(vec![1, 2, 3, 4, 5, 6]),
            PixelData::U32(vec![1, 2, 3, 4, 5, 6]),
            PixelData::I32(vec![1, 2, 3, 4, 5, 6]),
            PixelData::U64(vec![1, 2, 3, 4, 5, 6]),
            PixelData::I64(vec![1, 2, 3, 4, 5, 6]),
            PixelData::F32(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            PixelData::F64(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        ] {
            let mut img = Image::new(pd, shapes.clone());
            img.interleaved_to_planar();
            img.planar_to_interleaved();
            // After roundtrip the data should be unchanged.
            match &img.pixels {
                PixelData::U8(v) => assert_eq!(v, &vec![1u8, 2, 3, 4, 5, 6]),
                PixelData::I8(v) => assert_eq!(v, &vec![1i8, 2, 3, 4, 5, 6]),
                PixelData::U16(v) => assert_eq!(v, &vec![1u16, 2, 3, 4, 5, 6]),
                PixelData::I16(v) => assert_eq!(v, &vec![1i16, 2, 3, 4, 5, 6]),
                PixelData::U32(v) => assert_eq!(v, &vec![1u32, 2, 3, 4, 5, 6]),
                PixelData::U64(v) => assert_eq!(v, &vec![1u64, 2, 3, 4, 5, 6]),
                PixelData::I64(v) => assert_eq!(v, &vec![1i64, 2, 3, 4, 5, 6]),
                PixelData::I32(v) => assert_eq!(v, &vec![1i32, 2, 3, 4, 5, 6]),
                PixelData::F32(v) => assert_eq!(v, &vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0]),
                PixelData::F64(v) => assert_eq!(v, &vec![1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0]),
            }
        }
    }

    #[test]
    fn normalize_scales_u8_to_full_range() {
        let mut img = Image::new(PixelData::U8(vec![50, 100, 150, 200]), vec![2usize, 2]);
        img.normalize();
        assert_eq!(
            img.pixels.as_f64(),
            Some(&vec![0.0, 50.0 / 150.0, 100.0 / 150.0, 1.0])
        );
    }

    #[test]
    fn normalize_converts_integer_data_to_f64() {
        let mut img = Image::new(PixelData::U16(vec![0, 1000, 30000, 65535]), vec![4usize]);
        img.normalize();
        assert_eq!(
            img.pixels.as_f64(),
            Some(&vec![0.0, 1000.0 / 65535.0, 30000.0 / 65535.0, 1.0])
        );
    }

    fn with_blank(pixels: PixelData, blank: i64) -> Image {
        let mut img = Image::new(pixels, vec![4usize]);
        img.set_keyword(crate::fits::header::keywords::Keyword::integer(
            "BLANK", blank,
        ));
        img
    }

    #[test]
    fn blank_flags_use_the_last_blank_card() {
        let mut img = with_blank(PixelData::I16(vec![5, 6, 7, 8]), 5);
        img.keywords
            .push(crate::fits::header::keywords::Keyword::integer("BLANK", 6));
        assert_eq!(img.blank_flags(), Some(vec![false, true, false, false]));
    }

    #[test]
    fn normalize_turns_blank_pixels_into_nan() {
        let mut img = with_blank(PixelData::I16(vec![0, i16::MIN, 100, 50]), -32768);
        img.normalize();
        let bits: Vec<u64> = img
            .pixels
            .as_f64()
            .unwrap()
            .iter()
            .map(|v| v.to_bits())
            .collect();
        assert_eq!(bits, [0.0, f64::NAN, 1.0, 0.5].map(f64::to_bits));
        assert!(img.get_keyword("BLANK").is_none());
    }

    #[test]
    fn normalize_to_f32_turns_blank_pixels_into_nan() {
        let mut img = with_blank(PixelData::U16(vec![0, 7, 100, 50]), 7);
        img.normalize_to_f32();
        let bits: Vec<u32> = img
            .pixels
            .as_f32()
            .unwrap()
            .iter()
            .map(|v| v.to_bits())
            .collect();
        assert_eq!(bits, [0.0, f32::NAN, 1.0, 0.5].map(f32::to_bits));
        assert!(img.get_keyword("BLANK").is_none());
    }

    #[test]
    fn normalize_scales_f32_into_unit_range() {
        let mut img = Image::new(
            PixelData::F32(vec![10.0, 20.0, 30.0, 40.0]),
            vec![2usize, 2],
        );
        img.normalize();
        // min = 10 and range = 30, so the quotients are 0, 1/3, 2/3 and 1,
        // each rounded to f64 and then to f32. 1/3 is 1.0101...b x 2^-2 and
        // 2/3 is the same mantissa x 2^-1. The bits after the 23rd are
        // 1010..., so both roundings go up to mantissa 0x2AAAAB: 0x3EAAAAAB
        // and 0x3F2AAAAB. The f64 rounding is far from the f32 halfway
        // point, so it does not change the f32 result.
        let third = f32::from_bits(0x3EAA_AAAB);
        let two_thirds = f32::from_bits(0x3F2A_AAAB);
        assert_eq!(
            img.pixels.as_f32().unwrap(),
            &vec![0.0, third, two_thirds, 1.0]
        );
    }

    #[test]
    fn normalize_to_f32_replaces_variant() {
        let mut img = Image::new(
            PixelData::U16(vec![0, 32768, 65535, 16384]),
            vec![2usize, 2],
        );
        img.normalize_to_f32();
        assert!(matches!(img.pixels, PixelData::F32(_)));
        // min = 0 and range = 65535 = 2^16 - 1. 2^15 / (2^16 - 1) is
        // 2^-1 * (1 + 2^-16 + 2^-32 + ...), so its mantissa has ones at bits
        // 16, 32, 48 and so on. Rounding to f64 (52 bits) and then to f32
        // (23 bits) both drop bits 32 and up, which are less than half a
        // unit, leaving 0.5 + 2^-17. 2^14 / (2^16 - 1) gives 0.25 + 2^-18
        // the same way.
        let expected = vec![0.0f32, 0.5 + 1.0 / 131072.0, 1.0, 0.25 + 1.0 / 262144.0];
        assert_eq!(img.pixels.as_f32().unwrap(), &expected);
    }

    #[test]
    fn normalize_handles_constant_image_without_division_by_zero() {
        let mut img = Image::new(PixelData::U8(vec![7; 4]), vec![2usize, 2]);
        img.normalize();
        // All values equal, so min == max; division-by-zero is guarded by MIN_POSITIVE.
        assert_eq!(img.pixels.as_f64(), Some(&vec![0.0; 4]));
    }

    #[test]
    fn normalize_covers_every_variant() {
        for pd in [
            PixelData::U8(vec![0, 255]),
            PixelData::I8(vec![-128, 127]),
            PixelData::U16(vec![0, 65535]),
            PixelData::I16(vec![-100, 100]),
            PixelData::U32(vec![0, u32::MAX]),
            PixelData::U64(vec![0, u64::MAX]),
            PixelData::I64(vec![i64::MIN, i64::MAX]),
            PixelData::I32(vec![-1000, 1000]),
            PixelData::F32(vec![1.0, 10.0]),
            PixelData::F64(vec![-5.0, 5.0]),
        ] {
            let mut img = Image::new(pd, vec![2usize]);
            img.normalize();
        }
    }

    #[test]
    fn normalize_to_f32_covers_every_variant() {
        for pd in [
            PixelData::U8(vec![0, 255]),
            PixelData::I8(vec![-128, 127]),
            PixelData::U16(vec![0, 65535]),
            PixelData::I16(vec![-100, 100]),
            PixelData::U32(vec![0, u32::MAX]),
            PixelData::U64(vec![0, u64::MAX]),
            PixelData::I64(vec![i64::MIN, i64::MAX]),
            PixelData::I32(vec![-1000, 1000]),
            PixelData::F32(vec![1.0, 10.0]),
            PixelData::F64(vec![-5.0, 5.0]),
        ] {
            let mut img = Image::new(pd, vec![2usize]);
            img.normalize_to_f32();
            assert!(matches!(img.pixels, PixelData::F32(_)));
        }
    }

    #[test]
    fn debayer_skips_multi_channel_images() {
        let mut img = rgb_interleaved_u8();
        let before = img.pixels.as_u8().unwrap().clone();
        img.debayer(BayerPattern::Rggb).unwrap();
        assert_eq!(img.pixels.as_u8().unwrap(), &before);
        assert_eq!(img.dimensions, vec![2usize, 1, 3]);
    }

    #[test]
    fn debayer_expands_u8_mono_to_rgb() {
        let mut img = Image::new(PixelData::U8(vec![0; 16]), vec![4usize, 4]);
        img.debayer(BayerPattern::Rggb).unwrap();
        assert_eq!(img.channels(), 3);
        assert_eq!(img.dimensions, vec![4usize, 4, 3]);
        assert_eq!(img.pixels.as_u8().unwrap().len(), 4 * 4 * 3);
    }

    #[test]
    fn debayer_expands_u16_mono_to_rgb() {
        let mut img = Image::new(PixelData::U16(vec![0; 16]), vec![4usize, 4]);
        img.debayer(BayerPattern::Rggb).unwrap();
        assert_eq!(img.channels(), 3);
        assert_eq!(img.pixels.as_u16().unwrap().len(), 4 * 4 * 3);
    }

    #[test]
    fn debayer_rejects_pixels_that_do_not_fill_the_dimensions() {
        let mut img = Image::new(PixelData::U16(vec![0; 15]), vec![4usize, 4]);
        let result = img.debayer(BayerPattern::Rggb);
        assert!(matches!(
            result,
            Err(crate::core::errors::ImageError::InvalidDimensions(_))
        ));
        assert_eq!(img.dimensions, vec![4usize, 4]);
    }

    #[test]
    fn debayer_skips_unsupported_bit_depths() {
        let mut img = Image::new(PixelData::F32(vec![0.0; 16]), vec![4usize, 4]);
        img.debayer(BayerPattern::Rggb).unwrap();
        // Because we return early, dimensions shouldn't have gotten the third axis.
        assert_eq!(img.channels(), 1);
    }
}
