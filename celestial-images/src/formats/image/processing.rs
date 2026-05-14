use super::Image;
use crate::debayer::{debayer_bilinear_u16, debayer_bilinear_u8, BayerPattern};
use crate::formats::pixel_data::PixelData;

impl Image {
    pub fn interleaved_to_planar(&mut self) {
        if self.channels() != 3 {
            return;
        }
        let pixel_count = self.width() * self.height();

        macro_rules! convert {
            ($data:expr) => {{
                let mut planar = vec![Default::default(); $data.len()];
                for i in 0..pixel_count {
                    planar[i] = $data[i * 3];
                    planar[pixel_count + i] = $data[i * 3 + 1];
                    planar[pixel_count * 2 + i] = $data[i * 3 + 2];
                }
                *$data = planar;
            }};
        }

        match &mut self.pixels {
            PixelData::U8(data) => convert!(data),
            PixelData::U16(data) => convert!(data),
            PixelData::I16(data) => convert!(data),
            PixelData::I32(data) => convert!(data),
            PixelData::F32(data) => convert!(data),
            PixelData::F64(data) => convert!(data),
        }
    }

    pub fn planar_to_interleaved(&mut self) {
        if self.channels() != 3 {
            return;
        }
        let pixel_count = self.width() * self.height();

        macro_rules! convert {
            ($data:expr) => {{
                let mut interleaved = vec![Default::default(); $data.len()];
                for i in 0..pixel_count {
                    interleaved[i * 3] = $data[i];
                    interleaved[i * 3 + 1] = $data[pixel_count + i];
                    interleaved[i * 3 + 2] = $data[pixel_count * 2 + i];
                }
                *$data = interleaved;
            }};
        }

        match &mut self.pixels {
            PixelData::U8(data) => convert!(data),
            PixelData::U16(data) => convert!(data),
            PixelData::I16(data) => convert!(data),
            PixelData::I32(data) => convert!(data),
            PixelData::F32(data) => convert!(data),
            PixelData::F64(data) => convert!(data),
        }
    }

    pub fn normalize(&mut self) {
        let (min, max) = self.pixel_range();
        let range = (max - min).max(f64::MIN_POSITIVE);

        macro_rules! normalize {
            ($data:expr, $t:ty) => {{
                for v in $data.iter_mut() {
                    *v = ((*v as f64 - min) / range) as $t;
                }
            }};
        }

        match &mut self.pixels {
            PixelData::U8(data) => normalize!(data, u8),
            PixelData::U16(data) => normalize!(data, u16),
            PixelData::I16(data) => normalize!(data, i16),
            PixelData::I32(data) => normalize!(data, i32),
            PixelData::F32(data) => normalize!(data, f32),
            PixelData::F64(data) => normalize!(data, f64),
        }
    }

    pub fn normalize_to_f32(&mut self) {
        let (min, max) = self.pixel_range();
        let range = (max - min).max(f64::MIN_POSITIVE);

        let normalized: Vec<f32> = match &self.pixels {
            PixelData::U8(d) => d
                .iter()
                .map(|&v| ((v as f64 - min) / range) as f32)
                .collect(),
            PixelData::U16(d) => d
                .iter()
                .map(|&v| ((v as f64 - min) / range) as f32)
                .collect(),
            PixelData::I16(d) => d
                .iter()
                .map(|&v| ((v as f64 - min) / range) as f32)
                .collect(),
            PixelData::I32(d) => d
                .iter()
                .map(|&v| ((v as f64 - min) / range) as f32)
                .collect(),
            PixelData::F32(d) => d
                .iter()
                .map(|&v| ((v as f64 - min) / range) as f32)
                .collect(),
            PixelData::F64(d) => d.iter().map(|&v| ((v - min) / range) as f32).collect(),
        };
        self.pixels = PixelData::F32(normalized);
    }

    fn pixel_range(&self) -> (f64, f64) {
        match &self.pixels {
            PixelData::U8(d) => {
                let min = d.iter().copied().min().unwrap_or(0) as f64;
                let max = d.iter().copied().max().unwrap_or(0) as f64;
                (min, max)
            }
            PixelData::U16(d) => {
                let min = d.iter().copied().min().unwrap_or(0) as f64;
                let max = d.iter().copied().max().unwrap_or(0) as f64;
                (min, max)
            }
            PixelData::I16(d) => {
                let min = d.iter().copied().min().unwrap_or(0) as f64;
                let max = d.iter().copied().max().unwrap_or(0) as f64;
                (min, max)
            }
            PixelData::I32(d) => {
                let min = d.iter().copied().min().unwrap_or(0) as f64;
                let max = d.iter().copied().max().unwrap_or(0) as f64;
                (min, max)
            }
            PixelData::F32(d) => {
                let min = d.iter().copied().fold(f32::INFINITY, f32::min) as f64;
                let max = d.iter().copied().fold(f32::NEG_INFINITY, f32::max) as f64;
                (min, max)
            }
            PixelData::F64(d) => {
                let min = d.iter().copied().fold(f64::INFINITY, f64::min);
                let max = d.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                (min, max)
            }
        }
    }

    pub fn debayer(&mut self, pattern: BayerPattern) {
        if self.channels() != 1 {
            return;
        }

        let width = self.width();
        let height = self.height();

        match &self.pixels {
            PixelData::U8(data) => {
                let rgb = debayer_bilinear_u8(data, width, height, pattern);
                self.pixels = PixelData::U8(rgb);
            }
            PixelData::U16(data) => {
                let rgb = debayer_bilinear_u16(data, width, height, pattern);
                self.pixels = PixelData::U16(rgb);
            }
            _ => return,
        }

        self.dimensions = vec![width, height, 3];
    }
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
        assert_eq!(img.pixels.as_u8().unwrap(), original.pixels.as_u8().unwrap());
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
            PixelData::U16(vec![1, 2, 3, 4, 5, 6]),
            PixelData::I16(vec![1, 2, 3, 4, 5, 6]),
            PixelData::I32(vec![1, 2, 3, 4, 5, 6]),
            PixelData::F32(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            PixelData::F64(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        ] {
            let mut img = Image::new(pd, shapes.clone());
            img.interleaved_to_planar();
            img.planar_to_interleaved();
            // After roundtrip the data should be unchanged.
            match &img.pixels {
                PixelData::U8(v) => assert_eq!(v, &vec![1u8, 2, 3, 4, 5, 6]),
                PixelData::U16(v) => assert_eq!(v, &vec![1u16, 2, 3, 4, 5, 6]),
                PixelData::I16(v) => assert_eq!(v, &vec![1i16, 2, 3, 4, 5, 6]),
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
        let data = img.pixels.as_u8().unwrap();
        assert_eq!(data[0], 0);
        assert_eq!(data[3], 1);
    }

    #[test]
    fn normalize_scales_f32_into_unit_range() {
        let mut img = Image::new(
            PixelData::F32(vec![10.0, 20.0, 30.0, 40.0]),
            vec![2usize, 2],
        );
        img.normalize();
        let data = img.pixels.as_f32().unwrap();
        assert!((data[0] - 0.0).abs() < 1e-6);
        assert!((data[3] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn normalize_to_f32_replaces_variant() {
        let mut img = Image::new(
            PixelData::U16(vec![0, 32768, 65535, 16384]),
            vec![2usize, 2],
        );
        img.normalize_to_f32();
        let data = img.pixels.as_f32().unwrap();
        assert_eq!(data[0], 0.0);
        assert!((data[2] - 1.0).abs() < 1e-5);
        assert!(matches!(img.pixels, PixelData::F32(_)));
    }

    #[test]
    fn normalize_handles_constant_image_without_division_by_zero() {
        let mut img = Image::new(PixelData::U8(vec![7; 4]), vec![2usize, 2]);
        img.normalize();
        // All values equal, so min == max; division-by-zero is guarded by MIN_POSITIVE.
        // Output is undefined-but-finite — just ensure no panic occurred.
        let _data = img.pixels.as_u8().unwrap();
    }

    #[test]
    fn normalize_covers_every_variant() {
        for pd in [
            PixelData::U8(vec![0, 255]),
            PixelData::U16(vec![0, 65535]),
            PixelData::I16(vec![-100, 100]),
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
            PixelData::U16(vec![0, 65535]),
            PixelData::I16(vec![-100, 100]),
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
        img.debayer(BayerPattern::Rggb);
        assert_eq!(img.pixels.as_u8().unwrap(), &before);
        assert_eq!(img.dimensions, vec![2usize, 1, 3]);
    }

    #[test]
    fn debayer_expands_u8_mono_to_rgb() {
        let mut img = Image::new(PixelData::U8(vec![0; 16]), vec![4usize, 4]);
        img.debayer(BayerPattern::Rggb);
        assert_eq!(img.channels(), 3);
        assert_eq!(img.dimensions, vec![4usize, 4, 3]);
        assert_eq!(img.pixels.as_u8().unwrap().len(), 4 * 4 * 3);
    }

    #[test]
    fn debayer_expands_u16_mono_to_rgb() {
        let mut img = Image::new(PixelData::U16(vec![0; 16]), vec![4usize, 4]);
        img.debayer(BayerPattern::Rggb);
        assert_eq!(img.channels(), 3);
        assert_eq!(img.pixels.as_u16().unwrap().len(), 4 * 4 * 3);
    }

    #[test]
    fn debayer_skips_unsupported_bit_depths() {
        let mut img = Image::new(
            PixelData::F32(vec![0.0; 16]),
            vec![4usize, 4],
        );
        img.debayer(BayerPattern::Rggb);
        // Because we return early, dimensions shouldn't have gotten the third axis.
        assert_eq!(img.channels(), 1);
    }
}
