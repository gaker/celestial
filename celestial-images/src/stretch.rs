use std::cmp::Ordering;

use crate::core::{ImageError, Result};
use crate::formats::{Image, PixelData};

pub struct StfParams {
    pub shadow_clip: f64,
    pub target_background: f64,
}

impl Default for StfParams {
    fn default() -> Self {
        Self {
            shadow_clip: 2.8,
            target_background: 0.25,
        }
    }
}

pub trait Stretch {
    fn stf(&self, params: &StfParams) -> Result<Image>;
    fn mad_normalized(&self) -> Result<Image>;
    fn to_u8(&self) -> Result<Image>;
}

impl Stretch for Image {
    fn stf(&self, params: &StfParams) -> Result<Image> {
        let data = self.pixels.to_f64_normalized();
        if data.is_empty() {
            return Err(ImageError::UnsupportedFormat);
        }
        let channels = self.channels();
        let stretched = apply_per_channel_f64(&data, channels, |ch| {
            apply_stf(ch, params)
        });
        let mut img = Image::new(PixelData::F64(stretched), self.dimensions.clone());
        img.keywords = self.keywords.clone();
        Ok(img)
    }

    fn mad_normalized(&self) -> Result<Image> {
        let mut data = self.pixels.to_f32_normalized();
        if data.is_empty() {
            return Err(ImageError::UnsupportedFormat);
        }
        let channels = self.channels();
        apply_per_channel_f32_mut(&mut data, channels, mad_normalize);
        let mut img = Image::new(PixelData::F32(data), self.dimensions.clone());
        img.keywords = self.keywords.clone();
        Ok(img)
    }

    fn to_u8(&self) -> Result<Image> {
        let data = self.pixels.to_f64_normalized();
        let bytes = stf_to_u8(&data);
        let mut img = Image::new(PixelData::U8(bytes), self.dimensions.clone());
        img.keywords = self.keywords.clone();
        Ok(img)
    }
}

fn apply_per_channel_f64(
    data: &[f64],
    channels: usize,
    f: impl Fn(&[f64]) -> Vec<f64>,
) -> Vec<f64> {
    if channels <= 1 {
        return f(data);
    }
    let pixel_count = data.len() / channels;
    let mut channel_bufs: Vec<Vec<f64>> = (0..channels)
        .map(|ch| {
            (0..pixel_count)
                .map(|i| data[i * channels + ch])
                .collect()
        })
        .collect();
    for buf in &mut channel_bufs {
        *buf = f(buf);
    }
    let mut result = vec![0.0; data.len()];
    for ch in 0..channels {
        for i in 0..pixel_count {
            result[i * channels + ch] = channel_bufs[ch][i];
        }
    }
    result
}

fn apply_per_channel_f32_mut(
    data: &mut [f32],
    channels: usize,
    f: impl Fn(&mut [f32]),
) {
    if channels <= 1 {
        f(data);
        return;
    }
    let pixel_count = data.len() / channels;
    let mut channel_buf: Vec<f32> = vec![0.0; pixel_count];
    for ch in 0..channels {
        for i in 0..pixel_count {
            channel_buf[i] = data[i * channels + ch];
        }
        f(&mut channel_buf);
        for i in 0..pixel_count {
            data[i * channels + ch] = channel_buf[i];
        }
    }
}

struct StfStatistics {
    median: f64,
    normalized_mad: f64,
}

fn compute_stf_statistics(image: &[f64]) -> Option<StfStatistics> {
    if image.is_empty() {
        return None;
    }
    let mut sorted: Vec<f64> = image.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let median = sorted[sorted.len() / 2];
    let mut deviations: Vec<f64> = image.iter().map(|&x| libm::fabs(x - median)).collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mad = deviations[deviations.len() / 2];
    Some(StfStatistics {
        median,
        normalized_mad: mad * 1.4826,
    })
}

fn mtf(x: f64, m: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    if m <= 0.0 || m >= 1.0 {
        return x;
    }
    let numerator = (m - 1.0) * x;
    let denominator = (2.0 * m - 1.0) * x - m;
    numerator / denominator
}

fn compute_midtone_balance(target_bg: f64, median_normalized: f64) -> f64 {
    if median_normalized <= 0.0 || median_normalized >= 1.0 {
        return 0.5;
    }
    let m_num = (target_bg - 1.0) * median_normalized;
    let m_den = (2.0 * target_bg - 1.0) * median_normalized - target_bg;
    if libm::fabs(m_den) < 1e-15 {
        return 0.5;
    }
    let m = m_num / m_den;
    m.clamp(0.0001, 0.9999)
}

fn apply_stf(image: &[f64], params: &StfParams) -> Vec<f64> {
    let stats = match compute_stf_statistics(image) {
        Some(s) => s,
        None => return Vec::new(),
    };
    let shadow_clip = (stats.median - params.shadow_clip * stats.normalized_mad).max(0.0);
    let range = 1.0 - shadow_clip;
    if range <= 0.0 {
        return vec![0.5; image.len()];
    }
    let median_normalized = (stats.median - shadow_clip) / range;
    let midtone_balance = compute_midtone_balance(params.target_background, median_normalized);
    image
        .iter()
        .map(|&px| {
            let clamped = ((px - shadow_clip) / range).clamp(0.0, 1.0);
            mtf(clamped, midtone_balance)
        })
        .collect()
}

fn stf_to_u8(stretched: &[f64]) -> Vec<u8> {
    stretched
        .iter()
        .map(|&v| libm::round(v.clamp(0.0, 1.0) * 255.0) as u8)
        .collect()
}

pub fn mad_normalize(image: &mut [f32]) {
    let len = image.len();
    if len == 0 {
        return;
    }
    let mut sorted: Vec<f32> = image.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let median = sorted[len / 2];
    let mut deviations: Vec<f32> = image
        .iter()
        .map(|&x| libm::fabsf(x - median))
        .collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mad = deviations[len / 2];
    let sigma = (mad * 1.4826).max(1e-10);
    for pixel in image.iter_mut() {
        *pixel = ((*pixel - median) / sigma).clamp(-5.0, 100.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::Keyword;

    #[test]
    fn background_centered_near_zero() {
        let mut image = vec![100.0f32; 256 * 256];
        for (i, px) in image.iter_mut().enumerate() {
            *px += (i % 7) as f32 * 0.5;
        }
        image[1000] = 5000.0;
        image[2000] = 10000.0;

        mad_normalize(&mut image);

        let bg_mean: f32 = image[0..100].iter().sum::<f32>() / 100.0;
        assert!(
            bg_mean.abs() < 2.0,
            "background mean should be near 0, got {bg_mean}"
        );
        assert!(image[1000] > 5.0, "star pixel should be positive outlier");
        assert!(image[2000] > 5.0, "bright star should be positive outlier");
    }

    #[test]
    fn clamps_output_range() {
        let mut image = vec![0.0f32; 100];
        image[0] = 1e10;
        mad_normalize(&mut image);
        assert!(image.iter().all(|&v| v >= -5.0 && v <= 100.0));
    }

    #[test]
    fn empty_image_noop() {
        let mut image: Vec<f32> = vec![];
        mad_normalize(&mut image);
        assert!(image.is_empty());
    }

    #[test]
    fn stf_empty_image() {
        let image: Vec<f64> = vec![];
        let result = apply_stf(&image, &StfParams::default());
        assert!(result.is_empty());
    }

    #[test]
    fn stf_uniform_image_stays_uniform() {
        let image = vec![0.1; 1000];
        let result = apply_stf(&image, &StfParams::default());
        assert_eq!(result.len(), 1000);
        let first = result[0];
        for (i, &v) in result.iter().enumerate() {
            assert!(
                libm::fabs(v - first) < 1e-10,
                "pixel {} differs: {} vs {}",
                i,
                v,
                first
            );
        }
    }

    #[test]
    fn stf_dark_pixels_get_boosted() {
        let mut image = vec![0.05; 1000];
        for (i, px) in image.iter_mut().enumerate() {
            *px += (i % 10) as f64 * 0.001;
        }
        let params = StfParams::default();
        let result = apply_stf(&image, &params);
        let median_input = 0.05 + 0.0045;
        let median_output: f64 = {
            let mut sorted = result.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sorted[sorted.len() / 2]
        };
        assert!(
            median_output > median_input,
            "median should be boosted from {} to above, got {}",
            median_input,
            median_output
        );
    }

    #[test]
    fn stf_bright_pixels_dont_clip() {
        let mut image = vec![0.05; 1000];
        image[500] = 0.9;
        let result = apply_stf(&image, &StfParams::default());
        assert!(
            result[500] <= 1.0,
            "bright pixel should not exceed 1.0, got {}",
            result[500]
        );
        assert!(
            result[500] > 0.5,
            "bright pixel should remain bright, got {}",
            result[500]
        );
    }

    #[test]
    fn stf_output_in_valid_range() {
        let mut image = vec![0.02; 10000];
        for (i, px) in image.iter_mut().enumerate() {
            *px += (i as f64) * 0.0001;
        }
        image[9999] = 1.0;
        let result = apply_stf(&image, &StfParams::default());
        for (i, &v) in result.iter().enumerate() {
            assert!(
                v >= 0.0 && v <= 1.0,
                "pixel {} out of range: {}",
                i,
                v
            );
        }
    }

    #[test]
    fn stf_to_u8_conversion() {
        let stretched = vec![0.0, 0.5, 1.0, 0.25, 0.75];
        let bytes = stf_to_u8(&stretched);
        assert_eq!(bytes[0], 0);
        assert_eq!(bytes[1], 128);
        assert_eq!(bytes[2], 255);
        assert_eq!(bytes[3], 64);
        assert_eq!(bytes[4], 191);
    }

    #[test]
    fn stf_typical_astro_image() {
        let mut image = vec![0.02; 10000];
        for (i, px) in image.iter_mut().enumerate() {
            *px += ((i % 100) as f64) * 0.0001;
        }
        image[5000] = 0.3;
        image[5001] = 0.5;
        image[5002] = 0.8;
        let params = StfParams::default();
        let result = apply_stf(&image, &params);
        let median_result: f64 = {
            let mut sorted = result.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sorted[sorted.len() / 2]
        };
        assert!(
            libm::fabs(median_result - params.target_background) < 0.1,
            "median should be near target_background {}, got {}",
            params.target_background,
            median_result
        );
    }

    #[test]
    fn mtf_boundary_conditions() {
        assert_eq!(mtf(0.0, 0.5), 0.0);
        assert_eq!(mtf(1.0, 0.5), 1.0);
        assert_eq!(mtf(0.5, 0.5), 0.5);
    }

    #[test]
    fn mtf_preserves_midpoint() {
        let m = 0.25;
        let result = mtf(m, m);
        assert!(
            libm::fabs(result - 0.5) < 1e-10,
            "MTF(m, m) should equal 0.5, got {}",
            result
        );
    }

    #[test]
    fn trait_stf_mono() {
        let data: Vec<f32> = (0..100).map(|i| i as f32 / 100.0).collect();
        let img = Image::new(PixelData::F32(data), vec![10, 10]);
        let result = img.stf(&StfParams::default()).unwrap();
        assert!(matches!(result.pixels, PixelData::F64(_)));
        assert_eq!(result.dimensions, vec![10, 10]);
    }

    #[test]
    fn trait_stf_rgb() {
        let mut data: Vec<f64> = vec![0.0; 300];
        for i in 0..100 {
            data[i * 3] = 0.02 + (i % 10) as f64 * 0.001;
            data[i * 3 + 1] = 0.05 + (i % 10) as f64 * 0.002;
            data[i * 3 + 2] = 0.08 + (i % 10) as f64 * 0.003;
        }
        let img = Image::new(PixelData::F64(data), vec![10, 10, 3]);
        let result = img.stf(&StfParams::default()).unwrap();
        assert_eq!(result.dimensions, vec![10, 10, 3]);
        let pixels = result.pixels.as_f64().unwrap();
        assert_eq!(pixels.len(), 300);
        for &v in pixels {
            assert!(v >= 0.0 && v <= 1.0, "out of range: {v}");
        }
    }

    #[test]
    fn trait_stf_propagates_keywords() {
        let data = vec![0.1f64; 100];
        let mut img = Image::new(PixelData::F64(data), vec![10, 10]);
        img.keywords.push(Keyword::string("OBJECT", "M31"));
        let result = img.stf(&StfParams::default()).unwrap();
        assert_eq!(result.keywords.len(), 1);
        assert_eq!(result.keywords[0].name, "OBJECT");
    }

    #[test]
    fn trait_mad_normalized_mono() {
        let data: Vec<u16> = (0..100).map(|i| i * 100).collect();
        let img = Image::new(PixelData::U16(data), vec![10, 10]);
        let result = img.mad_normalized().unwrap();
        assert!(matches!(result.pixels, PixelData::F32(_)));
        assert_eq!(result.dimensions, vec![10, 10]);
    }

    #[test]
    fn trait_mad_normalized_rgb() {
        let mut data: Vec<f32> = vec![100.0; 300];
        for i in 0..100 {
            data[i * 3] = 100.0 + (i % 5) as f32;
            data[i * 3 + 1] = 200.0 + (i % 5) as f32;
            data[i * 3 + 2] = 50.0 + (i % 5) as f32;
        }
        data[150] = 5000.0;
        let img = Image::new(PixelData::F32(data), vec![10, 10, 3]);
        let result = img.mad_normalized().unwrap();
        assert_eq!(result.dimensions, vec![10, 10, 3]);
        let pixels = result.pixels.as_f32().unwrap();
        assert_eq!(pixels.len(), 300);
    }

    #[test]
    fn trait_to_u8_converts() {
        let data: Vec<f64> = vec![0.0, 0.5, 1.0, 0.25];
        let img = Image::new(PixelData::F64(data), vec![2, 2]);
        let result = img.to_u8().unwrap();
        let bytes = result.pixels.as_u8().unwrap();
        assert_eq!(bytes[0], 0);
        assert_eq!(bytes[1], 128);
        assert_eq!(bytes[2], 255);
        assert_eq!(bytes[3], 64);
    }

    #[test]
    fn trait_stf_empty_returns_error() {
        let img = Image::new(PixelData::F64(vec![]), vec![0, 0]);
        assert!(img.stf(&StfParams::default()).is_err());
    }

    #[test]
    fn per_channel_f64_mono_passthrough() {
        let data = vec![1.0, 2.0, 3.0];
        let result = apply_per_channel_f64(&data, 1, |ch| {
            ch.iter().map(|&v| v * 2.0).collect()
        });
        assert_eq!(result, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn per_channel_f64_rgb_splits_correctly() {
        let data = vec![1.0, 10.0, 100.0, 2.0, 20.0, 200.0];
        let result = apply_per_channel_f64(&data, 3, |ch| {
            ch.iter().map(|&v| v * 2.0).collect()
        });
        assert_eq!(result, vec![2.0, 20.0, 200.0, 4.0, 40.0, 400.0]);
    }

    #[test]
    fn per_channel_f32_mut_rgb() {
        let mut data = vec![1.0f32, 10.0, 100.0, 2.0, 20.0, 200.0];
        apply_per_channel_f32_mut(&mut data, 3, |ch| {
            for v in ch.iter_mut() {
                *v *= 2.0;
            }
        });
        assert_eq!(data, vec![2.0, 20.0, 200.0, 4.0, 40.0, 400.0]);
    }
}
