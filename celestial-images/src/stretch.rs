use std::cmp::Ordering;
use std::ops::{Add, Div, Sub};

use crate::core::errors::{ImageError, Result};
use crate::core::sample_order::{to_interleaved, to_planar};
use crate::formats::image::Image;
use crate::formats::pixel_data::PixelData;

// For normally distributed data the standard deviation is 1.4826 times the
// median absolute deviation (MAD).
const MAD_TO_SIGMA: f64 = 1.4826;

pub struct StfParams {
    pub shadow_clip: f64,
    pub target_background: f64,
}

impl Default for StfParams {
    fn default() -> Self {
        Self {
            shadow_clip: -2.8,
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
        params.validate()?;
        let mut data = self.pixels.to_f64_normalized();
        blank_to_nan(self, &mut data, f64::NAN);
        for_each_channel(&mut data, self.channels(), |ch| stretch_channel(ch, params))?;
        Ok(with_pixels(self, PixelData::F64(data)))
    }

    fn mad_normalized(&self) -> Result<Image> {
        let mut data = self.pixels.to_f32_normalized();
        blank_to_nan(self, &mut data, f32::NAN);
        for_each_channel(&mut data, self.channels(), mad_normalize)?;
        Ok(with_pixels(self, PixelData::F32(data)))
    }

    fn to_u8(&self) -> Result<Image> {
        let data = self.pixels.to_f64_normalized();
        Ok(with_pixels(self, PixelData::U8(stf_to_u8(&data))))
    }
}

impl StfParams {
    // A positive clip puts c0 above the median, which sends most of the
    // background to 0 and lets the inverted branch leave [0, 1].
    fn validate(&self) -> Result<()> {
        if !self.shadow_clip.is_finite() || self.shadow_clip > 0.0 {
            return Err(ImageError::InvalidInput(format!(
                "STF shadow_clip is a non-positive multiple of MADN (default -2.8), got {}",
                self.shadow_clip
            )));
        }
        if !(0.0..=1.0).contains(&self.target_background) {
            return Err(ImageError::InvalidInput(format!(
                "STF target_background must lie in [0, 1], got {}",
                self.target_background
            )));
        }
        Ok(())
    }
}

// The result's samples are no longer the stored values, so a BLANK card
// would mark the wrong pixels.
fn with_pixels(source: &Image, pixels: PixelData) -> Image {
    let mut img = Image::new(pixels, source.dimensions.clone());
    img.keywords = source.keywords.clone();
    img.remove_keyword("BLANK");
    img
}

// BLANK pixels hold no measurement, so they are left out of the statistics
// the same way NaN pixels are.
fn blank_to_nan<T: Copy>(source: &Image, normalized: &mut [T], nan: T) {
    let Some(flags) = source.blank_flags() else {
        return;
    };
    for (v, blank) in normalized.iter_mut().zip(flags) {
        if blank {
            *v = nan;
        }
    }
}

// Runs `f` on each channel of interleaved samples, one contiguous plane at a time.
fn for_each_channel<T: Copy>(
    data: &mut Vec<T>,
    channels: usize,
    mut f: impl FnMut(&mut [T]) -> Result<()>,
) -> Result<()> {
    if data.is_empty() {
        return Err(ImageError::UnsupportedFormat);
    }
    if channels == 0 || !data.len().is_multiple_of(channels) {
        return Err(ImageError::InvalidDimensions(format!(
            "{} samples do not divide into {channels} channels",
            data.len()
        )));
    }
    if channels == 1 {
        return f(data);
    }
    let mut planar = to_planar(data, channels);
    planar
        .chunks_exact_mut(data.len() / channels)
        .try_for_each(&mut f)?;
    *data = to_interleaved(&planar, channels);
    Ok(())
}

trait Real: Copy + PartialOrd + Add<Output = Self> + Sub<Output = Self> + Div<Output = Self> {
    const TWO: Self;
    fn finite(self) -> bool;
    fn magnitude(self) -> Self;
    fn order(&self, other: &Self) -> Ordering;
}

impl Real for f64 {
    const TWO: Self = 2.0;
    fn finite(self) -> bool {
        self.is_finite()
    }
    fn magnitude(self) -> Self {
        libm::fabs(self)
    }
    fn order(&self, other: &Self) -> Ordering {
        self.total_cmp(other)
    }
}

impl Real for f32 {
    const TWO: Self = 2.0;
    fn finite(self) -> bool {
        self.is_finite()
    }
    fn magnitude(self) -> Self {
        libm::fabsf(self)
    }
    fn order(&self, other: &Self) -> Ordering {
        self.total_cmp(other)
    }
}

// Median and median absolute deviation of the finite samples, or None when
// there are none. NaN marks a blank pixel. Infinities are left out too,
// because the normalization that runs first ignores them for its min and max.
fn median_mad<T: Real>(samples: &[T]) -> Option<(T, T)> {
    let mut buf: Vec<T> = samples.iter().copied().filter(|v| v.finite()).collect();
    let median = median_in_place(&mut buf)?;
    for v in buf.iter_mut() {
        *v = (*v - median).magnitude();
    }
    let mad = median_in_place(&mut buf)?;
    Some((median, mad))
}

// An even count averages the two middle values. The buffer holds no NaN, so
// `>` finds the same largest lower value as `total_cmp` up to the sign of a
// zero, which does not change the mean, and it runs about 20 times faster.
fn median_in_place<T: Real>(buf: &mut [T]) -> Option<T> {
    if buf.is_empty() {
        return None;
    }
    let mid = buf.len() / 2;
    let odd = buf.len() % 2 == 1;
    let (lower, &mut upper, _) = buf.select_nth_unstable_by(mid, T::order);
    if odd {
        return Some(upper);
    }
    let below = lower
        .iter()
        .copied()
        .reduce(|a, b| if b > a { b } else { a })?;
    Some((below + upper) / T::TWO)
}

fn no_finite_samples() -> ImageError {
    ImageError::InvalidInput("channel has no finite samples".to_string())
}

fn stretch_channel(channel: &mut [f64], params: &StfParams) -> Result<()> {
    let (median, mad) = median_mad(channel).ok_or_else(no_finite_samples)?;
    let stf = AutoStf::new(median, MAD_TO_SIGMA * mad, params);
    for v in channel.iter_mut() {
        *v = stf.apply(*v);
    }
    Ok(())
}

// Automatic screen transfer function in the PixInsight AutoSTF convention,
// with C = shadow_clip, B = target_background and MADN = 1.4826 * MAD:
//   median < 0.5:  c0 = clamp(median + C * MADN, 0, 1), c1 = 1, m = MTF(B, median - c0)
//   median >= 0.5, an inverted channel (the unlinked rule):
//                  c0 = 0, c1 = clamp(median - C * MADN, 0, 1), m = MTF(c1 - median, B)
// A pixel x maps to MTF(m, x') with x' = (x - c0) / (c1 - c0) clipped to [0, 1].
struct AutoStf {
    shadows: f64,
    highlights: f64,
    midtones: f64,
}

impl AutoStf {
    fn new(median: f64, madn: f64, params: &StfParams) -> Self {
        // The convention skips the clip when MADN is too small to change 1.0.
        let clip = (1.0 + madn != 1.0).then_some(params.shadow_clip * madn);
        let background = params.target_background;
        if median < 0.5 {
            let shadows = clip.map_or(0.0, |c| (median + c).clamp(0.0, 1.0));
            return Self {
                shadows,
                highlights: 1.0,
                midtones: mtf(background, median - shadows),
            };
        }
        let highlights = clip.map_or(1.0, |c| (median - c).clamp(0.0, 1.0));
        Self {
            shadows: 0.0,
            highlights,
            midtones: mtf(highlights - median, background),
        }
    }

    fn apply(&self, x: f64) -> f64 {
        if x.is_nan() {
            return x;
        }
        let clipped = if x <= self.shadows {
            0.0
        } else if x >= self.highlights {
            1.0
        } else {
            (x - self.shadows) / (self.highlights - self.shadows)
        };
        mtf(self.midtones, clipped)
    }
}

// Midtones transfer function with midtones balance m:
// MTF(m, 0) = 0, MTF(m, m) = 1/2, MTF(m, 1) = 1.
fn mtf(m: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    (m - 1.0) * x / ((2.0 * m - 1.0) * x - m)
}

fn stf_to_u8(stretched: &[f64]) -> Vec<u8> {
    stretched
        .iter()
        .map(|&v| libm::round(v.clamp(0.0, 1.0) * 255.0) as u8)
        .collect()
}

fn mad_normalize(channel: &mut [f32]) -> Result<()> {
    let (median, mad) = median_mad(channel).ok_or_else(no_finite_samples)?;
    let sigma = (mad * MAD_TO_SIGMA as f32).max(1e-10);
    for pixel in channel.iter_mut() {
        *pixel = ((*pixel - median) / sigma).clamp(-5.0, 100.0);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::Keyword;

    #[test]
    fn background_centered_near_zero() {
        let mut image = vec![100.0f32; 256 * 256];
        for (i, px) in image.iter_mut().enumerate() {
            *px += (i % 7) as f32 * 0.5;
        }
        image[1000] = 5000.0;
        image[2000] = 10000.0;

        mad_normalize(&mut image).unwrap();

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
        mad_normalize(&mut image).unwrap();
        assert!(image.iter().all(|&v| (-5.0..=100.0).contains(&v)));
    }

    #[test]
    fn mad_normalize_of_no_finite_samples_is_an_error() {
        assert!(matches!(
            mad_normalize(&mut []),
            Err(ImageError::InvalidInput(_))
        ));
        assert!(matches!(
            mad_normalize(&mut [f32::NAN, f32::INFINITY]),
            Err(ImageError::InvalidInput(_))
        ));
    }

    // Background near 1000 ADU with a few bright pixels, as a detector frame.
    fn adu_frame() -> Vec<i16> {
        (0..4096)
            .map(|i| 900 + ((i * 7919) % 1000) as i16 + if i % 97 == 0 { 20000 } else { 0 })
            .collect()
    }

    fn physical_image(scale: impl Fn(f64) -> f64) -> Image {
        let values = adu_frame().iter().map(|&v| scale(v as f64)).collect();
        Image::new(PixelData::F64(values), vec![64usize, 64])
    }

    #[test]
    fn stf_and_to_u8_of_scaled_integer_data_match_the_unscaled_values() {
        let plain = physical_image(|v| v);
        let half = physical_image(|v| 0.5 * v);
        let offset = physical_image(|v| 32767.0 + v);

        let stf = |img: &Image| img.stf(&StfParams::default()).unwrap().pixels;
        let bytes = |img: &Image| img.to_u8().unwrap().pixels;
        let reference = stf(&plain);
        let reference = reference.as_f64().unwrap();
        assert!(reference.iter().any(|&v| v != reference[0]), "flat stf");
        assert_eq!(stf(&half).as_f64(), Some(reference));
        assert_eq!(stf(&offset).as_f64(), Some(reference));

        let reference = bytes(&plain);
        let reference = reference.as_u8().unwrap();
        assert!(reference.iter().any(|&v| v != reference[0]), "flat to_u8");
        assert_eq!(bytes(&half).as_u8(), Some(reference));
        assert_eq!(bytes(&offset).as_u8(), Some(reference));
    }

    #[test]
    fn stf_of_64_bit_integer_data_matches_the_same_values_as_f64() {
        let reference = physical_image(|v| v)
            .stf(&StfParams::default())
            .unwrap()
            .pixels;
        let adu = adu_frame();
        let shifted_i64 = adu.iter().map(|&v| (1i64 << 62) + v as i64).collect();
        let shifted_u64 = adu.iter().map(|&v| (1u64 << 63) + v as u64).collect();
        for pixels in [PixelData::I64(shifted_i64), PixelData::U64(shifted_u64)] {
            let img = Image::new(pixels, vec![64usize, 64]);
            let stf = img.stf(&StfParams::default()).unwrap().pixels;
            assert_eq!(stf.as_f64(), reference.as_f64());
        }
    }

    fn sky_u16() -> Image {
        let pixels = (0..10000u32)
            .map(|i| {
                let star = if i % 97 == 0 { 30000 } else { 0 };
                (6000 + (i.wrapping_mul(2654435761) >> 20) % 1500 + star) as u16
            })
            .collect();
        Image::new(PixelData::U16(pixels), vec![100usize, 100])
    }

    fn stf_pixels(img: &Image) -> Vec<f64> {
        let out = img.stf(&StfParams::default()).unwrap();
        out.pixels.as_f64().unwrap().clone()
    }

    #[test]
    fn stf_of_a_sky_matches_the_auto_stf_reference() {
        let out = stf_pixels(&sky_u16());
        assert_eq!(out[0], 0.9395827656604909);
        assert_eq!(out[1], 0.3115140882046876);
        assert_eq!(out[50], 0.2668813744751756);
        let mut sorted = out.clone();
        sorted.sort_by(f64::total_cmp);
        assert_eq!(sorted[5000], 0.2664617492876044);
    }

    #[test]
    fn stf_rejects_a_positive_shadow_clip() {
        let sky = sky_u16();
        for shadow_clip in [2.8, 1e-9] {
            let params = StfParams {
                shadow_clip,
                target_background: 0.25,
            };
            let result = sky.stf(&params);
            assert!(
                matches!(result, Err(ImageError::InvalidInput(_))),
                "shadow_clip {shadow_clip}"
            );
        }
        let zero = StfParams {
            shadow_clip: 0.0,
            target_background: 0.25,
        };
        assert!(sky.stf(&zero).is_ok());
    }

    // Every 13th pixel holds the BLANK value 0. Those pixels carry no
    // measurement, so the result must equal the stretch of the same
    // normalized samples with NaN in their place.
    #[test]
    fn stf_treats_blank_pixels_as_nan() {
        let mut raw = sky_u16().pixels.as_u16().unwrap().clone();
        raw.iter_mut().step_by(13).for_each(|v| *v = 0);
        let masked: Vec<f64> = raw
            .iter()
            .map(|&v| if v == 0 { f64::NAN } else { v as f64 / 65535.0 })
            .collect();
        let mut img = Image::new(PixelData::U16(raw), vec![100usize, 100]);
        img.keywords.push(Keyword::integer("BLANK", 0));
        let got = stf_pixels(&img);
        let want = stf_pixels(&Image::new(PixelData::F64(masked), vec![100usize, 100]));
        assert_eq!(got.len(), want.len());
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(g.to_bits(), w.to_bits(), "pixel {i}: {g} vs {w}");
        }
        assert!(got[0].is_nan() && got[13].is_nan());
        let out = img.stf(&StfParams::default()).unwrap();
        assert!(out.get_keyword("BLANK").is_none());
    }

    #[test]
    fn mad_normalized_treats_blank_pixels_as_nan() {
        let raw: Vec<u16> = (0..100)
            .map(|i| if i % 9 == 0 { 7 } else { 100 + i })
            .collect();
        let masked: Vec<f32> = raw
            .iter()
            .map(|&v| if v == 7 { f32::NAN } else { v as f32 / 65535.0 })
            .collect();
        let mut img = Image::new(PixelData::U16(raw), vec![10usize, 10]);
        img.keywords.push(Keyword::integer("BLANK", 7));
        let got = img.mad_normalized().unwrap();
        let want = Image::new(PixelData::F32(masked), vec![10usize, 10])
            .mad_normalized()
            .unwrap();
        let bits = |img: &Image| -> Vec<u32> {
            img.pixels
                .as_f32()
                .unwrap()
                .iter()
                .map(|v| v.to_bits())
                .collect()
        };
        assert_eq!(bits(&got), bits(&want));
        assert!(got.pixels.as_f32().unwrap()[0].is_nan());
    }

    #[test]
    fn stf_of_adu_floats_rescales_by_the_data_range() {
        let adu = (0..1000).map(|i| 1000.0 + (i % 100) as f32).collect();
        let out = stf_pixels(&Image::new(PixelData::F32(adu), vec![50usize, 20]));
        assert_eq!(out[0], 0.0);
        assert_eq!(out[1], 0.029702970297029705);
        assert_eq!(out[2], 0.05825242718446602);
        assert_eq!(out[50], 0.7537688442211055);
        assert_eq!(out[99], 1.0);
    }

    #[test]
    fn stf_of_an_inverted_image_clips_the_highlights() {
        let sky = sky_u16();
        let inverted = sky.pixels.as_u16().unwrap().iter().map(|&v| 65535 - v);
        let img = Image::new(PixelData::U16(inverted.collect()), vec![100usize, 100]);
        let out = stf_pixels(&img);
        assert_eq!(out[0], 0.06041723433950893);
        assert_eq!(out[1], 0.6884859117953128);
        assert_eq!(out[50], 0.7331186255248239);
    }

    #[test]
    fn stf_of_a_constant_image_skips_the_clip() {
        let dark = Image::new(PixelData::F32(vec![0.1; 1000]), vec![40usize, 25]);
        assert!(stf_pixels(&dark).iter().all(|&v| v == 0.25000000000000006));
        let bright = Image::new(PixelData::U16(vec![40000; 100]), vec![10usize, 10]);
        assert!(stf_pixels(&bright).iter().all(|&v| v == 0.7499999999999998));
    }

    #[test]
    fn stf_leaves_nan_pixels_out_of_the_statistics() {
        let mut values: Vec<f32> = (0..4000)
            .map(|i| ((i * 7919) % 1000) as f32 / 4000.0)
            .collect();
        for v in values.iter_mut().step_by(7) {
            *v = f32::NAN;
        }
        let finite: Vec<f32> = values.iter().copied().filter(|v| !v.is_nan()).collect();
        let finite_len = finite.len();
        let with_nan = stf_pixels(&Image::new(PixelData::F32(values), vec![80usize, 50]));
        let without = stf_pixels(&Image::new(PixelData::F32(finite), vec![finite_len, 1]));
        let mut kept = without.iter();
        for (i, &v) in with_nan.iter().enumerate() {
            if i % 7 == 0 {
                assert!(v.is_nan(), "pixel {i} is {v}");
            } else {
                assert_eq!(Some(&v), kept.next(), "pixel {i}");
            }
        }
        assert_eq!(with_nan[1], 0.4100961708168686);
    }

    #[test]
    fn median_mad_of_even_and_odd_counts_skips_non_finite_samples() {
        assert_eq!(median_mad(&[4.0, 1.0, 3.0, 2.0]), Some((2.5, 1.0)));
        assert_eq!(
            median_mad(&[5.0f32, 1.0, f32::NAN, 3.0, f32::INFINITY]),
            Some((3.0, 2.0))
        );
        assert_eq!(median_mad::<f64>(&[f64::NAN, f64::NEG_INFINITY]), None);
    }

    #[test]
    fn mad_normalize_averages_the_two_middle_values() {
        let mut image = vec![1.0f32, 2.0, 3.0, 4.0];
        mad_normalize(&mut image).unwrap();
        assert_eq!(image, vec![-1.0117362, -0.33724537, 0.33724537, 1.0117362]);
    }

    fn apply_stf(image: &[f64], params: &StfParams) -> Vec<f64> {
        let mut channel = image.to_vec();
        stretch_channel(&mut channel, params).unwrap();
        channel
    }

    #[test]
    fn stf_of_no_finite_samples_is_an_error() {
        for channel in [vec![], vec![f64::NAN; 4]] {
            let mut channel = channel;
            let result = stretch_channel(&mut channel, &StfParams::default());
            assert!(matches!(result, Err(ImageError::InvalidInput(_))));
        }
        let img = Image::new(PixelData::F32(vec![f32::NAN; 4]), vec![2usize, 2]);
        assert!(img.stf(&StfParams::default()).is_err());
    }

    #[test]
    fn stf_uniform_image_stays_uniform() {
        let image = vec![0.1; 1000];
        let result = apply_stf(&image, &StfParams::default());
        assert_eq!(result.len(), 1000);
        let first = result[0];
        for (i, &v) in result.iter().enumerate() {
            assert_eq!(v, first, "pixel {i}");
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
            assert!((0.0..=1.0).contains(&v), "pixel {} out of range: {}", i, v);
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
        assert_eq!(mtf(0.5, 0.0), 0.0);
        assert_eq!(mtf(0.5, 1.0), 1.0);
        assert_eq!(mtf(0.5, 0.5), 0.5);
    }

    #[test]
    fn mtf_preserves_midpoint() {
        // (0.25 - 1) * 0.25 / ((2 * 0.25 - 1) * 0.25 - 0.25) = -0.1875 / -0.375,
        // and every term is a binary fraction, so the quotient is exactly 0.5.
        let m = 0.25;
        assert_eq!(mtf(m, m), 0.5);
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
            assert!((0.0..=1.0).contains(&v), "out of range: {v}");
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

    fn doubled(channel: &mut [f64]) -> Result<()> {
        channel.iter_mut().for_each(|v| *v *= 2.0);
        Ok(())
    }

    #[test]
    fn for_each_channel_mono_passthrough() {
        let mut data = vec![1.0, 2.0, 3.0];
        for_each_channel(&mut data, 1, doubled).unwrap();
        assert_eq!(data, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn for_each_channel_rgb_splits_correctly() {
        let mut data = vec![1.0, 10.0, 100.0, 2.0, 20.0, 200.0];
        let mut seen = Vec::new();
        for_each_channel(&mut data, 3, |ch| {
            seen.push(ch.to_vec());
            doubled(ch)
        })
        .unwrap();
        assert_eq!(
            seen,
            vec![vec![1.0, 2.0], vec![10.0, 20.0], vec![100.0, 200.0]]
        );
        assert_eq!(data, vec![2.0, 20.0, 200.0, 4.0, 40.0, 400.0]);
    }

    #[test]
    fn for_each_channel_rejects_a_partial_pixel() {
        let mut data = vec![1.0f32; 5];
        let result = for_each_channel(&mut data, 3, |_| Ok(()));
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
        let result = for_each_channel(&mut data, 0, |_| Ok(()));
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
    }
}
