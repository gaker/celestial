use super::Image;
use crate::core::{ImageError, Result};
use crate::formats::pixel_data::PixelData;
use std::path::Path;
use tiff::encoder::{colortype, TiffEncoder};
use tiff::ColorType;

impl Image {
    pub fn open_tiff<P: AsRef<Path>>(path: P) -> Result<Self> {
        use tiff::decoder::Decoder;

        let file = std::fs::File::open(path)?;
        let mut decoder = Decoder::new(file)
            .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF decode error: {}", e)))?;

        let (width, height) = decoder.dimensions().map_err(|e| {
            ImageError::FormatDetectionFailed(format!("TIFF dimensions error: {}", e))
        })?;

        // Samples-per-pixel comes from the TIFF colortype, not from the
        // decoded buffer (which is always flat). Without this, RGB TIFFs
        // get reported as 1-channel and downstream code treats the
        // interleaved R/G/B samples as 3× as many grayscale pixels.
        let colortype = decoder.colortype().map_err(|e| {
            ImageError::FormatDetectionFailed(format!("TIFF colortype error: {}", e))
        })?;

        let image = decoder
            .read_image()
            .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF read error: {}", e)))?;

        let pixels = decode_tiff_image(image)?;

        // `read_image` returns every stored sample, including extra samples the
        // colortype doesn't account for — most commonly an associated-alpha
        // channel on an RGB image (PixInsight writes RGBA with photometric=RGB,
        // so colortype() reports RGB / 3 samples while the buffer holds 4). The
        // true sample count is the buffer length per pixel; reconcile against
        // the colortype and drop the trailing extra samples so the in-memory
        // image is a clean 1- or 3-channel buffer.
        let pixel_count = (width as usize) * (height as usize);
        let actual = actual_samples(pixels.len(), pixel_count);
        let samples = reconcile_samples(pixels.len(), pixel_count, colortype)?;
        let pixels = pixels.keep_leading_samples(actual, samples);
        let dimensions = build_tiff_dimensions(width, height, samples);

        Ok(Self {
            pixels,
            dimensions,
            keywords: Vec::new(),
            xisf_properties: Vec::new(),
        })
    }

    pub(super) fn save_tiff(&self, path: &Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut encoder = TiffEncoder::new(file)
            .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF encoder error: {}", e)))?;

        let (width, height, channels) = self.extract_dimensions();
        write_tiff_image(
            &self.pixels,
            &mut encoder,
            width as u32,
            height as u32,
            channels,
        )
    }
}

/// Samples per pixel actually present in the decoded buffer.
fn actual_samples(buffer_len: usize, pixel_count: usize) -> usize {
    if pixel_count == 0 {
        return 0;
    }
    buffer_len / pixel_count
}

/// Resolve how many channels to keep, given the decoded buffer's true sample
/// count and the colortype. When the buffer carries more samples than the
/// colortype (extra/alpha samples), keep only the logical channels (1 for
/// gray-like, 3 for RGB-like) and drop the rest. When the two agree, keep the
/// colortype's channels minus any alpha it declares.
fn reconcile_samples(buffer_len: usize, pixel_count: usize, colortype: ColorType) -> Result<usize> {
    let actual = actual_samples(buffer_len, pixel_count);
    if pixel_count == 0 || actual == 0 || !buffer_len.is_multiple_of(pixel_count) {
        return Err(ImageError::FormatDetectionFailed(format!(
            "TIFF buffer ({buffer_len}) is not a whole number of samples per pixel ({pixel_count})"
        )));
    }
    let logical = samples_per_pixel(colortype);
    if actual <= logical {
        return Ok(kept_samples(colortype).min(actual));
    }
    // More samples than the colortype expects: keep the logical channels.
    // Gray-like (1–2) collapse to 1, RGB-like (3–4) collapse to 3.
    Ok(if logical >= 3 { 3 } else { 1 })
}

/// Samples per pixel the colortype itself declares.
fn samples_per_pixel(colortype: ColorType) -> usize {
    match colortype {
        ColorType::Gray(_) | ColorType::Palette(_) => 1,
        ColorType::GrayA(_) => 2,
        ColorType::RGB(_) | ColorType::YCbCr(_) | ColorType::Lab(_) => 3,
        ColorType::RGBA(_) | ColorType::CMYK(_) => 4,
        ColorType::CMYKA(_) => 5,
        ColorType::Multiband { num_samples, .. } => num_samples as usize,
        _ => 1,
    }
}

/// Channels kept in memory: an alpha channel is dropped, so GrayA lands as
/// gray and RGBA as RGB. Colortypes with no alpha keep every sample.
fn kept_samples(colortype: ColorType) -> usize {
    match colortype {
        ColorType::GrayA(_) => 1,
        ColorType::RGBA(_) => 3,
        ColorType::CMYKA(_) => 4,
        other => samples_per_pixel(other),
    }
}

fn decode_tiff_image(image: tiff::decoder::DecodingResult) -> Result<PixelData> {
    use tiff::decoder::DecodingResult;

    match image {
        DecodingResult::U8(data) => Ok(PixelData::U8(data)),
        DecodingResult::U16(data) => Ok(PixelData::U16(data)),
        DecodingResult::U32(data) => {
            let converted: Vec<i32> = data.iter().map(|&v| v as i32).collect();
            Ok(PixelData::I32(converted))
        }
        DecodingResult::F32(data) => Ok(PixelData::F32(data)),
        DecodingResult::F64(data) => Ok(PixelData::F64(data)),
        _ => Err(ImageError::UnsupportedFormat),
    }
}

fn build_tiff_dimensions(width: u32, height: u32, channels: usize) -> Vec<usize> {
    if channels == 1 {
        vec![width as usize, height as usize]
    } else {
        vec![width as usize, height as usize, channels]
    }
}

fn write_tiff_image(
    pixels: &PixelData,
    encoder: &mut TiffEncoder<std::fs::File>,
    width: u32,
    height: u32,
    channels: usize,
) -> Result<()> {
    match (pixels, channels) {
        (PixelData::U8(data), 1) => encoder.write_image::<colortype::Gray8>(width, height, data),
        (PixelData::U8(data), 3) => encoder.write_image::<colortype::RGB8>(width, height, data),
        (PixelData::U16(data), 1) => encoder.write_image::<colortype::Gray16>(width, height, data),
        (PixelData::U16(data), 3) => encoder.write_image::<colortype::RGB16>(width, height, data),
        (PixelData::I32(data), 1) => {
            let converted: Vec<u32> = data.iter().map(|&v| v as u32).collect();
            encoder.write_image::<colortype::Gray32>(width, height, &converted)
        }
        (PixelData::F32(data), 1) => {
            encoder.write_image::<colortype::Gray32Float>(width, height, data)
        }
        (PixelData::F32(data), 3) => {
            encoder.write_image::<colortype::RGB32Float>(width, height, data)
        }
        _ => return Err(ImageError::UnsupportedFormat),
    }
    .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF write error: {}", e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_tiff() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".tiff").tempfile().unwrap()
    }

    #[test]
    fn roundtrip_u8_gray() {
        let tmp = tmp_tiff();
        let original: Vec<u8> = (0u8..64).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u16_gray() {
        let tmp = tmp_tiff();
        let original: Vec<u16> = (0u16..64).map(|i| i * 100).collect();
        let img = Image::new(PixelData::U16(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f32_gray() {
        let tmp = tmp_tiff();
        let original: Vec<f32> = (0..64).map(|i| i as f32 * 0.25).collect();
        let img = Image::new(PixelData::F32(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
    }

    #[test]
    fn save_rgb_u8() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::U8(vec![0; 48]), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();
    }

    #[test]
    fn save_rgb_u16() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::U16(vec![0; 48]), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();
    }

    #[test]
    fn roundtrip_rgb_u8_reports_three_channels() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::U8((0u8..48).collect()), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert!(
            restored.is_rgb(),
            "RGB TIFF should round-trip with channels==3"
        );
        assert_eq!(restored.channels(), 3);
        assert_eq!(restored.width(), 4);
        assert_eq!(restored.height(), 4);
        assert_eq!(restored.pixels.as_u8().unwrap().len(), 48);
    }

    #[test]
    fn roundtrip_rgb_u16_reports_three_channels() {
        let tmp = tmp_tiff();
        let img = Image::new(
            PixelData::U16((0u16..48).map(|i| i * 1000).collect()),
            vec![4usize, 4, 3],
        );
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert!(restored.is_rgb());
        assert_eq!(restored.channels(), 3);
        assert_eq!(restored.pixels.as_u16().unwrap().len(), 48);
    }

    #[test]
    fn save_rejects_unsupported_combos() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::I16(vec![0; 16]), vec![4usize, 4]);
        assert!(img.save(tmp.path()).is_err());
    }

    #[test]
    fn reconcile_keeps_matching_samples() {
        // 3 samples/pixel buffer, colortype says RGB → keep 3.
        assert_eq!(reconcile_samples(4 * 3, 4, ColorType::RGB(8)).unwrap(), 3);
        // 1 sample/pixel gray.
        assert_eq!(reconcile_samples(4, 4, ColorType::Gray(8)).unwrap(), 1);
    }

    #[test]
    fn reconcile_strips_rgb_alpha() {
        // PixInsight RGBA: buffer has 4 samples/pixel but colortype reports
        // RGB (3). Keep 3, drop the assoc-alpha.
        assert_eq!(reconcile_samples(4 * 4, 4, ColorType::RGB(8)).unwrap(), 3);
    }

    #[test]
    fn reconcile_strips_declared_alpha() {
        // Colortype and buffer agree on 4 samples: it is a real RGBA TIFF, and
        // the alpha still has to go so the image lands as 3-channel RGB.
        assert_eq!(reconcile_samples(4 * 4, 4, ColorType::RGBA(8)).unwrap(), 3);
        assert_eq!(reconcile_samples(4 * 2, 4, ColorType::GrayA(8)).unwrap(), 1);
    }

    #[test]
    fn reconcile_strips_gray_alpha() {
        // GrayA stored against a Gray colortype: 2 samples/pixel, logical 1.
        assert_eq!(reconcile_samples(4 * 2, 4, ColorType::Gray(8)).unwrap(), 1);
    }

    #[test]
    fn reconcile_rejects_non_integer_samples() {
        // A buffer that isn't a whole number of samples per pixel.
        assert!(reconcile_samples(10, 4, ColorType::RGB(8)).is_err());
    }

    #[test]
    fn open_rgba_tiff_drops_alpha() {
        let tmp = tmp_tiff();
        // One row of two pixels: red then green, both opaque.
        let rgba: Vec<u8> = vec![10, 20, 30, 255, 40, 50, 60, 255];
        let file = std::fs::File::create(tmp.path()).unwrap();
        TiffEncoder::new(file)
            .unwrap()
            .write_image::<colortype::RGBA8>(2, 1, &rgba)
            .unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert!(restored.is_rgb());
        assert_eq!(restored.channels(), 3);
        assert_eq!(
            restored.pixels.as_u8().unwrap(),
            &vec![10, 20, 30, 40, 50, 60]
        );
    }

    #[test]
    fn open_returns_error_on_invalid_tiff() {
        let tmp = tempfile::Builder::new().suffix(".tiff").tempfile().unwrap();
        std::fs::write(tmp.path(), b"garbage").unwrap();
        assert!(Image::open_tiff(tmp.path()).is_err());
    }

    #[test]
    fn i32_saves_as_gray32() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::I32(vec![0i32, 1000, -500, 42]), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();
    }
}
