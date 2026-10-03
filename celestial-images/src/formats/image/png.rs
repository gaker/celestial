use super::Image;
use crate::core::{ImageError, Result};
use crate::formats::pixel_data::PixelData;
use std::io::BufReader;
use std::path::Path;

impl Image {
    pub fn open_png<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let buf_reader = BufReader::new(file);
        let decoder = png::Decoder::new(buf_reader);
        let mut reader = decoder
            .read_info()
            .map_err(|e| ImageError::FormatDetectionFailed(format!("PNG decode error: {}", e)))?;

        let buf_size = reader.output_buffer_size().ok_or_else(|| {
            ImageError::FormatDetectionFailed("Cannot determine PNG output buffer size".to_string())
        })?;
        let mut buf = vec![0; buf_size];
        let info = reader
            .next_frame(&mut buf)
            .map_err(|e| ImageError::FormatDetectionFailed(format!("PNG frame error: {}", e)))?;
        buf.truncate(info.buffer_size());

        let (pixels, dimensions) = parse_png_data(&buf, &info)?;
        Ok(Self {
            pixels,
            dimensions,
            keywords: Vec::new(),
            xisf_properties: Vec::new(),
        })
    }

    pub(super) fn save_png(&self, path: &Path) -> Result<()> {
        let (width, height, channels) = self.extract_dimensions();
        let color_type = channels_to_png_color_type(channels)?;
        let file = std::fs::File::create(path)?;

        match &self.pixels {
            PixelData::U8(data) => write_png_u8(file, width, height, color_type, data),
            PixelData::U16(data) => write_png_u16(file, width, height, color_type, data),
            _ => Err(ImageError::UnsupportedFormat),
        }
    }
}

fn parse_png_data(buf: &[u8], info: &png::OutputInfo) -> Result<(PixelData, Vec<usize>)> {
    let stored = stored_samples(info.color_type)?;
    let channels = kept_samples(stored);

    // An alpha channel has no meaning downstream, where `channels()` is read
    // as 1 (gray) or 3 (RGB). Drop it at load, the way the TIFF reader drops
    // extra samples, so a GrayA/RGBA PNG never reaches code that would treat
    // its samples as extra pixels.
    let pixels = decode_png_samples(buf, info.bit_depth)?.keep_leading_samples(stored, channels);

    Ok((pixels, build_png_dimensions(info, channels)))
}

/// Samples per pixel stored in the decoded buffer for a PNG color type.
fn stored_samples(color_type: png::ColorType) -> Result<usize> {
    use png::ColorType;

    match color_type {
        ColorType::Grayscale => Ok(1),
        ColorType::GrayscaleAlpha => Ok(2),
        ColorType::Rgb => Ok(3),
        ColorType::Rgba => Ok(4),
        ColorType::Indexed => Err(ImageError::UnsupportedFormat),
    }
}

/// Channels kept in memory: gray+alpha collapses to gray, RGBA to RGB.
fn kept_samples(stored: usize) -> usize {
    match stored {
        2 => 1,
        4 => 3,
        other => other,
    }
}

fn decode_png_samples(buf: &[u8], bit_depth: png::BitDepth) -> Result<PixelData> {
    match bit_depth {
        png::BitDepth::Eight => Ok(PixelData::U8(buf.to_vec())),
        png::BitDepth::Sixteen => Ok(PixelData::U16(
            buf.chunks_exact(2)
                .map(|b| u16::from_be_bytes([b[0], b[1]]))
                .collect(),
        )),
        _ => Err(ImageError::UnsupportedFormat),
    }
}

fn build_png_dimensions(info: &png::OutputInfo, channels: usize) -> Vec<usize> {
    let (width, height) = (info.width as usize, info.height as usize);
    if channels == 1 {
        vec![width, height]
    } else {
        vec![width, height, channels]
    }
}

fn channels_to_png_color_type(channels: usize) -> Result<png::ColorType> {
    match channels {
        1 => Ok(png::ColorType::Grayscale),
        2 => Ok(png::ColorType::GrayscaleAlpha),
        3 => Ok(png::ColorType::Rgb),
        4 => Ok(png::ColorType::Rgba),
        _ => Err(ImageError::UnsupportedFormat),
    }
}

fn write_png_u8(
    file: std::fs::File,
    width: usize,
    height: usize,
    color_type: png::ColorType,
    data: &[u8],
) -> Result<()> {
    let mut encoder = png::Encoder::new(file, width as u32, height as u32);
    encoder.set_color(color_type);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder
        .write_header()
        .map_err(|e| ImageError::FormatDetectionFailed(format!("PNG header error: {}", e)))?;
    writer
        .write_image_data(data)
        .map_err(|e| ImageError::FormatDetectionFailed(format!("PNG write error: {}", e)))
}

fn write_png_u16(
    file: std::fs::File,
    width: usize,
    height: usize,
    color_type: png::ColorType,
    data: &[u16],
) -> Result<()> {
    let mut encoder = png::Encoder::new(file, width as u32, height as u32);
    encoder.set_color(color_type);
    encoder.set_depth(png::BitDepth::Sixteen);
    let mut writer = encoder
        .write_header()
        .map_err(|e| ImageError::FormatDetectionFailed(format!("PNG header error: {}", e)))?;
    let bytes: Vec<u8> = data.iter().flat_map(|&v| v.to_be_bytes()).collect();
    writer
        .write_image_data(&bytes)
        .map_err(|e| ImageError::FormatDetectionFailed(format!("PNG write error: {}", e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_png() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".png").tempfile().unwrap()
    }

    #[test]
    fn roundtrip_u8_grayscale() {
        let tmp = tmp_png();
        let original: Vec<u8> = (0u8..64).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![8usize, 8]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u8_rgb() {
        let tmp = tmp_png();
        let original: Vec<u8> = (0u8..48).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 4, 3]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn open_rgba_drops_alpha_to_three_channels() {
        let tmp = tmp_png();
        // Two RGBA pixels; the alpha samples must not survive the load.
        let rgba: Vec<u8> = vec![10, 20, 30, 200, 40, 50, 60, 100];
        let img = Image::new(PixelData::U8(rgba), vec![2usize, 1, 4]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![2usize, 1, 3]);
        assert!(restored.is_rgb());
        assert_eq!(
            restored.pixels.as_u8().unwrap(),
            &vec![10, 20, 30, 40, 50, 60]
        );
    }

    #[test]
    fn open_rgba_u16_drops_alpha_to_three_channels() {
        let tmp = tmp_png();
        let rgba: Vec<u16> = vec![1000, 2000, 3000, 65535, 4000, 5000, 6000, 65535];
        let img = Image::new(PixelData::U16(rgba), vec![2usize, 1, 4]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![2usize, 1, 3]);
        assert_eq!(
            restored.pixels.as_u16().unwrap(),
            &vec![1000, 2000, 3000, 4000, 5000, 6000]
        );
    }

    #[test]
    fn open_gray_alpha_drops_alpha_to_one_channel() {
        let tmp = tmp_png();
        let gray_alpha: Vec<u8> = vec![10, 255, 20, 128, 30, 0, 40, 255];
        let img = Image::new(PixelData::U8(gray_alpha), vec![4usize, 1, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 1]);
        assert_eq!(restored.channels(), 1);
        assert_eq!(restored.pixels.as_u8().unwrap(), &vec![10, 20, 30, 40]);
    }

    #[test]
    fn roundtrip_u16_grayscale() {
        let tmp = tmp_png();
        let original: Vec<u16> = (0u16..64).map(|i| i * 1000).collect();
        let img = Image::new(PixelData::U16(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn save_rejects_unsupported_pixel_type() {
        let tmp = tmp_png();
        let img = Image::new(PixelData::F32(vec![0.0; 16]), vec![4usize, 4]);
        assert!(img.save(tmp.path()).is_err());
    }

    #[test]
    fn save_rejects_unsupported_channel_count() {
        let tmp = tmp_png();
        let img = Image::new(PixelData::U8(vec![0; 80]), vec![4usize, 4, 5]);
        assert!(img.save(tmp.path()).is_err());
    }

    #[test]
    fn open_returns_error_for_invalid_file() {
        let tmp = tempfile::Builder::new().suffix(".png").tempfile().unwrap();
        std::fs::write(tmp.path(), b"not a png").unwrap();
        assert!(Image::open_png(tmp.path()).is_err());
    }

    #[test]
    fn save_rgba_grayscale_alpha_and_rgba_channel_types() {
        // 2-channel (grayscale-alpha)
        let tmp = tmp_png();
        let img = Image::new(PixelData::U8(vec![0; 32]), vec![4usize, 4, 2]);
        img.save(tmp.path()).unwrap();

        // 4-channel (rgba)
        let tmp = tmp_png();
        let img = Image::new(PixelData::U8(vec![0; 64]), vec![4usize, 4, 4]);
        img.save(tmp.path()).unwrap();
    }
}
