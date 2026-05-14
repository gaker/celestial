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
    use png::ColorType;

    let channels = match info.color_type {
        ColorType::Grayscale => 1,
        ColorType::Rgb => 3,
        ColorType::GrayscaleAlpha => 2,
        ColorType::Rgba => 4,
        ColorType::Indexed => {
            return Err(ImageError::UnsupportedFormat);
        }
    };

    let dimensions = if channels == 1 {
        vec![info.width as usize, info.height as usize]
    } else {
        vec![info.width as usize, info.height as usize, channels]
    };

    let pixels = match info.bit_depth {
        png::BitDepth::Eight => PixelData::U8(buf.to_vec()),
        png::BitDepth::Sixteen => {
            let data: Vec<u16> = buf
                .chunks_exact(2)
                .map(|b| u16::from_be_bytes([b[0], b[1]]))
                .collect();
            PixelData::U16(data)
        }
        _ => return Err(ImageError::UnsupportedFormat),
    };

    Ok((pixels, dimensions))
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
        tempfile::Builder::new()
            .suffix(".png")
            .tempfile()
            .unwrap()
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
