use super::Image;
use crate::core::{ImageError, Result};
use crate::formats::pixel_data::PixelData;
use std::path::Path;
use tiff::encoder::{colortype, TiffEncoder};

impl Image {
    pub fn open_tiff<P: AsRef<Path>>(path: P) -> Result<Self> {
        use tiff::decoder::Decoder;

        let file = std::fs::File::open(path)?;
        let mut decoder = Decoder::new(file)
            .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF decode error: {}", e)))?;

        let (width, height) = decoder.dimensions().map_err(|e| {
            ImageError::FormatDetectionFailed(format!("TIFF dimensions error: {}", e))
        })?;

        let image = decoder
            .read_image()
            .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF read error: {}", e)))?;

        let (pixels, channels) = decode_tiff_image(image)?;
        let dimensions = build_tiff_dimensions(width, height, channels);

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
        write_tiff_image(&self.pixels, &mut encoder, width as u32, height as u32, channels)
    }
}

fn decode_tiff_image(image: tiff::decoder::DecodingResult) -> Result<(PixelData, usize)> {
    use tiff::decoder::DecodingResult;

    match image {
        DecodingResult::U8(data) => Ok((PixelData::U8(data), 1)),
        DecodingResult::U16(data) => Ok((PixelData::U16(data), 1)),
        DecodingResult::U32(data) => {
            let converted: Vec<i32> = data.iter().map(|&v| v as i32).collect();
            Ok((PixelData::I32(converted), 1))
        }
        DecodingResult::F32(data) => Ok((PixelData::F32(data), 1)),
        DecodingResult::F64(data) => Ok((PixelData::F64(data), 1)),
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
        _ => return Err(ImageError::UnsupportedFormat),
    }
    .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF write error: {}", e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_tiff() -> tempfile::NamedTempFile {
        tempfile::Builder::new()
            .suffix(".tiff")
            .tempfile()
            .unwrap()
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
    fn save_rejects_unsupported_combos() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::I16(vec![0; 16]), vec![4usize, 4]);
        assert!(img.save(tmp.path()).is_err());
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
        let img = Image::new(
            PixelData::I32(vec![0i32, 1000, -500, 42]),
            vec![2usize, 2],
        );
        img.save(tmp.path()).unwrap();
    }
}
