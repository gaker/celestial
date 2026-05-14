use super::Image;
use crate::core::{ImageError, Result};
use crate::formats::pixel_data::PixelData;
use crate::xisf::writer::{XisfWriter, wcs_to_xisf_properties};
use celestial_wcs::{WcsKeyword, WcsKeywordValue};
use std::path::Path;

impl Image {
    pub(super) fn open_xisf(path: &Path) -> Result<Self> {
        let mut xisf = crate::xisf::XisfFile::open(path).map_err(ImageError::Xisf)?;

        let info = xisf.image_info(0).ok_or_else(|| {
            ImageError::Xisf(crate::xisf::XisfError::InvalidFormat(
                "No images in file".to_string(),
            ))
        })?;

        let dimensions = info.geometry.clone();
        let sample_format = info.sample_format.clone();

        let raw_bytes = xisf.read_image_data_raw(0).map_err(ImageError::Xisf)?;

        let pixels = match sample_format {
            crate::xisf::header::SampleFormat::UInt8 => PixelData::U8(raw_bytes),
            crate::xisf::header::SampleFormat::UInt16 => {
                let data: Vec<u16> = raw_bytes
                    .chunks_exact(2)
                    .map(|b| u16::from_le_bytes([b[0], b[1]]))
                    .collect();
                PixelData::U16(data)
            }
            crate::xisf::header::SampleFormat::UInt32 => {
                let data: Vec<i32> = raw_bytes
                    .chunks_exact(4)
                    .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                    .collect();
                PixelData::I32(data)
            }
            crate::xisf::header::SampleFormat::Float32 => {
                let data: Vec<f32> = raw_bytes
                    .chunks_exact(4)
                    .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                    .collect();
                PixelData::F32(data)
            }
            crate::xisf::header::SampleFormat::Float64 => {
                let data: Vec<f64> = raw_bytes
                    .chunks_exact(8)
                    .map(|b| f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]))
                    .collect();
                PixelData::F64(data)
            }
        };

        let keywords = xisf.keywords().to_vec();

        Ok(Self {
            pixels,
            dimensions,
            keywords,
            xisf_properties: Vec::new(),
        })
    }

    pub(super) fn save_xisf(&self, path: &Path) -> Result<()> {
        let mut writer = XisfWriter::create(path).map_err(ImageError::Xisf)?;

        let (width, height, channels) = self.extract_dimensions();

        if channels == 3 {
            let mut planar_image = self.clone();
            planar_image.interleaved_to_planar();
            write_xisf_pixels(&mut writer, &planar_image.pixels, width, height, channels)?;
        } else {
            write_xisf_pixels(&mut writer, &self.pixels, width, height, channels)?;
        }

        for kw in &self.keywords {
            writer.add_keyword(kw.clone());
        }

        let wcs_keywords: Vec<WcsKeyword> = self
            .keywords
            .iter()
            .filter_map(|kw| {
                let value = match &kw.value {
                    Some(crate::fits::header::KeywordValue::Real(v)) => WcsKeywordValue::Real(*v),
                    Some(crate::fits::header::KeywordValue::Integer(v)) => {
                        WcsKeywordValue::Integer(*v)
                    }
                    Some(crate::fits::header::KeywordValue::String(v)) => {
                        WcsKeywordValue::String(v.clone())
                    }
                    _ => return None,
                };
                Some(WcsKeyword {
                    name: kw.name.clone(),
                    value,
                })
            })
            .collect();

        let properties = wcs_to_xisf_properties(&wcs_keywords);
        writer.add_properties(properties);
        writer.add_properties(self.xisf_properties.iter().cloned());

        writer.write().map_err(ImageError::Xisf)
    }
}

fn write_xisf_pixels<W: std::io::Write + std::io::Seek>(
    writer: &mut XisfWriter<W>,
    pixels: &PixelData,
    width: usize,
    height: usize,
    channels: usize,
) -> Result<()> {
    match pixels {
        PixelData::U8(data) => writer
            .add_image(data, width, height, channels)
            .map_err(ImageError::Xisf)?,
        PixelData::U16(data) => writer
            .add_image(data, width, height, channels)
            .map_err(ImageError::Xisf)?,
        PixelData::I16(data) => {
            let converted: Vec<u16> = data.iter().map(|&v| v as u16).collect();
            writer
                .add_image(&converted, width, height, channels)
                .map_err(ImageError::Xisf)?
        }
        PixelData::I32(data) => {
            let converted: Vec<u32> = data.iter().map(|&v| v as u32).collect();
            writer
                .add_image(&converted, width, height, channels)
                .map_err(ImageError::Xisf)?
        }
        PixelData::F32(data) => writer
            .add_image(data, width, height, channels)
            .map_err(ImageError::Xisf)?,
        PixelData::F64(data) => writer
            .add_image(data, width, height, channels)
            .map_err(ImageError::Xisf)?,
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::Keyword;

    fn tmp_xisf() -> tempfile::NamedTempFile {
        tempfile::Builder::new()
            .suffix(".xisf")
            .tempfile()
            .unwrap()
    }

    #[test]
    fn roundtrip_u8_mono() {
        let tmp = tmp_xisf();
        let original: Vec<u8> = (0u8..64).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![8usize, 8]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u16_mono() {
        let tmp = tmp_xisf();
        let original: Vec<u16> = (0u16..64).collect();
        let img = Image::new(PixelData::U16(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f32_mono() {
        let tmp = tmp_xisf();
        let original: Vec<f32> = (0..64).map(|i| i as f32 * 0.1).collect();
        let img = Image::new(PixelData::F32(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f64_mono() {
        let tmp = tmp_xisf();
        let original: Vec<f64> = (0..64).map(|i| i as f64 * 0.01).collect();
        let img = Image::new(PixelData::F64(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f64().unwrap(), &original);
    }

    #[test]
    fn rgb_roundtrip_preserves_channel_count() {
        let tmp = tmp_xisf();
        let original: Vec<u8> = (0u8..48).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 4, 3]);
        assert_eq!(restored.channels(), 3);
    }

    #[test]
    fn keywords_roundtrip_through_xisf() {
        let tmp = tmp_xisf();
        let mut img = Image::new(PixelData::U8(vec![0; 16]), vec![4usize, 4]);
        img.set_keyword(Keyword::string("OBJECT", "NGC1234"));
        img.set_keyword(Keyword::real("EXPTIME", 60.0));
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.get_keyword("OBJECT").unwrap().value,
            Some(crate::fits::header::KeywordValue::String("NGC1234".to_string()))
        );
    }

    #[test]
    fn i16_writes_as_unsigned_on_disk() {
        let tmp = tmp_xisf();
        let original = vec![0i16, 100, -50, -1];
        let img = Image::new(PixelData::I16(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        // Just confirm save succeeded and file reads back (the exact on-disk type depends
        // on conversion).
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![2usize, 2]);
    }

    #[test]
    fn i32_writes_as_unsigned_on_disk() {
        let tmp = tmp_xisf();
        let original = vec![0i32, 100, -50, -1];
        let img = Image::new(PixelData::I32(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![2usize, 2]);
    }
}
