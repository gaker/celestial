use super::Image;
use crate::core::{BitPix, ImageError, Result};
use crate::formats::pixel_data::PixelData;
use crate::fits::header::Keyword;
use crate::fits::io::writer::FitsWriter;
use std::path::Path;

impl Image {
    pub(super) fn open_fits(path: &Path) -> Result<Self> {
        let mut fits = crate::fits::FitsFile::open(path).map_err(ImageError::Fits)?;
        let (dimensions, bitpix) = fits.get_image_info(0).map_err(ImageError::Fits)?;
        let header = fits.get_header(0).map_err(ImageError::Fits)?;
        let is_unsigned_16 = is_u16_encoding(&header);

        let pixels = match bitpix {
            BitPix::U8 => {
                let (_, data): (_, Vec<u8>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::U8(data)
            }
            BitPix::I16 if is_unsigned_16 => {
                let (_, data): (_, Vec<i16>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                let u16_data: Vec<u16> =
                    data.iter().map(|&v| (v as i32 + 32768) as u16).collect();
                PixelData::U16(u16_data)
            }
            BitPix::I16 => {
                let (_, data): (_, Vec<i16>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::I16(data)
            }
            BitPix::I32 => {
                let (_, data): (_, Vec<i32>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::I32(data)
            }
            BitPix::I64 => {
                let (_, data): (_, Vec<i32>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::I32(data)
            }
            BitPix::F32 => {
                let (_, data): (_, Vec<f32>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::F32(data)
            }
            BitPix::F64 => {
                let (_, data): (_, Vec<f64>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::F64(data)
            }
        };

        let keywords = header.keywords().to_vec();

        Ok(Self {
            pixels,
            dimensions,
            keywords,
            xisf_properties: Vec::new(),
        })
    }

    pub(super) fn save_fits(&self, path: &Path) -> Result<()> {
        let mut writer = FitsWriter::create(path).map_err(ImageError::Fits)?;

        match &self.pixels {
            PixelData::U8(data) => writer
                .write_primary_image(data, &self.dimensions, &self.keywords)
                .map_err(ImageError::Fits)?,
            PixelData::U16(data) => {
                let shifted: Vec<i16> =
                    data.iter().map(|&v| (v as i32 - 32768) as i16).collect();
                let keywords = with_u16_scaling(&self.keywords);
                writer
                    .write_primary_image(&shifted, &self.dimensions, &keywords)
                    .map_err(ImageError::Fits)?
            }
            PixelData::I16(data) => writer
                .write_primary_image(data, &self.dimensions, &self.keywords)
                .map_err(ImageError::Fits)?,
            PixelData::I32(data) => writer
                .write_primary_image(data, &self.dimensions, &self.keywords)
                .map_err(ImageError::Fits)?,
            PixelData::F32(data) => writer
                .write_primary_image(data, &self.dimensions, &self.keywords)
                .map_err(ImageError::Fits)?,
            PixelData::F64(data) => writer
                .write_primary_image(data, &self.dimensions, &self.keywords)
                .map_err(ImageError::Fits)?,
        }

        Ok(())
    }
}

fn is_u16_encoding(header: &crate::fits::header::Header) -> bool {
    let bzero = header
        .get_keyword_value("BZERO")
        .and_then(|v| v.as_real())
        .unwrap_or(0.0);
    let bscale = header
        .get_keyword_value("BSCALE")
        .and_then(|v| v.as_real())
        .unwrap_or(1.0);
    bzero == 32768.0 && bscale == 1.0
}

fn with_u16_scaling(keywords: &[Keyword]) -> Vec<Keyword> {
    let mut out: Vec<Keyword> = keywords
        .iter()
        .filter(|k| k.name != "BZERO" && k.name != "BSCALE")
        .cloned()
        .collect();
    out.push(Keyword::real("BZERO", 32768.0).with_comment("offset data range to that of unsigned short"));
    out.push(Keyword::real("BSCALE", 1.0).with_comment("default scaling factor"));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::Keyword;

    fn tmp_fits() -> tempfile::NamedTempFile {
        tempfile::Builder::new()
            .suffix(".fits")
            .tempfile()
            .unwrap()
    }

    #[test]
    fn roundtrip_u8_preserves_data_and_dimensions() {
        let tmp = tmp_fits();
        let original: Vec<u8> = (0u8..100).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![10usize, 10]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_i16() {
        let tmp = tmp_fits();
        let original: Vec<i16> = (-50i16..50).collect();
        let img = Image::new(PixelData::I16(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i16().unwrap(), &original);
    }

    #[test]
    fn roundtrip_i32() {
        let tmp = tmp_fits();
        let original: Vec<i32> = (0..100).map(|i| i * 1000).collect();
        let img = Image::new(PixelData::I32(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i32().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f32() {
        let tmp = tmp_fits();
        let original: Vec<f32> = (0..100).map(|i| i as f32 * 0.1).collect();
        let img = Image::new(PixelData::F32(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f64() {
        let tmp = tmp_fits();
        let original: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
        let img = Image::new(PixelData::F64(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f64().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u16_uses_bzero_bscale() {
        let tmp = tmp_fits();
        let original: Vec<u16> = vec![0, 100, 32767, 32768, 40000, 65535];
        let img = Image::new(PixelData::U16(original.clone()), vec![3usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);

        let bzero = restored.get_keyword("BZERO").unwrap().value.clone().unwrap();
        let bscale = restored.get_keyword("BSCALE").unwrap().value.clone().unwrap();
        assert_eq!(bzero.as_real(), Some(32768.0));
        assert_eq!(bscale.as_real(), Some(1.0));
    }

    #[test]
    fn keywords_roundtrip() {
        let tmp = tmp_fits();
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.set_keyword(Keyword::string("OBJECT", "M31"));
        img.set_keyword(Keyword::real("EXPTIME", 30.5));
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.get_keyword("OBJECT").unwrap().value,
            Some(crate::fits::header::KeywordValue::String("M31".to_string()))
        );
        assert_eq!(
            restored.get_keyword("EXPTIME").unwrap().value,
            Some(crate::fits::header::KeywordValue::Real(30.5))
        );
    }
}
