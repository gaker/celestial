use crate::core::ImageError;
use crate::core::{BitPix, Result};
use crate::debayer::{debayer_bilinear_u16, debayer_bilinear_u8, BayerPattern};
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::DataArray;
use crate::fits::header::Keyword;
use crate::fits::io::writer::FitsWriter;
use crate::xisf::writer::{XisfDataType, XisfWriter};
use celestial_wcs::{Wcs, WcsKeyword, WcsKeywordValue};
use std::io::{Read, Seek};
use std::path::Path;

#[cfg(feature = "standard-formats")]
use {std::io::BufReader, tiff::encoder::colortype, tiff::encoder::TiffEncoder};

#[derive(Debug, Clone)]
pub struct TelescopeParams {
    pub exposure_time: f64,
    pub temperature: Option<f64>,
    pub gain: Option<f64>,
    pub binning: Option<(u32, u32)>,
    pub filter: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageFormat {
    Fits,
    Xisf,
    #[cfg(feature = "standard-formats")]
    Png,
    #[cfg(feature = "standard-formats")]
    Tiff,
}

impl ImageFormat {
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_lowercase().as_str() {
            "fits" | "fit" | "fts" => Some(Self::Fits),
            "xisf" => Some(Self::Xisf),
            #[cfg(feature = "standard-formats")]
            "png" => Some(Self::Png),
            #[cfg(feature = "standard-formats")]
            "tiff" | "tif" => Some(Self::Tiff),
            _ => None,
        }
    }

    pub fn from_magic_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.starts_with(b"SIMPLE  ") {
            Some(Self::Fits)
        } else if bytes.starts_with(b"<?xml") || bytes.starts_with(b"XISF") {
            Some(Self::Xisf)
        } else {
            #[cfg(feature = "standard-formats")]
            {
                if bytes.starts_with(b"\x89PNG") {
                    return Some(Self::Png);
                }
                if bytes.starts_with(b"II*\0") || bytes.starts_with(b"MM\0*") {
                    return Some(Self::Tiff);
                }
            }
            None
        }
    }

    pub fn detect<R: Read + Seek>(reader: &mut R) -> Result<Self> {
        use crate::core::ImageError;

        let mut magic_bytes = [0u8; 16];
        reader.read_exact(&mut magic_bytes)?;
        reader.seek(std::io::SeekFrom::Start(0))?;

        Self::from_magic_bytes(&magic_bytes)
            .ok_or_else(|| ImageError::FormatDetectionFailed("Unknown magic bytes".to_string()))
    }

    pub fn extension(&self) -> &'static str {
        match self {
            Self::Fits => "fits",
            Self::Xisf => "xisf",
            #[cfg(feature = "standard-formats")]
            Self::Png => "png",
            #[cfg(feature = "standard-formats")]
            Self::Tiff => "tiff",
        }
    }
}

const DEFAULT_TILE_SIZE: usize = 32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageKind {
    Mono,
    Rgb,
    Cube,
}

impl ImageKind {
    pub fn from_dimensions(dims: &[usize]) -> Self {
        match dims.len() {
            1 | 2 => Self::Mono,
            3 if dims[2] == 3 => Self::Rgb,
            _ => Self::Cube,
        }
    }
}

/// Pixel data storage for different bit depths.
#[derive(Debug, Clone)]
pub enum PixelData {
    U8(Vec<u8>),
    U16(Vec<u16>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    F32(Vec<f32>),
    F64(Vec<f64>),
}

impl PixelData {
    pub fn len(&self) -> usize {
        match self {
            Self::U8(v) => v.len(),
            Self::U16(v) => v.len(),
            Self::I16(v) => v.len(),
            Self::I32(v) => v.len(),
            Self::F32(v) => v.len(),
            Self::F64(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn bitpix(&self) -> BitPix {
        match self {
            Self::U8(_) => BitPix::U8,
            Self::U16(_) => BitPix::I16,
            Self::I16(_) => BitPix::I16,
            Self::I32(_) => BitPix::I32,
            Self::F32(_) => BitPix::F32,
            Self::F64(_) => BitPix::F64,
        }
    }

    pub fn as_u8(&self) -> Option<&Vec<u8>> {
        match self {
            Self::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u8_mut(&mut self) -> Option<&mut Vec<u8>> {
        match self {
            Self::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u16(&self) -> Option<&Vec<u16>> {
        match self {
            Self::U16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u16_mut(&mut self) -> Option<&mut Vec<u16>> {
        match self {
            Self::U16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i16(&self) -> Option<&Vec<i16>> {
        match self {
            Self::I16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i16_mut(&mut self) -> Option<&mut Vec<i16>> {
        match self {
            Self::I16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i32(&self) -> Option<&Vec<i32>> {
        match self {
            Self::I32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i32_mut(&mut self) -> Option<&mut Vec<i32>> {
        match self {
            Self::I32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f32(&self) -> Option<&Vec<f32>> {
        match self {
            Self::F32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f32_mut(&mut self) -> Option<&mut Vec<f32>> {
        match self {
            Self::F32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<&Vec<f64>> {
        match self {
            Self::F64(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f64_mut(&mut self) -> Option<&mut Vec<f64>> {
        match self {
            Self::F64(v) => Some(v),
            _ => None,
        }
    }

    /// Convert to f32, normalizing integer types to 0.0-1.0 range.
    pub fn to_f32_normalized(&self) -> Vec<f32> {
        match self {
            Self::U8(v) => v.iter().map(|&x| x as f32 / 255.0).collect(),
            Self::U16(v) => v.iter().map(|&x| x as f32 / 65535.0).collect(),
            Self::I16(v) => v.iter().map(|&x| (x as f32 + 32768.0) / 65535.0).collect(),
            Self::I32(v) => v
                .iter()
                .map(|&x| (x as f64 + 2147483648.0) as f32 / 4294967295.0)
                .collect(),
            Self::F32(v) => v.clone(),
            Self::F64(v) => v.iter().map(|&x| x as f32).collect(),
        }
    }

    /// Convert in-place to f32 (normalized for integer types).
    pub fn convert_to_f32(&mut self) {
        let converted = self.to_f32_normalized();
        *self = PixelData::F32(converted);
    }
}

impl From<&[u8]> for PixelData {
    fn from(data: &[u8]) -> Self {
        Self::U8(data.to_vec())
    }
}

impl From<&[u16]> for PixelData {
    fn from(data: &[u16]) -> Self {
        Self::U16(data.to_vec())
    }
}

impl From<&[i16]> for PixelData {
    fn from(data: &[i16]) -> Self {
        Self::I16(data.to_vec())
    }
}

impl From<&[i32]> for PixelData {
    fn from(data: &[i32]) -> Self {
        Self::I32(data.to_vec())
    }
}

impl From<&[f32]> for PixelData {
    fn from(data: &[f32]) -> Self {
        Self::F32(data.to_vec())
    }
}

impl From<&[f64]> for PixelData {
    fn from(data: &[f64]) -> Self {
        Self::F64(data.to_vec())
    }
}

impl From<Vec<u8>> for PixelData {
    fn from(data: Vec<u8>) -> Self {
        Self::U8(data)
    }
}

impl From<Vec<u16>> for PixelData {
    fn from(data: Vec<u16>) -> Self {
        Self::U16(data)
    }
}

impl From<Vec<i16>> for PixelData {
    fn from(data: Vec<i16>) -> Self {
        Self::I16(data)
    }
}

impl From<Vec<i32>> for PixelData {
    fn from(data: Vec<i32>) -> Self {
        Self::I32(data)
    }
}

impl From<Vec<f32>> for PixelData {
    fn from(data: Vec<f32>) -> Self {
        Self::F32(data)
    }
}

impl From<Vec<f64>> for PixelData {
    fn from(data: Vec<f64>) -> Self {
        Self::F64(data)
    }
}

/// A fully-loaded, mutable astronomical image.
///
/// `Image` owns all pixel data and metadata in memory, allowing full
/// read/write access. Load from any supported format and save to any format.
///
/// # Example
/// ```ignore
/// let mut img = Image::open("input.fits")?;
///
/// // Modify pixels
/// if let Some(pixels) = img.pixels.as_f32_mut() {
///     for p in pixels.iter_mut() {
///         *p = (*p * 1.5).min(1.0);
///     }
/// }
///
/// // Modify metadata
/// img.set_keyword(Keyword::string("OBJECT", "M31"));
///
/// // Save to different format
/// img.save("output.xisf")?;
/// ```
#[derive(Debug, Clone)]
pub struct Image {
    pub pixels: PixelData,
    pub dimensions: Vec<usize>,
    pub keywords: Vec<Keyword>,
}

impl Image {
    /// Create a new image from pixel data.
    pub fn new(pixels: PixelData, dimensions: impl Into<Vec<usize>>) -> Self {
        Self {
            pixels,
            dimensions: dimensions.into(),
            keywords: Vec::new(),
        }
    }

    /// Open an image from a file (FITS or XISF).
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let extension = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("");

        let format = match ImageFormat::from_extension(extension) {
            Some(f) => f,
            None => {
                let mut file = std::fs::File::open(path_ref)?;
                ImageFormat::detect(&mut file)?
            }
        };

        match format {
            ImageFormat::Fits => Self::open_fits(path_ref),
            ImageFormat::Xisf => Self::open_xisf(path_ref),
            #[cfg(feature = "standard-formats")]
            ImageFormat::Png => Self::open_png(path_ref),
            #[cfg(feature = "standard-formats")]
            ImageFormat::Tiff => Self::open_tiff(path_ref),
        }
    }

    fn open_fits(path: &Path) -> Result<Self> {
        let mut fits = crate::fits::FitsFile::open(path).map_err(ImageError::Fits)?;
        let (dimensions, bitpix) = fits.get_image_info(0).map_err(ImageError::Fits)?;

        let pixels = match bitpix {
            BitPix::U8 => {
                let (_, data): (_, Vec<u8>) =
                    fits.primary_hdu_with_data().map_err(ImageError::Fits)?;
                PixelData::U8(data)
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

        let header = fits.get_header(0).map_err(ImageError::Fits)?;
        let keywords = header.keywords().to_vec();

        Ok(Self {
            pixels,
            dimensions,
            keywords,
        })
    }

    fn open_xisf(path: &Path) -> Result<Self> {
        let mut xisf = crate::xisf::XisfFile::open(path).map_err(ImageError::Xisf)?;

        let info = xisf.image_info(0).ok_or_else(|| {
            ImageError::Xisf(crate::xisf::XisfError::InvalidFormat(
                "No images in file".to_string(),
            ))
        })?;

        let dimensions = info.geometry.clone();
        let sample_format = info.sample_format.clone();

        // Read raw bytes and convert based on sample format
        let raw_bytes = xisf.read_image_data_raw(0).map_err(ImageError::Xisf)?;

        let pixels = match sample_format {
            crate::xisf::SampleFormat::UInt8 => PixelData::U8(raw_bytes),
            crate::xisf::SampleFormat::UInt16 => {
                let data: Vec<u16> = raw_bytes
                    .chunks_exact(2)
                    .map(|b| u16::from_le_bytes([b[0], b[1]]))
                    .collect();
                PixelData::U16(data)
            }
            crate::xisf::SampleFormat::UInt32 => {
                let data: Vec<i32> = raw_bytes
                    .chunks_exact(4)
                    .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                    .collect();
                PixelData::I32(data)
            }
            crate::xisf::SampleFormat::Float32 => {
                let data: Vec<f32> = raw_bytes
                    .chunks_exact(4)
                    .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                    .collect();
                PixelData::F32(data)
            }
            crate::xisf::SampleFormat::Float64 => {
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
        })
    }

    #[cfg(feature = "standard-formats")]
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

        let (pixels, dimensions) = Self::parse_png_data(&buf, &info)?;
        Ok(Self {
            pixels,
            dimensions,
            keywords: Vec::new(),
        })
    }

    #[cfg(feature = "standard-formats")]
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

    #[cfg(feature = "standard-formats")]
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

        let (pixels, channels) = Self::decode_tiff_image(image)?;
        let dimensions = Self::build_tiff_dimensions(width, height, channels);

        Ok(Self {
            pixels,
            dimensions,
            keywords: Vec::new(),
        })
    }

    #[cfg(feature = "standard-formats")]
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

    #[cfg(feature = "standard-formats")]
    fn build_tiff_dimensions(width: u32, height: u32, channels: usize) -> Vec<usize> {
        if channels == 1 {
            vec![width as usize, height as usize]
        } else {
            vec![width as usize, height as usize, channels]
        }
    }

    /// Save the image to a file. Format is determined by extension.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path_ref = path.as_ref();
        let extension = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("");

        match ImageFormat::from_extension(extension) {
            Some(ImageFormat::Fits) => self.save_fits(path_ref),
            Some(ImageFormat::Xisf) => self.save_xisf(path_ref),
            #[cfg(feature = "standard-formats")]
            Some(ImageFormat::Png) => self.save_png(path_ref),
            #[cfg(feature = "standard-formats")]
            Some(ImageFormat::Tiff) => self.save_tiff(path_ref),
            None => Err(ImageError::UnsupportedFormat),
        }
    }

    fn save_fits(&self, path: &Path) -> Result<()> {
        let mut writer = FitsWriter::create(path).map_err(ImageError::Fits)?;

        match &self.pixels {
            PixelData::U8(data) => writer
                .write_primary_image(data, &self.dimensions, &self.keywords)
                .map_err(ImageError::Fits)?,
            PixelData::U16(data) => {
                // FITS doesn't have u16, convert to i16 with offset
                let converted: Vec<i16> = data.iter().map(|&v| v as i16).collect();
                writer
                    .write_primary_image(&converted, &self.dimensions, &self.keywords)
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

    fn save_xisf(&self, path: &Path) -> Result<()> {
        use crate::xisf::wcs_to_xisf_properties;

        let mut writer = XisfWriter::create(path).map_err(ImageError::Xisf)?;

        let (width, height, channels) = self.extract_dimensions();

        // XISF expects planar format for RGB images - convert from interleaved
        if channels == 3 {
            let mut planar_image = self.clone();
            planar_image.interleaved_to_planar();
            Self::write_xisf_pixels(&mut writer, &planar_image.pixels, width, height, channels)?;
        } else {
            Self::write_xisf_pixels(&mut writer, &self.pixels, width, height, channels)?;
        }

        for kw in &self.keywords {
            writer.add_keyword(kw.clone());
        }

        // Convert WCS keywords to native XISF properties for PixInsight compatibility
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

        writer.write().map_err(ImageError::Xisf)
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

    #[cfg(feature = "standard-formats")]
    fn save_png(&self, path: &Path) -> Result<()> {
        let (width, height, channels) = self.extract_dimensions();
        let color_type = Self::channels_to_png_color_type(channels)?;
        let file = std::fs::File::create(path)?;

        match &self.pixels {
            PixelData::U8(data) => Self::write_png_u8(file, width, height, color_type, data),
            PixelData::U16(data) => Self::write_png_u16(file, width, height, color_type, data),
            _ => Err(ImageError::UnsupportedFormat),
        }
    }

    #[cfg(feature = "standard-formats")]
    fn channels_to_png_color_type(channels: usize) -> Result<png::ColorType> {
        match channels {
            1 => Ok(png::ColorType::Grayscale),
            2 => Ok(png::ColorType::GrayscaleAlpha),
            3 => Ok(png::ColorType::Rgb),
            4 => Ok(png::ColorType::Rgba),
            _ => Err(ImageError::UnsupportedFormat),
        }
    }

    #[cfg(feature = "standard-formats")]
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

    #[cfg(feature = "standard-formats")]
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

    #[cfg(feature = "standard-formats")]
    fn save_tiff(&self, path: &Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut encoder = TiffEncoder::new(file)
            .map_err(|e| ImageError::FormatDetectionFailed(format!("TIFF encoder error: {}", e)))?;

        let (width, height, channels) = self.extract_dimensions();
        self.write_tiff_image(&mut encoder, width as u32, height as u32, channels)
    }

    #[cfg(feature = "standard-formats")]
    fn write_tiff_image(
        &self,
        encoder: &mut tiff::encoder::TiffEncoder<std::fs::File>,
        width: u32,
        height: u32,
        channels: usize,
    ) -> Result<()> {
        match (&self.pixels, channels) {
            (PixelData::U8(data), 1) => {
                encoder.write_image::<colortype::Gray8>(width, height, data)
            }
            (PixelData::U8(data), 3) => encoder.write_image::<colortype::RGB8>(width, height, data),
            (PixelData::U16(data), 1) => {
                encoder.write_image::<colortype::Gray16>(width, height, data)
            }
            (PixelData::U16(data), 3) => {
                encoder.write_image::<colortype::RGB16>(width, height, data)
            }
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

    fn extract_dimensions(&self) -> (usize, usize, usize) {
        let width = self.dimensions.first().copied().unwrap_or(1);
        let height = self.dimensions.get(1).copied().unwrap_or(1);
        let channels = self.dimensions.get(2).copied().unwrap_or(1);
        (width, height, channels)
    }

    // Dimension accessors
    pub fn width(&self) -> usize {
        self.dimensions.first().copied().unwrap_or(0)
    }

    pub fn height(&self) -> usize {
        self.dimensions.get(1).copied().unwrap_or(1)
    }

    pub fn channels(&self) -> usize {
        self.dimensions.get(2).copied().unwrap_or(1)
    }

    pub fn is_rgb(&self) -> bool {
        self.channels() == 3
    }

    pub fn kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    // Keyword helpers
    pub fn get_keyword(&self, name: &str) -> Option<&Keyword> {
        self.keywords.iter().find(|k| k.name == name)
    }

    pub fn set_keyword(&mut self, kw: Keyword) {
        if let Some(existing) = self.keywords.iter_mut().find(|k| k.name == kw.name) {
            *existing = kw;
        } else {
            self.keywords.push(kw);
        }
    }

    pub fn remove_keyword(&mut self, name: &str) {
        self.keywords.retain(|k| k.name != name);
    }

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

/// Unified astronomical image builder that can write to FITS or XISF format.
///
/// # Example
/// ```ignore
/// let data: Vec<u16> = capture_image();
/// AstroImage::new(&data, [1920, 1080])
///     .wcs(&wcs)
///     .keyword(Keyword::string("OBJECT", "M31"))
///     .keyword(Keyword::real("EXPTIME", 30.0))
///     .write_fits("output.fits")?;
/// ```
pub struct AstroImage<'a, T> {
    data: &'a [T],
    dimensions: Vec<usize>,
    keywords: Vec<Keyword>,
    wcs: Option<&'a Wcs>,
    compressed: bool,
    tile_size: Option<(usize, usize)>,
}

impl<'a, T> AstroImage<'a, T>
where
    T: DataArray + XisfDataType + Clone,
{
    pub fn new(data: &'a [T], dimensions: impl Into<Vec<usize>>) -> Self {
        Self {
            data,
            dimensions: dimensions.into(),
            keywords: Vec::new(),
            wcs: None,
            compressed: true,
            tile_size: None,
        }
    }

    pub fn wcs(mut self, wcs: &'a Wcs) -> Self {
        self.wcs = Some(wcs);
        self
    }

    pub fn keyword(mut self, kw: Keyword) -> Self {
        self.keywords.push(kw);
        self
    }

    pub fn keywords(mut self, keywords: impl IntoIterator<Item = Keyword>) -> Self {
        self.keywords.extend(keywords);
        self
    }

    pub fn compressed(mut self, compress: bool) -> Self {
        self.compressed = compress;
        self
    }

    pub fn tile_size(mut self, width: usize, height: usize) -> Self {
        self.tile_size = Some((width, height));
        self
    }

    pub fn image_kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    /// Write to FITS format.
    pub fn write_fits<P: AsRef<Path>>(&self, path: P) -> crate::fits::Result<()> {
        let all_keywords = self.build_keywords();
        let mut writer = FitsWriter::create(path)?;

        if self.compressed {
            let tile_size = self.compute_tile_size();
            writer.write_compressed_image(
                self.data,
                &self.dimensions,
                tile_size,
                CompressionAlgorithm::Rice,
                &all_keywords,
            )
        } else {
            writer.write_primary_image_with_checksum(self.data, &self.dimensions, &all_keywords)
        }
    }

    /// Write to XISF format.
    pub fn write_xisf<P: AsRef<Path>>(&self, path: P) -> crate::xisf::Result<()> {
        let mut writer = XisfWriter::create(path)?;

        let (width, height, channels) = self.extract_dimensions();
        writer.add_image(self.data, width, height, channels)?;

        for kw in self.build_keywords() {
            writer.add_keyword(kw);
        }

        writer.write()
    }

    /// Write to the format determined by file extension.
    pub fn write_to<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use crate::core::ImageError;

        let path_ref = path.as_ref();
        let ext = path_ref.extension().and_then(|e| e.to_str()).unwrap_or("");

        match ImageFormat::from_extension(ext) {
            Some(ImageFormat::Fits) => self.write_fits(path).map_err(ImageError::Fits),
            Some(ImageFormat::Xisf) => self.write_xisf(path).map_err(ImageError::Xisf),
            #[cfg(feature = "standard-formats")]
            Some(ImageFormat::Png) | Some(ImageFormat::Tiff) => Err(ImageError::UnsupportedFormat),
            None => Err(ImageError::FormatDetectionFailed(format!(
                "Unknown extension: {}",
                ext
            ))),
        }
    }

    fn extract_dimensions(&self) -> (usize, usize, usize) {
        let width = self.dimensions.first().copied().unwrap_or(1);
        let height = self.dimensions.get(1).copied().unwrap_or(1);
        let channels = self.dimensions.get(2).copied().unwrap_or(1);
        (width, height, channels)
    }

    fn build_keywords(&self) -> Vec<Keyword> {
        let mut keywords = Vec::new();

        if let Some(wcs) = self.wcs {
            keywords.extend(wcs_to_keywords(wcs));
        }

        keywords.extend(self.keywords.iter().cloned());

        keywords
    }

    fn compute_tile_size(&self) -> (usize, usize) {
        if let Some(size) = self.tile_size {
            return size;
        }

        let width = self.dimensions.first().copied().unwrap_or(1);
        let height = self.dimensions.get(1).copied().unwrap_or(1);

        let tile_w = width.min(DEFAULT_TILE_SIZE);
        let tile_h = height.min(DEFAULT_TILE_SIZE);

        (tile_w, tile_h)
    }
}

fn wcs_to_keywords(wcs: &Wcs) -> Vec<Keyword> {
    wcs.to_keywords()
        .into_iter()
        .map(wcs_keyword_to_keyword)
        .collect()
}

fn wcs_keyword_to_keyword(wk: WcsKeyword) -> Keyword {
    match wk.value {
        WcsKeywordValue::Real(v) => Keyword::real(wk.name, v),
        WcsKeywordValue::Integer(v) => Keyword::integer(wk.name, v),
        WcsKeywordValue::String(v) => Keyword::string(wk.name, v),
    }
}

#[derive(Debug, Clone)]
pub struct ImageInfo {
    pub dimensions: Vec<usize>,
    pub bitpix: BitPix,
    pub is_signed: bool,
    pub bytes_per_pixel: usize,
    pub total_pixels: usize,
}

impl ImageInfo {
    pub fn new(dimensions: Vec<usize>, bitpix: BitPix) -> Self {
        let total_pixels = if dimensions.is_empty() {
            0
        } else {
            dimensions.iter().product()
        };
        let bytes_per_pixel = bitpix.bytes_per_pixel();
        let is_signed = matches!(bitpix, BitPix::I16 | BitPix::I32 | BitPix::I64);

        Self {
            dimensions,
            bitpix,
            is_signed,
            bytes_per_pixel,
            total_pixels,
        }
    }

    pub fn data_size_bytes(&self) -> usize {
        self.total_pixels * self.bytes_per_pixel
    }

    pub fn is_2d(&self) -> bool {
        self.dimensions.len() == 2
    }

    pub fn width(&self) -> Option<usize> {
        self.dimensions.first().copied()
    }

    pub fn height(&self) -> Option<usize> {
        self.dimensions.get(1).copied()
    }
}

pub struct ImageWriter {
    format: ImageFormat,
    path: std::path::PathBuf,
}

impl ImageWriter {
    pub fn new<P: AsRef<Path>>(path: P, format: ImageFormat) -> Self {
        Self {
            format,
            path: path.as_ref().to_path_buf(),
        }
    }

    pub fn write_image<T>(&self, data: &[T], info: &ImageInfo, keywords: &[Keyword]) -> Result<()>
    where
        T: crate::fits::data::array::DataArray + Clone,
    {
        match self.format {
            ImageFormat::Fits => {
                let mut writer = crate::fits::FitsWriter::create(&self.path)
                    .map_err(crate::core::ImageError::Fits)?;

                writer
                    .write_primary_image(data, &info.dimensions, keywords)
                    .map_err(crate::core::ImageError::Fits)?;

                Ok(())
            }
            ImageFormat::Xisf => {
                use crate::core::ImageError;
                Err(ImageError::UnsupportedFormat)
            }
            #[cfg(feature = "standard-formats")]
            ImageFormat::Png | ImageFormat::Tiff => {
                use crate::core::ImageError;
                Err(ImageError::UnsupportedFormat)
            }
        }
    }

    pub fn write_telescope_image<T>(
        &self,
        data: &[T],
        info: &ImageInfo,
        params: &TelescopeParams,
    ) -> Result<()>
    where
        T: crate::fits::data::array::DataArray + Clone,
    {
        let keywords = Self::build_telescope_keywords(params);
        self.write_image(data, info, &keywords)
    }

    fn build_telescope_keywords(params: &TelescopeParams) -> Vec<Keyword> {
        let mut keywords = Vec::new();
        keywords.push(Keyword::real("EXPTIME", params.exposure_time));

        if let Some(temp) = params.temperature {
            keywords.push(Keyword::real("CCD-TEMP", temp));
        }
        if let Some(gain) = params.gain {
            keywords.push(Keyword::real("GAIN", gain));
        }
        if let Some((x, y)) = params.binning {
            keywords.push(Keyword::integer("XBINNING", x as i64));
            keywords.push(Keyword::integer("YBINNING", y as i64));
        }
        if let Some(ref filt) = params.filter {
            keywords.push(Keyword::string("FILTER", filt));
        }

        keywords
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::KeywordValue;
    use crate::test_utils::*;
    use std::io::Cursor;

    #[test]
    fn image_format_from_extension() {
        let fits_extensions = ["fits", "FITS", "fit", "FIT", "fts", "FTS"];
        let xisf_extensions = ["xisf", "XISF"];

        for ext in fits_extensions {
            assert_eq!(ImageFormat::from_extension(ext), Some(ImageFormat::Fits));
        }

        for ext in xisf_extensions {
            assert_eq!(ImageFormat::from_extension(ext), Some(ImageFormat::Xisf));
        }

        assert_eq!(ImageFormat::from_extension("jpg"), None);
        assert_eq!(ImageFormat::from_extension(""), None);
        assert_eq!(ImageFormat::from_extension("unknown"), None);
    }

    #[test]
    fn image_format_from_magic_bytes() {
        let fits_magic = b"SIMPLE  =                    T";
        let xisf_magic = b"XISF0100";

        assert_eq!(
            ImageFormat::from_magic_bytes(fits_magic),
            Some(ImageFormat::Fits)
        );
        assert_eq!(
            ImageFormat::from_magic_bytes(xisf_magic),
            Some(ImageFormat::Xisf)
        );

        assert_eq!(ImageFormat::from_magic_bytes(b"INVALID"), None);
        assert_eq!(ImageFormat::from_magic_bytes(b""), None);
        assert_eq!(ImageFormat::from_magic_bytes(&[0xFF, 0xFE]), None);
    }

    #[test]
    fn image_format_detect_from_reader() {
        let fits_data = create_minimal_fits();
        let mut cursor = Cursor::new(fits_data);
        let format = ImageFormat::detect(&mut cursor).unwrap();
        assert_eq!(format, ImageFormat::Fits);

        let invalid_data = vec![0x00, 0x01, 0x02];
        let mut invalid_cursor = Cursor::new(invalid_data);
        let result = ImageFormat::detect(&mut invalid_cursor);
        assert!(result.is_err());
    }

    #[test]
    fn image_info_creation_and_methods() {
        let dims = vec![1024, 1024];
        let bitpix = BitPix::I16;
        let info = ImageInfo::new(dims.clone(), bitpix);

        assert_eq!(info.dimensions, dims);
        assert_eq!(info.bitpix, bitpix);
        assert_eq!(info.data_size_bytes(), 1024 * 1024 * 2);
        assert!(info.is_2d());
        assert_eq!(info.width(), Some(1024));
        assert_eq!(info.height(), Some(1024));
    }

    #[test]
    fn image_info_edge_cases() {
        let info_1d = ImageInfo::new(vec![1000], BitPix::U8);
        assert!(!info_1d.is_2d());
        assert_eq!(info_1d.width(), Some(1000));
        assert_eq!(info_1d.height(), None);

        let info_3d = ImageInfo::new(vec![10, 10, 10], BitPix::F32);
        assert!(!info_3d.is_2d());
        assert_eq!(info_3d.width(), Some(10));
        assert_eq!(info_3d.height(), Some(10));

        let info_empty = ImageInfo::new(vec![], BitPix::I32);
        assert!(!info_empty.is_2d());
        assert_eq!(info_empty.width(), None);
        assert_eq!(info_empty.height(), None);
        assert_eq!(info_empty.data_size_bytes(), 0);
    }

    #[test]
    fn image_writer_build_telescope_keywords() {
        use crate::fits::header::KeywordValue;

        let params = TelescopeParams {
            exposure_time: 30.0,
            temperature: Some(-15.5),
            gain: Some(100.0),
            binning: Some((2, 2)),
            filter: Some("Ha".to_string()),
        };
        let keywords = ImageWriter::build_telescope_keywords(&params);

        assert_eq!(keywords.len(), 6);
        assert_eq!(keywords[0].name, "EXPTIME");
        assert_eq!(keywords[0].value, Some(KeywordValue::Real(30.0)));

        let params_none = TelescopeParams {
            exposure_time: 60.0,
            temperature: None,
            gain: None,
            binning: None,
            filter: None,
        };
        let keywords_none = ImageWriter::build_telescope_keywords(&params_none);

        assert_eq!(keywords_none.len(), 1);
        assert_eq!(keywords_none[0].name, "EXPTIME");
    }

    #[test]
    fn telescope_keywords_extreme_binning() {
        use crate::fits::header::KeywordValue;

        let extreme_binning_cases: [(u32, u32); 5] = [
            (0, 0),
            (1, 1),
            (u32::MAX, u32::MAX),
            (1, u32::MAX),
            (u32::MAX, 1),
        ];

        for (x, y) in extreme_binning_cases {
            let params = TelescopeParams {
                exposure_time: 1.0,
                temperature: None,
                gain: None,
                binning: Some((x, y)),
                filter: None,
            };
            let keywords = ImageWriter::build_telescope_keywords(&params);

            assert_eq!(keywords.len(), 3);
            assert_eq!(keywords[1].name, "XBINNING");
            assert_eq!(keywords[1].value, Some(KeywordValue::Integer(x as i64)));
            assert_eq!(keywords[2].name, "YBINNING");
            assert_eq!(keywords[2].value, Some(KeywordValue::Integer(y as i64)));
        }
    }

    #[test]
    fn telescope_keywords_filter_names() {
        use crate::fits::header::KeywordValue;

        let filter_names = ["Ha", "V", "R", "test"];

        for filter_name in filter_names {
            let params = TelescopeParams {
                exposure_time: 1.0,
                temperature: None,
                gain: None,
                binning: None,
                filter: Some(filter_name.to_string()),
            };
            let keywords = ImageWriter::build_telescope_keywords(&params);

            assert_eq!(keywords.len(), 2);
            assert_eq!(keywords[1].name, "FILTER");
            assert_eq!(
                keywords[1].value,
                Some(KeywordValue::String(filter_name.to_string()))
            );
        }
    }

    #[test]
    fn image_format_detect_error_cases() {
        let empty_data: Vec<u8> = vec![];
        let mut cursor = Cursor::new(empty_data);
        let result = ImageFormat::detect(&mut cursor);
        assert!(result.is_err());

        let short_data = vec![0x00, 0x01];
        let mut short_cursor = Cursor::new(short_data);
        let result = ImageFormat::detect(&mut short_cursor);
        assert!(result.is_err());
    }

    #[test]
    fn image_format_extension() {
        assert_eq!(ImageFormat::Fits.extension(), "fits");
        assert_eq!(ImageFormat::Xisf.extension(), "xisf");
    }

    #[test]
    fn image_format_from_magic_bytes_xml() {
        let xml_magic = b"<?xml version=\"1.0\"?>";
        assert_eq!(
            ImageFormat::from_magic_bytes(xml_magic),
            Some(ImageFormat::Xisf)
        );
    }

    #[test]
    fn image_writer_new() {
        use std::path::Path;
        let path = Path::new("/tmp/test.fits");
        let writer = ImageWriter::new(path, ImageFormat::Fits);
        assert_eq!(writer.format, ImageFormat::Fits);
        assert_eq!(writer.path, path);
    }

    #[test]
    fn image_writer_write_telescope_image() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let path = dir.path().join("telescope.fits");

        let writer = ImageWriter::new(&path, ImageFormat::Fits);
        let data = vec![100i16, 200, 300, 400];
        let info = ImageInfo::new(vec![2, 2], BitPix::I16);
        let params = TelescopeParams {
            exposure_time: 60.0,
            temperature: Some(-20.0),
            gain: Some(200.0),
            binning: Some((1, 1)),
            filter: Some("R".to_string()),
        };

        let result = writer.write_telescope_image(&data, &info, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn image_writer_write_image_xisf_unsupported() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.xisf");

        let writer = ImageWriter::new(&path, ImageFormat::Xisf);
        let data = vec![1u8, 2, 3, 4];
        let info = ImageInfo::new(vec![2, 2], BitPix::U8);

        let result = writer.write_image(&data, &info, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn telescope_params() {
        let params = TelescopeParams {
            exposure_time: 30.0,
            temperature: Some(-15.0),
            gain: None,
            binning: Some((2, 2)),
            filter: Some("Ha".to_string()),
        };

        let debug_str = format!("{:?}", params);
        assert!(debug_str.contains("TelescopeParams"));
        assert!(debug_str.contains("exposure_time: 30.0"));
    }

    #[test]
    fn image_format_and_partial_eq() {
        assert_eq!(ImageFormat::Fits, ImageFormat::Fits);
        assert_ne!(ImageFormat::Fits, ImageFormat::Xisf);

        let debug_str = format!("{:?}", ImageFormat::Fits);
        assert!(debug_str.contains("Fits"));
    }

    #[test]
    fn image_info_and_clone() {
        let info = ImageInfo::new(vec![100, 100], BitPix::F32);
        let cloned = info.clone();

        assert_eq!(info.dimensions, cloned.dimensions);
        assert_eq!(info.bitpix, cloned.bitpix);

        let debug_str = format!("{:?}", info);
        assert!(debug_str.contains("ImageInfo"));
    }

    #[test]
    fn image_info_signed_types() {
        let info_i16 = ImageInfo::new(vec![10], BitPix::I16);
        assert!(info_i16.is_signed);

        let info_i32 = ImageInfo::new(vec![10], BitPix::I32);
        assert!(info_i32.is_signed);

        let info_i64 = ImageInfo::new(vec![10], BitPix::I64);
        assert!(info_i64.is_signed);

        let info_u8 = ImageInfo::new(vec![10], BitPix::U8);
        assert!(!info_u8.is_signed);

        let info_f32 = ImageInfo::new(vec![10], BitPix::F32);
        assert!(!info_f32.is_signed);

        let info_f64 = ImageInfo::new(vec![10], BitPix::F64);
        assert!(!info_f64.is_signed);
    }

    #[test]
    fn image_writer_write_fits_with_keywords() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let path = dir.path().join("test_keywords.fits");

        let writer = ImageWriter::new(&path, ImageFormat::Fits);
        let data = vec![1u8, 2, 3, 4];
        let info = ImageInfo::new(vec![2, 2], BitPix::U8);
        let keywords = vec![
            Keyword::string("OBJECT", "M31"),
            Keyword::real("EXPTIME", 30.0).with_comment("Exposure time in seconds"),
        ];

        let result = writer.write_image(&data, &info, &keywords);
        assert!(result.is_ok());
    }

    // ==================== AstroImage tests ====================

    // Note: AstroImage requires types that implement both DataArray (FITS) and XisfDataType (XISF).
    // Common types: u8, f32, f64. FITS-only: i16, i32, i64. XISF-only: u16, u32.

    #[test]
    fn astro_image_write_fits() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".fits").unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        let result = AstroImage::new(&data, [10, 10])
            .compressed(false)
            .write_fits(temp_file.path());

        assert!(result.is_ok());
    }

    #[test]
    fn astro_image_write_xisf() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".xisf").unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        let result = AstroImage::new(&data, [10, 10]).write_xisf(temp_file.path());

        assert!(result.is_ok());
    }

    #[test]
    fn astro_image_write_to_by_extension() {
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        // Write to FITS by extension
        let fits_path = dir.path().join("test.fits");
        let result = AstroImage::new(&data, [10, 10])
            .compressed(false)
            .write_to(&fits_path);
        assert!(result.is_ok());

        // Write to XISF by extension
        let xisf_path = dir.path().join("test.xisf");
        let result = AstroImage::new(&data, [10, 10]).write_to(&xisf_path);
        assert!(result.is_ok());
    }

    #[test]
    fn astro_image_with_keywords() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".fits").unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        let result = AstroImage::new(&data, [10, 10])
            .keyword(Keyword::string("OBJECT", "M31"))
            .keyword(Keyword::real("EXPTIME", 30.0).with_comment("Exposure time"))
            .keyword(Keyword::integer("GAIN", 100))
            .keyword(Keyword::logical("PREVIEW", true))
            .compressed(false)
            .write_fits(temp_file.path());

        assert!(result.is_ok());

        // Verify the keywords were written
        let mut fits = crate::fits::FitsFile::open(temp_file.path()).unwrap();
        let header = fits.get_header(0).unwrap();

        assert!(header.get_keyword_value("OBJECT").is_some());
        assert!(header.get_keyword_value("EXPTIME").is_some());
        assert!(header.get_keyword_value("GAIN").is_some());
        assert!(header.get_keyword_value("PREVIEW").is_some());
    }

    #[test]
    fn astro_image_with_wcs() {
        use celestial_wcs::{Projection, WcsBuilder};
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".fits").unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        let wcs = WcsBuilder::new()
            .crpix(5.0, 5.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build()
            .unwrap();

        let result = AstroImage::new(&data, [10, 10])
            .wcs(&wcs)
            .compressed(false)
            .write_fits(temp_file.path());

        assert!(result.is_ok());

        // Verify WCS keywords were written
        let mut fits = crate::fits::FitsFile::open(temp_file.path()).unwrap();
        let header = fits.get_header(0).unwrap();

        assert!(header.get_keyword_value("CTYPE1").is_some());
        assert!(header.get_keyword_value("CRPIX1").is_some());
        assert!(header.get_keyword_value("CRVAL1").is_some());
        assert!(header.get_keyword_value("CD1_1").is_some());
    }

    #[test]
    fn astro_image_xisf_with_keywords() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".xisf").unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        let result = AstroImage::new(&data, [10, 10])
            .keyword(Keyword::string("OBJECT", "M42"))
            .keyword(Keyword::real("EXPTIME", 60.0))
            .write_xisf(temp_file.path());

        assert!(result.is_ok());

        // Verify the keywords were written
        let reader = crate::xisf::XisfFile::open(temp_file.path()).unwrap();
        let keywords = reader.keywords();

        assert!(keywords.iter().any(
            |k| k.name == "OBJECT" && k.value == Some(KeywordValue::String("M42".to_string()))
        ));
        assert!(keywords.iter().any(|k| k.name == "EXPTIME"));
    }

    #[test]
    fn astro_image_rgb() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".xisf").unwrap();
        let data: Vec<f32> = (0..300).map(|i| i as f32).collect();

        let image = AstroImage::new(&data, [10, 10, 3]);
        assert_eq!(image.image_kind(), ImageKind::Rgb);

        let result = image.write_xisf(temp_file.path());
        assert!(result.is_ok());
    }

    #[test]
    fn astro_image_u8() {
        use tempfile::NamedTempFile;

        // u8 works for both FITS and XISF
        let temp_file = NamedTempFile::with_suffix(".fits").unwrap();
        let data: Vec<u8> = (0..100).collect();

        let result = AstroImage::new(&data, [10, 10])
            .compressed(false)
            .write_fits(temp_file.path());
        assert!(result.is_ok());

        let temp_file = NamedTempFile::with_suffix(".xisf").unwrap();
        let result = AstroImage::new(&data, [10, 10]).write_xisf(temp_file.path());
        assert!(result.is_ok());
    }

    #[test]
    fn image_kind_detection() {
        assert_eq!(ImageKind::from_dimensions(&[100]), ImageKind::Mono);
        assert_eq!(ImageKind::from_dimensions(&[100, 100]), ImageKind::Mono);
        assert_eq!(ImageKind::from_dimensions(&[100, 100, 3]), ImageKind::Rgb);
        assert_eq!(ImageKind::from_dimensions(&[100, 100, 4]), ImageKind::Cube);
        assert_eq!(ImageKind::from_dimensions(&[100, 100, 10]), ImageKind::Cube);
    }

    #[test]
    fn keyword_types() {
        let real = Keyword::real("EXPTIME", 30.0);
        assert!(matches!(real.value, Some(KeywordValue::Real(_))));

        let int = Keyword::integer("GAIN", 100);
        assert!(matches!(int.value, Some(KeywordValue::Integer(100))));

        let string = Keyword::string("OBJECT", "M31");
        if let Some(KeywordValue::String(s)) = &string.value {
            assert_eq!(s, "M31");
        } else {
            panic!("Expected String variant");
        }

        let boolean = Keyword::logical("PREVIEW", true);
        assert!(matches!(boolean.value, Some(KeywordValue::Logical(true))));

        let with_comment = Keyword::real("TEMP", -15.0).with_comment("CCD temperature");
        assert_eq!(with_comment.comment, Some("CCD temperature".to_string()));
    }

    #[test]
    fn astro_image_unknown_extension_error() {
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();

        let bad_path = dir.path().join("test.unknown");
        let result = AstroImage::new(&data, [10, 10]).write_to(&bad_path);

        assert!(result.is_err());
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn png_roundtrip_u8_grayscale() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".png").unwrap();
        let original_data: Vec<u8> = (0..100).collect();
        let img = Image::new(PixelData::U8(original_data.clone()), vec![10, 10]);

        img.save(temp_file.path()).unwrap();
        let loaded = Image::open(temp_file.path()).unwrap();

        assert_eq!(loaded.dimensions, vec![10, 10]);
        assert_eq!(loaded.pixels.as_u8().unwrap(), &original_data);
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn png_roundtrip_u16_grayscale() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".png").unwrap();
        let original_data: Vec<u16> = (0..100).map(|i| i * 100).collect();
        let img = Image::new(PixelData::U16(original_data.clone()), vec![10, 10]);

        img.save(temp_file.path()).unwrap();
        let loaded = Image::open(temp_file.path()).unwrap();

        assert_eq!(loaded.dimensions, vec![10, 10]);
        assert_eq!(loaded.pixels.as_u16().unwrap(), &original_data);
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn png_roundtrip_u8_rgb() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".png").unwrap();
        let original_data: Vec<u8> = (0..75).collect();
        let img = Image::new(PixelData::U8(original_data.clone()), vec![5, 5, 3]);

        img.save(temp_file.path()).unwrap();
        let loaded = Image::open(temp_file.path()).unwrap();

        assert_eq!(loaded.dimensions, vec![5, 5, 3]);
        assert_eq!(loaded.pixels.as_u8().unwrap(), &original_data);
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn tiff_roundtrip_u8_grayscale() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".tiff").unwrap();
        let original_data: Vec<u8> = (0..100).collect();
        let img = Image::new(PixelData::U8(original_data.clone()), vec![10, 10]);

        img.save(temp_file.path()).unwrap();
        let loaded = Image::open(temp_file.path()).unwrap();

        assert_eq!(loaded.dimensions, vec![10, 10]);
        assert_eq!(loaded.pixels.as_u8().unwrap(), &original_data);
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn tiff_roundtrip_u16_grayscale() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".tiff").unwrap();
        let original_data: Vec<u16> = (0..100).map(|i| i * 100).collect();
        let img = Image::new(PixelData::U16(original_data.clone()), vec![10, 10]);

        img.save(temp_file.path()).unwrap();
        let loaded = Image::open(temp_file.path()).unwrap();

        assert_eq!(loaded.dimensions, vec![10, 10]);
        assert_eq!(loaded.pixels.as_u16().unwrap(), &original_data);
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn tiff_roundtrip_f32_grayscale() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::with_suffix(".tiff").unwrap();
        let original_data: Vec<f32> = (0..100).map(|i| i as f32 * 0.01).collect();
        let img = Image::new(PixelData::F32(original_data.clone()), vec![10, 10]);

        img.save(temp_file.path()).unwrap();
        let loaded = Image::open(temp_file.path()).unwrap();

        assert_eq!(loaded.dimensions, vec![10, 10]);
        let loaded_data = loaded.pixels.as_f32().unwrap();
        for (a, b) in original_data.iter().zip(loaded_data.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn image_format_png_tiff_extension() {
        assert_eq!(ImageFormat::from_extension("png"), Some(ImageFormat::Png));
        assert_eq!(ImageFormat::from_extension("PNG"), Some(ImageFormat::Png));
        assert_eq!(ImageFormat::from_extension("tiff"), Some(ImageFormat::Tiff));
        assert_eq!(ImageFormat::from_extension("tif"), Some(ImageFormat::Tiff));
        assert_eq!(ImageFormat::from_extension("TIFF"), Some(ImageFormat::Tiff));
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn image_format_png_tiff_magic_bytes() {
        let png_magic = b"\x89PNG\r\n\x1a\n";
        let tiff_le_magic = b"II*\0";
        let tiff_be_magic = b"MM\0*";

        assert_eq!(
            ImageFormat::from_magic_bytes(png_magic),
            Some(ImageFormat::Png)
        );
        assert_eq!(
            ImageFormat::from_magic_bytes(tiff_le_magic),
            Some(ImageFormat::Tiff)
        );
        assert_eq!(
            ImageFormat::from_magic_bytes(tiff_be_magic),
            Some(ImageFormat::Tiff)
        );
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn image_format_png_tiff_extension_method() {
        assert_eq!(ImageFormat::Png.extension(), "png");
        assert_eq!(ImageFormat::Tiff.extension(), "tiff");
    }

    #[test]
    fn interleaved_to_planar_u8() {
        // 2x2 RGB image: [R0,G0,B0, R1,G1,B1, R2,G2,B2, R3,G3,B3]
        let interleaved: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        let mut img = Image::new(PixelData::U8(interleaved), vec![2, 2, 3]);

        img.interleaved_to_planar();

        // Expected planar: [R0,R1,R2,R3, G0,G1,G2,G3, B0,B1,B2,B3]
        let expected: Vec<u8> = vec![1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12];
        assert_eq!(img.pixels.as_u8().unwrap(), &expected);
    }

    #[test]
    fn planar_to_interleaved_u8() {
        // 2x2 RGB planar: [R0,R1,R2,R3, G0,G1,G2,G3, B0,B1,B2,B3]
        let planar: Vec<u8> = vec![1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12];
        let mut img = Image::new(PixelData::U8(planar), vec![2, 2, 3]);

        img.planar_to_interleaved();

        // Expected interleaved: [R0,G0,B0, R1,G1,B1, R2,G2,B2, R3,G3,B3]
        let expected: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        assert_eq!(img.pixels.as_u8().unwrap(), &expected);
    }

    #[test]
    fn interleaved_planar_roundtrip() {
        let original: Vec<u16> = vec![100, 200, 300, 400, 500, 600, 700, 800, 900];
        let mut img = Image::new(PixelData::U16(original.clone()), vec![3, 1, 3]);

        img.interleaved_to_planar();
        img.planar_to_interleaved();

        assert_eq!(img.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn interleaved_to_planar_non_rgb_unchanged() {
        let mono: Vec<u8> = vec![1, 2, 3, 4];
        let mut img = Image::new(PixelData::U8(mono.clone()), vec![2, 2]);

        img.interleaved_to_planar();

        assert_eq!(img.pixels.as_u8().unwrap(), &mono);
    }

    #[test]
    fn planar_to_interleaved_non_rgb_unchanged() {
        let mono: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0];
        let mut img = Image::new(PixelData::F32(mono.clone()), vec![2, 2]);

        img.planar_to_interleaved();

        assert_eq!(img.pixels.as_f32().unwrap(), &mono);
    }

    #[test]
    fn interleaved_to_planar_f32() {
        let interleaved: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut img = Image::new(PixelData::F32(interleaved), vec![2, 1, 3]);

        img.interleaved_to_planar();

        let expected: Vec<f32> = vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0];
        assert_eq!(img.pixels.as_f32().unwrap(), &expected);
    }

    #[test]
    fn interleaved_to_planar_i16() {
        let interleaved: Vec<i16> = vec![-1, -2, -3, 4, 5, 6];
        let mut img = Image::new(PixelData::I16(interleaved), vec![2, 1, 3]);

        img.interleaved_to_planar();

        let expected: Vec<i16> = vec![-1, 4, -2, 5, -3, 6];
        assert_eq!(img.pixels.as_i16().unwrap(), &expected);
    }

    #[test]
    fn normalize_f32() {
        let data: Vec<f32> = vec![10.0, 20.0, 30.0, 40.0];
        let mut img = Image::new(PixelData::F32(data), vec![2, 2]);

        img.normalize();

        let result = img.pixels.as_f32().unwrap();
        assert_eq!(result[0], 0.0);
        assert_eq!(result[3], 1.0);
        assert!((result[1] - 1.0 / 3.0).abs() < 1e-6);
        assert!((result[2] - 2.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn normalize_to_f32_converts_type() {
        let data: Vec<u16> = vec![0, 100, 200, 300];
        let mut img = Image::new(PixelData::U16(data), vec![2, 2]);

        img.normalize_to_f32();

        let result = img.pixels.as_f32().unwrap();
        assert_eq!(result[0], 0.0);
        assert_eq!(result[3], 1.0);
    }

    #[test]
    fn normalize_constant_image() {
        let data: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0];
        let mut img = Image::new(PixelData::F32(data), vec![2, 2]);

        img.normalize();

        let result = img.pixels.as_f32().unwrap();
        for &v in result {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn pixel_range_f32() {
        let data: Vec<f32> = vec![-5.0, 0.0, 10.0, 100.0];
        let img = Image::new(PixelData::F32(data), vec![2, 2]);

        let (min, max) = img.pixel_range();
        assert_eq!(min, -5.0);
        assert_eq!(max, 100.0);
    }

    #[test]
    fn image_debayer_u8() {
        use crate::debayer::BayerPattern;

        let raw: Vec<u8> = vec![100, 50, 60, 200];
        let mut img = Image::new(PixelData::U8(raw), vec![2, 2]);

        img.debayer(BayerPattern::Rggb);

        assert_eq!(img.dimensions, vec![2, 2, 3]);
        assert_eq!(img.channels(), 3);
        assert!(img.is_rgb());

        let pixels = img.pixels.as_u8().unwrap();
        assert_eq!(pixels.len(), 12);
        assert_eq!(pixels[0], 100); // R at (0,0) preserved
        assert_eq!(pixels[11], 200); // B at (1,1) preserved
    }

    #[test]
    fn image_debayer_u16() {
        use crate::debayer::BayerPattern;

        let raw: Vec<u16> = vec![1000, 500, 600, 2000];
        let mut img = Image::new(PixelData::U16(raw), vec![2, 2]);

        img.debayer(BayerPattern::Rggb);

        assert_eq!(img.dimensions, vec![2, 2, 3]);
        let pixels = img.pixels.as_u16().unwrap();
        assert_eq!(pixels[0], 1000); // R at (0,0)
        assert_eq!(pixels[11], 2000); // B at (1,1)
    }

    #[test]
    fn image_debayer_ignores_rgb() {
        use crate::debayer::BayerPattern;

        let rgb: Vec<u8> = vec![1, 2, 3, 4, 5, 6];
        let mut img = Image::new(PixelData::U8(rgb.clone()), vec![2, 1, 3]);

        img.debayer(BayerPattern::Rggb);

        // Should be unchanged
        assert_eq!(img.dimensions, vec![2, 1, 3]);
        assert_eq!(img.pixels.as_u8().unwrap(), &rgb);
    }
}
