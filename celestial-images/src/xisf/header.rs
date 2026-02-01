use crate::core::BitPix;
use crate::fits::header::Keyword;
use crate::xisf::{Result, XisfError};

#[derive(Debug, Clone, PartialEq)]
pub enum XisfPropertyValue {
    String(String),
    Float64(f64),
    Int32(i32),
    F64Vector(Vec<f64>),
    F64Matrix {
        rows: usize,
        cols: usize,
        data: Vec<f64>,
    },
}

impl XisfPropertyValue {
    pub fn type_name(&self) -> &'static str {
        match self {
            Self::String(_) => "String",
            Self::Float64(_) => "Float64",
            Self::Int32(_) => "Int32",
            Self::F64Vector(_) => "F64Vector",
            Self::F64Matrix { .. } => "F64Matrix",
        }
    }

    pub fn format_value(&self) -> String {
        match self {
            Self::String(s) => s.clone(),
            Self::Float64(v) => format!("{}", v),
            Self::Int32(v) => format!("{}", v),
            Self::F64Vector(v) => v
                .iter()
                .map(|x| format!("{}", x))
                .collect::<Vec<_>>()
                .join(" "),
            Self::F64Matrix { data, .. } => data
                .iter()
                .map(|x| format!("{}", x))
                .collect::<Vec<_>>()
                .join(" "),
        }
    }

    pub fn is_scalar(&self) -> bool {
        matches!(self, Self::Float64(_) | Self::Int32(_))
    }

    pub fn needs_data_block(&self) -> bool {
        matches!(self, Self::F64Vector(_) | Self::F64Matrix { .. })
    }

    pub fn to_le_bytes(&self) -> Option<Vec<u8>> {
        match self {
            Self::F64Vector(v) => {
                let mut bytes = Vec::with_capacity(v.len() * 8);
                for &val in v {
                    bytes.extend_from_slice(&val.to_le_bytes());
                }
                Some(bytes)
            }
            Self::F64Matrix { data, .. } => {
                let mut bytes = Vec::with_capacity(data.len() * 8);
                for &val in data {
                    bytes.extend_from_slice(&val.to_le_bytes());
                }
                Some(bytes)
            }
            _ => None,
        }
    }

    pub fn data_size(&self) -> usize {
        match self {
            Self::F64Vector(v) => v.len() * 8,
            Self::F64Matrix { data, .. } => data.len() * 8,
            _ => 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct XisfProperty {
    pub id: String,
    pub value: XisfPropertyValue,
}

impl XisfProperty {
    pub fn new(id: impl Into<String>, value: XisfPropertyValue) -> Self {
        Self {
            id: id.into(),
            value,
        }
    }

    pub fn string(id: impl Into<String>, value: impl Into<String>) -> Self {
        Self::new(id, XisfPropertyValue::String(value.into()))
    }

    pub fn float64(id: impl Into<String>, value: f64) -> Self {
        Self::new(id, XisfPropertyValue::Float64(value))
    }

    pub fn int32(id: impl Into<String>, value: i32) -> Self {
        Self::new(id, XisfPropertyValue::Int32(value))
    }

    pub fn f64_vector(id: impl Into<String>, values: Vec<f64>) -> Self {
        Self::new(id, XisfPropertyValue::F64Vector(values))
    }

    pub fn f64_matrix(id: impl Into<String>, rows: usize, cols: usize, data: Vec<f64>) -> Self {
        Self::new(id, XisfPropertyValue::F64Matrix { rows, cols, data })
    }
}

#[derive(Debug, Clone)]
pub struct XisfHeader {
    pub version: String,
    pub images: Vec<ImageInfo>,
    pub keywords: Vec<Keyword>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum XisfCompression {
    #[default]
    None,
    Lz4,
    Lz4Hc,
    Zlib,
    Zstd,
}

impl XisfCompression {
    pub fn parse(s: &str) -> Self {
        let lower = s.to_lowercase();
        if lower.starts_with("lz4+hc") || lower.starts_with("lz4-hc") {
            Self::Lz4Hc
        } else if lower.starts_with("lz4") {
            Self::Lz4
        } else if lower.starts_with("zlib") {
            Self::Zlib
        } else if lower.starts_with("zstd") || lower.starts_with("zstandard") {
            Self::Zstd
        } else {
            Self::None
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::None => "",
            Self::Lz4 => "lz4",
            Self::Lz4Hc => "lz4-hc",
            Self::Zlib => "zlib",
            Self::Zstd => "zstd",
        }
    }

    pub fn is_compressed(&self) -> bool {
        !matches!(self, Self::None)
    }
}

#[derive(Debug, Clone)]
pub struct ImageInfo {
    pub geometry: Vec<usize>,
    pub sample_format: SampleFormat,
    pub bounds: (f64, f64),
    pub color_space: ColorSpace,
    pub pixel_storage: PixelStorage,
    pub location: DataLocation,
    pub compression: XisfCompression,
    pub uncompressed_size: Option<u64>,
}

#[derive(Debug, Clone)]
pub enum SampleFormat {
    UInt8,
    UInt16,
    UInt32,
    Float32,
    Float64,
}

impl SampleFormat {
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "UInt8" => Ok(Self::UInt8),
            "UInt16" => Ok(Self::UInt16),
            "UInt32" => Ok(Self::UInt32),
            "Float32" => Ok(Self::Float32),
            "Float64" => Ok(Self::Float64),
            _ => Err(XisfError::UnsupportedFormat(s.to_string())),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::UInt8 => "UInt8",
            Self::UInt16 => "UInt16",
            Self::UInt32 => "UInt32",
            Self::Float32 => "Float32",
            Self::Float64 => "Float64",
        }
    }

    pub fn to_bitpix(&self) -> BitPix {
        match self {
            Self::UInt8 => BitPix::U8,
            Self::UInt16 => BitPix::I16,
            Self::UInt32 => BitPix::I32,
            Self::Float32 => BitPix::F32,
            Self::Float64 => BitPix::F64,
        }
    }

    pub fn bytes_per_sample(&self) -> usize {
        match self {
            Self::UInt8 => 1,
            Self::UInt16 => 2,
            Self::UInt32 => 4,
            Self::Float32 => 4,
            Self::Float64 => 8,
        }
    }

    pub fn is_floating_point(&self) -> bool {
        matches!(self, Self::Float32 | Self::Float64)
    }
}

#[derive(Debug, Clone)]
pub enum ColorSpace {
    Gray,
    RGB,
    Unknown(String),
}

impl ColorSpace {
    pub fn parse(s: &str) -> Self {
        match s {
            "Gray" => Self::Gray,
            "RGB" => Self::RGB,
            other => Self::Unknown(other.to_string()),
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            Self::Gray => "Gray",
            Self::RGB => "RGB",
            Self::Unknown(s) => s,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PixelStorage {
    #[default]
    Planar,
    Normal,
}

impl PixelStorage {
    pub fn parse(s: &str) -> Self {
        match s {
            "Normal" => Self::Normal,
            _ => Self::Planar,
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            Self::Planar => "Planar",
            Self::Normal => "Normal",
        }
    }
}

#[derive(Debug, Clone)]
pub struct DataLocation {
    pub offset: u64,
    pub size: u64,
}

impl DataLocation {
    pub fn new(offset: u64, size: u64) -> Self {
        Self { offset, size }
    }

    pub fn parse(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split(':').collect();
        if parts.len() != 3 || parts[0] != "attachment" {
            return Err(XisfError::InvalidFormat(format!(
                "Invalid location format: {}",
                s
            )));
        }

        let offset = parts[1]
            .parse::<u64>()
            .map_err(|_| XisfError::InvalidFormat(format!("Invalid offset: {}", parts[1])))?;

        let size = parts[2]
            .parse::<u64>()
            .map_err(|_| XisfError::InvalidFormat(format!("Invalid size: {}", parts[2])))?;

        Ok(Self { offset, size })
    }

    pub fn format(&self) -> String {
        format!("attachment:{}:{}", self.offset, self.size)
    }
}

pub fn parse_geometry(geometry_str: &str) -> Result<Vec<usize>> {
    let parts: Vec<&str> = geometry_str.split(':').collect();
    if parts.len() < 2 || parts.len() > 3 {
        return Err(XisfError::InvalidGeometry(geometry_str.to_string()));
    }

    let mut dimensions = Vec::new();

    let width = parts[0]
        .parse::<usize>()
        .map_err(|_| XisfError::InvalidGeometry(format!("Invalid width: {}", parts[0])))?;
    let height = parts[1]
        .parse::<usize>()
        .map_err(|_| XisfError::InvalidGeometry(format!("Invalid height: {}", parts[1])))?;

    if width == 0 {
        return Err(XisfError::InvalidGeometry(
            "Width cannot be zero".to_string(),
        ));
    }
    if height == 0 {
        return Err(XisfError::InvalidGeometry(
            "Height cannot be zero".to_string(),
        ));
    }

    dimensions.push(width);
    dimensions.push(height);

    if parts.len() == 3 {
        let channels = parts[2]
            .parse::<usize>()
            .map_err(|_| XisfError::InvalidGeometry(format!("Invalid channels: {}", parts[2])))?;
        if channels == 0 {
            return Err(XisfError::InvalidGeometry(
                "Channels cannot be zero".to_string(),
            ));
        }
        if channels > 1 {
            dimensions.push(channels);
        }
    }

    Ok(dimensions)
}

pub fn format_geometry(geometry: &[usize]) -> String {
    geometry
        .iter()
        .map(|d| d.to_string())
        .collect::<Vec<_>>()
        .join(":")
}

pub fn format_geometry_with_channels(geometry: &[usize]) -> String {
    if geometry.len() == 2 {
        format!("{}:{}:1", geometry[0], geometry[1])
    } else {
        format_geometry(geometry)
    }
}

impl ImageInfo {
    pub fn data_size(&self) -> u64 {
        let total_pixels: u64 = self.geometry.iter().map(|&d| d as u64).product();
        total_pixels * self.sample_format.bytes_per_sample() as u64
    }

    pub fn format_bounds(&self) -> String {
        format!("{}:{}", self.bounds.0, self.bounds.1)
    }

    pub fn num_channels(&self) -> usize {
        if self.geometry.len() >= 3 {
            self.geometry[2]
        } else {
            1
        }
    }

    pub fn pixels_per_channel(&self) -> usize {
        if self.geometry.len() >= 2 {
            self.geometry[0] * self.geometry[1]
        } else {
            0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_format_from_str_valid() {
        let valid_formats = [
            ("UInt8", SampleFormat::UInt8),
            ("UInt16", SampleFormat::UInt16),
            ("UInt32", SampleFormat::UInt32),
            ("Float32", SampleFormat::Float32),
            ("Float64", SampleFormat::Float64),
        ];

        for (input, _expected) in valid_formats {
            let result = SampleFormat::parse(input).unwrap();
            assert!(matches!(
                result,
                SampleFormat::UInt8
                    | SampleFormat::UInt16
                    | SampleFormat::UInt32
                    | SampleFormat::Float32
                    | SampleFormat::Float64
            ));
        }
    }

    #[test]
    fn sample_format_from_str_invalid() {
        let invalid_formats = [
            "uint8",
            "UINT8",
            "Int8",
            "Complex64",
            "Float16",
            "",
            "unknown",
            "42",
        ];

        for input in invalid_formats {
            assert!(SampleFormat::parse(input).is_err());
        }
    }

    #[test]
    fn sample_format_to_bitpix_conversion() {
        let conversions = [
            (SampleFormat::UInt8, BitPix::U8),
            (SampleFormat::UInt16, BitPix::I16),
            (SampleFormat::UInt32, BitPix::I32),
            (SampleFormat::Float32, BitPix::F32),
            (SampleFormat::Float64, BitPix::F64),
        ];

        for (sample_format, expected_bitpix) in conversions {
            assert_eq!(sample_format.to_bitpix(), expected_bitpix);
        }
    }

    #[test]
    fn sample_format_bytes_per_sample() {
        let byte_sizes = [
            (SampleFormat::UInt8, 1),
            (SampleFormat::UInt16, 2),
            (SampleFormat::UInt32, 4),
            (SampleFormat::Float32, 4),
            (SampleFormat::Float64, 8),
        ];

        for (sample_format, expected_bytes) in byte_sizes {
            assert_eq!(sample_format.bytes_per_sample(), expected_bytes);
        }
    }

    #[test]
    fn color_space_from_str() {
        assert!(matches!(ColorSpace::parse("Gray"), ColorSpace::Gray));
        assert!(matches!(ColorSpace::parse("RGB"), ColorSpace::RGB));
        assert!(matches!(ColorSpace::parse("CMYK"), ColorSpace::Unknown(_)));
        assert!(matches!(ColorSpace::parse(""), ColorSpace::Unknown(_)));

        if let ColorSpace::Unknown(s) = ColorSpace::parse("HSV") {
            assert_eq!(s, "HSV");
        } else {
            panic!("Expected Unknown variant");
        }
    }

    #[test]
    fn pixel_storage_parse_and_as_str() {
        assert_eq!(PixelStorage::parse("Planar"), PixelStorage::Planar);
        assert_eq!(PixelStorage::parse("Normal"), PixelStorage::Normal);
        assert_eq!(PixelStorage::parse("unknown"), PixelStorage::Planar);
        assert_eq!(PixelStorage::parse(""), PixelStorage::Planar);

        assert_eq!(PixelStorage::Planar.as_str(), "Planar");
        assert_eq!(PixelStorage::Normal.as_str(), "Normal");
    }

    #[test]
    fn pixel_storage_default() {
        assert_eq!(PixelStorage::default(), PixelStorage::Planar);
    }

    #[test]
    fn data_location_from_str_valid() {
        let location = DataLocation::parse("attachment:1024:2048").unwrap();
        assert_eq!(location.offset, 1024);
        assert_eq!(location.size, 2048);

        let location = DataLocation::parse("attachment:0:999999").unwrap();
        assert_eq!(location.offset, 0);
        assert_eq!(location.size, 999999);
    }

    #[test]
    fn data_location_from_str_invalid() {
        let invalid_locations = [
            "embedded:1024:2048",
            "attachment:abc:2048",
            "attachment:1024:def",
            "attachment:1024",
            "attachment:1024:2048:extra",
            "",
            "attachment::2048",
            "1024:2048",
        ];

        for invalid in invalid_locations {
            assert!(DataLocation::parse(invalid).is_err());
        }
    }

    #[test]
    fn parse_geometry_valid() {
        let geometry = parse_geometry("1920:1080").unwrap();
        assert_eq!(geometry, vec![1920, 1080]);

        let geometry = parse_geometry("1920:1080:3").unwrap();
        assert_eq!(geometry, vec![1920, 1080, 3]);

        let geometry = parse_geometry("1920:1080:1").unwrap();
        assert_eq!(geometry, vec![1920, 1080]);

        let geometry = parse_geometry("1:1").unwrap();
        assert_eq!(geometry, vec![1, 1]);

        let geometry = parse_geometry("65535:65535:4").unwrap();
        assert_eq!(geometry, vec![65535, 65535, 4]);
    }

    #[test]
    fn parse_geometry_invalid() {
        let invalid_geometries = [
            "",
            "1920",
            "1920:1080:3:4",
            "abc:1080",
            "1920:def",
            "1920:1080:xyz",
            "0:1080",
            "1920:0",
            "-1920:1080",
            "1920:-1080",
            "1920:1080:0",
            "1920:1080:-1",
        ];

        for invalid in invalid_geometries {
            assert!(parse_geometry(invalid).is_err());
        }
    }

    #[test]
    fn xisf_header_creation_and_access() {
        let mut header = XisfHeader {
            version: "1.0".to_string(),
            images: Vec::new(),
            keywords: Vec::new(),
        };

        header
            .keywords
            .push(Keyword::string("TELESCOP", "Hubble").with_comment("Telescope name"));

        assert_eq!(header.version, "1.0");
        assert_eq!(header.images.len(), 0);
        assert_eq!(header.keywords.len(), 1);
        assert_eq!(header.keywords[0].name, "TELESCOP");
        assert_eq!(
            header.keywords[0].value,
            Some(crate::fits::header::KeywordValue::String(
                "Hubble".to_string()
            ))
        );
        assert_eq!(
            header.keywords[0].comment,
            Some("Telescope name".to_string())
        );
    }

    #[test]
    fn image_info_creation_and_access() {
        let image_info = ImageInfo {
            geometry: vec![1920, 1080, 3],
            sample_format: SampleFormat::UInt16,
            bounds: (0.0, 65535.0),
            color_space: ColorSpace::RGB,
            pixel_storage: PixelStorage::Normal,
            location: DataLocation {
                offset: 1024,
                size: 12441600,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };

        assert_eq!(image_info.geometry, vec![1920, 1080, 3]);
        assert!(matches!(image_info.sample_format, SampleFormat::UInt16));
        assert_eq!(image_info.bounds, (0.0, 65535.0));
        assert!(matches!(image_info.color_space, ColorSpace::RGB));
        assert_eq!(image_info.pixel_storage, PixelStorage::Normal);
        assert_eq!(image_info.location.offset, 1024);
        assert_eq!(image_info.location.size, 12441600);
    }

    #[test]
    fn data_size_calculations() {
        let formats_and_sizes = [
            (SampleFormat::UInt8, vec![1920, 1080], 1920 * 1080),
            (SampleFormat::UInt16, vec![1920, 1080], 1920 * 1080 * 2),
            (
                SampleFormat::UInt32,
                vec![1920, 1080, 3],
                1920 * 1080 * 3 * 4,
            ),
            (SampleFormat::Float32, vec![512, 512], 512 * 512 * 4),
            (SampleFormat::Float64, vec![100, 100, 4], 100 * 100 * 4 * 8),
        ];

        for (format, geometry, expected_size) in formats_and_sizes {
            let total_pixels: usize = geometry.iter().product();
            let actual_size = total_pixels * format.bytes_per_sample();
            assert_eq!(actual_size, expected_size);
        }
    }

    #[test]
    fn extreme_geometry_values() {
        let large_geometry = parse_geometry("65535:65535").unwrap();
        assert_eq!(large_geometry, vec![65535, 65535]);

        let total_pixels: u64 = large_geometry.iter().map(|&x| x as u64).product();
        assert_eq!(total_pixels, 65535u64 * 65535u64);
    }

    #[test]
    fn roundtrip_data_location() {
        let original_locations = [
            "attachment:0:1",
            "attachment:1024:2048",
            "attachment:18446744073709551615:18446744073709551615",
        ];

        for original in original_locations {
            let _parsed = DataLocation::parse(original).unwrap();
        }
    }

    #[test]
    fn image_info_num_channels() {
        let gray = ImageInfo {
            geometry: vec![100, 100],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::Gray,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation {
                offset: 0,
                size: 10000,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };
        assert_eq!(gray.num_channels(), 1);

        let rgb = ImageInfo {
            geometry: vec![100, 100, 3],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::RGB,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation {
                offset: 0,
                size: 30000,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };
        assert_eq!(rgb.num_channels(), 3);
    }

    #[test]
    fn image_info_pixels_per_channel() {
        let gray = ImageInfo {
            geometry: vec![100, 50],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::Gray,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation {
                offset: 0,
                size: 5000,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };
        assert_eq!(gray.pixels_per_channel(), 5000);

        let rgb = ImageInfo {
            geometry: vec![100, 50, 3],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::RGB,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation {
                offset: 0,
                size: 15000,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };
        assert_eq!(rgb.pixels_per_channel(), 5000);
    }

    #[test]
    fn xisf_property_value_type_names() {
        assert_eq!(
            XisfPropertyValue::String("test".into()).type_name(),
            "String"
        );
        assert_eq!(XisfPropertyValue::Float64(1.0).type_name(), "Float64");
        assert_eq!(XisfPropertyValue::Int32(42).type_name(), "Int32");
        assert_eq!(
            XisfPropertyValue::F64Vector(vec![1.0]).type_name(),
            "F64Vector"
        );
        assert_eq!(
            XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 2,
                data: vec![1.0, 2.0, 3.0, 4.0]
            }
            .type_name(),
            "F64Matrix"
        );
    }

    #[test]
    fn xisf_property_value_formatting() {
        assert_eq!(
            XisfPropertyValue::String("TAN".into()).format_value(),
            "TAN"
        );
        assert_eq!(XisfPropertyValue::Float64(3.14).format_value(), "3.14");
        assert_eq!(XisfPropertyValue::Int32(-42).format_value(), "-42");
        assert_eq!(
            XisfPropertyValue::F64Vector(vec![1.5, 2.5, 3.5]).format_value(),
            "1.5 2.5 3.5"
        );
        assert_eq!(
            XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 2,
                data: vec![1.0, 0.0, 0.0, 1.0]
            }
            .format_value(),
            "1 0 0 1"
        );
    }

    #[test]
    fn xisf_property_constructors() {
        let prop = XisfProperty::string("Test:Id", "value");
        assert_eq!(prop.id, "Test:Id");
        assert_eq!(prop.value, XisfPropertyValue::String("value".into()));

        let prop = XisfProperty::float64("Test:Float", 3.14);
        assert_eq!(prop.id, "Test:Float");
        assert_eq!(prop.value, XisfPropertyValue::Float64(3.14));

        let prop = XisfProperty::int32("Test:Int", 42);
        assert_eq!(prop.id, "Test:Int");
        assert_eq!(prop.value, XisfPropertyValue::Int32(42));

        let prop = XisfProperty::f64_vector("Test:Vec", vec![1.0, 2.0]);
        assert_eq!(prop.id, "Test:Vec");
        assert_eq!(prop.value, XisfPropertyValue::F64Vector(vec![1.0, 2.0]));

        let prop = XisfProperty::f64_matrix("Test:Mat", 2, 2, vec![1.0, 0.0, 0.0, 1.0]);
        assert_eq!(prop.id, "Test:Mat");
        assert_eq!(
            prop.value,
            XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 2,
                data: vec![1.0, 0.0, 0.0, 1.0]
            }
        );
    }
}
