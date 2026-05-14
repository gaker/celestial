mod data_type;
mod layout;
mod wcs;
mod xml;

pub use data_type::XisfDataType;
pub use wcs::wcs_to_xisf_properties;

use crate::fits::header::Keyword;
use crate::xisf::header::{
    ColorSpace, DataLocation, ImageInfo, PixelStorage, XisfCompression, XisfProperty,
};
use crate::xisf::Result;
use std::fs::File;
use std::io::{BufWriter, Cursor, Seek, Write};
use std::path::Path;

const XISF_SIGNATURE: &[u8] = b"XISF0100";
const HEADER_ALIGNMENT: usize = 16;

pub struct XisfWriter<W> {
    pub(super) writer: W,
    pub(super) images: Vec<ImageData>,
    pub(super) keywords: Vec<Keyword>,
    pub(super) properties: Vec<XisfProperty>,
    pub(super) compression: XisfCompression,
    pub(super) property_blocks: Vec<PropertyDataBlock>,
}

pub(super) struct ImageData {
    pub(super) info: ImageInfo,
    pub(super) data: Vec<u8>,
}

pub(super) struct PropertyDataBlock {
    pub(super) property_index: usize,
    pub(super) data: Vec<u8>,
    pub(super) location: DataLocation,
}

impl XisfWriter<BufWriter<File>> {
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        Ok(Self::new(BufWriter::new(file)))
    }
}

impl XisfWriter<Cursor<Vec<u8>>> {
    pub fn write_to_vec(mut self) -> Result<Vec<u8>> {
        self.write_internal()?;
        Ok(self.writer.into_inner())
    }
}

impl<W: Write + Seek> XisfWriter<W> {
    pub fn new(writer: W) -> Self {
        Self {
            writer,
            images: Vec::new(),
            keywords: Vec::new(),
            properties: Vec::new(),
            compression: XisfCompression::None,
            property_blocks: Vec::new(),
        }
    }

    pub fn compression(mut self, compression: XisfCompression) -> Self {
        self.compression = compression;
        self
    }

    pub fn set_compression(&mut self, compression: XisfCompression) {
        self.compression = compression;
    }

    pub fn add_image<T: XisfDataType>(
        &mut self,
        data: &[T],
        width: usize,
        height: usize,
        channels: usize,
    ) -> Result<()> {
        let bounds = T::calculate_bounds(data);
        self.add_image_with_bounds(data, width, height, channels, bounds)
    }

    pub fn add_image_with_bounds<T: XisfDataType>(
        &mut self,
        data: &[T],
        width: usize,
        height: usize,
        channels: usize,
        bounds: (f64, f64),
    ) -> Result<()> {
        let geometry = layout::build_geometry(width, height, channels);
        layout::validate_data_length(data.len(), &geometry)?;

        let raw_bytes = T::to_le_bytes(data);
        let uncompressed_size = raw_bytes.len();
        let (compressed_bytes, compression_used) =
            layout::compress_data(self.compression, &raw_bytes);

        let info = ImageInfo {
            geometry,
            sample_format: T::SAMPLE_FORMAT,
            bounds,
            color_space: if channels == 3 {
                ColorSpace::Rgb
            } else {
                ColorSpace::Gray
            },
            pixel_storage: PixelStorage::Planar,
            location: DataLocation::new(0, compressed_bytes.len() as u64),
            compression: compression_used,
            uncompressed_size: if compression_used.is_compressed() {
                Some(uncompressed_size as u64)
            } else {
                None
            },
        };

        self.images.push(ImageData {
            info,
            data: compressed_bytes,
        });
        Ok(())
    }

    pub fn set_keyword(&mut self, name: &str, value: &str) {
        self.keywords.push(Keyword::string(name, value));
    }

    pub fn add_keyword(&mut self, keyword: Keyword) {
        self.keywords.push(keyword);
    }

    #[deprecated(note = "Use set_keyword() instead")]
    pub fn set_fits_keyword(&mut self, name: &str, value: &str) {
        self.set_keyword(name, value);
    }

    #[deprecated(note = "Use add_keyword() instead")]
    pub fn add_fits_keyword(&mut self, keyword: Keyword) {
        self.add_keyword(keyword);
    }

    pub fn add_property(&mut self, property: XisfProperty) {
        self.properties.push(property);
    }

    pub fn add_properties(&mut self, properties: impl IntoIterator<Item = XisfProperty>) {
        self.properties.extend(properties);
    }

    pub fn write(mut self) -> Result<()> {
        self.write_internal()?;
        self.writer.flush()?;
        Ok(())
    }

    fn write_internal(&mut self) -> Result<()> {
        self.property_blocks = layout::build_property_blocks(&self.properties);
        self.writer.write_all(XISF_SIGNATURE)?;

        let padded_xml = layout::generate_final_xml(self)?;
        layout::write_header_length(&mut self.writer, padded_xml.len() as u32)?;
        self.writer.write_all(&padded_xml)?;
        layout::write_property_data(&mut self.writer, &self.property_blocks)?;
        layout::write_image_data(&mut self.writer, &self.images)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Keyword, KeywordValue};
    use crate::xisf::header::{SampleFormat, XisfProperty};
    use crate::xisf::{XisfError, XisfFile};
    use std::io::Cursor;

    fn roundtrip<T: XisfDataType>(
        data: &[T],
        width: usize,
        height: usize,
        channels: usize,
        compression: XisfCompression,
    ) -> XisfFile<Cursor<Vec<u8>>> {
        let mut writer = XisfWriter::new(Cursor::new(Vec::new())).compression(compression);
        writer.add_image(data, width, height, channels).unwrap();
        let bytes = writer.write_to_vec().unwrap();
        XisfFile::new(Cursor::new(bytes)).unwrap()
    }

    #[test]
    fn new_writer_starts_empty() {
        let w = XisfWriter::new(Cursor::new(Vec::new()));
        assert!(w.images.is_empty());
        assert!(w.keywords.is_empty());
        assert!(w.properties.is_empty());
        assert_eq!(w.compression, XisfCompression::None);
        assert!(w.property_blocks.is_empty());
    }

    #[test]
    fn compression_builder_sets_field() {
        let w = XisfWriter::new(Cursor::new(Vec::new())).compression(XisfCompression::Lz4);
        assert_eq!(w.compression, XisfCompression::Lz4);
    }

    #[test]
    fn set_compression_mutates_in_place() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.set_compression(XisfCompression::Zlib);
        assert_eq!(w.compression, XisfCompression::Zlib);
    }

    #[test]
    fn add_image_length_mismatch_is_rejected() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        let err = w.add_image::<u8>(&[0, 1, 2], 2, 2, 1).unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)));
        assert!(err.to_string().contains("Data length"));
    }

    #[test]
    fn roundtrip_u8_gray() {
        let data: Vec<u8> = (0..100).collect();
        let mut reader = roundtrip(&data, 10, 10, 1, XisfCompression::None);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![10, 10]);
        assert!(matches!(info.sample_format, SampleFormat::UInt8));

        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<u8 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_u16_gray() {
        let data: Vec<u16> = (0u16..64).collect();
        let mut reader = roundtrip(&data, 8, 8, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().sample_format,
            SampleFormat::UInt16
        ));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<u16 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_u32_gray() {
        let data: Vec<u32> = (0u32..36).collect();
        let mut reader = roundtrip(&data, 6, 6, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().sample_format,
            SampleFormat::UInt32
        ));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<u32 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_f32_rgb_preserves_bounds() {
        let data: Vec<f32> = (0..300).map(|i| i as f32 / 300.0).collect();
        let mut reader = roundtrip(&data, 10, 10, 3, XisfCompression::None);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![10, 10, 3]);
        assert!(matches!(info.sample_format, SampleFormat::Float32));
        assert!(matches!(info.color_space, ColorSpace::Rgb));
        let (min, max) = info.bounds;
        assert!(min >= 0.0 && max <= 1.0);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<f32 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_f64_gray() {
        let data: Vec<f64> = (0..16).map(|i| i as f64 * 0.25).collect();
        let mut reader = roundtrip(&data, 4, 4, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().sample_format,
            SampleFormat::Float64
        ));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<f64 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn single_channel_gets_gray_colorspace() {
        let data = vec![128u8; 25];
        let reader = roundtrip(&data, 5, 5, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().color_space,
            ColorSpace::Gray
        ));
    }

    #[test]
    fn three_channel_gets_rgb_colorspace() {
        let data = vec![128u8; 75];
        let reader = roundtrip(&data, 5, 5, 3, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().color_space,
            ColorSpace::Rgb
        ));
    }

    #[test]
    fn lz4_compressed_roundtrip_preserves_data() {
        let data = vec![42u8; 4096];
        let mut reader = roundtrip(&data, 64, 64, 1, XisfCompression::Lz4);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Lz4);
        assert_eq!(info.uncompressed_size, Some(4096));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn lz4_hc_compressed_roundtrip_preserves_data() {
        let data = vec![7u8; 2048];
        let mut reader = roundtrip(&data, 32, 64, 1, XisfCompression::Lz4Hc);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Lz4Hc);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn zlib_compressed_roundtrip_preserves_data() {
        let data = vec![3u8; 4096];
        let mut reader = roundtrip(&data, 64, 64, 1, XisfCompression::Zlib);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Zlib);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn zstd_falls_back_to_uncompressed() {
        let data = vec![5u8; 512];
        let mut reader = roundtrip(&data, 16, 32, 1, XisfCompression::Zstd);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::None);
        assert_eq!(info.uncompressed_size, None);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn incompressible_lz4_falls_back_to_uncompressed() {
        let data: Vec<u8> = (0..255).collect();
        let mut reader = roundtrip(&data, 255, 1, 1, XisfCompression::Lz4);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::None);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn add_image_with_bounds_uses_explicit_bounds() {
        let data: Vec<f32> = vec![0.5; 4];
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image_with_bounds(&data, 2, 2, 1, (-10.0, 10.0)).unwrap();
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.image_info(0).unwrap().bounds, (-10.0, 10.0));
    }

    #[test]
    fn set_keyword_appears_in_roundtrip() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.set_keyword("TELESCOP", "Hubble");
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let kw = reader.get_keyword("TELESCOP").unwrap();
        assert_eq!(kw.value, Some(KeywordValue::String("Hubble".to_string())));
    }

    #[test]
    fn add_keyword_appears_in_roundtrip() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_keyword(Keyword::integer("EXPTIME", 30).with_comment("seconds"));
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let kw = reader.get_keyword("EXPTIME").unwrap();
        assert_eq!(kw.value, Some(KeywordValue::Integer(30)));
    }

    #[test]
    #[allow(deprecated)]
    fn deprecated_set_fits_keyword_still_works() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.set_fits_keyword("OBJECT", "M31");
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let kw = reader.get_keyword("OBJECT").unwrap();
        assert_eq!(kw.value, Some(KeywordValue::String("M31".to_string())));
    }

    #[test]
    #[allow(deprecated)]
    fn deprecated_add_fits_keyword_still_works() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_fits_keyword(Keyword::real("GAIN", 1.5));
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert!(reader.get_keyword("GAIN").is_some());
    }

    #[test]
    fn add_property_and_add_properties_accumulate() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_property(XisfProperty::int32("A", 1));
        w.add_properties(vec![
            XisfProperty::int32("B", 2),
            XisfProperty::int32("C", 3),
        ]);
        assert_eq!(w.properties.len(), 3);
        assert_eq!(w.properties[0].id, "A");
        assert_eq!(w.properties[2].id, "C");
    }

    #[test]
    fn roundtrip_property_data_blocks_for_vector_and_matrix() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_property(XisfProperty::f64_vector("Vec", vec![1.0, 2.0, 3.0]));
        w.add_property(XisfProperty::f64_matrix(
            "Mat",
            2,
            2,
            vec![1.0, 0.0, 0.0, 1.0],
        ));
        w.add_property(XisfProperty::float64("Scalar", 0.5));
        w.add_property(XisfProperty::boolean("Flag", true));
        w.add_property(XisfProperty::int32("Count", 42));
        w.add_property(XisfProperty::string("Name", "value"));

        // Just verify write succeeds and produces a readable file.
        let bytes = w.write_to_vec().unwrap();
        assert!(bytes.starts_with(b"XISF0100"));
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 1);
    }

    #[test]
    fn write_to_vec_starts_with_xisf_signature() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        let bytes = w.write_to_vec().unwrap();
        assert!(bytes.starts_with(b"XISF0100"));
    }

    #[test]
    fn create_writes_to_filesystem() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.xisf");
        let mut w = XisfWriter::create(&path).unwrap();
        w.add_image::<u8>(&[0u8; 16], 4, 4, 1).unwrap();
        w.write().unwrap();

        let mut reader = XisfFile::new(std::fs::File::open(&path).unwrap()).unwrap();
        assert_eq!(reader.num_images(), 1);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, vec![0u8; 16]);
    }

    #[test]
    fn write_flushes_underlying_writer() {
        // write() takes self by value and calls flush. Using Cursor keeps bytes in memory —
        // we just confirm the call path compiles and succeeds.
        let w = XisfWriter::new(Cursor::new(Vec::new()));
        assert!(w.write().is_ok());
    }

    #[test]
    fn multiple_images_get_consecutive_locations() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[1u8; 16], 4, 4, 1).unwrap();
        w.add_image::<u8>(&[2u8; 16], 4, 4, 1).unwrap();
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 2);
        let first = reader.image_info(0).unwrap().location.offset;
        let second = reader.image_info(1).unwrap().location.offset;
        assert!(second > first);
    }

    #[test]
    fn writer_without_images_still_produces_valid_file() {
        let w = XisfWriter::new(Cursor::new(Vec::new()));
        let bytes = w.write_to_vec().unwrap();
        assert!(bytes.starts_with(b"XISF0100"));
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 0);
    }
}
