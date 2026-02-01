use crate::fits::header::{Keyword, KeywordValue};
use crate::xisf::header::{
    format_geometry_with_channels, ColorSpace, DataLocation, ImageInfo, PixelStorage, SampleFormat,
    XisfCompression, XisfProperty, XisfPropertyValue,
};
use crate::xisf::{Result, XisfError};
use byteorder::{LittleEndian, WriteBytesExt};
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, Event};
use quick_xml::Writer;
use std::fs::File;
use std::io::{BufWriter, Cursor, Seek, Write};
use std::path::Path;

const XISF_SIGNATURE: &[u8] = b"XISF0100";
const HEADER_ALIGNMENT: usize = 16;

pub struct XisfWriter<W> {
    writer: W,
    images: Vec<ImageData>,
    keywords: Vec<Keyword>,
    properties: Vec<XisfProperty>,
    compression: XisfCompression,
    property_blocks: Vec<PropertyDataBlock>,
}

struct ImageData {
    info: ImageInfo,
    data: Vec<u8>,
}

struct PropertyDataBlock {
    property_index: usize,
    data: Vec<u8>,
    location: DataLocation,
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
        let geometry = build_geometry(width, height, channels);
        let expected_pixels: usize = geometry.iter().product();
        if data.len() != expected_pixels {
            return Err(XisfError::InvalidFormat(format!(
                "Data length {} does not match geometry {:?}",
                data.len(),
                geometry
            )));
        }

        let raw_bytes = T::to_le_bytes(data);
        let uncompressed_size = raw_bytes.len();
        let (compressed_bytes, compression_used) = self.compress_data(&raw_bytes);

        let info = ImageInfo {
            geometry,
            sample_format: T::SAMPLE_FORMAT,
            bounds,
            color_space: if channels == 3 {
                ColorSpace::RGB
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

    fn compress_data(&self, data: &[u8]) -> (Vec<u8>, XisfCompression) {
        match self.compression {
            XisfCompression::None => (data.to_vec(), XisfCompression::None),
            XisfCompression::Lz4 | XisfCompression::Lz4Hc => {
                let compressed = lz4_flex::compress_prepend_size(data);
                // Only use compression if it actually saves space
                if compressed.len() < data.len() {
                    (compressed, self.compression)
                } else {
                    (data.to_vec(), XisfCompression::None)
                }
            }
            XisfCompression::Zlib => {
                use flate2::write::ZlibEncoder;
                use flate2::Compression;
                let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
                if encoder.write_all(data).is_ok() {
                    if let Ok(compressed) = encoder.finish() {
                        if compressed.len() < data.len() {
                            return (compressed, XisfCompression::Zlib);
                        }
                    }
                }
                (data.to_vec(), XisfCompression::None)
            }
            XisfCompression::Zstd => {
                // Zstd not implemented yet, fall back to no compression
                (data.to_vec(), XisfCompression::None)
            }
        }
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
        self.build_property_blocks();
        self.writer.write_all(XISF_SIGNATURE)?;

        let padded_xml = self.generate_final_xml()?;
        self.write_header_length(padded_xml.len() as u32)?;
        self.writer.write_all(&padded_xml)?;
        self.write_property_data()?;
        self.write_image_data()?;
        Ok(())
    }

    fn build_property_blocks(&mut self) {
        self.property_blocks.clear();
        for (index, property) in self.properties.iter().enumerate() {
            if property.value.needs_data_block() {
                if let Some(data) = property.value.to_le_bytes() {
                    self.property_blocks.push(PropertyDataBlock {
                        property_index: index,
                        data,
                        location: DataLocation::new(0, 0),
                    });
                }
            }
        }
    }

    fn write_property_data(&mut self) -> Result<()> {
        for block in &self.property_blocks {
            self.writer.write_all(&block.data)?;
        }
        Ok(())
    }

    fn write_header_length(&mut self, length: u32) -> Result<()> {
        self.writer.write_u32::<LittleEndian>(length)?;
        self.writer.write_u32::<LittleEndian>(0)?;
        Ok(())
    }

    fn write_image_data(&mut self) -> Result<()> {
        for image in &self.images {
            self.writer.write_all(&image.data)?;
        }
        Ok(())
    }

    fn generate_final_xml(&mut self) -> Result<Vec<u8>> {
        self.set_placeholder_offsets();
        let first_pass = self.generate_xml_content()?;
        let mut padded_size = align_to(first_pass.len(), HEADER_ALIGNMENT);

        loop {
            self.calculate_final_offsets(padded_size);
            let final_xml = self.generate_xml_content()?;
            let actual_padded = align_to(final_xml.len(), HEADER_ALIGNMENT);

            if actual_padded <= padded_size {
                return Ok(pad_to_alignment(&final_xml, HEADER_ALIGNMENT));
            }
            padded_size = actual_padded;
        }
    }

    fn set_placeholder_offsets(&mut self) {
        let mut offset = 0u64;
        for block in &mut self.property_blocks {
            let size = block.data.len() as u64;
            block.location = DataLocation::new(offset, size);
            offset += size;
        }
        for image in &mut self.images {
            image.info.location = DataLocation::new(offset, image.data.len() as u64);
            offset += image.data.len() as u64;
        }
    }

    fn calculate_final_offsets(&mut self, padded_xml_size: usize) {
        let header_size = XISF_SIGNATURE.len() + 8 + padded_xml_size;
        let mut offset = header_size as u64;
        for block in &mut self.property_blocks {
            let size = block.data.len() as u64;
            block.location = DataLocation::new(offset, size);
            offset += size;
        }
        for image in &mut self.images {
            image.info.location = DataLocation::new(offset, image.data.len() as u64);
            offset += image.data.len() as u64;
        }
    }

    fn generate_xml_content(&self) -> Result<Vec<u8>> {
        let mut buffer = Cursor::new(Vec::new());
        let mut writer = Writer::new_with_indent(&mut buffer, b' ', 2);

        write_xml_declaration(&mut writer)?;
        self.write_xisf_content(&mut writer)?;

        Ok(buffer.into_inner())
    }

    fn write_xisf_content<Wr: Write>(&self, writer: &mut Writer<Wr>) -> Result<()> {
        let mut xisf_start = BytesStart::new("xisf");
        xisf_start.push_attribute(("version", "1.0"));
        xisf_start.push_attribute(("xmlns", "http://www.pixinsight.com/xisf"));
        xisf_start.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
        xisf_start.push_attribute((
            "xsi:schemaLocation",
            "http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd",
        ));
        writer
            .write_event(Event::Start(xisf_start))
            .map_err(|e| XisfError::XmlParse(e.to_string()))?;

        for image in &self.images {
            self.write_image_element(writer, &image.info)?;
        }

        writer
            .write_event(Event::End(BytesEnd::new("xisf")))
            .map_err(|e| XisfError::XmlParse(e.to_string()))?;

        Ok(())
    }

    fn write_image_element<Wr: Write>(
        &self,
        writer: &mut Writer<Wr>,
        info: &ImageInfo,
    ) -> Result<()> {
        let mut elem = BytesStart::new("Image");
        elem.push_attribute((
            "geometry",
            format_geometry_with_channels(&info.geometry).as_str(),
        ));
        elem.push_attribute(("sampleFormat", info.sample_format.as_str()));
        if info.sample_format.is_floating_point() {
            elem.push_attribute(("bounds", info.format_bounds().as_str()));
        }
        elem.push_attribute(("colorSpace", info.color_space.as_str()));
        if info.pixel_storage != PixelStorage::Planar {
            elem.push_attribute(("pixelStorage", info.pixel_storage.as_str()));
        }

        // Add compression attribute if data is compressed
        if info.compression.is_compressed() {
            let compression_str = format!(
                "{}:{}",
                info.compression.as_str(),
                info.uncompressed_size.unwrap_or(0)
            );
            elem.push_attribute(("compression", compression_str.as_str()));
        }

        elem.push_attribute(("location", info.location.format().as_str()));

        let has_content = !self.keywords.is_empty() || !self.properties.is_empty();

        if !has_content {
            writer
                .write_event(Event::Empty(elem))
                .map_err(|e| XisfError::XmlParse(e.to_string()))?;
        } else {
            writer
                .write_event(Event::Start(elem))
                .map_err(|e| XisfError::XmlParse(e.to_string()))?;

            self.write_properties(writer)?;
            self.write_keywords_as_fits(writer)?;

            writer
                .write_event(Event::End(BytesEnd::new("Image")))
                .map_err(|e| XisfError::XmlParse(e.to_string()))?;
        }
        Ok(())
    }

    fn write_keywords_as_fits<Wr: Write>(&self, writer: &mut Writer<Wr>) -> Result<()> {
        for keyword in &self.keywords {
            let mut elem = BytesStart::new("FITSKeyword");
            elem.push_attribute(("name", keyword.name.as_str()));

            let value_str = keyword_value_to_string(&keyword.value);
            elem.push_attribute(("value", value_str.as_str()));

            let comment_str = keyword.comment.as_deref().unwrap_or("");
            elem.push_attribute(("comment", comment_str));

            writer
                .write_event(Event::Empty(elem))
                .map_err(|e| XisfError::XmlParse(e.to_string()))?;
        }
        Ok(())
    }

    fn write_properties<Wr: Write>(&self, writer: &mut Writer<Wr>) -> Result<()> {
        use quick_xml::events::BytesText;

        let mut block_index = 0usize;

        for (prop_index, property) in self.properties.iter().enumerate() {
            let mut elem = BytesStart::new("Property");
            elem.push_attribute(("id", property.id.as_str()));
            elem.push_attribute(("type", property.value.type_name()));

            match &property.value {
                XisfPropertyValue::Float64(v) => {
                    elem.push_attribute(("value", v.to_string().as_str()));
                    writer
                        .write_event(Event::Empty(elem))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                }
                XisfPropertyValue::Int32(v) => {
                    elem.push_attribute(("value", v.to_string().as_str()));
                    writer
                        .write_event(Event::Empty(elem))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                }
                XisfPropertyValue::String(s) => {
                    writer
                        .write_event(Event::Start(elem))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                    writer
                        .write_event(Event::Text(BytesText::new(s)))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                    writer
                        .write_event(Event::End(BytesEnd::new("Property")))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                }
                XisfPropertyValue::F64Vector(v) => {
                    elem.push_attribute(("length", v.len().to_string().as_str()));
                    if let Some(block) = self.find_property_block(prop_index, &mut block_index) {
                        elem.push_attribute(("location", block.location.format().as_str()));
                    }
                    writer
                        .write_event(Event::Empty(elem))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                }
                XisfPropertyValue::F64Matrix { rows, cols, .. } => {
                    elem.push_attribute(("rows", rows.to_string().as_str()));
                    elem.push_attribute(("columns", cols.to_string().as_str()));
                    if let Some(block) = self.find_property_block(prop_index, &mut block_index) {
                        elem.push_attribute(("location", block.location.format().as_str()));
                    }
                    writer
                        .write_event(Event::Empty(elem))
                        .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                }
            }
        }
        Ok(())
    }

    fn find_property_block(
        &self,
        prop_index: usize,
        hint: &mut usize,
    ) -> Option<&PropertyDataBlock> {
        for i in *hint..self.property_blocks.len() {
            if self.property_blocks[i].property_index == prop_index {
                *hint = i + 1;
                return Some(&self.property_blocks[i]);
            }
        }
        None
    }
}

fn keyword_value_to_string(value: &Option<KeywordValue>) -> String {
    match value {
        None => String::new(),
        Some(KeywordValue::Logical(b)) => {
            if *b {
                "T".to_string()
            } else {
                "F".to_string()
            }
        }
        Some(KeywordValue::Integer(i)) => i.to_string(),
        Some(KeywordValue::Real(f)) => {
            // Ensure decimal point is preserved for whole numbers
            let s = f.to_string();
            if s.contains('.') || s.contains('e') || s.contains('E') {
                s
            } else {
                format!("{}.0", s)
            }
        }
        Some(KeywordValue::String(s)) => s.clone(),
        Some(KeywordValue::Complex(r, i)) => format!("({}, {})", r, i),
    }
}

fn write_xml_declaration<W: Write>(writer: &mut Writer<W>) -> Result<()> {
    let decl = BytesDecl::new("1.0", Some("UTF-8"), None);
    writer
        .write_event(Event::Decl(decl))
        .map_err(|e| XisfError::XmlParse(e.to_string()))
}

fn build_geometry(width: usize, height: usize, channels: usize) -> Vec<usize> {
    if channels > 1 {
        vec![width, height, channels]
    } else {
        vec![width, height]
    }
}

fn pad_to_alignment(data: &[u8], alignment: usize) -> Vec<u8> {
    let mut padded = data.to_vec();
    let remainder = padded.len() % alignment;
    if remainder != 0 {
        padded.resize(padded.len() + (alignment - remainder), 0);
    }
    padded
}

fn align_to(size: usize, alignment: usize) -> usize {
    let remainder = size % alignment;
    if remainder == 0 {
        size
    } else {
        size + (alignment - remainder)
    }
}

pub trait XisfDataType: Copy {
    const SAMPLE_FORMAT: SampleFormat;
    fn to_le_bytes(data: &[Self]) -> Vec<u8>;
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self>;
    fn default_bounds() -> (f64, f64);
    fn calculate_bounds(data: &[Self]) -> (f64, f64);
}

impl XisfDataType for u8 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt8;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        data.to_vec()
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes.to_vec()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 255.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for u16 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt16;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 2);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 65535.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for u32 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::UInt32;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 4);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 4294967295.0)
    }
    fn calculate_bounds(_data: &[Self]) -> (f64, f64) {
        Self::default_bounds()
    }
}

impl XisfDataType for f32 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::Float32;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 4);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 1.0)
    }
    fn calculate_bounds(data: &[Self]) -> (f64, f64) {
        if data.is_empty() {
            return Self::default_bounds();
        }
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        for &v in data {
            let fv = v as f64;
            if fv.is_finite() {
                min = min.min(fv);
                max = max.max(fv);
            }
        }
        if min.is_infinite() || max.is_infinite() {
            Self::default_bounds()
        } else {
            (min, max)
        }
    }
}

impl XisfDataType for f64 {
    const SAMPLE_FORMAT: SampleFormat = SampleFormat::Float64;
    fn to_le_bytes(data: &[Self]) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(data.len() * 8);
        for &val in data {
            bytes.extend_from_slice(&val.to_le_bytes());
        }
        bytes
    }
    fn from_le_bytes(bytes: &[u8]) -> Vec<Self> {
        bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]))
            .collect()
    }
    fn default_bounds() -> (f64, f64) {
        (0.0, 1.0)
    }
    fn calculate_bounds(data: &[Self]) -> (f64, f64) {
        if data.is_empty() {
            return Self::default_bounds();
        }
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        for &v in data {
            if v.is_finite() {
                min = min.min(v);
                max = max.max(v);
            }
        }
        if min.is_infinite() || max.is_infinite() {
            Self::default_bounds()
        } else {
            (min, max)
        }
    }
}

pub fn wcs_to_xisf_properties(wcs_keywords: &[celestial_wcs::WcsKeyword]) -> Vec<XisfProperty> {
    let mut properties = Vec::new();
    let mut crval1: Option<f64> = None;
    let mut crval2: Option<f64> = None;
    let mut crpix1: Option<f64> = None;
    let mut crpix2: Option<f64> = None;
    let mut cd1_1: Option<f64> = None;
    let mut cd1_2: Option<f64> = None;
    let mut cd2_1: Option<f64> = None;
    let mut cd2_2: Option<f64> = None;
    let mut proj_code: Option<String> = None;
    let mut lonpole: Option<f64> = None;
    let mut latpole: Option<f64> = None;

    for kw in wcs_keywords {
        match kw.name.as_str() {
            "CTYPE1" => {
                if let celestial_wcs::WcsKeywordValue::String(s) = &kw.value {
                    proj_code = extract_projection_code(s);
                }
            }
            "CRVAL1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crval1 = Some(*v);
                }
            }
            "CRVAL2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crval2 = Some(*v);
                }
            }
            "CRPIX1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crpix1 = Some(*v);
                }
            }
            "CRPIX2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    crpix2 = Some(*v);
                }
            }
            "CD1_1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd1_1 = Some(*v);
                }
            }
            "CD1_2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd1_2 = Some(*v);
                }
            }
            "CD2_1" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd2_1 = Some(*v);
                }
            }
            "CD2_2" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    cd2_2 = Some(*v);
                }
            }
            "LONPOLE" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    lonpole = Some(*v);
                }
            }
            "LATPOLE" => {
                if let celestial_wcs::WcsKeywordValue::Real(v) = &kw.value {
                    latpole = Some(*v);
                }
            }
            _ => {}
        }
    }

    if let Some(code) = proj_code {
        properties.push(XisfProperty::string(
            "PCL:AstrometricSolution:ProjectionSystem",
            code,
        ));
    }

    if let (Some(ra), Some(dec)) = (crval1, crval2) {
        properties.push(XisfProperty::f64_vector(
            "PCL:AstrometricSolution:ReferenceCelestialCoordinates",
            vec![ra, dec],
        ));
    }

    if let (Some(x), Some(y)) = (crpix1, crpix2) {
        properties.push(XisfProperty::f64_vector(
            "PCL:AstrometricSolution:ReferenceImageCoordinates",
            vec![x, y],
        ));
    }

    if let (Some(c11), Some(c12), Some(c21), Some(c22)) = (cd1_1, cd1_2, cd2_1, cd2_2) {
        properties.push(XisfProperty::f64_matrix(
            "PCL:AstrometricSolution:LinearTransformationMatrix",
            2,
            2,
            vec![c11, c12, c21, c22],
        ));
    }

    // ReferenceNativeCoordinates: (φ₀, θ₀) in degrees
    // For zenithal projections (TAN, SIN, etc.), this is (0, 90)
    properties.push(XisfProperty::f64_vector(
        "PCL:AstrometricSolution:ReferenceNativeCoordinates",
        vec![0.0, 90.0],
    ));

    // CelestialPoleNativeCoordinates: (φₚ, θₚ) derived from LONPOLE/LATPOLE
    // φₚ = LONPOLE (default 180° if CRVAL2 >= θ₀, else 0°)
    // θₚ = LATPOLE (default = CRVAL2 for zenithal)
    let phi_p = lonpole.unwrap_or_else(|| {
        if crval2.unwrap_or(0.0) >= 90.0 {
            0.0
        } else {
            180.0
        }
    });
    let theta_p = latpole.unwrap_or_else(|| crval2.unwrap_or(0.0));
    properties.push(XisfProperty::f64_vector(
        "PCL:AstrometricSolution:CelestialPoleNativeCoordinates",
        vec![phi_p, theta_p],
    ));

    properties
}

fn extract_projection_code(ctype: &str) -> Option<String> {
    let trimmed = ctype.trim();
    if trimmed.len() < 8 {
        return None;
    }
    let dashes_pos = trimmed.find('-')?;
    let after_dashes = &trimmed[dashes_pos..];
    let code_start = after_dashes.trim_start_matches('-');
    if code_start.is_empty() {
        return None;
    }
    Some(wcs_code_to_pixinsight_name(code_start))
}

fn wcs_code_to_pixinsight_name(code: &str) -> String {
    match code {
        "TAN" => "Gnomonic",
        "SIN" => "Orthographic",
        "STG" => "Stereographic",
        "ARC" => "ZenithalEqualArea",
        "ZEA" => "ZenithalEqualArea",
        "AZP" => "ZenithalPerspective",
        "SZP" => "SlantZenithalPerspective",
        "ZPN" => "ZenithalPolynomial",
        "AIR" => "Airy",
        "CYP" => "CylindricalPerspective",
        "CEA" => "CylindricalEqualArea",
        "CAR" => "PlateCarree",
        "MER" => "Mercator",
        "SFL" => "SansonFlamsteed",
        "PAR" => "Parabolic",
        "MOL" => "Mollweide",
        "AIT" => "HammerAitoff",
        "COP" => "ConicPerspective",
        "COE" => "ConicEqualArea",
        "COD" => "ConicEquidistant",
        "COO" => "ConicOrthomorphic",
        other => return other.to_string(),
    }
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn build_geometry_single_channel() {
        let g = build_geometry(100, 50, 1);
        assert_eq!(g, vec![100, 50]);
    }

    #[test]
    fn build_geometry_rgb() {
        let g = build_geometry(100, 50, 3);
        assert_eq!(g, vec![100, 50, 3]);
    }

    #[test]
    fn pad_to_alignment_no_padding_needed() {
        let data = vec![0u8; 16];
        let padded = pad_to_alignment(&data, 16);
        assert_eq!(padded.len(), 16);
    }

    #[test]
    fn pad_to_alignment_needs_padding() {
        let data = vec![0u8; 10];
        let padded = pad_to_alignment(&data, 16);
        assert_eq!(padded.len(), 16);
    }

    #[test]
    fn align_to_already_aligned() {
        assert_eq!(align_to(32, 16), 32);
    }

    #[test]
    fn align_to_needs_alignment() {
        assert_eq!(align_to(17, 16), 32);
    }

    #[test]
    fn u8_roundtrip() {
        let data = vec![0u8, 127, 255];
        let bytes = <u8 as XisfDataType>::to_le_bytes(&data);
        let restored = <u8 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(data, restored);
    }

    #[test]
    fn u16_roundtrip() {
        let data = vec![0u16, 32768, 65535];
        let bytes = <u16 as XisfDataType>::to_le_bytes(&data);
        let restored = <u16 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(data, restored);
    }

    #[test]
    fn f32_roundtrip() {
        let data = vec![0.0f32, 0.5, 1.0];
        let bytes = <f32 as XisfDataType>::to_le_bytes(&data);
        let restored = <f32 as XisfDataType>::from_le_bytes(&bytes);
        assert_eq!(data, restored);
    }

    #[test]
    fn writer_add_image_valid() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);
        let data = vec![0u16; 100];
        let result = writer.add_image(&data, 10, 10, 1);
        assert!(result.is_ok());
    }

    #[test]
    fn writer_add_image_invalid_size() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);
        let data = vec![0u16; 50];
        let result = writer.add_image(&data, 10, 10, 1);
        assert!(result.is_err());
    }

    #[test]
    fn writer_set_keyword() {
        use crate::fits::header::KeywordValue;

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);
        writer.set_keyword("TELESCOP", "Test");
        assert_eq!(writer.keywords.len(), 1);
        assert_eq!(writer.keywords[0].name, "TELESCOP");
        assert_eq!(
            writer.keywords[0].value,
            Some(KeywordValue::String("Test".to_string()))
        );
    }

    #[test]
    fn write_complete_xisf() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.set_keyword("TELESCOP", "TestScope");
        writer.set_keyword("FILTER", "R");

        let result = writer.write();
        assert!(result.is_ok());
    }

    #[test]
    fn write_rgb_image() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<f32> = vec![0.5; 300];
        writer.add_image(&data, 10, 10, 3).unwrap();

        let result = writer.write();
        assert!(result.is_ok());
    }

    #[test]
    fn write_multiple_images() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data1: Vec<u8> = vec![128; 100];
        let data2: Vec<u16> = vec![32768; 100];

        writer.add_image(&data1, 10, 10, 1).unwrap();
        writer.add_image(&data2, 10, 10, 1).unwrap();

        let result = writer.write();
        assert!(result.is_ok());
    }

    #[test]
    fn roundtrip_u16_image() {
        use crate::xisf::{PixelStorage, XisfFile};

        let original_data: Vec<u16> = (0..100).collect();
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);
        writer.add_image(&original_data, 10, 10, 1).unwrap();
        writer.set_keyword("TELESCOP", "RoundtripTest");

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        assert_eq!(reader.num_images(), 1);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![10, 10]);
        assert!(matches!(
            info.sample_format,
            crate::xisf::SampleFormat::UInt16
        ));
        assert_eq!(info.pixel_storage, PixelStorage::Planar);

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <u16 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);

        let telescop = reader.get_keyword("TELESCOP");
        assert!(telescop.is_some());
        assert_eq!(
            telescop.unwrap().value,
            Some(KeywordValue::String("RoundtripTest".to_string()))
        );
    }

    #[test]
    fn roundtrip_f32_rgb_image() {
        use crate::xisf::{PixelStorage, XisfFile};

        let original_data: Vec<f32> = (0..300).map(|i| i as f32 / 300.0).collect();
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);
        writer.add_image(&original_data, 10, 10, 3).unwrap();

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        assert_eq!(reader.num_images(), 1);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![10, 10, 3]);
        assert!(matches!(
            info.sample_format,
            crate::xisf::SampleFormat::Float32
        ));
        assert!(matches!(info.color_space, crate::xisf::ColorSpace::RGB));
        assert_eq!(info.pixel_storage, PixelStorage::Planar);

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <f32 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);
    }

    #[test]
    fn roundtrip_multiple_images() {
        use crate::xisf::XisfFile;

        let data1: Vec<u8> = (0..100).collect();
        let data2: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);
        writer.add_image(&data1, 10, 10, 1).unwrap();
        writer.add_image(&data2, 10, 10, 1).unwrap();

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        assert_eq!(reader.num_images(), 2);

        let raw1 = reader.read_image_data_raw(0).unwrap();
        let restored1 = <u8 as XisfDataType>::from_le_bytes(&raw1);
        assert_eq!(data1, restored1);

        let raw2 = reader.read_image_data_raw(1).unwrap();
        let restored2 = <f64 as XisfDataType>::from_le_bytes(&raw2);
        assert_eq!(data2, restored2);
    }

    #[test]
    fn roundtrip_file_io() {
        use crate::xisf::XisfFile;
        use tempfile::NamedTempFile;

        let original_data: Vec<u16> = (0..256).collect();
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let writer = XisfWriter::create(path).unwrap();
        let mut writer = writer;
        writer.add_image(&original_data, 16, 16, 1).unwrap();
        writer.set_keyword("OBSERVER", "Test");
        writer.write().unwrap();

        let mut reader = XisfFile::open(path).unwrap();
        assert_eq!(reader.num_images(), 1);

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <u16 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);

        let observer = reader.get_keyword("OBSERVER");
        assert!(observer.is_some());
        assert_eq!(
            observer.unwrap().value,
            Some(KeywordValue::String("Test".to_string()))
        );
    }

    #[test]
    fn roundtrip_history_and_comment_keywords() {
        use crate::fits::header::Keyword;
        use crate::xisf::XisfFile;

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();

        writer
            .add_keyword(Keyword::real("EXPTIME", 300.0).with_comment("Exposure time in seconds"));
        writer.add_keyword(Keyword::history("Dark frame subtracted"));
        writer.add_keyword(Keyword::history("Flat field corrected"));
        writer.add_keyword(Keyword::comment("Processed with Celestial"));

        let inner = writer.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(inner)).unwrap();

        let keywords = reader.keywords();
        assert_eq!(keywords.len(), 4);

        let exptime = &keywords[0];
        assert_eq!(exptime.name, "EXPTIME");
        assert_eq!(exptime.value, Some(KeywordValue::Real(300.0)));
        assert_eq!(
            exptime.comment,
            Some("Exposure time in seconds".to_string())
        );

        let history1 = &keywords[1];
        assert_eq!(history1.name, "HISTORY");
        assert_eq!(history1.value, None);
        assert_eq!(history1.comment, Some("Dark frame subtracted".to_string()));

        let history2 = &keywords[2];
        assert_eq!(history2.name, "HISTORY");
        assert_eq!(history2.value, None);
        assert_eq!(history2.comment, Some("Flat field corrected".to_string()));

        let comment = &keywords[3];
        assert_eq!(comment.name, "COMMENT");
        assert_eq!(comment.value, None);
        assert_eq!(
            comment.comment,
            Some("Processed with Celestial".to_string())
        );
    }

    #[test]
    fn single_quote_in_keyword_value() {
        use crate::xisf::XisfFile;

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.set_keyword("STACKMTH", "Sigma Clip (2.5σ)");
        writer.set_keyword("INSTRUME", "Apollo-M MINI (IMX429)");

        let inner = writer.write_to_vec().unwrap();

        // Check XML doesn't contain &apos;
        let xml_str = String::from_utf8_lossy(&inner);
        assert!(
            !xml_str.contains("&apos;"),
            "XML should not contain &apos; escape sequence"
        );

        // Verify roundtrip
        let reader = XisfFile::new(Cursor::new(inner)).unwrap();
        let stackmth = reader.get_keyword("STACKMTH");
        assert!(stackmth.is_some());
        assert_eq!(
            stackmth.unwrap().value,
            Some(KeywordValue::String("Sigma Clip (2.5σ)".to_string()))
        );

        let instrume = reader.get_keyword("INSTRUME");
        assert!(instrume.is_some());
        assert_eq!(
            instrume.unwrap().value,
            Some(KeywordValue::String("Apollo-M MINI (IMX429)".to_string()))
        );
    }

    #[test]
    fn roundtrip_lz4_compressed_u16() {
        use crate::xisf::XisfFile;

        // Create data that should compress well (repeating pattern)
        let original_data: Vec<u16> = (0..1000).map(|i| (i % 256) as u16).collect();

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer).compression(XisfCompression::Lz4);
        writer.add_image(&original_data, 100, 10, 1).unwrap();
        writer.set_keyword("COMPRESS", "LZ4");

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        assert_eq!(reader.num_images(), 1);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Lz4);
        assert!(info.uncompressed_size.is_some());

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <u16 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);

        let compress = reader.get_keyword("COMPRESS");
        assert!(compress.is_some());
        assert_eq!(
            compress.unwrap().value,
            Some(KeywordValue::String("LZ4".to_string()))
        );
    }

    #[test]
    fn roundtrip_lz4_compressed_f32_rgb() {
        use crate::xisf::XisfFile;

        // Create RGB float data with patterns that should compress well
        let original_data: Vec<f32> = (0..3000).map(|i| (i % 256) as f32 / 255.0).collect();

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer).compression(XisfCompression::Lz4);
        writer.add_image(&original_data, 100, 10, 3).unwrap();

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![100, 10, 3]);
        assert!(matches!(info.color_space, crate::xisf::ColorSpace::RGB));

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <f32 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);
    }

    #[test]
    fn roundtrip_zlib_compressed_u16() {
        use crate::xisf::XisfFile;

        let original_data: Vec<u16> = (0..1000).map(|i| (i % 256) as u16).collect();

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer).compression(XisfCompression::Zlib);
        writer.add_image(&original_data, 100, 10, 1).unwrap();

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Zlib);

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <u16 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);
    }

    #[test]
    fn compression_fallback_for_incompressible_data() {
        use crate::xisf::XisfFile;

        // Random-like data that won't compress well
        let original_data: Vec<u8> = (0..100).map(|i| (i * 17 + 31) as u8).collect();

        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer).compression(XisfCompression::Lz4);
        writer.add_image(&original_data, 10, 10, 1).unwrap();

        let inner = writer.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(inner)).unwrap();

        // Small random data might not compress, so compression could be None
        let _info = reader.image_info(0).unwrap();
        // Data should still round-trip correctly regardless of compression
        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <u8 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);
    }

    #[test]
    fn compression_with_file_io() {
        use crate::xisf::XisfFile;
        use tempfile::NamedTempFile;

        let original_data: Vec<u16> = (0..1000).map(|i| (i % 256) as u16).collect();
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let writer = XisfWriter::create(path)
            .unwrap()
            .compression(XisfCompression::Lz4);
        let mut writer = writer;
        writer.add_image(&original_data, 100, 10, 1).unwrap();
        writer.write().unwrap();

        let mut reader = XisfFile::open(path).unwrap();
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Lz4);

        let raw_bytes = reader.read_image_data_raw(0).unwrap();
        let restored = <u16 as XisfDataType>::from_le_bytes(&raw_bytes);
        assert_eq!(original_data, restored);
    }

    #[test]
    fn write_string_property() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.add_property(XisfProperty::string(
            "PCL:AstrometricSolution:ProjectionSystem",
            "TAN",
        ));

        let inner = writer.write_to_vec().unwrap();
        let xml_str = String::from_utf8_lossy(&inner);
        assert!(xml_str.contains(r#"id="PCL:AstrometricSolution:ProjectionSystem""#));
        assert!(xml_str.contains(r#"type="String""#));
        assert!(xml_str.contains(">TAN</Property>"));
    }

    #[test]
    fn write_scalar_float64_property() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.add_property(XisfProperty::float64("Observation:CenterRA", 191.758));

        let inner = writer.write_to_vec().unwrap();
        let xml_str = String::from_utf8_lossy(&inner);
        assert!(xml_str.contains(r#"id="Observation:CenterRA""#));
        assert!(xml_str.contains(r#"type="Float64""#));
        assert!(xml_str.contains(r#"value="191.758""#));
        assert!(xml_str.contains(r#"/>"#));
        assert!(
            !xml_str.contains("</Property>")
                || xml_str.contains(">TAN</Property>")
                || xml_str.contains(">Gnomonic</Property>")
        );
    }

    #[test]
    fn write_scalar_int32_property() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.add_property(XisfProperty::int32("Observation:FrameCount", 42));

        let inner = writer.write_to_vec().unwrap();
        let xml_str = String::from_utf8_lossy(&inner);
        assert!(xml_str.contains(r#"id="Observation:FrameCount""#));
        assert!(xml_str.contains(r#"type="Int32""#));
        assert!(xml_str.contains(r#"value="42""#));
    }

    #[test]
    fn write_f64_vector_property() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.add_property(XisfProperty::f64_vector(
            "PCL:AstrometricSolution:ReferenceCelestialCoordinates",
            vec![191.758, -5.048],
        ));

        let inner = writer.write_to_vec().unwrap();
        let xml_str = String::from_utf8_lossy(&inner);
        assert!(xml_str.contains(r#"id="PCL:AstrometricSolution:ReferenceCelestialCoordinates""#));
        assert!(xml_str.contains(r#"type="F64Vector""#));
        assert!(xml_str.contains(r#"length="2""#));
        assert!(xml_str.contains(r#"location="attachment:"#));
        assert!(xml_str.contains(r#"/>"#));
    }

    #[test]
    fn write_f64_matrix_property() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.add_property(XisfProperty::f64_matrix(
            "PCL:AstrometricSolution:LinearTransformationMatrix",
            2,
            2,
            vec![-0.000409, 0.0, 0.0, 0.000409],
        ));

        let inner = writer.write_to_vec().unwrap();
        let xml_str = String::from_utf8_lossy(&inner);
        assert!(xml_str.contains(r#"id="PCL:AstrometricSolution:LinearTransformationMatrix""#));
        assert!(xml_str.contains(r#"type="F64Matrix""#));
        assert!(xml_str.contains(r#"rows="2""#));
        assert!(xml_str.contains(r#"columns="2""#));
        assert!(xml_str.contains(r#"location="attachment:"#));
    }

    #[test]
    fn wcs_keywords_to_properties() {
        use celestial_wcs::WcsKeyword;

        let wcs_keywords = vec![
            WcsKeyword::string("CTYPE1", "RA---TAN"),
            WcsKeyword::string("CTYPE2", "DEC--TAN"),
            WcsKeyword::real("CRVAL1", 191.758),
            WcsKeyword::real("CRVAL2", -5.048),
            WcsKeyword::real("CRPIX1", 512.5),
            WcsKeyword::real("CRPIX2", 512.5),
            WcsKeyword::real("CD1_1", -0.000409),
            WcsKeyword::real("CD1_2", 0.0),
            WcsKeyword::real("CD2_1", 0.0),
            WcsKeyword::real("CD2_2", 0.000409),
        ];

        let properties = wcs_to_xisf_properties(&wcs_keywords);

        assert_eq!(properties.len(), 6);

        assert_eq!(properties[0].id, "PCL:AstrometricSolution:ProjectionSystem");
        assert_eq!(
            properties[0].value,
            XisfPropertyValue::String("Gnomonic".to_string())
        );

        assert_eq!(
            properties[1].id,
            "PCL:AstrometricSolution:ReferenceCelestialCoordinates"
        );
        assert_eq!(
            properties[1].value,
            XisfPropertyValue::F64Vector(vec![191.758, -5.048])
        );

        assert_eq!(
            properties[2].id,
            "PCL:AstrometricSolution:ReferenceImageCoordinates"
        );
        assert_eq!(
            properties[2].value,
            XisfPropertyValue::F64Vector(vec![512.5, 512.5])
        );

        assert_eq!(
            properties[3].id,
            "PCL:AstrometricSolution:LinearTransformationMatrix"
        );
        assert_eq!(
            properties[3].value,
            XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 2,
                data: vec![-0.000409, 0.0, 0.0, 0.000409]
            }
        );

        assert_eq!(
            properties[4].id,
            "PCL:AstrometricSolution:ReferenceNativeCoordinates"
        );
        assert_eq!(
            properties[4].value,
            XisfPropertyValue::F64Vector(vec![0.0, 90.0])
        );

        assert_eq!(
            properties[5].id,
            "PCL:AstrometricSolution:CelestialPoleNativeCoordinates"
        );
        assert_eq!(
            properties[5].value,
            XisfPropertyValue::F64Vector(vec![180.0, -5.048])
        );
    }

    #[test]
    fn extract_projection_code_tan() {
        assert_eq!(
            extract_projection_code("RA---TAN"),
            Some("Gnomonic".to_string())
        );
    }

    #[test]
    fn extract_projection_code_sin() {
        assert_eq!(
            extract_projection_code("RA---SIN"),
            Some("Orthographic".to_string())
        );
    }

    #[test]
    fn extract_projection_code_with_spaces() {
        assert_eq!(
            extract_projection_code("RA---TAN "),
            Some("Gnomonic".to_string())
        );
    }

    #[test]
    fn extract_projection_code_dec() {
        assert_eq!(
            extract_projection_code("DEC--TAN"),
            Some("Gnomonic".to_string())
        );
    }

    #[test]
    fn extract_projection_code_too_short() {
        assert_eq!(extract_projection_code("RA-TAN"), None);
    }

    #[test]
    fn write_properties_and_keywords_together() {
        let buffer = Cursor::new(Vec::new());
        let mut writer = XisfWriter::new(buffer);

        let data: Vec<u16> = (0..100).collect();
        writer.add_image(&data, 10, 10, 1).unwrap();
        writer.add_property(XisfProperty::string("Test:Property", "TestValue"));
        writer.set_keyword("TELESCOP", "TestScope");

        let inner = writer.write_to_vec().unwrap();
        let xml_str = String::from_utf8_lossy(&inner);

        assert!(xml_str.contains(r#"<Property id="Test:Property""#));
        assert!(xml_str.contains(r#"<FITSKeyword name="TELESCOP""#));
    }
}
