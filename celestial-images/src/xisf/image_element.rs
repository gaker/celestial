use crate::core::types::ByteOrder;
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::block_attrs::BlockAttrs;
use crate::xisf::header::color_space::ColorSpace;
use crate::xisf::header::compression::XisfCompression;
use crate::xisf::header::geometry::parse_geometry;
use crate::xisf::header::image_info::ImageInfo;
use crate::xisf::header::location::{DataLocation, LocationAttr};
use crate::xisf::header::pixel_storage::PixelStorage;
use crate::xisf::header::sample_format::SampleFormat;
use crate::xisf::reader::attribute_text;
use quick_xml::events::BytesStart;

// The values an Image element has for the attributes it leaves out.
pub(super) fn image_defaults() -> ImageInfo {
    ImageInfo {
        geometry: Vec::new(),
        sample_format: SampleFormat::Float32,
        bounds: (0.0, 1.0),
        color_space: ColorSpace::Gray,
        pixel_storage: PixelStorage::Planar,
        location: DataLocation::new(0, 0),
        compression: XisfCompression::None,
        uncompressed_size: None,
        shuffle_item_size: None,
        subblocks: Vec::new(),
        byte_order: ByteOrder::LittleEndian,
        checksum: None,
        block_defect: None,
    }
}

// Reads the attributes of an Image element. An element without a location
// is listed with a defect, so the image indices match the property
// scanner's (XISF 1.0 section 7).
pub(super) fn parse_image_element(element: &BytesStart) -> Result<ImageInfo> {
    let mut info = image_defaults();
    let mut block = BlockAttrs::default();
    let mut bounds = None;
    let mut has_location = false;

    for attr in element.attributes() {
        let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
        let value = attribute_text(&attr)?;
        match attr.key.as_ref() {
            "geometry" => info.geometry = parse_geometry(&value)?,
            "sampleFormat" => info.sample_format = SampleFormat::parse(&value)?,
            "bounds" => bounds = Some(parse_bounds(&value)),
            "colorSpace" => info.color_space = ColorSpace::parse(&value),
            "pixelStorage" => info.pixel_storage = PixelStorage::parse(&value),
            "location" => {
                has_location = true;
                apply_location(&mut info, &mut block, &value);
            }
            key => block.apply(key, &value),
        }
    }

    if !has_location {
        record_missing_location(&mut info);
    }
    apply_bounds(&mut info, &mut block, bounds);
    Ok(block.into_image_info(info))
}

fn record_missing_location(info: &mut ImageInfo) {
    let mut missing = BlockAttrs::default();
    missing.record(XisfError::InvalidFormat(
        "The Image element has no location attribute".to_string(),
    ));
    missing.apply_to(info);
}

// Sets the location given by a location attribute. A header block starts
// empty and gets its bytes when its text is read.
fn apply_location(info: &mut ImageInfo, block: &mut BlockAttrs, value: &str) {
    match LocationAttr::parse(value) {
        Ok(LocationAttr::Attachment(location)) => info.location = location,
        Ok(LocationAttr::Inline(encoding)) => {
            info.location = DataLocation::inline(encoding, Vec::new())
        }
        Ok(LocationAttr::Embedded) => info.location = DataLocation::embedded(Vec::new()),
        Err(error) => block.record(error),
    }
}

// XISF 1.0 section 11.5.1 requires bounds for floating point images, whose
// samples have no default range, and makes them optional for integer
// images. Malformed bounds make a floating point image unreadable. On an
// integer image no sample decoding uses them, so they are ignored as if the
// attribute were absent, and the image has the range of its sample type.
fn apply_bounds(info: &mut ImageInfo, block: &mut BlockAttrs, bounds: Option<Result<(f64, f64)>>) {
    match (bounds, integer_range(&info.sample_format)) {
        (Some(Ok(bounds)), _) => info.bounds = bounds,
        (_, Some(range)) => info.bounds = range,
        (Some(Err(error)), None) => block.record(error),
        (None, None) => {}
    }
}

// [0, 2^k - 1] for a k-bit unsigned sample. u64::MAX rounds to 2^64 in f64.
fn integer_range(format: &SampleFormat) -> Option<(f64, f64)> {
    let high = match format {
        SampleFormat::UInt8 => f64::from(u8::MAX),
        SampleFormat::UInt16 => f64::from(u16::MAX),
        SampleFormat::UInt32 => f64::from(u32::MAX),
        SampleFormat::UInt64 => u64::MAX as f64,
        SampleFormat::Float32 | SampleFormat::Float64 => return None,
    };
    Some((0.0, high))
}

// The black point and the white point.
pub(super) fn parse_bounds(value: &str) -> Result<(f64, f64)> {
    let invalid = || {
        XisfError::InvalidFormat(format!(
            "Invalid bounds attribute \"{value}\": \
             it is not two finite numbers lower:upper with lower <= upper"
        ))
    };
    let (lower, upper) = value.split_once(':').ok_or_else(invalid)?;
    let finite = |text: &str| text.parse::<f64>().ok().filter(|v| v.is_finite());
    match (finite(lower), finite(upper)) {
        (Some(lower), Some(upper)) if lower <= upper => Ok((lower, upper)),
        _ => Err(invalid()),
    }
}
