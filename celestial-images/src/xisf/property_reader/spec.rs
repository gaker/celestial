use super::scalar::scalar_value;
use crate::xisf::data_block::read_stored;
use crate::xisf::header::checksum::Checksum;
use crate::xisf::header::location::{DataLocation, LocationAttr};
use crate::xisf::header::property::{XisfProperty, XisfPropertyValue};
use crate::xisf::header::text_encoding::TextEncoding;
use crate::xisf::reader::attribute_text;
use quick_xml::events::BytesStart;
use std::collections::HashMap;
use std::io::{Read, Seek};

// A Property child of an Image element. Values in an attached data block
// stay there until the caller reads them.
#[derive(Debug, Clone)]
pub(in crate::xisf) struct PropertySpec {
    id: String,
    value: SpecValue,
}

#[derive(Debug, Clone)]
enum SpecValue {
    Ready(XisfPropertyValue),
    Block {
        shape: Shape,
        location: DataLocation,
        big_endian: bool,
        checksum: Option<Checksum>,
    },
}

#[derive(Debug, Clone, Copy)]
enum Shape {
    Vector(usize),
    Matrix(usize, usize),
    Utf8,
}

// What the content of a Property element holds (XISF 1.0 sections 10.3 and
// 11.1.6).
pub(super) enum Content {
    Text,
    Inline(TextEncoding),
    Embedded,
    Nothing,
}

#[derive(Default)]
pub(super) struct PropertyAttrs {
    id: Option<String>,
    type_name: Option<String>,
    value: Option<String>,
    length: Option<String>,
    rows: Option<String>,
    columns: Option<String>,
    location: Option<String>,
    compressed: bool,
    byte_order: Option<String>,
    checksum: Option<String>,
    // The decoded block of an inline or embedded location.
    pub(super) header_block: Option<DataLocation>,
}

impl PropertyAttrs {
    // A String without a location holds its value as text. A malformed
    // location leaves the content unused, and the property unavailable.
    pub(super) fn content(&self) -> Content {
        let Some(location) = self.location.as_deref() else {
            let is_string = self.type_name.as_deref() == Some("String");
            return if is_string {
                Content::Text
            } else {
                Content::Nothing
            };
        };
        match LocationAttr::parse(location) {
            Ok(LocationAttr::Inline(encoding)) => Content::Inline(encoding),
            Ok(LocationAttr::Embedded) => Content::Embedded,
            _ => Content::Nothing,
        }
    }

    pub(super) fn mark_compressed(&mut self) {
        self.compressed = true;
    }
}

// Attributes that do not parse leave the property unavailable.
pub(super) fn property_attrs(element: &BytesStart) -> Option<PropertyAttrs> {
    let mut attrs = PropertyAttrs::default();
    for attr in element.attributes() {
        let attr = attr.ok()?;
        let text = Some(attribute_text(&attr).ok()?);
        match attr.key.as_ref() {
            "id" => attrs.id = text,
            "type" => attrs.type_name = text,
            "value" => attrs.value = text,
            "length" => attrs.length = text,
            "rows" => attrs.rows = text,
            "columns" => attrs.columns = text,
            "location" => attrs.location = text,
            "compression" => attrs.compressed = true,
            "byteOrder" => attrs.byte_order = text,
            "checksum" => attrs.checksum = text,
            _ => {}
        }
    }
    Some(attrs)
}

// Returns None for a property this reader skips: a malformed one, a type
// this crate has no value for, a data block that is compressed or in
// another file, or a checksum this reader cannot check. XISF 1.0 section 7
// makes such a property unavailable while the rest of the unit stays
// readable. `text` is None when a String property's content is not plain
// text.
pub(super) fn property_spec(attrs: &PropertyAttrs, text: Option<String>) -> Option<PropertySpec> {
    let id = attrs.id.clone()?;
    let value = match attrs.type_name.as_deref()? {
        "F64Vector" => block_spec(attrs, Shape::Vector(count(&attrs.length)?))?,
        "F64Matrix" => {
            let shape = Shape::Matrix(count(&attrs.rows)?, count(&attrs.columns)?);
            block_spec(attrs, shape)?
        }
        "String" if attrs.location.is_some() => block_spec(attrs, Shape::Utf8)?,
        type_name => {
            let value = attrs.value.as_deref().unwrap_or("");
            SpecValue::Ready(scalar_value(type_name, value, text)?)
        }
    };
    Some(PropertySpec { id, value })
}

// Dimensions are plain decimal digits. Rust's integer parser would also
// take a leading '+', which the data block sizes refuse too.
fn count(text: &Option<String>) -> Option<usize> {
    let digits = text.as_deref()?.trim();
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    digits.parse().ok()
}

fn block_spec(attrs: &PropertyAttrs, shape: Shape) -> Option<SpecValue> {
    if attrs.compressed {
        return None;
    }
    let location = match LocationAttr::parse(attrs.location.as_deref()?).ok()? {
        LocationAttr::Attachment(location) => location,
        _ => attrs.header_block.clone()?,
    };
    let big_endian = match attrs.byte_order.as_deref() {
        None | Some("little") => false,
        Some("big") => true,
        Some(_) => return None,
    };
    let checksum = match attrs.checksum.as_deref() {
        Some(value) => Some(Checksum::parse(value).ok()?),
        None => None,
    };
    shape.fits(location.size).then_some(SpecValue::Block {
        shape,
        location,
        big_endian,
        checksum,
    })
}

impl Shape {
    // A string fits a block of any length. Vectors and matrices hold
    // Float64 values.
    fn fits(self, size: u64) -> bool {
        let count = match self {
            Self::Vector(length) => Some(length),
            Self::Matrix(rows, columns) => rows.checked_mul(columns),
            Self::Utf8 => return true,
        };
        let bytes = count.and_then(|n| n.checked_mul(8));
        bytes.and_then(|b| u64::try_from(b).ok()) == Some(size)
    }
}

// A data block that cannot be read or fails its checksum leaves its
// property unavailable. Each distinct attached block is read once, but every
// property decodes its own copy, so one block shared by many properties could
// decode to many times the file length. The copies of attached blocks are
// limited to `stream_len` bytes in all, and a property past that budget is
// unavailable too. Inline and embedded blocks are not counted: each one is
// decoded from its own header text, so their total is bounded by the header.
pub(in crate::xisf) fn read_properties<R: Read + Seek>(
    reader: &mut R,
    stream_len: u64,
    specs: &[PropertySpec],
) -> Vec<XisfProperty> {
    let mut blocks: Blocks = HashMap::new();
    let mut budget = stream_len;
    let mut property = |spec: &PropertySpec| {
        let value = match &spec.value {
            SpecValue::Ready(value) => value.clone(),
            SpecValue::Block {
                shape,
                location,
                big_endian,
                checksum,
            } => {
                let cost = match location.header_bytes() {
                    Some(_) => 0,
                    None => location.size,
                };
                let remaining = budget.checked_sub(cost)?;
                let bytes = stored_bytes(&mut blocks, reader, stream_len, location)?;
                if let Some(checksum) = checksum {
                    checksum.verify(bytes).ok()?;
                }
                let value = block_value(bytes, *shape, *big_endian)?;
                budget = remaining;
                value
            }
        };
        Some(XisfProperty::new(spec.id.clone(), value))
    };
    specs.iter().filter_map(&mut property).collect()
}

// Attached blocks read so far, by offset and size. None marks a block that
// could not be read.
type Blocks = HashMap<(u64, u64), Option<Vec<u8>>>;

// The bytes of an inline or embedded block are already in the location.
fn stored_bytes<'a, R: Read + Seek>(
    blocks: &'a mut Blocks,
    reader: &mut R,
    stream_len: u64,
    location: &'a DataLocation,
) -> Option<&'a [u8]> {
    if let Some(bytes) = location.header_bytes() {
        return Some(bytes);
    }
    blocks
        .entry((location.offset, location.size))
        .or_insert_with(|| read_stored(reader, stream_len, location).ok())
        .as_deref()
}

fn block_value(bytes: &[u8], shape: Shape, big_endian: bool) -> Option<XisfPropertyValue> {
    match shape {
        Shape::Vector(_) => Some(XisfPropertyValue::F64Vector(f64_values(bytes, big_endian))),
        Shape::Matrix(rows, cols) => Some(XisfPropertyValue::F64Matrix {
            rows,
            cols,
            data: f64_values(bytes, big_endian),
        }),
        Shape::Utf8 => utf8_value(bytes),
    }
}

fn f64_values(bytes: &[u8], big_endian: bool) -> Vec<f64> {
    let decode = if big_endian {
        f64::from_be_bytes
    } else {
        f64::from_le_bytes
    };
    bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|&b| decode(b))
        .collect()
}

// A string is UTF-8 without U+0000 (XISF 1.0 section 8.4.4.3).
fn utf8_value(bytes: &[u8]) -> Option<XisfPropertyValue> {
    let text = String::from_utf8(bytes.to_vec()).ok()?;
    (!text.contains('\0')).then_some(XisfPropertyValue::String(text))
}
