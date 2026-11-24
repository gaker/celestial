mod scalar;
pub(super) mod spec;

use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::location::DataLocation;
use crate::xisf::header_text::DataElement;
use quick_xml::escape::resolve_predefined_entity;
use quick_xml::events::{BytesRef, BytesStart, Event};
use quick_xml::{Reader, XmlVersion};
use spec::{property_attrs, property_spec, Content, PropertyAttrs, PropertySpec};

// The open Image element and its properties.
struct ImageScope {
    depth: usize,
    specs: Vec<PropertySpec>,
}

#[derive(Default)]
struct Scan {
    depth: usize,
    image: Option<ImageScope>,
    images: Vec<Vec<PropertySpec>>,
}

// Returns the properties of each Image element, in the order the header
// parser lists the images. Properties elsewhere, such as in Metadata, are
// not image properties and are skipped.
pub(super) fn image_property_specs(xml: &str) -> Result<Vec<Vec<PropertySpec>>> {
    let mut reader = Reader::from_str(xml);
    let mut scan = Scan::default();
    loop {
        match reader.read_event().map_err(xml_error)? {
            Event::Start(e) => scan.start(&e, &mut reader)?,
            Event::Empty(e) => scan.empty(&e),
            Event::End(e) => scan.end(e.name().as_ref() == "Image"),
            Event::Eof => return Ok(scan.images),
            _ => {}
        }
    }
}

impl Scan {
    fn start(&mut self, element: &BytesStart, reader: &mut Reader<&[u8]>) -> Result<()> {
        if element.name().as_ref() == "Image" {
            self.image = Some(ImageScope {
                depth: self.depth,
                specs: Vec::new(),
            });
        } else if self.is_image_property(element) {
            if let Some(attrs) = property_attrs(element) {
                return self.property(attrs, reader);
            }
        }
        self.depth += 1;
        Ok(())
    }

    // The content readers consume the closing tag, so depth is unchanged
    // after them.
    fn property(&mut self, mut attrs: PropertyAttrs, reader: &mut Reader<&[u8]>) -> Result<()> {
        let text = match attrs.content() {
            Content::Text => element_text(reader)?,
            Content::Inline(encoding) => {
                let text = element_text(reader)?;
                let bytes = text.and_then(|t| encoding.decode(&t).ok());
                attrs.header_block = bytes.map(|b| DataLocation::inline(encoding, b));
                None
            }
            Content::Embedded => {
                let bytes = embedded_data(reader, &mut attrs)?;
                attrs.header_block = bytes.map(DataLocation::embedded);
                None
            }
            Content::Nothing => {
                self.depth += 1;
                None
            }
        };
        self.push(attrs, text);
        Ok(())
    }

    // An empty element has no text and no Data child.
    fn empty(&mut self, element: &BytesStart) {
        if element.name().as_ref() == "Image" {
            self.images.push(Vec::new());
        } else if self.is_image_property(element) {
            let Some(mut attrs) = property_attrs(element) else {
                return;
            };
            let text = match attrs.content() {
                Content::Text => Some(String::new()),
                Content::Inline(encoding) => {
                    let bytes = encoding.decode("").ok();
                    attrs.header_block = bytes.map(|b| DataLocation::inline(encoding, b));
                    None
                }
                Content::Embedded | Content::Nothing => None,
            };
            self.push(attrs, text);
        }
    }

    fn end(&mut self, is_image: bool) {
        self.depth = self.depth.saturating_sub(1);
        let closes_image = self.image.as_ref().is_some_and(|i| i.depth == self.depth);
        if is_image && closes_image {
            if let Some(scope) = self.image.take() {
                self.images.push(scope.specs);
            }
        }
    }

    fn is_image_property(&self, element: &BytesStart) -> bool {
        element.name().as_ref() == "Property"
            && self
                .image
                .as_ref()
                .is_some_and(|i| i.depth + 1 == self.depth)
    }

    fn push(&mut self, attrs: PropertyAttrs, text: Option<String>) {
        let spec = property_spec(&attrs, text);
        if let (Some(spec), Some(scope)) = (spec, self.image.as_mut()) {
            scope.specs.push(spec);
        }
    }
}

// Reads the children of a Property element with an embedded block, with
// its closing tag, and returns the decoded text of its one Data child (XISF
// 1.0 section 10.3). A missing, repeated or undecodable Data element gives
// None. A compressed block marks the property compressed, which makes it
// unavailable.
fn embedded_data(reader: &mut Reader<&[u8]>, attrs: &mut PropertyAttrs) -> Result<Option<Vec<u8>>> {
    let mut data = None;
    let mut data_elements = 0usize;
    let mut depth = 0usize;
    loop {
        match reader.read_event().map_err(xml_error)? {
            Event::Start(e) if depth == 0 && e.name().as_ref() == "Data" => {
                let text = element_text(reader)?;
                data_elements += 1;
                data = decode_data(&e, text, attrs);
            }
            Event::Empty(e) if depth == 0 && e.name().as_ref() == "Data" => {
                data_elements += 1;
                data = decode_data(&e, Some(String::new()), attrs);
            }
            Event::Start(_) => depth += 1,
            Event::End(_) if depth == 0 => return Ok(data.filter(|_| data_elements == 1)),
            Event::End(_) => depth -= 1,
            Event::Eof => return Err(xml_error("a Property element is not closed")),
            _ => {}
        }
    }
}

// The Data element is read as for an image block, so a property ignores its
// checksum and byteOrder attributes too.
fn decode_data(
    data: &BytesStart,
    text: Option<String>,
    attrs: &mut PropertyAttrs,
) -> Option<Vec<u8>> {
    let element = DataElement::parse(data).ok()?;
    if !element.block.is_plain() {
        attrs.mark_compressed();
    }
    element.decode(text).ok()
}

// XISF 1.0 section 11.1: the character content of a String property is its
// value, with whitespace kept. Returns None when the content holds markup or
// a reference to no character. The closing tag is consumed either way.
pub(super) fn element_text(reader: &mut Reader<&[u8]>) -> Result<Option<String>> {
    let mut text = Some(String::new());
    let mut depth = 0usize;
    loop {
        match reader.read_event().map_err(xml_error)? {
            Event::Text(t) => append(&mut text, Some(t.xml_content(XmlVersion::Implicit1_0))),
            Event::CData(c) => append(&mut text, Some(c.xml_content(XmlVersion::Implicit1_0))),
            Event::GeneralRef(r) => append(&mut text, resolve_reference(&r)),
            Event::Start(_) => {
                depth += 1;
                text = None;
            }
            Event::Empty(_) => text = None,
            Event::End(_) if depth == 0 => return Ok(text),
            Event::End(_) => depth -= 1,
            Event::Eof => return Err(xml_error("a String property is not closed")),
            _ => {}
        }
    }
}

fn append(text: &mut Option<String>, part: Option<impl AsRef<str>>) {
    match (text.as_mut(), part) {
        (Some(text), Some(part)) => text.push_str(part.as_ref()),
        _ => *text = None,
    }
}

fn resolve_reference(reference: &BytesRef) -> Option<String> {
    if let Some(ch) = reference.resolve_char_ref().ok()? {
        return Some(ch.to_string());
    }
    resolve_predefined_entity(reference).map(str::to_string)
}

fn xml_error(error: impl std::fmt::Display) -> XisfError {
    XisfError::XmlParse(error.to_string())
}

#[cfg(test)]
mod tests {
    use super::spec::read_properties;
    use super::*;
    use crate::xisf::header::property::{XisfProperty, XisfPropertyValue};
    use std::io::Cursor;

    fn read(xml: &str, file: &[u8]) -> Result<Vec<Vec<XisfProperty>>> {
        let mut cursor = Cursor::new(file.to_vec());
        let specs = image_property_specs(xml)?;
        Ok(specs
            .iter()
            .map(|specs| read_properties(&mut cursor, file.len() as u64, specs))
            .collect())
    }

    fn image(children: &str) -> String {
        format!(
            r#"<xisf><Image geometry="1:1:1" sampleFormat="UInt8" location="attachment:0:1">{children}</Image></xisf>"#
        )
    }

    #[test]
    fn scalar_forms_of_section_8_3_are_read() {
        let xml = image(concat!(
            r#"<Property id="a" type="Boolean" value="true"/>"#,
            r#"<Property id="b" type="Boolean" value="0"/>"#,
            r#"<Property id="c" type="Int32" value=" -123 "/>"#,
            r#"<Property id="d" type="Int32" value="0x80E950AB"/>"#,
            r#"<Property id="e" type="Int32" value="0b101"/>"#,
            r#"<Property id="f" type="Int32" value="0o17"/>"#,
            r#"<Property id="g" type="Float64" value="+.123"/>"#,
            r#"<Property id="h" type="Float64" value="-Inf"/>"#,
        ));
        let props = read(&xml, &[0]).unwrap();
        let values: Vec<XisfPropertyValue> = props[0].iter().map(|p| p.value.clone()).collect();
        assert_eq!(
            values,
            vec![
                XisfPropertyValue::Boolean(true),
                XisfPropertyValue::Boolean(false),
                XisfPropertyValue::Int32(-123),
                XisfPropertyValue::Int32(-2132193109),
                XisfPropertyValue::Int32(5),
                XisfPropertyValue::Int32(15),
                XisfPropertyValue::Float64(0.123),
                XisfPropertyValue::Float64(f64::NEG_INFINITY),
            ]
        );
    }

    // The forms of XISF 1.0 section 8.3.3 that this writer wrote before it
    // wrote +Inf and -Inf are read too.
    #[test]
    fn float64_special_values_are_read_in_every_section_8_3_3_form() {
        let forms = ["+Inf", "-Inf", "inf", "-inf", "NaN", "nan", "-nan"];
        let children: String = forms
            .iter()
            .map(|f| format!(r#"<Property id="{f}" type="Float64" value="{f}"/>"#))
            .collect();
        let props = read(&image(&children), &[0]).unwrap();
        let bits: Vec<u64> = props[0]
            .iter()
            .map(|p| match p.value {
                XisfPropertyValue::Float64(v) if v.is_nan() => f64::NAN.to_bits(),
                XisfPropertyValue::Float64(v) => v.to_bits(),
                _ => 0,
            })
            .collect();
        let inf = f64::INFINITY.to_bits();
        let neg = f64::NEG_INFINITY.to_bits();
        let nan = f64::NAN.to_bits();
        assert_eq!(bits, [inf, neg, inf, neg, nan, nan, nan]);
    }

    #[test]
    fn string_content_keeps_whitespace_and_decodes_references() {
        let xml = image(
            r#"<Property id="s" type="String">  a &lt; b &amp; &#x41;&#66; <![CDATA[<c>]]> </Property>"#,
        );
        let props = read(&xml, &[0]).unwrap();
        assert_eq!(
            props[0][0].value,
            XisfPropertyValue::String("  a < b & AB <c> ".to_string())
        );
    }

    #[test]
    fn big_endian_matrix_block_is_read() {
        let mut file = vec![0u8];
        for v in [1.5f64, -2.0, 0.25, 8.0] {
            file.extend_from_slice(&v.to_be_bytes());
        }
        let xml = image(
            r#"<Property id="m" type="F64Matrix" rows="2" columns="2" byteOrder="big" location="attachment:1:32"/>"#,
        );
        let props = read(&xml, &file).unwrap();
        assert_eq!(
            props[0][0].value,
            XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 2,
                data: vec![1.5, -2.0, 0.25, 8.0]
            }
        );
    }

    #[test]
    fn unsupported_types_and_blocks_are_skipped() {
        let xml = image(concat!(
            r#"<Property id="a" type="UInt16" value="3"/>"#,
            r#"<Property id="b" type="F64Vector" length="1" location="inline:base64">AAAAAAAA8D8=</Property>"#,
            r#"<Property id="c" type="F64Vector" length="1" compression="zlib:8" location="attachment:0:1"/>"#,
            r#"<Property id="d" type="Int32" value="7"/>"#,
        ));
        let props = read(&xml, &[0]).unwrap();
        assert_eq!(
            props[0],
            vec![
                XisfProperty::new("b", XisfPropertyValue::F64Vector(vec![1.0])),
                XisfProperty::int32("d", 7)
            ]
        );
    }

    fn values(xml: &str, file: &[u8]) -> Vec<(String, XisfPropertyValue)> {
        read(xml, file).unwrap()[0]
            .iter()
            .map(|p| (p.id.clone(), p.value.clone()))
            .collect()
    }

    // XISF 1.0 sections 10.3, 11.1.6 and 11.1.8. Encodings from Python's
    // base64 and struct modules.
    #[test]
    fn inline_and_embedded_property_blocks_are_read() {
        let xml = image(concat!(
            r#"<Property id="b64" type="F64Vector" length="1" location="inline:base64">AAAAAAAA8D8=</Property>"#,
            r#"<Property id="hex" type="F64Vector" length="1" location="inline:hex"> 000000000000f03f </Property>"#,
            r#"<Property id="m" type="F64Matrix" rows="2" columns="1" location="inline:base64">AAAAAAAA+D8AAAAAAAAAwA==</Property>"#,
            r#"<Property id="empty" type="F64Vector" length="0" location="inline:base64"></Property>"#,
            r#"<Property id="empty2" type="F64Vector" length="0" location="inline:hex"/>"#,
            r#"<Property id="emb" type="F64Vector" length="1" location="embedded"><Data encoding="hex">000000000000f03f</Data></Property>"#,
            r#"<Property id="s" type="String" location="inline:base64">VGhlIHF1aWNrIGJyb3duIGZveA==</Property>"#,
            r#"<Property id="s2" type="String" location="embedded"><Data encoding="base64">w4VuZ3N0csO2bQ==</Data></Property>"#,
            r#"<Property id="sum" type="F64Vector" length="1" location="inline:base64" checksum="sha1:3eafb06982d26cb8d0847b681f31c9055684eee2">AAAAAAAA8D8=</Property>"#,
        ));
        let string = |s: &str| XisfPropertyValue::String(s.to_string());
        let expected = vec![
            ("b64", XisfPropertyValue::F64Vector(vec![1.0])),
            ("hex", XisfPropertyValue::F64Vector(vec![1.0])),
            (
                "m",
                XisfPropertyValue::F64Matrix {
                    rows: 2,
                    cols: 1,
                    data: vec![1.5, -2.0],
                },
            ),
            ("empty", XisfPropertyValue::F64Vector(vec![])),
            ("empty2", XisfPropertyValue::F64Vector(vec![])),
            ("emb", XisfPropertyValue::F64Vector(vec![1.0])),
            ("s", string("The quick brown fox")),
            ("s2", string("Ångström")),
            ("sum", XisfPropertyValue::F64Vector(vec![1.0])),
        ];
        let expected: Vec<(String, XisfPropertyValue)> = expected
            .into_iter()
            .map(|(id, v)| (id.to_string(), v))
            .collect();
        assert_eq!(values(&xml, &[0]), expected);
    }

    // XISF 1.0 section 11.1.6: a String value may be a data block, which
    // holds UTF-8 without U+0000 (section 8.4.4.3).
    #[test]
    fn a_string_with_an_attached_block_reads_the_block() {
        let mut file = vec![0u8];
        file.extend_from_slice("Ångström".as_bytes());
        file.extend_from_slice(&[0xFF, 0xFE, b'a', 0, b'b']);
        // Room in the decoding budget for both copies of the shared block.
        file.resize(64, 0);
        let xml = image(concat!(
            r#"<Property id="s" type="String" location="attachment:1:10"/>"#,
            r#"<Property id="t" type="String" location="attachment:1:10">ignored</Property>"#,
            r#"<Property id="bad" type="String" location="attachment:11:2"/>"#,
            r#"<Property id="nul" type="String" location="attachment:13:3"/>"#,
        ));
        let string = XisfPropertyValue::String("Ångström".to_string());
        assert_eq!(
            values(&xml, &file),
            vec![("s".to_string(), string.clone()), ("t".to_string(), string)]
        );
    }

    #[test]
    fn malformed_inline_and_embedded_property_blocks_are_skipped() {
        let cases = [
            r#"<Property id="a" type="F64Vector" length="1" location="inline:base64">AAAAAAAA8D8</Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="inline:base64">AAAA</Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="inline:hex">000000000000f03g</Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="inline:base32">AAAAAAAA8D8=</Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="inline:base64">AAAA<b/>AAAA8D8=</Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="embedded"/>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="embedded"></Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="embedded"><Data>AAAAAAAA8D8=</Data></Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="embedded"><Data encoding="base64" compression="zlib:8">AAAAAAAA8D8=</Data></Property>"#,
            r#"<Property id="a" type="F64Vector" length="1" location="inline:base64" checksum="sha1:00">AAAAAAAA8D8=</Property>"#,
            r#"<Property id="a" type="String" location="inline:base64">//4=</Property>"#,
            r#"<Property id="a" type="String" location="inline:base64">YQBi</Property>"#,
            r#"<Property id="a" type="String" location="url(http://example.com/a.txt)"/>"#,
        ];
        for case in cases {
            let good = r#"<Property id="first" type="Int32" value="1"/>"#;
            let last = r#"<Property id="last" type="Int32" value="2"/>"#;
            let xml = image(&format!("{good}{case}{last}"));
            let ids: Vec<String> = values(&xml, &[0; 16])
                .into_iter()
                .map(|(id, _)| id)
                .collect();
            assert_eq!(ids, vec!["first", "last"], "{case}");
        }
    }

    // XISF 1.0 section 10.6 moves only the compression and subblocks
    // attributes to the Data element of an embedded block. Its checksum and
    // byteOrder stay on the Property element (section 10.5).
    #[test]
    fn checksum_and_byte_order_on_a_data_element_are_ignored() {
        let xml = image(concat!(
            r#"<Property id="ck" type="F64Vector" length="1" location="embedded">"#,
            r#"<Data encoding="base64" checksum="sha1:0000000000000000000000000000000000000000">AAAAAAAA8D8=</Data></Property>"#,
            r#"<Property id="bo" type="F64Vector" length="1" location="embedded">"#,
            r#"<Data encoding="base64" byteOrder="big">AAAAAAAA8D8=</Data></Property>"#,
            r#"<Property id="parent" type="F64Vector" length="1" location="embedded" checksum="sha1:0000000000000000000000000000000000000000">"#,
            r#"<Data encoding="base64">AAAAAAAA8D8=</Data></Property>"#,
        ));
        let one = XisfPropertyValue::F64Vector(vec![1.0]);
        assert_eq!(
            values(&xml, &[0]),
            vec![("ck".to_string(), one.clone()), ("bo".to_string(), one)]
        );
    }

    #[test]
    fn a_property_with_two_data_elements_is_skipped() {
        let xml = image(concat!(
            r#"<Property id="a" type="F64Vector" length="1" location="embedded">"#,
            r#"<Data encoding="base64">AAAAAAAA8D8=</Data><Data encoding="base64">AAAAAAAA8D8=</Data></Property>"#,
            r#"<Property id="b" type="Int32" value="2"/>"#,
        ));
        assert_eq!(
            values(&xml, &[0]),
            vec![("b".to_string(), XisfPropertyValue::Int32(2))]
        );
    }

    // Indices follow the images the header parser lists, which is every
    // Image element: one without a location is listed as unreadable.
    #[test]
    fn properties_outside_images_are_not_image_properties() {
        let xml = concat!(
            r#"<xisf><Metadata><Property id="m" type="Int32" value="1"/></Metadata>"#,
            r#"<Property id="root" type="Int32" value="2"/>"#,
            r#"<Image geometry="1:1:1"><Property id="x" type="Int32" value="3"/></Image>"#,
            r#"<Image geometry="1:1:1" location="attachment:0:1"><Property id="y" type="Int32" value="4"/></Image>"#,
            r#"<Image geometry="1:1:1" location="attachment:0:1"/></xisf>"#,
        );
        let props = read(xml, &[0]).unwrap();
        assert_eq!(
            props,
            vec![
                vec![XisfProperty::int32("x", 3)],
                vec![XisfProperty::int32("y", 4)],
                vec![]
            ]
        );
    }

    // XISF 1.0 section 7: an object a decoder cannot use is unavailable, and
    // the rest of the unit stays readable. Each case sits between two good
    // properties, and both of them are still read.
    #[test]
    fn malformed_properties_are_skipped() {
        let mut file = vec![0u8];
        file.extend_from_slice(&1.5f64.to_le_bytes());
        file.extend_from_slice(&(-2.5f64).to_le_bytes());
        let cases = [
            r#"<Property type="Int32" value="1"/>"#,
            r#"<Property id="a" value="1"/>"#,
            r#"<Property id="a" type="Boolean" value="True"/>"#,
            r#"<Property id="a" type="Boolean" value="yes"/>"#,
            r#"<Property id="a" type="Int32" value="2147483648"/>"#,
            r#"<Property id="a" type="Int32" value="0x1FFFFFFFF"/>"#,
            r#"<Property id="a" type="Int32" value="0x-1"/>"#,
            r#"<Property id="a" type="Float64"/>"#,
            r#"<Property id="a" type="Float64" value="1,5"/>"#,
            r#"<Property id="a" type="Float64" value="0x1p3"/>"#,
            r#"<Property id="a" type="F64Vector" length="2" location="attachment:0:8"/>"#,
            r#"<Property id="a" type="F64Vector" length="1" byteOrder="Big" location="attachment:1:8"/>"#,
            r#"<Property id="a" type="F64Vector" length="-2" location="attachment:1:16"/>"#,
            r#"<Property id="a" type="F64Vector" length="+2" location="attachment:1:16"/>"#,
            r#"<Property id="a" type="F64Vector" length="18446744073709551615" location="attachment:1:16"/>"#,
            r#"<Property id="a" type="F64Vector" length="2305843009213693954" location="attachment:1:16"/>"#,
            r#"<Property id="a" type="F64Vector" length="2" location="attachment:x:16"/>"#,
            r#"<Property id="a" type="F64Matrix" rows="4294967296" columns="4294967296" location="attachment:1:16"/>"#,
            r#"<Property id="a" type="F64Matrix" rows="2305843009213693952" columns="8" location="attachment:1:0"/>"#,
            r#"<Property id="a" type="String">&nbsp;</Property>"#,
            r#"<Property id="a" type="String">&#0;</Property>"#,
            r#"<Property id="a" type="String">&#xD800;</Property>"#,
            r#"<Property id="a" type="String">&#x110000;</Property>"#,
            r#"<Property id="a" type="String">&#99999999999999999999999;</Property>"#,
            r#"<Property id="a" type="String">&#x;</Property>"#,
            r#"<Property id="a" type="String">a<b/>c</Property>"#,
            r#"<Property id="a" type="String">a<b>x</b>c</Property>"#,
            r#"<Property id="a" id="b" type="Int32" value="1"/>"#,
            r#"<Property id="a" type="F64Vector" length="2" location="attachment:1:16"/>"#,
        ];
        for (i, case) in cases.iter().enumerate() {
            let good = r#"<Property id="first" type="Int32" value="1"/>"#;
            let last = r#"<Property id="last" type="Int32" value="2"/>"#;
            let xml = image(&format!("{good}{case}{last}"));
            let props = read(&xml, &file).unwrap();
            let ids: Vec<&str> = props[0].iter().map(|p| p.id.as_str()).collect();
            // The last case is well formed and reads.
            let expected = if i + 1 == cases.len() {
                vec!["first", "a", "last"]
            } else {
                vec!["first", "last"]
            };
            assert_eq!(ids, expected, "{case}");
        }
    }

    #[test]
    fn blocks_past_the_end_of_the_file_are_skipped() {
        for location in [
            "attachment:1:8",
            "attachment:18446744073709551615:8",
            "attachment:1:9223372036854775808",
        ] {
            let length = if location.ends_with(":8") {
                "1"
            } else {
                "1152921504606846976"
            };
            let xml = image(&format!(
                r#"<Property id="v" type="F64Vector" length="{length}" location="{location}"/><Property id="n" type="Int32" value="3"/>"#
            ));
            assert_eq!(
                read(&xml, &[0; 4]).unwrap(),
                vec![vec![XisfProperty::int32("n", 3)]]
            );
        }
    }

    // Each property decodes its own copy of a block, so the bound counts a
    // shared block once per property.
    #[test]
    fn properties_past_the_file_length_budget_are_skipped() {
        let mut file = vec![0u8];
        file.extend_from_slice(&1.5f64.to_le_bytes());
        file.extend_from_slice(&[0; 7]);
        let one = r#"<Property id="v" type="F64Vector" length="1" location="attachment:1:8"/>"#;
        let scalar = r#"<Property id="s" type="Int32" value="4"/>"#;
        let xml = image(&format!("{}{scalar}", one.repeat(3)));
        let ids: Vec<String> = read(&xml, &file).unwrap()[0]
            .iter()
            .map(|p| p.id.clone())
            .collect();
        assert_eq!(ids, vec!["v", "v", "s"]);
    }

    // A block past the end of the file is not read, so it uses none of the
    // budget.
    #[test]
    fn unreadable_blocks_use_none_of_the_budget() {
        let mut file = vec![0u8];
        file.extend_from_slice(&1.5f64.to_le_bytes());
        let past = r#"<Property id="p" type="F64Vector" length="1" location="attachment:9:8"/>"#;
        let one = r#"<Property id="v" type="F64Vector" length="1" location="attachment:1:8"/>"#;
        let xml = image(&format!("{past}{one}"));
        let ids: Vec<String> = read(&xml, &file).unwrap()[0]
            .iter()
            .map(|p| p.id.clone())
            .collect();
        assert_eq!(ids, vec!["v"]);
    }

    struct CountingReader {
        inner: Cursor<Vec<u8>>,
        bytes_read: std::rc::Rc<std::cell::Cell<usize>>,
    }

    impl std::io::Read for CountingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let n = self.inner.read(buf)?;
            self.bytes_read.set(self.bytes_read.get() + n);
            Ok(n)
        }
    }

    impl std::io::Seek for CountingReader {
        fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
            self.inner.seek(pos)
        }
    }

    #[test]
    fn a_shared_block_is_read_once() {
        let mut file = vec![0u8; 64];
        file[1..9].copy_from_slice(&1.5f64.to_le_bytes());
        file[9..17].copy_from_slice(&(-2.0f64).to_le_bytes());
        let xml = image(concat!(
            r#"<Property id="a" type="F64Vector" length="2" location="attachment:1:16"/>"#,
            r#"<Property id="b" type="F64Matrix" rows="1" columns="2" location="attachment:1:16"/>"#,
            r#"<Property id="c" type="F64Vector" length="2" location="attachment:1:16"/>"#,
        ));
        let bytes_read = std::rc::Rc::new(std::cell::Cell::new(0));
        let mut reader = CountingReader {
            inner: Cursor::new(file),
            bytes_read: bytes_read.clone(),
        };
        let specs = image_property_specs(&xml).unwrap();
        let props = read_properties(&mut reader, 64, &specs[0]);
        assert_eq!(bytes_read.get(), 16);
        let values: Vec<XisfPropertyValue> = props.into_iter().map(|p| p.value).collect();
        assert_eq!(
            values,
            vec![
                XisfPropertyValue::F64Vector(vec![1.5, -2.0]),
                XisfPropertyValue::F64Matrix {
                    rows: 1,
                    cols: 2,
                    data: vec![1.5, -2.0]
                },
                XisfPropertyValue::F64Vector(vec![1.5, -2.0]),
            ]
        );
    }

    // XISF 1.0 section 10.5. A property whose block fails its checksum, or
    // whose checksum cannot be checked, is unavailable.
    #[test]
    fn property_block_checksums_are_verified() {
        let mut file = vec![0u8];
        file.extend_from_slice(&1.5f64.to_le_bytes());
        // Room in the decoding budget for every property.
        file.resize(64, 0);
        let sha1 = "dd58a0f3fcf1993ddf7bacbaaa2507c430793478";
        let sha256 = "e163f8cb0f7067a7fc78ca859a77f849aea3214f38fb75b884e4a16be725c905";
        let vector = |id: &str, checksum: &str| {
            format!(
                r#"<Property id="{id}" type="F64Vector" length="1" location="attachment:1:8" checksum="{checksum}"/>"#
            )
        };
        let xml = image(
            &[
                vector("sha1", &format!("sha1:{sha1}")),
                vector("sha256", &format!("sha-256:{sha256}")),
                vector("mismatch", &format!("sha-256:{}", &sha1.repeat(2)[..64])),
                vector("short", "sha1:dd58"),
                vector("sha3", &format!("sha3-256:{sha256}")),
                vector("other", &format!("md5:{sha1}")),
            ]
            .concat(),
        );
        let ids: Vec<String> = read(&xml, &file).unwrap()[0]
            .iter()
            .map(|p| p.id.clone())
            .collect();
        assert_eq!(ids, vec!["sha1", "sha256"]);
    }

    // Markup that is not well formed makes the whole header unreadable.
    #[test]
    fn ill_formed_xml_is_an_error() {
        for child in [
            r#"<Property id="s" type="String"><![CDATA[abc</Property>"#,
            r#"<Property id="s" type="String">a & b</Property>"#,
            r#"<Property id="x" type="Int32" value="1">"#,
        ] {
            assert!(
                matches!(
                    image_property_specs(&image(child)),
                    Err(XisfError::XmlParse(_))
                ),
                "{child}"
            );
        }
    }
}
