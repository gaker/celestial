use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::block_attrs::BlockAttrs;
use crate::xisf::header::image_info::{ImageInfo, XisfHeader};
use crate::xisf::header::location::{DataLocation, LocationAttr};
use crate::xisf::header::text_encoding::TextEncoding;
use crate::xisf::keyword_reader::KeywordScan;
use crate::xisf::property_reader::element_text;
use crate::xisf::reader::attribute_text;
use quick_xml::events::BytesStart;
use quick_xml::Reader;

// State of the header parser. Image data blocks serialized in the header
// (XISF 1.0 section 10.3) are read here: the text of an Image element with
// an inline block, or of the Data child of one with an embedded block.
pub(super) struct HeaderScan {
    pub(super) header: XisfHeader,
    pub(super) image: Option<ImageInfo>,
    pub(super) keywords: KeywordScan,
    depth: usize,
    image_depth: usize,
    embedded: bool,
    data_elements: usize,
}

impl HeaderScan {
    pub(super) fn new() -> Self {
        let header = XisfHeader {
            version: String::new(),
            images: Vec::new(),
            keywords: Vec::new(),
            root_keywords: Vec::new(),
        };
        Self {
            header,
            image: None,
            keywords: KeywordScan::default(),
            depth: 0,
            image_depth: 0,
            embedded: false,
            data_elements: 0,
        }
    }

    pub(super) fn targets(&mut self) -> (&mut String, &mut Option<ImageInfo>) {
        (&mut self.header.version, &mut self.image)
    }

    // Whether an element read now, before `enter`, is a child of the open
    // Image element.
    pub(super) fn in_image(&self) -> bool {
        self.image.is_some() && self.depth == self.image_depth
    }

    // Whether an element read now, before `enter`, is a child of the root
    // element.
    pub(super) fn at_root(&self) -> bool {
        self.depth == 1
    }

    // Lists the open image with its keywords.
    pub(super) fn finish_image(&mut self) {
        if let Some(image) = self.image.take() {
            self.header.images.push(image);
            self.keywords.finish_image();
        }
    }

    // `header_len` is the length of the XML header in bytes.
    pub(super) fn finish(mut self, header_len: usize) -> Result<XisfHeader> {
        let (images, root) = self.keywords.finish(header_len)?;
        self.header.keywords = images;
        self.header.root_keywords = root;
        Ok(self.header)
    }

    pub(super) fn enter(&mut self) {
        self.depth += 1;
    }

    pub(super) fn leave(&mut self) {
        self.depth = self.depth.saturating_sub(1);
    }

    // Only a Data child of the Image element holds its block. A Data
    // element deeper down belongs to a property or another child.
    pub(super) fn is_image_data(&self, element: &BytesStart) -> bool {
        self.embedded && self.depth == self.image_depth && element.name().as_ref() == "Data"
    }

    // Called after the start tag of an Image element is parsed. The text of
    // an inline block is read here, with the end tag.
    pub(super) fn open_image(
        &mut self,
        element: &BytesStart,
        reader: &mut Reader<&[u8]>,
    ) -> Result<()> {
        self.image_depth = self.depth;
        self.embedded = false;
        self.data_elements = 0;
        self.keywords.open_image();
        match text_block(element)? {
            Some(LocationAttr::Inline(encoding)) => {
                let text = element_text(reader)?;
                self.leave();
                if let Some(image) = self.image.as_mut() {
                    set_inline(image, encoding, text);
                }
                self.finish_image();
            }
            Some(_) => self.embedded = true,
            None => {}
        }
        Ok(())
    }

    // Reads the Data element's text, with its end tag.
    pub(super) fn read_data(
        &mut self,
        data: &BytesStart,
        reader: &mut Reader<&[u8]>,
    ) -> Result<()> {
        let text = element_text(reader)?;
        self.set_data(data, text)
    }

    // XISF 1.0 section 10.3 gives an embedded block one Data child.
    pub(super) fn set_data(&mut self, data: &BytesStart, text: Option<String>) -> Result<()> {
        self.data_elements += 1;
        let Some(image) = self.image.as_mut() else {
            return Ok(());
        };
        if self.data_elements > 1 {
            record(
                image,
                invalid("An embedded data block has more than one Data element"),
            );
            return Ok(());
        }
        set_embedded(image, data, text)
    }

    // Called at the end tag of an Image element, before it is listed.
    pub(super) fn close_image(&mut self) {
        let missing = std::mem::take(&mut self.embedded) && self.data_elements == 0;
        if let (true, Some(image)) = (missing, self.image.as_mut()) {
            record(image, missing_data());
        }
    }

    // Called for an empty Image element, before it is listed. An inline
    // block then has no text, and an embedded one has no Data element.
    pub(super) fn close_empty_image(&mut self, element: &BytesStart) -> Result<()> {
        let Some(image) = self.image.as_mut() else {
            return Ok(());
        };
        match text_block(element)? {
            Some(LocationAttr::Inline(encoding)) => {
                set_inline(image, encoding, Some(String::new()))
            }
            Some(_) => record(image, missing_data()),
            None => {}
        }
        Ok(())
    }
}

// The Data child of an embedded block. XISF 1.0 section 10.6 moves the
// compression and subblocks attributes of the block to it. The checksum
// and byteOrder attributes stay on the element that has the block (section
// 10.5), so they are ignored here.
pub(super) struct DataElement {
    encoding: Option<String>,
    pub(super) block: BlockAttrs,
}

impl DataElement {
    pub(super) fn parse(data: &BytesStart) -> Result<Self> {
        let mut element = Self {
            encoding: None,
            block: BlockAttrs::default(),
        };
        for attr in data.attributes() {
            let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
            let value = attribute_text(&attr)?;
            match attr.key.as_ref() {
                "encoding" => element.encoding = Some(value),
                key @ ("compression" | "subblocks") => element.block.apply(key, &value),
                _ => {}
            }
        }
        Ok(element)
    }

    // `text` is None when the element holds markup.
    pub(super) fn decode(&self, text: Option<String>) -> Result<Vec<u8>> {
        let encoding = self.encoding.as_deref().ok_or_else(|| {
            invalid("The Data element of an embedded data block has no encoding attribute")
        })?;
        let encoding = TextEncoding::parse(encoding)?;
        let text =
            text.ok_or_else(|| invalid("The text of an embedded data block holds XML markup"))?;
        encoding.decode(&text)
    }
}

// The location of an Image element whose block is in the header. A
// malformed location is reported when the block is read.
fn text_block(element: &BytesStart) -> Result<Option<LocationAttr>> {
    for attr in element.attributes() {
        let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
        if attr.key.as_ref() == "location" {
            let value = attribute_text(&attr)?;
            return Ok(LocationAttr::parse(&value)
                .ok()
                .filter(|l| !matches!(l, LocationAttr::Attachment(_))));
        }
    }
    Ok(None)
}

fn set_inline(image: &mut ImageInfo, encoding: TextEncoding, text: Option<String>) {
    let markup = || invalid("The text of an inline data block holds XML markup");
    match text.ok_or_else(markup).and_then(|t| encoding.decode(&t)) {
        Ok(bytes) => image.location = DataLocation::inline(encoding, bytes),
        Err(error) => record(image, error),
    }
}

fn set_embedded(image: &mut ImageInfo, data: &BytesStart, text: Option<String>) -> Result<()> {
    let mut element = DataElement::parse(data)?;
    match element.decode(text) {
        Ok(bytes) => image.location = DataLocation::embedded(bytes),
        Err(error) => element.block.record(error),
    }
    element.block.apply_to(image);
    Ok(())
}

fn record(image: &mut ImageInfo, error: XisfError) {
    let mut attrs = BlockAttrs::default();
    attrs.record(error);
    attrs.apply_to(image);
}

fn missing_data() -> XisfError {
    invalid("An embedded data block has no Data element")
}

fn invalid(message: &str) -> XisfError {
    XisfError::InvalidFormat(message.to_string())
}
