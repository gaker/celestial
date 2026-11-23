use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::text_encoding::TextEncoding;
use std::fmt;

#[derive(Debug, Clone)]
pub(crate) struct DataLocation {
    pub offset: u64,
    pub size: u64,
    pub(crate) header_block: Option<HeaderBlock>,
}

// The decoded bytes of an inline or embedded block (XISF 1.0 section 10.3),
// with the location attribute that describes the block. The `offset` and
// `size` of its DataLocation address these bytes, not the file.
#[derive(Clone)]
pub(crate) struct HeaderBlock {
    location: &'static str,
    bytes: Vec<u8>,
}

// The bytes can be as long as the header, so only their count is printed.
impl fmt::Debug for HeaderBlock {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({} bytes)", self.location, self.bytes.len())
    }
}

// The value of a location attribute. Blocks in other files (`url(...)` and
// `path(...)`) belong to distributed units and are not supported.
#[derive(Debug, Clone)]
pub(crate) enum LocationAttr {
    Attachment(DataLocation),
    Inline(TextEncoding),
    Embedded,
}

impl LocationAttr {
    pub(crate) fn parse(s: &str) -> Result<Self> {
        if let Some(encoding) = s.strip_prefix("inline:") {
            return TextEncoding::parse(encoding).map(Self::Inline);
        }
        if s == "embedded" {
            return Ok(Self::Embedded);
        }
        DataLocation::parse(s).map(Self::Attachment)
    }
}

impl DataLocation {
    pub(crate) fn new(offset: u64, size: u64) -> Self {
        Self {
            offset,
            size,
            header_block: None,
        }
    }

    pub(crate) fn inline(encoding: TextEncoding, bytes: Vec<u8>) -> Self {
        Self::in_header(encoding.inline_location(), bytes)
    }

    pub(crate) fn embedded(bytes: Vec<u8>) -> Self {
        Self::in_header("embedded", bytes)
    }

    fn in_header(location: &'static str, bytes: Vec<u8>) -> Self {
        Self {
            offset: 0,
            size: bytes.len() as u64,
            header_block: Some(HeaderBlock { location, bytes }),
        }
    }

    // The decoded bytes of an inline or embedded block, None for an
    // attachment.
    pub(crate) fn header_bytes(&self) -> Option<&[u8]> {
        self.header_block
            .as_ref()
            .map(|block| block.bytes.as_slice())
    }

    pub(crate) fn parse(s: &str) -> Result<Self> {
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

        Ok(Self::new(offset, size))
    }

    pub(crate) fn format(&self) -> String {
        match &self.header_block {
            Some(block) => block.location.to_string(),
            None => format!("attachment:{}:{}", self.offset, self.size),
        }
    }

    pub(crate) fn ends_within(&self, stream_len: u64) -> bool {
        self.offset
            .checked_add(self.size)
            .is_some_and(|end| end <= stream_len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_stores_fields_as_given() {
        let loc = DataLocation::new(1024, 2048);
        assert_eq!(loc.offset, 1024);
        assert_eq!(loc.size, 2048);
    }

    #[test]
    fn parse_valid_attachment() {
        let loc = DataLocation::parse("attachment:1024:2048").unwrap();
        assert_eq!(loc.offset, 1024);
        assert_eq!(loc.size, 2048);
    }

    #[test]
    fn parse_accepts_zero_offset_and_size() {
        let loc = DataLocation::parse("attachment:0:0").unwrap();
        assert_eq!(loc.offset, 0);
        assert_eq!(loc.size, 0);
    }

    #[test]
    fn parse_accepts_u64_max() {
        let max = u64::MAX;
        let s = format!("attachment:{max}:{max}");
        let loc = DataLocation::parse(&s).unwrap();
        assert_eq!(loc.offset, max);
        assert_eq!(loc.size, max);
    }

    #[test]
    fn parse_rejects_non_attachment_scheme() {
        assert!(DataLocation::parse("embedded:100:200").is_err());
        assert!(DataLocation::parse("inline:100:200").is_err());
        assert!(DataLocation::parse("ATTACHMENT:100:200").is_err());
    }

    #[test]
    fn parse_rejects_wrong_segment_count() {
        assert!(DataLocation::parse("").is_err());
        assert!(DataLocation::parse("attachment").is_err());
        assert!(DataLocation::parse("attachment:100").is_err());
        assert!(DataLocation::parse("attachment:100:200:300").is_err());
        assert!(DataLocation::parse("100:200").is_err());
    }

    #[test]
    fn parse_rejects_non_numeric_segments() {
        assert!(DataLocation::parse("attachment:abc:200").is_err());
        assert!(DataLocation::parse("attachment:100:xyz").is_err());
        assert!(DataLocation::parse("attachment::200").is_err());
        assert!(DataLocation::parse("attachment:100:").is_err());
    }

    #[test]
    fn parse_rejects_negative_numbers() {
        assert!(DataLocation::parse("attachment:-1:200").is_err());
        assert!(DataLocation::parse("attachment:100:-1").is_err());
    }

    #[test]
    fn parse_rejects_overflow() {
        let overflow = "18446744073709551616";
        let s = format!("attachment:{overflow}:0");
        assert!(DataLocation::parse(&s).is_err());
    }

    #[test]
    fn format_produces_parseable_string() {
        let original = DataLocation::new(512, 65536);
        let formatted = original.format();
        assert_eq!(formatted, "attachment:512:65536");
        let reparsed = DataLocation::parse(&formatted).unwrap();
        assert_eq!(reparsed.offset, original.offset);
        assert_eq!(reparsed.size, original.size);
    }

    #[test]
    fn format_parse_roundtrip_at_boundaries() {
        for (offset, size) in [(0u64, 0u64), (0, u64::MAX), (u64::MAX, 0), (1, 1)] {
            let loc = DataLocation::new(offset, size);
            let reparsed = DataLocation::parse(&loc.format()).unwrap();
            assert_eq!(reparsed.offset, offset);
            assert_eq!(reparsed.size, size);
        }
    }

    #[test]
    fn location_attribute_forms_are_parsed() {
        assert!(matches!(
            LocationAttr::parse("inline:base64"),
            Ok(LocationAttr::Inline(TextEncoding::Base64))
        ));
        assert!(matches!(
            LocationAttr::parse("inline:hex"),
            Ok(LocationAttr::Inline(TextEncoding::Hex))
        ));
        assert!(matches!(
            LocationAttr::parse("embedded"),
            Ok(LocationAttr::Embedded)
        ));
        assert!(matches!(
            LocationAttr::parse("attachment:16:4"),
            Ok(LocationAttr::Attachment(DataLocation {
                offset: 16,
                size: 4,
                header_block: None
            }))
        ));
        for rejected in [
            "inline:",
            "inline:Base64",
            "embedded:1:2",
            "url(a.xisb)",
            "path(/a)",
        ] {
            assert!(LocationAttr::parse(rejected).is_err(), "{rejected}");
        }
    }

    #[test]
    fn header_blocks_address_the_decoded_bytes() {
        let location = DataLocation::inline(TextEncoding::Hex, vec![1, 2, 3]);
        assert_eq!((location.offset, location.size), (0, 3));
        assert_eq!(location.header_bytes(), Some(&[1u8, 2, 3][..]));
        assert_eq!(location.format(), "inline:hex");
        assert_eq!(DataLocation::embedded(Vec::new()).format(), "embedded");
        assert_eq!(DataLocation::new(4, 2).header_bytes(), None);
    }

    #[test]
    fn parse_error_message_mentions_context() {
        let err = DataLocation::parse("embedded:1:2").unwrap_err();
        assert!(err.to_string().contains("embedded:1:2"));

        let err = DataLocation::parse("attachment:bad:2").unwrap_err();
        assert!(err.to_string().contains("bad"));

        let err = DataLocation::parse("attachment:1:bad").unwrap_err();
        assert!(err.to_string().contains("bad"));
    }
}
