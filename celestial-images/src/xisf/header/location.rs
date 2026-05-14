use crate::xisf::{Result, XisfError};

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
    fn parse_error_message_mentions_context() {
        let err = DataLocation::parse("embedded:1:2").unwrap_err();
        assert!(err.to_string().contains("embedded:1:2"));

        let err = DataLocation::parse("attachment:bad:2").unwrap_err();
        assert!(err.to_string().contains("bad"));

        let err = DataLocation::parse("attachment:1:bad").unwrap_err();
        assert!(err.to_string().contains("bad"));
    }
}
