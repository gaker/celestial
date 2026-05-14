use crate::xisf::{Result, XisfError};

pub(crate) fn parse_geometry(geometry_str: &str) -> Result<Vec<usize>> {
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

pub(crate) fn format_geometry(geometry: &[usize]) -> String {
    geometry
        .iter()
        .map(|d| d.to_string())
        .collect::<Vec<_>>()
        .join(":")
}

pub(crate) fn format_geometry_with_channels(geometry: &[usize]) -> String {
    if geometry.len() == 2 {
        format!("{}:{}:1", geometry[0], geometry[1])
    } else {
        format_geometry(geometry)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_two_dim_geometry() {
        assert_eq!(parse_geometry("1920:1080").unwrap(), vec![1920, 1080]);
        assert_eq!(parse_geometry("1:1").unwrap(), vec![1, 1]);
    }

    #[test]
    fn parse_three_dim_with_multiple_channels() {
        assert_eq!(parse_geometry("1920:1080:3").unwrap(), vec![1920, 1080, 3]);
        assert_eq!(parse_geometry("100:100:4").unwrap(), vec![100, 100, 4]);
    }

    #[test]
    fn parse_three_dim_with_one_channel_drops_channel() {
        assert_eq!(parse_geometry("1920:1080:1").unwrap(), vec![1920, 1080]);
    }

    #[test]
    fn parse_rejects_wrong_segment_count() {
        let err = parse_geometry("1920").unwrap_err();
        assert!(err.to_string().contains("1920"));

        let err = parse_geometry("1920:1080:3:4").unwrap_err();
        assert!(err.to_string().contains("1920:1080:3:4"));

        let err = parse_geometry("").unwrap_err();
        assert!(matches!(err, XisfError::InvalidGeometry(_)));
    }

    #[test]
    fn parse_rejects_zero_dimensions() {
        let err = parse_geometry("0:1080").unwrap_err();
        assert!(err.to_string().to_lowercase().contains("width"));

        let err = parse_geometry("1920:0").unwrap_err();
        assert!(err.to_string().to_lowercase().contains("height"));

        let err = parse_geometry("1920:1080:0").unwrap_err();
        assert!(err.to_string().to_lowercase().contains("channels"));
    }

    #[test]
    fn parse_rejects_non_numeric() {
        assert!(parse_geometry("abc:1080").is_err());
        assert!(parse_geometry("1920:def").is_err());
        assert!(parse_geometry("1920:1080:xyz").is_err());
    }

    #[test]
    fn parse_rejects_negative_numbers() {
        assert!(parse_geometry("-1920:1080").is_err());
        assert!(parse_geometry("1920:-1080").is_err());
        assert!(parse_geometry("1920:1080:-1").is_err());
    }

    #[test]
    fn parse_error_messages_mention_offending_segment() {
        let err = parse_geometry("abc:1080").unwrap_err();
        assert!(err.to_string().contains("abc"));

        let err = parse_geometry("1920:def").unwrap_err();
        assert!(err.to_string().contains("def"));

        let err = parse_geometry("1920:1080:xyz").unwrap_err();
        assert!(err.to_string().contains("xyz"));
    }

    #[test]
    fn parse_accepts_large_values() {
        let g = parse_geometry("65535:65535:4").unwrap();
        assert_eq!(g, vec![65535, 65535, 4]);
    }

    #[test]
    fn format_geometry_joins_with_colons() {
        assert_eq!(format_geometry(&[1920, 1080]), "1920:1080");
        assert_eq!(format_geometry(&[1920, 1080, 3]), "1920:1080:3");
    }

    #[test]
    fn format_geometry_handles_empty_and_single() {
        assert_eq!(format_geometry(&[]), "");
        assert_eq!(format_geometry(&[42]), "42");
    }

    #[test]
    fn format_with_channels_appends_one_for_two_dim() {
        assert_eq!(format_geometry_with_channels(&[1920, 1080]), "1920:1080:1");
    }

    #[test]
    fn format_with_channels_passes_through_multi_channel() {
        assert_eq!(
            format_geometry_with_channels(&[1920, 1080, 3]),
            "1920:1080:3"
        );
    }

    #[test]
    fn format_with_channels_passes_through_non_two_dim() {
        assert_eq!(format_geometry_with_channels(&[]), "");
        assert_eq!(format_geometry_with_channels(&[100]), "100");
    }

    #[test]
    fn parse_format_roundtrip_two_dim() {
        let original = vec![1920, 1080];
        let reparsed = parse_geometry(&format_geometry(&original)).unwrap();
        assert_eq!(reparsed, original);
    }

    #[test]
    fn parse_format_roundtrip_three_dim() {
        let original = vec![1920, 1080, 3];
        let reparsed = parse_geometry(&format_geometry(&original)).unwrap();
        assert_eq!(reparsed, original);
    }

    #[test]
    fn parse_then_format_with_channels_normalizes_single_channel() {
        let parsed = parse_geometry("1920:1080:1").unwrap();
        assert_eq!(parsed, vec![1920, 1080]);
        assert_eq!(format_geometry_with_channels(&parsed), "1920:1080:1");
    }
}
