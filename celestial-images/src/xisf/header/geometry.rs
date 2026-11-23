use crate::xisf::errors::{Result, XisfError};

// XISF 1.0 section 11.5.1: geometry is dim1:...:dimN:channel-count. One or
// two image dimensions are supported. A one-dimensional image reads as [n]
// when it has one channel and as [n, 1, channels] otherwise, so the channel
// count is always the third dimension.
pub(crate) fn parse_geometry(geometry_str: &str) -> Result<Vec<usize>> {
    let parts: Vec<&str> = geometry_str.split(':').collect();
    let Some((channels, axes)) = parts.split_last().filter(|(_, axes)| !axes.is_empty()) else {
        return Err(XisfError::InvalidGeometry(geometry_str.to_string()));
    };
    if axes.len() > 2 {
        return Err(XisfError::InvalidGeometry(format!(
            "{geometry_str} has {} dimensions; images with more than two are not supported",
            axes.len()
        )));
    }
    let mut dimensions = Vec::with_capacity(3);
    for (part, name) in axes.iter().zip(["Width", "Height"]) {
        dimensions.push(parse_dimension(part, name)?);
    }
    let channels = parse_dimension(channels, "Channels")?;
    if channels > 1 {
        dimensions.resize(2, 1);
        dimensions.push(channels);
    }
    Ok(dimensions)
}

fn parse_dimension(text: &str, name: &str) -> Result<usize> {
    let value = text.parse::<usize>().map_err(|_| {
        XisfError::InvalidGeometry(format!("Invalid {}: {text}", name.to_lowercase()))
    })?;
    if value == 0 {
        return Err(XisfError::InvalidGeometry(format!("{name} cannot be zero")));
    }
    Ok(value)
}

pub(crate) fn format_geometry(geometry: &[usize]) -> String {
    geometry
        .iter()
        .map(|d| d.to_string())
        .collect::<Vec<_>>()
        .join(":")
}

// [w, h] is written with an explicit channel count of 1, since without it
// the height would read as the channel count of a one-dimensional image.
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
    fn parse_two_field_geometry_as_one_dimensional() {
        assert_eq!(parse_geometry("1920:1080").unwrap(), vec![1920, 1, 1080]);
        assert_eq!(parse_geometry("1:1").unwrap(), vec![1]);
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
        assert!(err.to_string().to_lowercase().contains("channels"));

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
        let reparsed = parse_geometry(&format_geometry_with_channels(&original)).unwrap();
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
