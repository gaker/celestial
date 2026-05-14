use thiserror::Error;

/// Error returned when metadata can't be extracted from an image.
///
/// Returned from [`crate::metadata::metadata_from_image`], [`crate::metadata::hint_from_image`],
/// [`crate::metadata::jd_from_image`], and [`crate::metadata::plate_scale_from_image`].
#[derive(Error, Debug)]
pub enum MetadataError {
    /// A required header is missing and no override was supplied.
    #[error("image is missing required header: {0}")]
    MissingHeader(&'static str),

    /// A header is present but its value isn't the expected type (e.g. string where a
    /// number is required).
    #[error("header {keyword} has no {kind} value")]
    InvalidHeaderType {
        /// Header keyword that was inspected.
        keyword: &'static str,
        /// Expected value kind ("numeric", "string", etc).
        kind: &'static str,
    },

    /// None of the supported position-hint header pairs were found: `CRVAL1/2`,
    /// `RA/DEC`, or `OBJCTRA/OBJCTDEC`.
    #[error("no position hint: need CRVAL1/2, RA/DEC, or OBJCTRA/OBJCTDEC")]
    NoPositionHint,

    /// `DATE-OBS` is a date-only string and `TIME-OBS` is missing.
    #[error("DATE-OBS has no time component and TIME-OBS is missing")]
    MissingTimeOfDay,

    /// `DATE-OBS` failed to parse as an ISO 8601 timestamp or a plain date.
    #[error("failed to parse DATE-OBS {value:?}: {source}")]
    DateParse {
        /// Offending header value.
        value: String,
        /// Underlying chrono parse error.
        #[source]
        source: chrono::ParseError,
    },

    /// `TIME-OBS` failed to parse as a time-of-day string.
    #[error("failed to parse TIME-OBS {value:?}: {source}")]
    TimeParse {
        /// Offending header value.
        value: String,
        /// Underlying chrono parse error.
        #[source]
        source: chrono::ParseError,
    },

    /// A sexagesimal angle (hours or degrees-minutes-seconds) failed to parse.
    #[error("failed to parse sexagesimal angle {value:?}: {source}")]
    AngleParse {
        /// Offending header value.
        value: String,
        /// Underlying parse error from [`celestial_core::AstroError`].
        #[source]
        source: celestial_core::AstroError,
    },

    /// The computed position is invalid (out-of-range declination, non-finite value).
    #[error("invalid position: {0}")]
    InvalidPosition(#[from] celestial_coords::CoordError),

    /// A numeric header that must be positive was zero or negative.
    #[error("{field} must be positive, got {value}")]
    NonPositiveValue {
        /// Name of the offending field ("focal length", "pixel size").
        field: &'static str,
        /// Rejected value.
        value: f64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_messages_include_context() {
        let cases: Vec<(MetadataError, &str)> = vec![
            (MetadataError::MissingHeader("FOCALLEN"), "FOCALLEN"),
            (
                MetadataError::InvalidHeaderType {
                    keyword: "CRVAL1",
                    kind: "numeric",
                },
                "CRVAL1",
            ),
            (MetadataError::NoPositionHint, "CRVAL"),
            (MetadataError::MissingTimeOfDay, "TIME-OBS"),
            (
                MetadataError::NonPositiveValue {
                    field: "focal length",
                    value: -1.0,
                },
                "focal length",
            ),
        ];
        for (err, needle) in cases {
            let msg = err.to_string();
            assert!(msg.contains(needle), "{msg:?} missing {needle:?}");
        }
    }

    #[test]
    fn date_parse_display_mentions_value() {
        let chrono_err = chrono::NaiveDate::parse_from_str("not-a-date", "%Y-%m-%d").unwrap_err();
        let err = MetadataError::DateParse {
            value: "not-a-date".to_string(),
            source: chrono_err,
        };
        assert!(err.to_string().contains("not-a-date"));
    }

    #[test]
    fn time_parse_display_mentions_value() {
        let chrono_err = chrono::NaiveTime::parse_from_str("nope", "%H:%M:%S").unwrap_err();
        let err = MetadataError::TimeParse {
            value: "nope".to_string(),
            source: chrono_err,
        };
        assert!(err.to_string().contains("nope"));
    }

    #[test]
    fn angle_parse_display_mentions_value() {
        let astro_err = celestial_core::angle::parse_hms("nope").unwrap_err();
        let err = MetadataError::AngleParse {
            value: "nope".to_string(),
            source: astro_err,
        };
        assert!(err.to_string().contains("nope"));
    }

    #[test]
    fn coord_error_converts_via_from_impl() {
        let coord_err =
            celestial_coords::ICRSPosition::from_degrees(f64::NAN, 0.0).unwrap_err();
        let err: MetadataError = coord_err.into();
        assert!(matches!(err, MetadataError::InvalidPosition(_)));
    }
}
