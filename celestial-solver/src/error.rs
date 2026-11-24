use celestial_catalog::error::CatalogError;
use celestial_images::core::errors::ImageError;
use celestial_wcs::error::WcsError;
use thiserror::Error;

use crate::metadata::error::MetadataError;

// The error of every fallible step of the solver. A caller tells the outcomes apart by
// variant: `NoMatch` when no match passed verification, `InvalidParameter`,
// `InvalidDimensions`, `BufferLength` and `Metadata` for bad input, `FitFailed` when the
// matched pairs do not determine a WCS, and the wrapped errors of the catalog, WCS and
// image crates for failures there.
#[derive(Debug, Error)]
pub enum SolveError {
    #[error("invalid image dimensions {width}x{height}")]
    InvalidDimensions { width: usize, height: usize },

    #[error("pixel buffer holds {actual} values but a {width}x{height} image needs {expected}")]
    BufferLength {
        width: usize,
        height: usize,
        expected: usize,
        actual: usize,
    },

    #[error("invalid parameter {name}: {reason}")]
    InvalidParameter { name: &'static str, reason: String },

    #[error("image has no finite pixels")]
    NoFinitePixels,

    #[error("{name} is not finite: the pixel values overflow f64 arithmetic")]
    NonFiniteStatistic { name: &'static str },

    #[error(transparent)]
    Metadata(#[from] MetadataError),

    #[error("no stars detected")]
    NoStars,

    // No quad match, or no solution fitted from one, passed verification. `reason` gives
    // the log-odds or the check each candidate failed.
    #[error("{reason}")]
    NoMatch { reason: String },

    // The pairs do not determine a WCS: too few of them, a singular system, or no
    // consistent inlier set.
    #[error("{reason}")]
    FitFailed { reason: String },

    // A failure while verifying the refined or initial solution, other than a rejection,
    // such as a catalog error. `solution` is "refined" or "initial". The message leaves
    // out the source, as anyhow context does, so a reporter that prints the chain prints
    // the source once.
    #[error("{solution} solution")]
    Verification {
        solution: &'static str,
        source: Box<SolveError>,
    },

    #[error(transparent)]
    Catalog(#[from] CatalogError),

    #[error(transparent)]
    Wcs(#[from] WcsError),

    #[error(transparent)]
    Image(#[from] ImageError),
}

impl SolveError {
    pub(crate) fn no_match(reason: impl Into<String>) -> Self {
        Self::NoMatch {
            reason: reason.into(),
        }
    }

    pub(crate) fn fit_failed(reason: impl Into<String>) -> Self {
        Self::FitFailed {
            reason: reason.into(),
        }
    }

    // The message and then each source's message, as `{:#}` of anyhow prints them. Some
    // errors repeat their source's message in their own and some do not, so a source
    // message already in the text is not added again.
    pub(crate) fn chain(&self) -> String {
        let mut text = self.to_string();
        let mut source = std::error::Error::source(self);
        while let Some(error) = source {
            let message = error.to_string();
            if !text.contains(&message) {
                text = format!("{text}: {message}");
            }
            source = error.source();
        }
        text
    }
}

#[cfg(test)]
impl SolveError {
    // The variant name and the message with its sources, for tests that check both.
    pub(crate) fn describe(&self) -> String {
        let variant = match self {
            Self::InvalidDimensions { .. } => "InvalidDimensions",
            Self::BufferLength { .. } => "BufferLength",
            Self::InvalidParameter { .. } => "InvalidParameter",
            Self::NoFinitePixels => "NoFinitePixels",
            Self::NonFiniteStatistic { .. } => "NonFiniteStatistic",
            Self::Metadata(_) => "Metadata",
            Self::NoStars => "NoStars",
            Self::NoMatch { .. } => "NoMatch",
            Self::FitFailed { .. } => "FitFailed",
            Self::Verification { .. } => "Verification",
            Self::Catalog(_) => "Catalog",
            Self::Wcs(_) => "Wcs",
            Self::Image(_) => "Image",
        };
        format!("{variant}: {}", self.chain())
    }
}
