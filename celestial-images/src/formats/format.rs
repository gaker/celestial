use crate::core::{ImageError, Result};
use std::io::{Read, Seek};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageFormat {
    Fits,
    Xisf,
    #[cfg(feature = "standard-formats")]
    Png,
    #[cfg(feature = "standard-formats")]
    Tiff,
}

impl ImageFormat {
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_lowercase().as_str() {
            "fits" | "fit" | "fts" => Some(Self::Fits),
            "xisf" => Some(Self::Xisf),
            #[cfg(feature = "standard-formats")]
            "png" => Some(Self::Png),
            #[cfg(feature = "standard-formats")]
            "tiff" | "tif" => Some(Self::Tiff),
            _ => None,
        }
    }

    pub fn from_magic_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.starts_with(b"SIMPLE  ") {
            Some(Self::Fits)
        } else if bytes.starts_with(b"<?xml") || bytes.starts_with(b"XISF") {
            Some(Self::Xisf)
        } else {
            #[cfg(feature = "standard-formats")]
            {
                if bytes.starts_with(b"\x89PNG") {
                    return Some(Self::Png);
                }
                if bytes.starts_with(b"II*\0") || bytes.starts_with(b"MM\0*") {
                    return Some(Self::Tiff);
                }
            }
            None
        }
    }

    pub fn detect<R: Read + Seek>(reader: &mut R) -> Result<Self> {
        let mut magic_bytes = [0u8; 16];
        reader.read_exact(&mut magic_bytes)?;
        reader.seek(std::io::SeekFrom::Start(0))?;

        Self::from_magic_bytes(&magic_bytes)
            .ok_or_else(|| ImageError::FormatDetectionFailed("Unknown magic bytes".to_string()))
    }

    pub fn extension(&self) -> &'static str {
        match self {
            Self::Fits => "fits",
            Self::Xisf => "xisf",
            #[cfg(feature = "standard-formats")]
            Self::Png => "png",
            #[cfg(feature = "standard-formats")]
            Self::Tiff => "tiff",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn from_extension_fits_variants() {
        for ext in ["fits", "fit", "fts"] {
            assert_eq!(ImageFormat::from_extension(ext), Some(ImageFormat::Fits));
        }
    }

    #[test]
    fn from_extension_is_case_insensitive() {
        assert_eq!(ImageFormat::from_extension("FITS"), Some(ImageFormat::Fits));
        assert_eq!(ImageFormat::from_extension("Xisf"), Some(ImageFormat::Xisf));
        assert_eq!(ImageFormat::from_extension("FiT"), Some(ImageFormat::Fits));
    }

    #[test]
    fn from_extension_xisf() {
        assert_eq!(ImageFormat::from_extension("xisf"), Some(ImageFormat::Xisf));
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn from_extension_png_and_tiff() {
        assert_eq!(ImageFormat::from_extension("png"), Some(ImageFormat::Png));
        assert_eq!(ImageFormat::from_extension("tiff"), Some(ImageFormat::Tiff));
        assert_eq!(ImageFormat::from_extension("tif"), Some(ImageFormat::Tiff));
    }

    #[test]
    fn from_extension_unknown_returns_none() {
        assert_eq!(ImageFormat::from_extension(""), None);
        assert_eq!(ImageFormat::from_extension("jpg"), None);
        assert_eq!(ImageFormat::from_extension("bmp"), None);
    }

    #[test]
    fn from_magic_bytes_fits_detects_simple_header() {
        assert_eq!(
            ImageFormat::from_magic_bytes(b"SIMPLE  = T"),
            Some(ImageFormat::Fits)
        );
    }

    #[test]
    fn from_magic_bytes_xisf_detects_xml_or_signature() {
        assert_eq!(
            ImageFormat::from_magic_bytes(b"<?xml version=\"1.0\""),
            Some(ImageFormat::Xisf)
        );
        assert_eq!(
            ImageFormat::from_magic_bytes(b"XISF0100"),
            Some(ImageFormat::Xisf)
        );
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn from_magic_bytes_png_and_tiff() {
        assert_eq!(
            ImageFormat::from_magic_bytes(b"\x89PNG\r\n\x1a\n"),
            Some(ImageFormat::Png)
        );
        assert_eq!(
            ImageFormat::from_magic_bytes(b"II*\0"),
            Some(ImageFormat::Tiff)
        );
        assert_eq!(
            ImageFormat::from_magic_bytes(b"MM\0*"),
            Some(ImageFormat::Tiff)
        );
    }

    #[test]
    fn from_magic_bytes_unknown_returns_none() {
        assert_eq!(ImageFormat::from_magic_bytes(b"random data"), None);
        assert_eq!(ImageFormat::from_magic_bytes(b""), None);
    }

    #[test]
    fn detect_reads_and_rewinds_cursor() {
        let mut data = b"SIMPLE  = T".to_vec();
        data.resize(100, b' ');
        let mut cursor = Cursor::new(data);
        let format = ImageFormat::detect(&mut cursor).unwrap();
        assert_eq!(format, ImageFormat::Fits);
        assert_eq!(cursor.position(), 0);
    }

    #[test]
    fn detect_errors_on_unrecognized_magic() {
        let mut cursor = Cursor::new(vec![0u8; 100]);
        assert!(ImageFormat::detect(&mut cursor).is_err());
    }

    #[test]
    fn detect_errors_on_short_input() {
        let mut cursor = Cursor::new(vec![0u8; 4]);
        assert!(ImageFormat::detect(&mut cursor).is_err());
    }

    #[test]
    fn extension_matches_canonical_name() {
        assert_eq!(ImageFormat::Fits.extension(), "fits");
        assert_eq!(ImageFormat::Xisf.extension(), "xisf");
    }

    #[cfg(feature = "standard-formats")]
    #[test]
    fn extension_png_and_tiff() {
        assert_eq!(ImageFormat::Png.extension(), "png");
        assert_eq!(ImageFormat::Tiff.extension(), "tiff");
    }

    #[test]
    fn equality_and_copy_semantics() {
        let a = ImageFormat::Fits;
        let b = a;
        assert_eq!(a, b);
        assert_ne!(ImageFormat::Fits, ImageFormat::Xisf);
    }
}
