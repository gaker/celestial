use super::fits_scaling::Form::{Offset, Physical, Plain};
use super::fits_scaling::{keywords_for_save, loaded_keywords, Scaling};
use super::Image;
use crate::core::errors::{ImageError, Result};
use crate::core::types::BitPix;
use crate::fits::compression::decode::{decode_image, DecodedImage};
use crate::fits::compression::spec::has_pixels;
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::{decode, DataArray};
use crate::fits::errors::FitsError;
use crate::fits::hdu::{Hdu, HduTrait};
use crate::fits::header::keywords::Keyword;
use crate::fits::header::parser::Header;
use crate::fits::io::reader::{FitsFile, FitsReader};
use crate::fits::io::writer::{compressed_image_bytes, primary_image_bytes};
use crate::formats::atomic_write::write_atomically;
use crate::formats::pixel_data::PixelData;
use std::path::Path;

impl Image {
    pub(super) fn open_fits(path: &Path) -> Result<Self> {
        let mut fits = FitsFile::open(path).map_err(ImageError::Fits)?;
        let DecodedImage {
            header,
            dimensions,
            bitpix,
            data,
        } = match first_image_hdu(&mut fits)? {
            Found::Plain(index) => plain_image(&mut fits, index)?,
            Found::Compressed(index) => {
                let header = fits.get_header(index)?;
                decode_image(&header, &fits.read_data_unit(index)?)?
            }
        };
        let scaling = Scaling::from_header(&header, bitpix)?;
        let pixels = read_pixels(&data, bitpix, &scaling)?;
        let keywords = loaded_keywords(header.keywords(), &pixels, scaling.blank);

        // FITS stores multi-channel images planar (all R, then all G, then
        // all B).
        Ok(Self::from_parts(
            pixels,
            dimensions,
            keywords,
            Vec::new(),
            false,
        ))
    }

    pub(super) fn save_fits(&self, path: &Path) -> Result<()> {
        self.write_fits(path, Layout::Primary)
    }

    // Writes a FITS file whose image is tile-compressed (FITS 4.0 §10) in the
    // first extension. Tiles are rows by default, or `tile_size` pixels in
    // the first two axes. Unsigned 16- and 32-bit and signed 8-bit pixels
    // carry the same BZERO offset as an uncompressed save.
    pub fn save_fits_compressed<P: AsRef<Path>>(
        &self,
        path: P,
        algorithm: CompressionAlgorithm,
        tile_size: Option<(usize, usize)>,
    ) -> Result<()> {
        // U64 tiles are stored as ZBITPIX 64 with BZERO 2^63, which other
        // FITS readers do not decode to the same values.
        if matches!(self.pixels, PixelData::U64(_)) {
            return Err(ImageError::Fits(FitsError::UnsupportedCompression(
                "other FITS readers cannot decode tiles of 64-bit unsigned pixels; \
                 save the image uncompressed"
                    .to_string(),
            )));
        }
        self.write_fits(path.as_ref(), Layout::Tiled(algorithm, tile_size))
    }

    fn write_fits(&self, path: &Path, layout: Layout) -> Result<()> {
        // FITS expects multi-channel data in planar layout (NAXIS3
        // chunks each holding one full channel). In-memory we keep
        // it interleaved, so flip for the write.
        let source = if self.has_rgb_planes() {
            let mut planar = self.clone();
            planar.interleaved_to_planar();
            std::borrow::Cow::Owned(planar)
        } else {
            std::borrow::Cow::Borrowed(self)
        };
        let keywords = keywords_for_save(&source.keywords, &source.pixels);
        let target = Target {
            dims: &source.dimensions,
            keywords: &keywords,
            layout,
        };
        let bytes = encode_pixels(target, &source.pixels).map_err(ImageError::Fits)?;
        write_atomically(path, &bytes).map_err(|e| ImageError::Fits(FitsError::Io(e)))
    }
}

#[derive(Clone, Copy)]
enum Layout {
    Primary,
    Tiled(CompressionAlgorithm, Option<(usize, usize)>),
}

struct Target<'a> {
    dims: &'a [usize],
    keywords: &'a [Keyword],
    layout: Layout,
}

impl Target<'_> {
    fn encode<T: DataArray>(self, data: &[T]) -> crate::fits::errors::Result<Vec<u8>> {
        match self.layout {
            Layout::Primary => primary_image_bytes(data, self.dims, self.keywords, false),
            Layout::Tiled(algorithm, tile_size) => {
                compressed_image_bytes(data, self.dims, tile_size, algorithm, self.keywords)
            }
        }
    }
}

enum Found {
    Plain(usize),
    Compressed(usize),
}

fn plain_image(
    fits: &mut FitsFile<FitsReader>,
    index: usize,
) -> crate::fits::errors::Result<DecodedImage> {
    let (dimensions, bitpix) = fits.get_image_info(index)?;
    Ok(DecodedImage {
        header: fits.get_header(index)?,
        dimensions,
        bitpix,
        data: fits.read_image_array(index)?,
    })
}

// The first HDU that holds image data, the HDU astropy's getdata returns. A
// BINTABLE with ZIMAGE = T counts as the image its tiles hold. Other tables,
// HDUs without pixels and extensions of unknown type are skipped, as FITS 4.0
// §3.4.1 asks of readers. With no image before the HDU the scan stopped at,
// the scan error is returned. Random groups (§6) are not an image.
fn first_image_hdu(fits: &mut FitsFile<FitsReader>) -> crate::fits::errors::Result<Found> {
    for index in 0..fits.num_hdus() {
        if index > 0 && !may_hold_image(&fits.get_header(index)?) {
            continue;
        }
        let dimensions = match fits.read_hdu(index)? {
            Hdu::Primary(hdu) => hdu.data_dimensions()?,
            Hdu::Image(hdu) => hdu.data_dimensions()?,
            Hdu::BinaryTable(hdu) if hdu.is_compressed_image() => {
                if has_pixels(hdu.header()) {
                    return Ok(Found::Compressed(index));
                }
                continue;
            }
            _ => continue,
        };
        // GCOUNT = 0 leaves an empty data unit even when every axis is set.
        let has_data = fits.hdu_info(index).is_some_and(|info| info.data_size > 0);
        if has_data && !dimensions.is_empty() && dimensions.iter().all(|&n| n > 0) {
            return Ok(Found::Plain(index));
        }
    }
    Err(no_image(fits)?)
}

// Why no HDU was found: the scan error, random groups in the primary HDU,
// or an image of no pixels in every HDU.
fn no_image(fits: &mut FitsFile<FitsReader>) -> crate::fits::errors::Result<FitsError> {
    if let Some(error) = fits.take_scan_error() {
        return Ok(error);
    }
    let groups = fits.num_hdus() > 0 && matches!(fits.read_hdu(0)?, Hdu::RandomGroups(_));
    let message = match groups {
        true => "the primary HDU holds random groups, not an image",
        false => "no HDU holds image data",
    };
    Ok(FitsError::InvalidFormat(message.to_string()))
}

// Tables are checked when read, so only compressed images are read as
// tables. A malformed table does not hide a later image.
fn may_hold_image(header: &Header) -> bool {
    match header.get("XTENSION").as_str() {
        Some("IMAGE") => true,
        Some("BINTABLE") => header.get("ZIMAGE").as_bool() == Some(true),
        _ => false,
    }
}

// Flipping the top bit adds the Table 11 offset modulo 2^n, which maps the
// signed stored range onto the unsigned in-memory range and back.
fn encode_pixels(target: Target, pixels: &PixelData) -> crate::fits::errors::Result<Vec<u8>> {
    match pixels {
        PixelData::U8(d) => target.encode(d),
        PixelData::I8(d) => {
            let stored: Vec<u8> = d.iter().map(|&v| (v as u8) ^ 0x80).collect();
            target.encode(&stored)
        }
        PixelData::U16(d) => {
            let stored: Vec<i16> = d.iter().map(|&v| (v ^ 0x8000) as i16).collect();
            target.encode(&stored)
        }
        PixelData::I16(d) => target.encode(d),
        PixelData::U32(d) => {
            let stored: Vec<i32> = d.iter().map(|&v| (v ^ (1 << 31)) as i32).collect();
            target.encode(&stored)
        }
        PixelData::I32(d) => target.encode(d),
        PixelData::U64(d) => {
            let stored: Vec<i64> = d.iter().map(|&v| (v ^ (1 << 63)) as i64).collect();
            target.encode(&stored)
        }
        PixelData::I64(d) => target.encode(d),
        PixelData::F32(d) => target.encode(d),
        PixelData::F64(d) => target.encode(d),
    }
}

// Integer data keeps its width when the scaling is the identity or the
// Table 11 offset for that width. Flipping the top bit adds that offset
// modulo 2^n. Any other scaling returns F64, which holds BZERO + BSCALE *
// stored as the f64 arithmetic gives it.
fn read_pixels(
    bytes: &[u8],
    bitpix: BitPix,
    scaling: &Scaling,
) -> crate::fits::errors::Result<PixelData> {
    let integer = |stored: i64| scaling.integer(stored);
    let real = |stored: f64| scaling.real(stored);
    Ok(match (bitpix, scaling.form(bitpix)) {
        (BitPix::U8, Plain) => PixelData::U8(decode(bytes, u8::from_be_bytes)?),
        (BitPix::U8, Offset) => PixelData::I8(decode(bytes, |b| i8::from_be_bytes(b) ^ i8::MIN)?),
        (BitPix::U8, Physical) => {
            PixelData::F64(decode(bytes, |b| integer(u8::from_be_bytes(b).into()))?)
        }
        (BitPix::I16, Plain) => PixelData::I16(decode(bytes, i16::from_be_bytes)?),
        (BitPix::I16, Offset) => PixelData::U16(decode(bytes, |b| u16::from_be_bytes(b) ^ 0x8000)?),
        (BitPix::I16, Physical) => {
            PixelData::F64(decode(bytes, |b| integer(i16::from_be_bytes(b).into()))?)
        }
        (BitPix::I32, Plain) => PixelData::I32(decode(bytes, i32::from_be_bytes)?),
        (BitPix::I32, Offset) => {
            PixelData::U32(decode(bytes, |b| u32::from_be_bytes(b) ^ (1 << 31))?)
        }
        (BitPix::I32, Physical) => {
            PixelData::F64(decode(bytes, |b| integer(i32::from_be_bytes(b).into()))?)
        }
        (BitPix::I64, Plain) => PixelData::I64(decode(bytes, i64::from_be_bytes)?),
        (BitPix::I64, Offset) => {
            PixelData::U64(decode(bytes, |b| u64::from_be_bytes(b) ^ (1 << 63))?)
        }
        (BitPix::I64, Physical) => {
            PixelData::F64(decode(bytes, |b| integer(i64::from_be_bytes(b)))?)
        }
        (BitPix::F32, Plain) => PixelData::F32(decode(bytes, f32::from_be_bytes)?),
        (BitPix::F32, _) => PixelData::F64(decode(bytes, |b| real(f32::from_be_bytes(b).into()))?),
        (BitPix::F64, Plain) => PixelData::F64(decode(bytes, f64::from_be_bytes)?),
        (BitPix::F64, _) => PixelData::F64(decode(bytes, |b| real(f64::from_be_bytes(b)))?),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::{Keyword, KeywordValue};
    use crate::fits::io::test_support::FitsWriter;

    fn tmp_fits() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".fits").tempfile().unwrap()
    }

    // One primary HDU with NAXIS = 1. `extra` cards follow NAXIS1 and
    // `data` is the big-endian data unit before padding.
    fn raw_fits(
        bitpix: i64,
        len: usize,
        extra: &[(&str, &str)],
        data: &[u8],
    ) -> tempfile::NamedTempFile {
        let bitpix = bitpix.to_string();
        let len = len.to_string();
        let mut cards = vec![
            ("SIMPLE", "T"),
            ("BITPIX", bitpix.as_str()),
            ("NAXIS", "1"),
            ("NAXIS1", len.as_str()),
        ];
        cards.extend_from_slice(extra);
        let mut bytes: Vec<u8> = cards
            .iter()
            .flat_map(|(name, value)| {
                format!("{:<80}", format!("{name:<8}= {value:>20}")).into_bytes()
            })
            .collect();
        bytes.extend_from_slice(format!("{:<80}", "END").as_bytes());
        bytes.resize(bytes.len().div_ceil(2880) * 2880, b' ');
        bytes.extend_from_slice(data);
        bytes.resize(bytes.len().div_ceil(2880) * 2880, 0);
        let tmp = tmp_fits();
        std::fs::write(tmp.path(), bytes).unwrap();
        tmp
    }

    fn stored_header(path: &Path) -> crate::fits::header::parser::Header {
        crate::fits::io::reader::FitsFile::open(path)
            .unwrap()
            .get_header(0)
            .unwrap()
    }

    #[test]
    fn open_reads_bitpix8_with_bzero_minus128_as_i8() {
        let tmp = raw_fits(8, 6, &[("BZERO", "-128")], &[0, 127, 128, 129, 228, 255]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i8(), Some(&vec![-128, -1, 0, 1, 100, 127]));
        assert!(img.get_keyword("BZERO").is_none());
    }

    #[test]
    fn open_reads_bitpix32_with_bzero_2_pow_31_as_u32() {
        let stored: Vec<u8> = [i32::MIN, -2147483647, -1, 0, 852516352, i32::MAX]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(
            32,
            6,
            &[("BZERO", "2147483648"), ("BSCALE", "1.0")],
            &stored,
        );
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(
            img.pixels.as_u32(),
            Some(&vec![0, 1, 2147483647, 2147483648, 3000000000, 4294967295])
        );
        assert!(img.get_keyword("BZERO").is_none());
        assert!(img.get_keyword("BSCALE").is_none());
    }

    #[test]
    fn open_applies_bscale_and_bzero_to_integer_data_as_f64() {
        let stored: Vec<u8> = [0i16, 1, -2, 10, i16::MIN, i16::MAX]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(16, 6, &[("BSCALE", "0.5"), ("BZERO", "100.0")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(
            img.pixels.as_f64(),
            Some(&vec![100.0, 100.5, 99.0, 105.0, -16284.0, 16483.5])
        );
        assert!(img.get_keyword("BZERO").is_none());
        assert!(img.get_keyword("BSCALE").is_none());
    }

    #[test]
    fn open_turns_blank_pixels_of_scaled_data_into_nan() {
        let stored: Vec<u8> = [0i16, i16::MIN, 7]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(
            16,
            3,
            &[("BSCALE", "0.5"), ("BZERO", "100.0"), ("BLANK", "-32768")],
            &stored,
        );
        let img = Image::open(tmp.path()).unwrap();
        let bits: Vec<u64> = img
            .pixels
            .as_f64()
            .unwrap()
            .iter()
            .map(|v| v.to_bits())
            .collect();
        assert_eq!(
            bits,
            vec![100f64.to_bits(), f64::NAN.to_bits(), 103.5f64.to_bits()]
        );
        assert!(img.get_keyword("BLANK").is_none());
    }

    #[test]
    fn unscaled_integer_blank_stays_raw_until_converted_to_float() {
        let stored: Vec<u8> = [5i16, i16::MIN, 7, 9]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(16, 4, &[("BLANK", "-32768")], &stored);
        let mut img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![5, i16::MIN, 7, 9]));
        assert_eq!(img.get("BLANK").as_i64(), Some(-32768));

        img.normalize();
        let bits: Vec<u64> = img
            .pixels
            .as_f64()
            .unwrap()
            .iter()
            .map(|v| v.to_bits())
            .collect();
        assert_eq!(bits, [0.0, f64::NAN, 0.5, 1.0].map(f64::to_bits));
        assert!(img.get_keyword("BLANK").is_none());
    }

    #[test]
    fn open_reads_a_real_blank_that_holds_an_integer() {
        let stored: Vec<u8> = [5i16, i16::MIN]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(16, 2, &[("BLANK", "-32768.0")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![5, i16::MIN]));
        assert_eq!(
            img.get_keyword("BLANK").unwrap().value,
            Some(KeywordValue::Integer(-32768))
        );
    }

    #[test]
    fn open_drops_blank_cards_it_cannot_use() {
        let stored: Vec<u8> = [5i16, 7].iter().flat_map(|v| v.to_be_bytes()).collect();
        for blank in ["'none'", "1.5", "70000"] {
            let tmp = raw_fits(16, 2, &[("BLANK", blank)], &stored);
            let img = Image::open(tmp.path()).unwrap();
            assert_eq!(img.pixels.as_i16(), Some(&vec![5, 7]), "{blank}");
            assert!(img.get_keyword("BLANK").is_none(), "{blank}");
        }
        let tmp = raw_fits(16, 2, &[("BSCALE", "2.0"), ("BLANK", "70000")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_f64(), Some(&vec![10.0, 14.0]));
    }

    #[test]
    fn open_stores_blank_of_offset_data_in_physical_units() {
        let stored: Vec<u8> = [i16::MIN, 0, 5]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(16, 3, &[("BZERO", "32768"), ("BLANK", "-32768")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_u16(), Some(&vec![0, 32768, 32773]));
        assert_eq!(img.get("BLANK").as_i64(), Some(0));
    }

    #[test]
    fn save_writes_blank_of_offset_data_in_stored_units() {
        let cases = [
            (PixelData::U16(vec![0, 9]), 0, -32768),
            (PixelData::I8(vec![-128, 9]), -128, 0),
            (PixelData::U32(vec![5, 9]), 5, -2147483643),
        ];
        for (pixels, physical, stored) in cases {
            let tmp = tmp_fits();
            let mut img = Image::new(pixels, vec![2usize]);
            img.set_keyword(Keyword::integer("BLANK", physical));
            img.save(tmp.path()).unwrap();
            assert_eq!(
                stored_header(tmp.path()).get("BLANK").as_i64(),
                Some(stored)
            );
            let restored = Image::open(tmp.path()).unwrap();
            assert_eq!(restored.get("BLANK").as_i64(), Some(physical));
        }
    }

    #[test]
    fn float_data_is_scaled_once_across_open_and_save() {
        let stored: Vec<u8> = [0f32, 1.0, 2.5, -1.0, 3.0, 4.0]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(-32, 6, &[("BSCALE", "2.0"), ("BZERO", "1.0")], &stored);
        let expected = vec![1.0f32, 3.0, 6.0, -1.0, 7.0, 9.0];

        let first = Image::open(tmp.path()).unwrap();
        assert_eq!(first.pixels.to_f32(), expected);
        let resaved = tmp_fits();
        first.save(resaved.path()).unwrap();
        let second = Image::open(resaved.path()).unwrap();
        assert_eq!(second.pixels.to_f32(), expected);
        assert!(matches!(second.pixels, PixelData::F64(_)));
        assert!(!stored_header(resaved.path()).get("BZERO").exists());
        assert!(!stored_header(resaved.path()).get("BSCALE").exists());
    }

    #[test]
    fn open_ignores_blank_on_float_data() {
        let stored: Vec<u8> = [1.5f32, -2.0]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let tmp = raw_fits(-32, 2, &[("BLANK", "'none'")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_f32(), Some(&vec![1.5, -2.0]));
        assert!(img.get_keyword("BLANK").is_none());
    }

    const I64_VALUES: [i64; 6] = [-(1 << 62), -1, 0, 1, (1 << 53) + 1, i64::MAX];
    const U64_VALUES: [u64; 5] = [0, 1, (1 << 53) + 1, 1 << 63, u64::MAX];

    fn u64_stored() -> Vec<i64> {
        U64_VALUES.iter().map(|&v| (v ^ (1 << 63)) as i64).collect()
    }

    #[test]
    fn open_reads_bitpix64_as_i64() {
        let stored: Vec<u8> = I64_VALUES.iter().flat_map(|v| v.to_be_bytes()).collect();
        let tmp = raw_fits(64, 6, &[], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i64(), Some(&I64_VALUES.to_vec()));
    }

    #[test]
    fn open_reads_bitpix64_with_bzero_2_pow_63_as_u64() {
        let stored: Vec<u8> = u64_stored().iter().flat_map(|v| v.to_be_bytes()).collect();
        let tmp = raw_fits(64, 5, &[("BZERO", "9223372036854775808")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_u64(), Some(&U64_VALUES.to_vec()));
        assert!(img.get_keyword("BZERO").is_none());
    }

    #[test]
    fn open_reads_bitpix64_as_u64_for_every_spelling_of_2_pow_63() {
        let stored: Vec<u8> = u64_stored().iter().flat_map(|v| v.to_be_bytes()).collect();
        for bzero in [
            "9223372036854775808",
            "9.223372036854775808E18",
            "9.223372036854776E+18",
            "9223372036854776000.0",
            "9223372036854775808.0",
        ] {
            let tmp = raw_fits(64, 5, &[("BZERO", bzero)], &stored);
            let img = Image::open(tmp.path()).unwrap();
            assert_eq!(img.pixels.as_u64(), Some(&U64_VALUES.to_vec()), "{bzero}");
        }
    }

    #[test]
    fn open_reads_bitpix64_with_bzero_i64_max_as_f64() {
        let stored: Vec<u8> = [0i64, 1, -1].iter().flat_map(|v| v.to_be_bytes()).collect();
        let tmp = raw_fits(64, 3, &[("BZERO", "9223372036854775807")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        let bzero = i64::MAX as f64;
        assert_eq!(
            img.pixels.as_f64(),
            Some(&vec![bzero, bzero + 1.0, bzero - 1.0])
        );
    }

    #[test]
    fn open_applies_bscale_to_bitpix64_as_f64() {
        let stored: Vec<u8> = [3i64, -4].iter().flat_map(|v| v.to_be_bytes()).collect();
        let tmp = raw_fits(64, 2, &[("BSCALE", "2.0")], &stored);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_f64(), Some(&vec![6.0, -8.0]));
    }

    #[test]
    fn roundtrip_i64() {
        let tmp = tmp_fits();
        Image::new(PixelData::I64(I64_VALUES.to_vec()), vec![6usize])
            .save(tmp.path())
            .unwrap();
        let header = stored_header(tmp.path());
        assert_eq!(header.get("BITPIX").as_i64(), Some(64));
        assert!(!header.get("BZERO").exists());
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i64(), Some(&I64_VALUES.to_vec()));
    }

    #[test]
    fn roundtrip_u64_uses_bzero_2_pow_63() {
        let tmp = tmp_fits();
        Image::new(PixelData::U64(U64_VALUES.to_vec()), vec![5usize])
            .save(tmp.path())
            .unwrap();
        let header = stored_header(tmp.path());
        assert_eq!(header.get("BITPIX").as_i64(), Some(64));
        assert_eq!(header.get("BZERO").as_f64(), Some(9223372036854775808.0));
        assert_eq!(header.get("BSCALE").as_f64(), Some(1.0));
        let (_, stored): (_, Vec<i64>) = crate::fits::io::reader::FitsFile::open(tmp.path())
            .unwrap()
            .primary_hdu_with_data()
            .unwrap();
        assert_eq!(stored, u64_stored());
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u64(), Some(&U64_VALUES.to_vec()));
    }

    #[test]
    fn open_rejects_non_finite_scaling() {
        let tmp = raw_fits(16, 1, &[("BZERO", "NaN")], &[0, 1]);
        assert!(matches!(
            Image::open(tmp.path()),
            Err(ImageError::Fits(crate::fits::errors::FitsError::InvalidKeywordValue { keyword, .. }))
                if keyword == "BZERO"
        ));
    }

    #[test]
    fn roundtrip_i8_uses_bzero_minus128() {
        let tmp = tmp_fits();
        let original = vec![-128i8, -1, 0, 1, 100, 127];
        Image::new(PixelData::I8(original.clone()), vec![6usize])
            .save(tmp.path())
            .unwrap();

        let header = stored_header(tmp.path());
        assert_eq!(header.get("BITPIX").as_i64(), Some(8));
        assert_eq!(header.get("BZERO").as_f64(), Some(-128.0));
        assert_eq!(header.get("BSCALE").as_f64(), Some(1.0));
        let (_, stored): (_, Vec<u8>) = crate::fits::io::reader::FitsFile::open(tmp.path())
            .unwrap()
            .primary_hdu_with_data()
            .unwrap();
        assert_eq!(stored, vec![0, 127, 128, 129, 228, 255]);
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i8(), Some(&original));
    }

    #[test]
    fn roundtrip_u32_uses_bzero_2_pow_31() {
        let tmp = tmp_fits();
        let original = vec![0u32, 1, 2147483647, 2147483648, 3000000000, 4294967295];
        Image::new(PixelData::U32(original.clone()), vec![6usize])
            .save(tmp.path())
            .unwrap();

        let header = stored_header(tmp.path());
        assert_eq!(header.get("BITPIX").as_i64(), Some(32));
        assert_eq!(header.get("BZERO").as_f64(), Some(2147483648.0));
        assert_eq!(header.get("BSCALE").as_f64(), Some(1.0));
        let (_, stored): (_, Vec<i32>) = crate::fits::io::reader::FitsFile::open(tmp.path())
            .unwrap()
            .primary_hdu_with_data()
            .unwrap();
        assert_eq!(
            stored,
            vec![i32::MIN, -2147483647, -1, 0, 852516352, i32::MAX]
        );
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u32(), Some(&original));
    }

    #[test]
    fn save_replaces_scaling_keywords_set_on_an_integer_image() {
        let tmp = tmp_fits();
        let mut img = Image::new(PixelData::I16(vec![-3, 4]), vec![2usize]);
        img.set_keyword(Keyword::real("BZERO", 32768.0));
        img.set_keyword(Keyword::real("BSCALE", 1.0));
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i16(), Some(&vec![-3, 4]));
    }

    #[test]
    fn roundtrip_u8_preserves_data_and_dimensions() {
        let tmp = tmp_fits();
        let original: Vec<u8> = (0u8..100).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![10usize, 10]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_i16() {
        let tmp = tmp_fits();
        let original: Vec<i16> = (-50i16..50).collect();
        let img = Image::new(PixelData::I16(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i16().unwrap(), &original);
    }

    #[test]
    fn roundtrip_i32() {
        let tmp = tmp_fits();
        let original: Vec<i32> = (0..100).map(|i| i * 1000).collect();
        let img = Image::new(PixelData::I32(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i32().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f32() {
        let tmp = tmp_fits();
        let original: Vec<f32> = (0..100).map(|i| i as f32 * 0.1).collect();
        let img = Image::new(PixelData::F32(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f64() {
        let tmp = tmp_fits();
        let original: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
        let img = Image::new(PixelData::F64(original.clone()), vec![10usize, 10]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f64().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u16_uses_bzero_bscale() {
        let tmp = tmp_fits();
        let original: Vec<u16> = vec![0, 100, 32767, 32768, 40000, 65535];
        let img = Image::new(PixelData::U16(original.clone()), vec![3usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);

        let header = stored_header(tmp.path());
        assert_eq!(header.get("BZERO").as_f64(), Some(32768.0));
        assert_eq!(header.get("BSCALE").as_f64(), Some(1.0));
        assert!(restored.get_keyword("BZERO").is_none());
        assert!(restored.get_keyword("BSCALE").is_none());
    }

    #[test]
    fn roundtrip_rgb_f32_stays_interleaved_in_memory() {
        // In-memory layout must be interleaved (RGB RGB RGB) both before
        // save and after load — the on-disk FITS layout is planar but
        // that conversion is handled transparently.
        let tmp = tmp_fits();
        let interleaved: Vec<f32> = vec![
            1.0, 2.0, 3.0, // pixel (0,0)
            4.0, 5.0, 6.0, // pixel (1,0)
            7.0, 8.0, 9.0, // pixel (0,1)
            10.0, 11.0, 12.0, // pixel (1,1)
        ];
        let img = Image::new(PixelData::F32(interleaved.clone()), vec![2usize, 2, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert!(restored.is_rgb());
        assert_eq!(restored.dimensions, vec![2usize, 2, 3]);
        assert_eq!(restored.pixels.as_f32().unwrap(), &interleaved);
    }

    #[test]
    fn keywords_roundtrip() {
        let tmp = tmp_fits();
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.set_keyword(Keyword::string("OBJECT", "M31"));
        img.set_keyword(Keyword::real("EXPTIME", 30.5));
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.get_keyword("OBJECT").unwrap().value,
            Some(KeywordValue::String("M31".to_string()))
        );
        assert_eq!(
            restored.get_keyword("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(30.5))
        );
    }

    #[test]
    fn string_keywords_survive_repeated_save_and_open() {
        let long_quoted: String = (0..150)
            .map(|i| {
                if i % 7 == 0 {
                    '\''
                } else {
                    (b'a' + (i % 26) as u8) as char
                }
            })
            .collect();
        let values = [
            ("OBSERVER", "O'Brien".to_string()),
            ("QUOTES", "''''".to_string()),
            ("OBJECT", "M31 / Andromeda".to_string()),
            ("LEADSP", "  lead".to_string()),
            ("LONGQ", long_quoted),
            ("LONGSL", format!("{} / {}", "x".repeat(70), "y".repeat(30))),
        ];
        let mut img = Image::new(PixelData::I16(vec![0; 4]), vec![2usize, 2]);
        for (name, value) in &values {
            img.set_keyword(Keyword::string(*name, value.as_str()).with_comment("c"));
        }

        let first = tmp_fits();
        let second = tmp_fits();
        img.save(first.path()).unwrap();
        Image::open(first.path())
            .unwrap()
            .save(second.path())
            .unwrap();

        for path in [first.path(), second.path()] {
            let restored = Image::open(path).unwrap();
            for (name, value) in &values {
                assert_eq!(
                    restored.get(name).as_str(),
                    Some(value.as_str()),
                    "{}",
                    name
                );
                assert_eq!(
                    restored.get_keyword(name).unwrap().comment.as_deref(),
                    Some("c"),
                    "{}",
                    name
                );
            }
        }
        let cards = |path: &std::path::Path| -> Vec<Vec<u8>> {
            let bytes = std::fs::read(path).unwrap();
            bytes[..2880]
                .chunks(80)
                .filter(|c| c.iter().any(|&b| b != b' '))
                .map(<[u8]>::to_vec)
                .collect()
        };
        assert_eq!(
            String::from_utf8(cards(first.path()).concat()).unwrap(),
            String::from_utf8(cards(second.path()).concat()).unwrap()
        );
    }

    fn hdu_bytes(cards: &[&str], data: &[u8]) -> Vec<u8> {
        let mut bytes: Vec<u8> = cards
            .iter()
            .chain(std::iter::once(&"END"))
            .flat_map(|c| format!("{:<80}", c).into_bytes())
            .collect();
        bytes.resize(bytes.len().div_ceil(2880) * 2880, b' ');
        bytes.extend_from_slice(data);
        bytes.resize(bytes.len().div_ceil(2880) * 2880, 0);
        bytes
    }

    fn file_with(parts: &[Vec<u8>]) -> tempfile::NamedTempFile {
        let tmp = tmp_fits();
        std::fs::write(tmp.path(), parts.concat()).unwrap();
        tmp
    }

    const EMPTY_PRIMARY: &[&str] = &["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 0", "EXTEND  = T"];

    const I16_PAIR_EXTENSION: &[&str] = &[
        "XTENSION= 'IMAGE   '",
        "BITPIX  = 16",
        "NAXIS   = 1",
        "NAXIS1  = 2",
        "PCOUNT  = 0",
        "GCOUNT  = 1",
        "EXTNAME = 'SCI     '",
    ];

    fn table_extension() -> Vec<u8> {
        hdu_bytes(
            &[
                "XTENSION= 'BINTABLE'",
                "BITPIX  = 8",
                "NAXIS   = 2",
                "NAXIS1  = 4",
                "NAXIS2  = 1",
                "PCOUNT  = 0",
                "GCOUNT  = 1",
                "TFIELDS = 1",
                "TFORM1  = '1J      '",
            ],
            &[0, 0, 0, 9],
        )
    }

    #[test]
    fn open_reads_first_image_extension_after_empty_primary() {
        let tmp = file_with(&[
            hdu_bytes(EMPTY_PRIMARY, &[]),
            hdu_bytes(I16_PAIR_EXTENSION, &[0, 5, 0, 6]),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.dimensions, vec![2usize]);
        assert_eq!(img.pixels.as_i16(), Some(&vec![5, 6]));
        assert_eq!(img.get("EXTNAME").as_str(), Some("SCI"));
    }

    #[test]
    fn open_skips_table_before_image_extension() {
        let tmp = file_with(&[
            hdu_bytes(EMPTY_PRIMARY, &[]),
            table_extension(),
            hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![7, 8]));
    }

    #[test]
    fn open_skips_malformed_and_ascii_tables_before_image_extension() {
        let wide_columns = hdu_bytes(
            &[
                "XTENSION= 'BINTABLE'",
                "BITPIX  = 8",
                "NAXIS   = 2",
                "NAXIS1  = 4",
                "NAXIS2  = 1",
                "PCOUNT  = 0",
                "GCOUNT  = 1",
                "TFIELDS = 1",
                "TFORM1  = '2J      '",
            ],
            &[0; 4],
        );
        let ascii = hdu_bytes(
            &[
                "XTENSION= 'TABLE   '",
                "BITPIX  = 8",
                "NAXIS   = 2",
                "NAXIS1  = 3",
                "NAXIS2  = 1",
                "PCOUNT  = 0",
                "GCOUNT  = 1",
                "TFIELDS = 1",
                "TBCOL1  = 1",
                "TFORM1  = 'I3      '",
            ],
            b" 42",
        );
        let tmp = file_with(&[
            hdu_bytes(EMPTY_PRIMARY, &[]),
            wide_columns,
            ascii,
            hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![7, 8]));
    }

    #[test]
    fn open_errors_when_no_hdu_holds_image_data() {
        let tmp = file_with(&[hdu_bytes(EMPTY_PRIMARY, &[]), table_extension()]);
        assert!(matches!(
            Image::open(tmp.path()),
            Err(ImageError::Fits(FitsError::InvalidFormat(_)))
        ));
    }

    fn i16_primary(axis_length: &str, extra: &[&str], data: &[u8]) -> Vec<u8> {
        let mut cards = vec!["SIMPLE  = T", "BITPIX  = 16", "NAXIS   = 1", axis_length];
        cards.extend_from_slice(extra);
        hdu_bytes(&cards, data)
    }

    #[test]
    fn open_reads_standard_keywords_beside_hierarch_namesakes() {
        let tmp = file_with(&[i16_primary(
            "NAXIS1  = 4",
            &[
                "HIERARCH NAXIS1 = 1",
                "HIERARCH BITPIX = 8",
                "HIERARCH BSCALE = 3",
                "HIERARCH foo = 1",
                "HIERARCH a.b = 'x'",
            ],
            &[0, 5, 0, 6, 0, 7, 0, 8],
        )]);
        let resaved = tmp_fits();
        for path in [tmp.path(), resaved.path()] {
            let img = Image::open(path).unwrap();
            assert_eq!(img.dimensions, vec![4]);
            assert_eq!(img.pixels.as_i16(), Some(&vec![5, 6, 7, 8]));
            let names: Vec<&str> = img
                .keywords
                .iter()
                .map(|k| k.name.as_str())
                .filter(|n| n.starts_with("HIERARCH"))
                .collect();
            assert_eq!(
                names,
                [
                    "HIERARCH NAXIS1",
                    "HIERARCH BITPIX",
                    "HIERARCH BSCALE",
                    "HIERARCH foo",
                    "HIERARCH a.b"
                ]
            );
            img.save(resaved.path()).unwrap();
        }
        let bytes = std::fs::read(resaved.path()).unwrap();
        assert!(bytes
            .chunks(80)
            .any(|card| card.starts_with(b"HIERARCH foo = 1 ")));
    }

    #[test]
    fn open_returns_primary_before_malformed_extension() {
        let tmp = file_with(&[
            i16_primary("NAXIS1  = 2", &[], &[0, 5, 0, 6]),
            hdu_bytes(
                &[
                    "XTENSION= 'IMAGE   '",
                    "BITPIX  = 16",
                    "NAXIS   = 1",
                    "NAXIS1  = 2000",
                ],
                &[],
            ),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![5, 6]));

        let garbage = file_with(&[
            i16_primary("NAXIS1  = 2", &[], &[0, 5, 0, 6]),
            b"garbage!".repeat(360),
        ]);
        let img = Image::open(garbage.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![5, 6]));
    }

    #[test]
    fn open_returns_scan_error_when_no_image_precedes_it() {
        let tmp = file_with(&[hdu_bytes(EMPTY_PRIMARY, &[]), b"garbage!".repeat(360)]);
        match Image::open(tmp.path()) {
            Err(ImageError::Fits(FitsError::InvalidFormat(message))) => {
                assert!(message.contains("END"), "{}", message)
            }
            other => panic!("expected the scan error, got {:?}", other.map(|_| ())),
        }
    }

    #[test]
    fn open_skips_image_extension_with_gcount_zero() {
        let tmp = file_with(&[
            hdu_bytes(EMPTY_PRIMARY, &[]),
            hdu_bytes(
                &[
                    "XTENSION= 'IMAGE   '",
                    "BITPIX  = 16",
                    "NAXIS   = 1",
                    "NAXIS1  = 2",
                    "PCOUNT  = 0",
                    "GCOUNT  = 0",
                ],
                &[],
            ),
            hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![7, 8]));
    }

    #[test]
    fn open_skips_unknown_extension_type() {
        let tmp = file_with(&[
            hdu_bytes(EMPTY_PRIMARY, &[]),
            hdu_bytes(
                &[
                    "XTENSION= 'FOREIGN '",
                    "BITPIX  = 16",
                    "NAXIS   = 1",
                    "NAXIS1  = 2",
                    "PCOUNT  = 0",
                    "GCOUNT  = 1",
                ],
                &[0, 1, 0, 2],
            ),
            hdu_bytes(I16_PAIR_EXTENSION, &[0, 7, 0, 8]),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i16(), Some(&vec![7, 8]));
    }

    #[test]
    fn open_decodes_tile_compressed_image() {
        let tmp = tmp_fits();
        let data: Vec<i32> = (0..70).map(|i| 1000 + (i % 9) - 4).collect();
        FitsWriter::create(tmp.path())
            .unwrap()
            .write_compressed_image(
                &data,
                &[10, 7],
                Some((4, 3)),
                crate::fits::compression::CompressionAlgorithm::Rice,
                &[],
            )
            .unwrap();
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.dimensions, vec![10, 7]);
        assert_eq!(img.pixels.as_i32(), Some(&data));
    }

    #[test]
    fn open_returns_fits_error_for_malformed_fits_without_fits_extension() {
        let tmp = tempfile::Builder::new().suffix(".dat").tempfile().unwrap();
        std::fs::write(
            tmp.path(),
            hdu_bytes(
                &["SIMPLE  = T", "BITPIX  = 8", "NAXIS   = 2", "NAXIS1  = 4"],
                &[1, 2, 3, 4],
            ),
        )
        .unwrap();
        assert!(matches!(
            Image::open(tmp.path()),
            Err(ImageError::Fits(FitsError::KeywordNotFound { .. }))
        ));
    }

    #[test]
    fn open_rejects_data_size_beyond_file_length() {
        let tmp = file_with(&[hdu_bytes(
            &[
                "SIMPLE  = T",
                "BITPIX  = 16",
                "NAXIS   = 2",
                "NAXIS1  = 1000000000",
                "NAXIS2  = 1000000000",
            ],
            &[0; 2880],
        )]);
        assert!(matches!(Image::open(tmp.path()), Err(ImageError::Fits(_))));
    }

    #[test]
    fn image_from_extension_saves_as_valid_primary() {
        let tmp = file_with(&[
            hdu_bytes(EMPTY_PRIMARY, &[]),
            hdu_bytes(I16_PAIR_EXTENSION, &[0, 5, 0, 6]),
        ]);
        let img = Image::open(tmp.path()).unwrap();
        let resaved = tmp_fits();
        img.save(resaved.path()).unwrap();
        let header = stored_header(resaved.path());
        assert!(header.is_primary());
        for name in ["XTENSION", "PCOUNT", "GCOUNT"] {
            assert!(header.get_keyword(name).is_none(), "{}", name);
        }
        assert_eq!(
            Image::open(resaved.path()).unwrap().pixels.as_i16(),
            Some(&vec![5, 6])
        );
    }

    #[test]
    fn repeated_save_keeps_one_copy_of_standard_comments() {
        let first = tmp_fits();
        Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2])
            .save(first.path())
            .unwrap();
        let second = tmp_fits();
        Image::open(first.path())
            .unwrap()
            .save(second.path())
            .unwrap();
        let third = tmp_fits();
        Image::open(second.path())
            .unwrap()
            .save(third.path())
            .unwrap();
        let comments = Image::open(third.path())
            .unwrap()
            .keywords
            .iter()
            .filter(|k| k.name == "COMMENT")
            .count();
        assert_eq!(comments, 2);
    }

    #[test]
    fn open_uses_last_bzero_and_bscale_when_they_repeat() {
        let tmp = file_with(&[i16_primary(
            "NAXIS1  = 2",
            &["BZERO   = 0", "BZERO   = 100", "BSCALE  = 1", "BSCALE  = 2"],
            &[0, 5, 0, 6],
        )]);
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_f64(), Some(&vec![110.0, 112.0]));
        assert!(img.get_keyword("BZERO").is_none());
    }

    fn cards_named(path: &Path, name: &str) -> Vec<Option<KeywordValue>> {
        let header = stored_header(path);
        let named = header.keywords().iter().filter(|k| k.name == name);
        named.map(|k| k.value.clone()).collect()
    }

    #[test]
    fn set_keyword_replaces_every_copy_of_a_repeated_wcs_card() {
        let tmp = file_with(&[i16_primary(
            "NAXIS1  = 2",
            &["CRVAL1  = 83.6", "CRVAL1  = 84.0"],
            &[0, 5, 0, 6],
        )]);
        let mut img = Image::open(tmp.path()).unwrap();
        img.set_keyword(Keyword::real("CRVAL1", 50.0));
        img.save(tmp.path()).unwrap();
        let crval1 = cards_named(tmp.path(), "CRVAL1");
        assert_eq!(crval1, vec![Some(KeywordValue::Real(50.0))]);
    }

    #[test]
    fn set_keyword_replaces_every_copy_of_a_repeated_blank() {
        let tmp = file_with(&[i16_primary(
            "NAXIS1  = 3",
            &["BLANK   = 5", "BLANK   = 6"],
            &[0, 5, 0, 6, 0, 7],
        )]);
        let mut img = Image::open(tmp.path()).unwrap();
        img.set_keyword(Keyword::integer("BLANK", 7));
        assert_eq!(img.blank_flags(), Some(vec![false, false, true]));
        img.save(tmp.path()).unwrap();
        let blank = cards_named(tmp.path(), "BLANK");
        assert_eq!(blank, vec![Some(KeywordValue::Integer(7))]);
    }

    #[test]
    fn open_uses_last_blank_when_it_repeats() {
        let tmp = file_with(&[i16_primary(
            "NAXIS1  = 2",
            &["BLANK   = 5", "BLANK   = 6"],
            &[0, 5, 0, 6],
        )]);
        let img = Image::open(tmp.path()).unwrap();
        let blanks: Vec<_> = img
            .keywords
            .iter()
            .filter(|k| k.name == "BLANK")
            .map(|k| k.value.clone())
            .collect();
        assert_eq!(blanks, vec![Some(KeywordValue::Integer(6)); 2]);
    }

    // NAXIS1 = 7, NAXIS2 = 5, NAXIS3 = 3, NAXIS4 = 2.
    fn four_axis_cube() -> (Vec<i32>, Vec<usize>) {
        ((0..210).collect(), vec![7, 5, 3, 2])
    }

    #[test]
    fn open_keeps_the_planes_of_a_four_axis_image_as_stored() {
        let tmp = tmp_fits();
        let (data, dims) = four_axis_cube();
        FitsWriter::create(tmp.path())
            .unwrap()
            .write_primary_image(&data, &dims, &[])
            .unwrap();

        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.dimensions, dims);
        assert_eq!(img.pixels.as_i32(), Some(&data));
    }

    #[test]
    fn open_interleaves_the_planes_of_a_three_axis_rgb_image() {
        let tmp = tmp_fits();
        FitsWriter::create(tmp.path())
            .unwrap()
            .write_primary_image(&[1i32, 4, 2, 5, 3, 6], &[2, 1, 3], &[])
            .unwrap();

        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_i32(), Some(&vec![1, 2, 3, 4, 5, 6]));
    }

    #[test]
    fn four_axis_image_is_saved_as_stored_and_round_trips() {
        let tmp = tmp_fits();
        let (data, dims) = four_axis_cube();
        let img = Image::new(PixelData::I32(data.clone()), dims.clone());
        img.save(tmp.path()).unwrap();

        let stored = FitsFile::open(tmp.path())
            .unwrap()
            .read_data_unit(0)
            .unwrap();
        let expected: Vec<u8> = data.iter().flat_map(|v| v.to_be_bytes()).collect();
        assert_eq!(&stored[..expected.len()], expected.as_slice());
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, dims);
        assert_eq!(restored.pixels.as_i32(), Some(&data));

        for algorithm in [
            CompressionAlgorithm::Rice,
            CompressionAlgorithm::Gzip,
            CompressionAlgorithm::Gzip2,
        ] {
            for tile_size in [None, Some((3, 2))] {
                img.save_fits_compressed(tmp.path(), algorithm, tile_size)
                    .unwrap();
                let restored = Image::open(tmp.path()).unwrap();
                assert_eq!(restored.dimensions, dims, "{algorithm:?} {tile_size:?}");
                assert_eq!(restored.pixels.as_i32(), Some(&data), "{algorithm:?}");
            }
        }
    }
}
