use super::Image;
use crate::core::errors::{ImageError, Result};
use crate::core::shape::{plane_dimensions, plane_shape_with_axes};
use crate::formats::atomic_write::write_atomically;
use crate::formats::pixel_data::PixelData;
use std::io::Cursor;
use std::path::Path;
use tiff::encoder::{colortype, TiffEncoder};
use tiff::ColorType;

impl Image {
    pub fn open_tiff<P: AsRef<Path>>(path: P) -> Result<Self> {
        use tiff::decoder::Decoder;

        let file = std::fs::File::open(path)?;
        let limits = decoding_limits(file.metadata()?.len());
        let mut decoder = Decoder::new(file)
            .map_err(decode_error)?
            .with_limits(limits);
        let (width, height) = decoder.dimensions().map_err(decode_error)?;

        // Samples-per-pixel comes from the TIFF colortype, not from the
        // decoded buffer (which is always flat). Without this, RGB TIFFs
        // get reported as 1-channel and downstream code treats the
        // interleaved R/G/B samples as 3× as many grayscale pixels.
        let colortype = decoder.colortype().map_err(decode_error)?;
        check_sample_depth(colortype)?;
        let image = decoder.read_image().map_err(decode_error)?;

        let pixels = decode_tiff_image(image)?;

        // `read_image` returns every stored sample, including extra samples the
        // colortype doesn't account for — most commonly an associated-alpha
        // channel on an RGB image (PixInsight writes RGBA with photometric=RGB,
        // so colortype() reports RGB / 3 samples while the buffer holds 4). The
        // true sample count is the buffer length per pixel; reconcile against
        // the colortype and drop the trailing extra samples so the in-memory
        // image is a clean 1- or 3-channel buffer.
        let pixel_count = (width as usize) * (height as usize);
        let actual = actual_samples(pixels.len(), pixel_count);
        let samples = reconcile_samples(pixels.len(), pixel_count, colortype)?;
        let pixels = pixels.keep_leading_samples(actual, samples);
        let dimensions = plane_dimensions(width as usize, height as usize, samples);
        Ok(Self::new(pixels, dimensions))
    }

    // The file is encoded in memory first, so a pixel type TIFF cannot hold
    // returns an error without truncating an existing file at `path`.
    pub(super) fn save_tiff(&self, path: &Path) -> Result<()> {
        let (width, height, channels) =
            plane_shape_with_axes(&self.dimensions, ImageError::InvalidDimensions)?;
        let mut buffer = Cursor::new(Vec::new());
        let mut encoder = TiffEncoder::new(&mut buffer).map_err(encode_error)?;

        write_tiff_image(
            &self.pixels,
            &mut encoder,
            width as u32,
            height as u32,
            channels,
        )?;
        Ok(write_atomically(path, &buffer.into_inner())?)
    }
}

// The decoder refuses an image whose declared decoded size passes
// `decoding_buffer_size`. Its 256 MiB default rejects ordinary full-frame
// RGB16 files. LZW, the TIFF compression with the largest ratio, emits at
// most a 4096-byte string for each 12-bit code, 2731 times the bytes it
// reads, so a file of `file_len` bytes decodes to at most 2731 times that.
fn decoding_limits(file_len: u64) -> tiff::decoder::Limits {
    const FLOOR: usize = 256 * 1024 * 1024;
    let mut limits = tiff::decoder::Limits::default();
    limits.decoding_buffer_size = file_len
        .checked_mul(2731)
        .and_then(|bytes| usize::try_from(bytes).ok())
        .map_or(usize::MAX, |bytes| bytes.max(FLOOR));
    limits
}

// The decoder returns samples of other widths bit-packed, which the
// reconciliation below would misread as whole samples.
fn check_sample_depth(colortype: ColorType) -> Result<()> {
    let bits = match colortype {
        ColorType::Gray(bits)
        | ColorType::GrayA(bits)
        | ColorType::RGB(bits)
        | ColorType::RGBA(bits)
        | ColorType::CMYK(bits)
        | ColorType::CMYKA(bits)
        | ColorType::YCbCr(bits)
        | ColorType::Lab(bits)
        | ColorType::Palette(bits) => bits,
        ColorType::Multiband { bit_depth, .. } => bit_depth,
        _ => return Ok(()),
    };
    if matches!(bits, 8 | 16 | 32 | 64) {
        return Ok(());
    }
    Err(ImageError::Tiff(format!(
        "{bits}-bit samples are not supported; samples must be 8, 16, 32 or 64 bits"
    )))
}

/// Samples per pixel actually present in the decoded buffer.
fn actual_samples(buffer_len: usize, pixel_count: usize) -> usize {
    if pixel_count == 0 {
        return 0;
    }
    buffer_len / pixel_count
}

/// Resolve how many channels to keep, given the decoded buffer's true sample
/// count and the colortype. When the buffer carries more samples than the
/// colortype (extra/alpha samples), keep only the logical channels (1 for
/// gray-like, 3 for RGB-like) and drop the rest. When the two agree, keep the
/// colortype's channels minus any alpha it declares.
fn reconcile_samples(buffer_len: usize, pixel_count: usize, colortype: ColorType) -> Result<usize> {
    let actual = actual_samples(buffer_len, pixel_count);
    if pixel_count == 0 || actual == 0 || !buffer_len.is_multiple_of(pixel_count) {
        return Err(ImageError::Tiff(format!(
            "{buffer_len} decoded samples do not divide evenly among {pixel_count} pixels"
        )));
    }
    let logical = samples_per_pixel(colortype);
    if actual <= logical {
        return Ok(kept_samples(colortype).min(actual));
    }
    // More samples than the colortype expects: keep the logical channels.
    // Gray-like (1–2) collapse to 1, RGB-like (3–4) collapse to 3.
    Ok(if logical >= 3 { 3 } else { 1 })
}

/// Samples per pixel the colortype itself declares.
fn samples_per_pixel(colortype: ColorType) -> usize {
    match colortype {
        ColorType::Gray(_) | ColorType::Palette(_) => 1,
        ColorType::GrayA(_) => 2,
        ColorType::RGB(_) | ColorType::YCbCr(_) | ColorType::Lab(_) => 3,
        ColorType::RGBA(_) | ColorType::CMYK(_) => 4,
        ColorType::CMYKA(_) => 5,
        ColorType::Multiband { num_samples, .. } => num_samples as usize,
        _ => 1,
    }
}

/// Channels kept in memory: an alpha channel is dropped, so GrayA lands as
/// gray and RGBA as RGB. Colortypes with no alpha keep every sample.
fn kept_samples(colortype: ColorType) -> usize {
    match colortype {
        ColorType::GrayA(_) => 1,
        ColorType::RGBA(_) => 3,
        ColorType::CMYKA(_) => 4,
        other => samples_per_pixel(other),
    }
}

fn decode_tiff_image(image: tiff::decoder::DecodingResult) -> Result<PixelData> {
    use tiff::decoder::DecodingResult;

    match image {
        DecodingResult::U8(data) => Ok(PixelData::U8(data)),
        DecodingResult::U16(data) => Ok(PixelData::U16(data)),
        DecodingResult::U32(data) => Ok(PixelData::U32(data)),
        DecodingResult::U64(data) => Ok(PixelData::U64(data)),
        DecodingResult::I8(data) => Ok(PixelData::I8(data)),
        DecodingResult::I32(data) => Ok(PixelData::I32(data)),
        DecodingResult::I64(data) => Ok(PixelData::I64(data)),
        DecodingResult::F32(data) => Ok(PixelData::F32(data)),
        DecodingResult::F64(data) => Ok(PixelData::F64(data)),
        _ => Err(ImageError::UnsupportedFormat),
    }
}

fn write_tiff_image(
    pixels: &PixelData,
    encoder: &mut TiffEncoder<&mut Cursor<Vec<u8>>>,
    width: u32,
    height: u32,
    channels: usize,
) -> Result<()> {
    match (pixels, channels) {
        (PixelData::U8(data), 1) => encoder.write_image::<colortype::Gray8>(width, height, data),
        (PixelData::U8(data), 3) => encoder.write_image::<colortype::RGB8>(width, height, data),
        (PixelData::U16(data), 1) => encoder.write_image::<colortype::Gray16>(width, height, data),
        (PixelData::U16(data), 3) => encoder.write_image::<colortype::RGB16>(width, height, data),
        (PixelData::I8(data), 1) => encoder.write_image::<colortype::GrayI8>(width, height, data),
        (PixelData::U32(data), 1) => encoder.write_image::<colortype::Gray32>(width, height, data),
        (PixelData::U32(data), 3) => encoder.write_image::<colortype::RGB32>(width, height, data),
        (PixelData::I32(data), 1) => encoder.write_image::<colortype::GrayI32>(width, height, data),
        (PixelData::U64(data), 1) => encoder.write_image::<colortype::Gray64>(width, height, data),
        (PixelData::U64(data), 3) => encoder.write_image::<colortype::RGB64>(width, height, data),
        (PixelData::I64(data), 1) => encoder.write_image::<colortype::GrayI64>(width, height, data),
        (PixelData::F32(data), 1) => {
            encoder.write_image::<colortype::Gray32Float>(width, height, data)
        }
        (PixelData::F32(data), 3) => {
            encoder.write_image::<colortype::RGB32Float>(width, height, data)
        }
        _ => return Err(ImageError::UnsupportedFormat),
    }
    .map_err(encode_error)
}

fn decode_error(error: tiff::TiffError) -> ImageError {
    ImageError::Tiff(format!("decode: {error}"))
}

fn encode_error(error: tiff::TiffError) -> ImageError {
    ImageError::Tiff(format!("encode: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_tiff() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".tiff").tempfile().unwrap()
    }

    #[test]
    fn roundtrip_u8_gray() {
        let tmp = tmp_tiff();
        let original: Vec<u8> = (0u8..64).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u16_gray() {
        let tmp = tmp_tiff();
        let original: Vec<u16> = (0u16..64).map(|i| i * 100).collect();
        let img = Image::new(PixelData::U16(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f32_gray() {
        let tmp = tmp_tiff();
        let original: Vec<f32> = (0..64).map(|i| i as f32 * 0.25).collect();
        let img = Image::new(PixelData::F32(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
    }

    #[test]
    fn open_reads_gray32_samples_as_u32() {
        let tmp = tmp_tiff();
        let file = std::fs::File::create(tmp.path()).unwrap();
        TiffEncoder::new(file)
            .unwrap()
            .write_image::<colortype::Gray32>(2, 1, &[3000000000, 7])
            .unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u32(), Some(&vec![3000000000, 7]));
    }

    #[test]
    fn roundtrip_u32_gray() {
        let tmp = tmp_tiff();
        let original = vec![0u32, 7, 2147483648, 4294967295];
        let img = Image::new(PixelData::U32(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u32(), Some(&original));
    }

    #[test]
    fn roundtrip_i32_gray_keeps_sign() {
        let tmp = tmp_tiff();
        let original = vec![0i32, 1000, -500, i32::MIN];
        let img = Image::new(PixelData::I32(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        // TIFF 6.0 SampleFormat: 1 is unsigned integer, 2 is signed integer.
        assert_eq!(sample_format_tag(tmp.path()), 2);
        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i32(), Some(&original));
    }

    #[test]
    fn roundtrip_i8_gray() {
        let tmp = tmp_tiff();
        let original = vec![-128i8, -1, 0, 127];
        let img = Image::new(PixelData::I8(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_i8(), Some(&original));
    }

    #[test]
    fn failed_save_leaves_existing_file_untouched() {
        let tmp = tmp_tiff();
        std::fs::write(tmp.path(), b"keep").unwrap();
        for pixels in [PixelData::I8(vec![0; 12]), PixelData::I64(vec![0; 12])] {
            let img = Image::new(pixels, vec![2usize, 2, 3]);
            assert!(img.save(tmp.path()).is_err());
            assert_eq!(std::fs::read(tmp.path()).unwrap(), b"keep");
        }
    }

    fn sample_format_tag(path: &Path) -> u32 {
        let file = std::fs::File::open(path).unwrap();
        tiff::decoder::Decoder::new(file)
            .unwrap()
            .get_tag_u32(tiff::tags::Tag::SampleFormat)
            .unwrap()
    }

    #[test]
    fn roundtrip_64_bit_integers_gray() {
        let tmp = tmp_tiff();
        let unsigned = vec![0u64, (1 << 53) + 1, 1 << 63, u64::MAX];
        Image::new(PixelData::U64(unsigned.clone()), vec![2usize, 2])
            .save(tmp.path())
            .unwrap();
        assert_eq!(sample_format_tag(tmp.path()), 1);
        assert_eq!(
            Image::open_tiff(tmp.path()).unwrap().pixels.as_u64(),
            Some(&unsigned)
        );

        let signed = vec![i64::MIN, -1, (1 << 53) + 1, i64::MAX];
        Image::new(PixelData::I64(signed.clone()), vec![2usize, 2])
            .save(tmp.path())
            .unwrap();
        assert_eq!(sample_format_tag(tmp.path()), 2);
        assert_eq!(
            Image::open_tiff(tmp.path()).unwrap().pixels.as_i64(),
            Some(&signed)
        );
    }

    #[test]
    fn roundtrip_rgb_u32_and_u64() {
        let tmp = tmp_tiff();
        let u32s: Vec<u32> = (0..12).map(|i| i * 350_000_000).collect();
        Image::new(PixelData::U32(u32s.clone()), vec![2usize, 2, 3])
            .save(tmp.path())
            .unwrap();
        assert_eq!(
            Image::open_tiff(tmp.path()).unwrap().pixels.as_u32(),
            Some(&u32s)
        );

        let u64s: Vec<u64> = (0..12).map(|i| i * 1_500_000_000_000_000_000).collect();
        Image::new(PixelData::U64(u64s.clone()), vec![2usize, 2, 3])
            .save(tmp.path())
            .unwrap();
        assert_eq!(
            Image::open_tiff(tmp.path()).unwrap().pixels.as_u64(),
            Some(&u64s)
        );
    }

    #[test]
    fn save_rgb_u8() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::U8(vec![0; 48]), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();
    }

    #[test]
    fn save_rgb_u16() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::U16(vec![0; 48]), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();
    }

    #[test]
    fn roundtrip_rgb_u8_reports_three_channels() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::U8((0u8..48).collect()), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert!(
            restored.is_rgb(),
            "RGB TIFF should round-trip with channels==3"
        );
        assert_eq!(restored.channels(), 3);
        assert_eq!(restored.width(), 4);
        assert_eq!(restored.height(), 4);
        assert_eq!(restored.pixels.as_u8().unwrap().len(), 48);
    }

    #[test]
    fn roundtrip_rgb_u16_reports_three_channels() {
        let tmp = tmp_tiff();
        let img = Image::new(
            PixelData::U16((0u16..48).map(|i| i * 1000).collect()),
            vec![4usize, 4, 3],
        );
        img.save(tmp.path()).unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert!(restored.is_rgb());
        assert_eq!(restored.channels(), 3);
        assert_eq!(restored.pixels.as_u16().unwrap().len(), 48);
    }

    #[test]
    fn save_rejects_unsupported_combos() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::I16(vec![0; 16]), vec![4usize, 4]);
        assert!(img.save(tmp.path()).is_err());
    }

    #[test]
    fn reconcile_keeps_matching_samples() {
        // 3 samples/pixel buffer, colortype says RGB → keep 3.
        assert_eq!(reconcile_samples(4 * 3, 4, ColorType::RGB(8)).unwrap(), 3);
        // 1 sample/pixel gray.
        assert_eq!(reconcile_samples(4, 4, ColorType::Gray(8)).unwrap(), 1);
    }

    #[test]
    fn reconcile_strips_rgb_alpha() {
        // PixInsight RGBA: buffer has 4 samples/pixel but colortype reports
        // RGB (3). Keep 3, drop the assoc-alpha.
        assert_eq!(reconcile_samples(4 * 4, 4, ColorType::RGB(8)).unwrap(), 3);
    }

    #[test]
    fn reconcile_strips_declared_alpha() {
        // Colortype and buffer agree on 4 samples: it is a real RGBA TIFF, and
        // the alpha still has to go so the image lands as 3-channel RGB.
        assert_eq!(reconcile_samples(4 * 4, 4, ColorType::RGBA(8)).unwrap(), 3);
        assert_eq!(reconcile_samples(4 * 2, 4, ColorType::GrayA(8)).unwrap(), 1);
    }

    #[test]
    fn reconcile_strips_gray_alpha() {
        // GrayA stored against a Gray colortype: 2 samples/pixel, logical 1.
        assert_eq!(reconcile_samples(4 * 2, 4, ColorType::Gray(8)).unwrap(), 1);
    }

    #[test]
    fn reconcile_rejects_non_integer_samples() {
        // A buffer that isn't a whole number of samples per pixel.
        assert!(reconcile_samples(10, 4, ColorType::RGB(8)).is_err());
    }

    #[test]
    fn open_rgba_tiff_drops_alpha() {
        let tmp = tmp_tiff();
        // One row of two pixels: red then green, both opaque.
        let rgba: Vec<u8> = vec![10, 20, 30, 255, 40, 50, 60, 255];
        let file = std::fs::File::create(tmp.path()).unwrap();
        TiffEncoder::new(file)
            .unwrap()
            .write_image::<colortype::RGBA8>(2, 1, &rgba)
            .unwrap();

        let restored = Image::open_tiff(tmp.path()).unwrap();
        assert!(restored.is_rgb());
        assert_eq!(restored.channels(), 3);
        assert_eq!(
            restored.pixels.as_u8().unwrap(),
            &vec![10, 20, 30, 40, 50, 60]
        );
    }

    #[test]
    fn open_returns_error_on_invalid_tiff() {
        let tmp = tempfile::Builder::new().suffix(".tiff").tempfile().unwrap();
        std::fs::write(tmp.path(), b"garbage").unwrap();
        assert!(Image::open_tiff(tmp.path()).is_err());
    }

    // A little-endian, uncompressed gray TIFF with `bits` per sample. Every
    // strip of `rows_per_strip` rows points at the same `strip` bytes, so the
    // file can declare an image far larger than itself. The file is padded
    // with zeros to at least `min_len` bytes.
    fn shared_strip_tiff(
        (width, height): (u32, u32),
        bits: u32,
        rows_per_strip: u32,
        strip: &[u8],
        min_len: usize,
    ) -> Vec<u8> {
        let strips = height.div_ceil(rows_per_strip);
        let offsets_at = 8 + 2 + 8 * 12 + 4;
        let counts_at = offsets_at + 4 * strips;
        let strip_at = counts_at + 4 * strips;
        let entries: [(u16, u16, u32, u32); 8] = [
            (256, 4, 1, width),
            (257, 4, 1, height),
            (258, 3, 1, bits),
            (259, 3, 1, 1),
            (262, 3, 1, 1),
            (273, 4, strips, offsets_at),
            (278, 4, 1, rows_per_strip),
            (279, 4, strips, counts_at),
        ];
        let mut file = b"II\x2a\x00\x08\x00\x00\x00\x08\x00".to_vec();
        for (tag, kind, count, value) in entries {
            file.extend_from_slice(&tag.to_le_bytes());
            file.extend_from_slice(&kind.to_le_bytes());
            file.extend_from_slice(&count.to_le_bytes());
            file.extend_from_slice(&value.to_le_bytes());
        }
        file.extend_from_slice(&[0; 4]);
        for value in [strip_at, strip.len() as u32] {
            for _ in 0..strips {
                file.extend_from_slice(&value.to_le_bytes());
            }
        }
        file.extend_from_slice(strip);
        file.resize(file.len().max(min_len), 0);
        file
    }

    fn open_bytes(bytes: &[u8]) -> Result<Image> {
        let tmp = tmp_tiff();
        std::fs::write(tmp.path(), bytes).unwrap();
        Image::open_tiff(tmp.path())
    }

    #[test]
    fn decoding_limit_scales_with_file_length_above_256_mib() {
        let limit = |len| decoding_limits(len).decoding_buffer_size;
        assert_eq!(limit(0), 256 * 1024 * 1024);
        assert_eq!(limit(98_000), 256 * 1024 * 1024);
        assert_eq!(limit(1_000_000), 2_731_000_000);
        assert_eq!(limit(u64::MAX), usize::MAX);
    }

    // 16384 x 16640 8-bit samples are 260 MiB, past the decoder's 256 MiB
    // default. The 1 MiB file allows up to 2.7 GiB.
    #[test]
    fn open_reads_an_image_past_256_mib_that_the_file_length_allows() {
        let strip: Vec<u8> = (0..16384 * 64).map(|i| (i % 251) as u8).collect();
        let image = open_bytes(&shared_strip_tiff((16384, 16640), 8, 64, &strip, 0)).unwrap();
        assert_eq!(image.dimensions, vec![16384usize, 16640]);
        let pixels = image.pixels.as_u8().unwrap();
        assert_eq!(pixels.len(), 16384 * 16640);
        assert_eq!(&pixels[..strip.len()], &strip[..]);
        assert_eq!(&pixels[pixels.len() - strip.len()..], &strip[..]);
    }

    // A 1,000,000-byte file allows 2,731,000,000 decoded bytes. The header
    // declares 2731 more.
    #[test]
    fn open_rejects_a_declared_size_past_the_file_length_bound() {
        let file = shared_strip_tiff((2731, 1_000_001), 8, 500_001, &[0; 4], 1_000_000);
        assert_eq!(file.len(), 1_000_000);
        let error = open_bytes(&file).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("TIFF error: decode: decoder limits exceeded")
        );
    }

    #[test]
    fn open_names_an_unsupported_sample_depth() {
        let error = open_bytes(&shared_strip_tiff((23, 37), 1, 1, &[0; 3], 0))
            .err()
            .map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some(
                "TIFF error: 1-bit samples are not supported; samples must be 8, 16, 32 or 64 bits"
            )
        );
    }

    #[test]
    fn decode_errors_name_tiff() {
        let tmp = tmp_tiff();
        std::fs::write(tmp.path(), b"garbage").unwrap();
        let error = Image::open_tiff(tmp.path()).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("TIFF error: decode: format error: TIFF signature not found")
        );
    }

    #[test]
    fn i32_saves_as_gray32() {
        let tmp = tmp_tiff();
        let img = Image::new(PixelData::I32(vec![0i32, 1000, -500, 42]), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();
    }
}
