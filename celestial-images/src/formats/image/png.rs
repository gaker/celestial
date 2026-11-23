use super::Image;
use crate::core::errors::{ImageError, Result};
use crate::core::shape::{plane_dimensions, plane_shape_with_axes};
use crate::formats::atomic_write::write_atomically;
use crate::formats::pixel_data::PixelData;
use std::borrow::Cow;
use std::io::BufReader;
use std::path::Path;

impl Image {
    pub fn open_png<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let file_len = file.metadata()?.len();
        let decoder = png::Decoder::new(BufReader::new(file));
        let mut reader = decoder.read_info().map_err(decode_error)?;
        check_inflated_size(reader.info(), file_len)?;

        // Palette and low-bit samples are expanded after decoding. The
        // decoder's EXPAND would size this buffer up to 32 times the
        // inflated data before any of that data is read.
        let buf_size = reader
            .output_buffer_size()
            .ok_or_else(|| ImageError::Png("the decoded image does not fit in memory".into()))?;
        let mut buf = zeroed(buf_size)?;
        let info = reader.next_frame(&mut buf).map_err(decode_error)?;
        buf.truncate(info.buffer_size());

        let palette = reader.info().palette.as_deref();
        let (pixels, dimensions) = parse_png_data(buf, &info, palette)?;
        Ok(Self::new(pixels, dimensions))
    }

    pub(super) fn save_png(&self, path: &Path) -> Result<()> {
        let (width, height, channels) =
            plane_shape_with_axes(&self.dimensions, ImageError::InvalidDimensions)?;
        let color_type = channels_to_png_color_type(channels)?;
        let (depth, samples) = match &self.pixels {
            PixelData::U8(data) => (png::BitDepth::Eight, Cow::Borrowed(data.as_slice())),
            PixelData::U16(data) => (
                png::BitDepth::Sixteen,
                Cow::Owned(data.iter().flat_map(|&v| v.to_be_bytes()).collect()),
            ),
            _ => return Err(ImageError::UnsupportedFormat),
        };
        let bytes = encode_png((width as u32, height as u32), color_type, depth, &samples)?;
        Ok(write_atomically(path, &bytes)?)
    }
}

// `samples` holds 8-bit samples, or 16-bit samples as big-endian byte pairs.
fn encode_png(
    (width, height): (u32, u32),
    color_type: png::ColorType,
    depth: png::BitDepth,
    samples: &[u8],
) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    let mut encoder = png::Encoder::new(&mut bytes, width, height);
    encoder.set_color(color_type);
    encoder.set_depth(depth);
    let mut writer = encoder.write_header().map_err(encode_error)?;
    writer.write_image_data(samples).map_err(encode_error)?;
    writer.finish().map_err(encode_error)?;
    Ok(bytes)
}

fn decode_error(error: png::DecodingError) -> ImageError {
    ImageError::Png(format!("decode: {error}"))
}

fn encode_error(error: png::EncodingError) -> ImageError {
    ImageError::Png(format!("encode: {error}"))
}

// The IDAT stream inflates to one filter byte plus the packed samples of
// each row. An interlaced image inflates to at least this many bytes.
fn check_inflated_size(info: &png::Info, file_len: u64) -> Result<()> {
    let (width, height) = (info.width, info.height);
    let inflated = u64::try_from(info.bits_per_pixel())
        .ok()
        .and_then(|bits| bits.checked_mul(u64::from(width)))
        .and_then(|bits| bits.div_ceil(8).checked_add(1))
        .and_then(|row| row.checked_mul(u64::from(height)));
    // 1032 to 1 is the largest ratio DEFLATE can compress to (zlib's technical notes).
    let limit = file_len.saturating_mul(1032);
    match inflated {
        Some(inflated) if inflated <= limit => Ok(()),
        _ => Err(ImageError::Png(format!(
            "{width}x{height} pixels inflate to {} bytes, more than a {file_len}-byte file can hold",
            inflated.map_or_else(|| format!("more than {}", u64::MAX), |n| n.to_string())
        ))),
    }
}

// `vec![0; len]` aborts the process when the allocation fails.
fn zeroed(len: usize) -> Result<Vec<u8>> {
    let mut buf = reserved(len)?;
    buf.resize(len, 0);
    Ok(buf)
}

// The whole output is reserved before it is filled. A vector that grows as
// it is filled aborts the process when a reallocation fails, and its
// capacity can reach twice the final length.
fn reserved<T>(len: usize) -> Result<Vec<T>> {
    let mut samples = Vec::new();
    samples
        .try_reserve_exact(len)
        .map_err(|_| ImageError::Png(format!("cannot allocate {len} samples for the image")))?;
    Ok(samples)
}

// Samples in the expanded image: width x height x `channels`.
fn expanded_len(info: &png::OutputInfo, channels: usize) -> Result<usize> {
    usize::try_from(u64::from(info.width) * u64::from(info.height))
        .ok()
        .and_then(|pixels| pixels.checked_mul(channels))
        .ok_or_else(|| {
            ImageError::Png(format!(
                "{}x{} pixels of {channels} samples do not fit in memory",
                info.width, info.height
            ))
        })
}

fn parse_png_data(
    buf: Vec<u8>,
    info: &png::OutputInfo,
    palette: Option<&[u8]>,
) -> Result<(PixelData, Vec<usize>)> {
    use png::BitDepth::{Eight, Four, One, Two};
    use png::ColorType::{Grayscale, Indexed};

    let (pixels, stored) = match (info.color_type, info.bit_depth) {
        (Indexed, One | Two | Four | Eight) => (expand_palette(&buf, info, palette)?, 3),
        (Grayscale, One | Two | Four) => (expand_gray(&buf, info)?, 1),
        (color_type, depth) => (decode_png_samples(buf, depth)?, stored_samples(color_type)?),
    };
    let channels = kept_samples(stored);

    // An alpha channel has no meaning downstream, where `channels()` is read
    // as 1 (gray) or 3 (RGB). Drop it at load, the way the TIFF reader drops
    // extra samples, so a GrayA/RGBA PNG never reaches code that would treat
    // its samples as extra pixels.
    let pixels = pixels.keep_leading_samples(stored, channels);

    let (width, height) = (info.width as usize, info.height as usize);
    Ok((pixels, plane_dimensions(width, height, channels)))
}

// Gray samples of 1, 2 or 4 bits scale to 8 bits, so the largest value maps
// to 255.
fn expand_gray(buf: &[u8], info: &png::OutputInfo) -> Result<PixelData> {
    let depth = info.bit_depth as u8;
    let scale = u8::MAX / (u8::MAX >> (8 - depth));
    let mut gray = reserved(expanded_len(info, 1)?)?;
    for row in rows(buf, info)? {
        gray.extend(packed_samples(row, depth, info.width).map(|v| v * scale));
    }
    Ok(PixelData::U8(gray))
}

// PNG 2nd edition (ISO/IEC 15948:2004) section 11.2.3 makes an index past
// the last palette entry an error. tRNS is not read, because the alpha it
// adds is dropped at load.
fn expand_palette(buf: &[u8], info: &png::OutputInfo, palette: Option<&[u8]>) -> Result<PixelData> {
    let palette = palette.ok_or_else(|| ImageError::Png("the image has no palette".into()))?;
    if !palette.len().is_multiple_of(3) {
        return Err(ImageError::Png(format!(
            "the PLTE chunk is {} bytes, not a whole number of 3-byte entries",
            palette.len()
        )));
    }
    let mut rgb = reserved(expanded_len(info, 3)?)?;
    for row in rows(buf, info)? {
        for index in packed_samples(row, info.bit_depth as u8, info.width) {
            let start = usize::from(index) * 3;
            let entry = palette.get(start..start + 3).ok_or_else(|| {
                ImageError::Png(format!(
                    "palette index {index} is past the last of {} entries",
                    palette.len() / 3
                ))
            })?;
            rgb.extend_from_slice(entry);
        }
    }
    Ok(PixelData::U8(rgb))
}

fn rows<'a>(buf: &'a [u8], info: &png::OutputInfo) -> Result<std::slice::ChunksExact<'a, u8>> {
    if info.line_size == 0 {
        return Err(ImageError::Png("the decoded rows are empty".into()));
    }
    Ok(buf.chunks_exact(info.line_size))
}

// Samples of `depth` bits (1, 2, 4 or 8), packed from the most significant
// bit. Each row starts on a byte boundary.
fn packed_samples(row: &[u8], depth: u8, width: u32) -> impl Iterator<Item = u8> + '_ {
    let mask = u8::MAX >> (8 - depth);
    let per_byte = 8 / depth;
    row.iter()
        .flat_map(move |&byte| {
            (0..per_byte)
                .rev()
                .map(move |k| (byte >> (k * depth)) & mask)
        })
        .take(width as usize)
}

/// Samples per pixel stored in the decoded buffer for a PNG color type.
fn stored_samples(color_type: png::ColorType) -> Result<usize> {
    use png::ColorType;

    match color_type {
        ColorType::Grayscale => Ok(1),
        ColorType::GrayscaleAlpha => Ok(2),
        ColorType::Rgb => Ok(3),
        ColorType::Rgba => Ok(4),
        ColorType::Indexed => Err(ImageError::UnsupportedFormat),
    }
}

/// Channels kept in memory: gray+alpha collapses to gray, RGBA to RGB.
fn kept_samples(stored: usize) -> usize {
    match stored {
        2 => 1,
        4 => 3,
        other => other,
    }
}

fn decode_png_samples(buf: Vec<u8>, bit_depth: png::BitDepth) -> Result<PixelData> {
    match bit_depth {
        png::BitDepth::Eight => Ok(PixelData::U8(buf)),
        png::BitDepth::Sixteen => {
            let (pairs, _) = buf.as_chunks::<2>();
            let mut samples = reserved(pairs.len())?;
            samples.extend(pairs.iter().map(|&b| u16::from_be_bytes(b)));
            Ok(PixelData::U16(samples))
        }
        _ => Err(ImageError::UnsupportedFormat),
    }
}

fn channels_to_png_color_type(channels: usize) -> Result<png::ColorType> {
    match channels {
        1 => Ok(png::ColorType::Grayscale),
        2 => Ok(png::ColorType::GrayscaleAlpha),
        3 => Ok(png::ColorType::Rgb),
        4 => Ok(png::ColorType::Rgba),
        _ => Err(ImageError::UnsupportedFormat),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_png() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".png").tempfile().unwrap()
    }

    #[test]
    fn roundtrip_u8_grayscale() {
        let tmp = tmp_png();
        let original: Vec<u8> = (0u8..64).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![8usize, 8]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn roundtrip_u8_rgb() {
        let tmp = tmp_png();
        let original: Vec<u8> = (0u8..48).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 4, 3]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn open_rgba_drops_alpha_to_three_channels() {
        let tmp = tmp_png();
        // Two RGBA pixels; the alpha samples must not survive the load.
        let rgba: Vec<u8> = vec![10, 20, 30, 200, 40, 50, 60, 100];
        let img = Image::new(PixelData::U8(rgba), vec![2usize, 1, 4]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![2usize, 1, 3]);
        assert!(restored.is_rgb());
        assert_eq!(
            restored.pixels.as_u8().unwrap(),
            &vec![10, 20, 30, 40, 50, 60]
        );
    }

    #[test]
    fn open_rgba_u16_drops_alpha_to_three_channels() {
        let tmp = tmp_png();
        let rgba: Vec<u16> = vec![1000, 2000, 3000, 65535, 4000, 5000, 6000, 65535];
        let img = Image::new(PixelData::U16(rgba), vec![2usize, 1, 4]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![2usize, 1, 3]);
        assert_eq!(
            restored.pixels.as_u16().unwrap(),
            &vec![1000, 2000, 3000, 4000, 5000, 6000]
        );
    }

    #[test]
    fn open_gray_alpha_drops_alpha_to_one_channel() {
        let tmp = tmp_png();
        let gray_alpha: Vec<u8> = vec![10, 255, 20, 128, 30, 0, 40, 255];
        let img = Image::new(PixelData::U8(gray_alpha), vec![4usize, 1, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 1]);
        assert_eq!(restored.channels(), 1);
        assert_eq!(restored.pixels.as_u8().unwrap(), &vec![10, 20, 30, 40]);
    }

    #[test]
    fn roundtrip_u16_grayscale() {
        let tmp = tmp_png();
        let original: Vec<u16> = (0u16..64).map(|i| i * 1000).collect();
        let img = Image::new(PixelData::U16(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open_png(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn save_rejects_unsupported_pixel_type() {
        let tmp = tmp_png();
        let img = Image::new(PixelData::F32(vec![0.0; 16]), vec![4usize, 4]);
        assert!(img.save(tmp.path()).is_err());
    }

    #[test]
    fn save_rejects_unsupported_channel_count() {
        let tmp = tmp_png();
        let img = Image::new(PixelData::U8(vec![0; 80]), vec![4usize, 4, 5]);
        assert!(img.save(tmp.path()).is_err());
    }

    #[test]
    fn open_returns_error_for_invalid_file() {
        let tmp = tempfile::Builder::new().suffix(".png").tempfile().unwrap();
        std::fs::write(tmp.path(), b"not a png").unwrap();
        assert!(Image::open_png(tmp.path()).is_err());
    }

    fn png_chunk(out: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
        let mut crc = flate2::Crc::new();
        crc.update(kind);
        crc.update(data);
        out.extend_from_slice(&(data.len() as u32).to_be_bytes());
        out.extend_from_slice(kind);
        out.extend_from_slice(data);
        out.extend_from_slice(&crc.sum().to_be_bytes());
    }

    // A PNG whose IHDR declares `width` x `height` RGBA16 and whose IDAT
    // holds 64 zero bytes.
    fn declared_rgba16(width: u32, height: u32) -> Vec<u8> {
        use std::io::Write;
        let mut out = b"\x89PNG\r\n\x1a\n".to_vec();
        let mut ihdr = Vec::new();
        ihdr.extend_from_slice(&width.to_be_bytes());
        ihdr.extend_from_slice(&height.to_be_bytes());
        ihdr.extend_from_slice(&[16, 6, 0, 0, 0]);
        png_chunk(&mut out, b"IHDR", &ihdr);
        let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        z.write_all(&[0u8; 64]).unwrap();
        png_chunk(&mut out, b"IDAT", &z.finish().unwrap());
        png_chunk(&mut out, b"IEND", &[]);
        out
    }

    fn encoded(
        (width, height): (u32, u32),
        color: png::ColorType,
        depth: png::BitDepth,
        palette: Option<(&[u8], &[u8])>,
        data: &[u8],
    ) -> tempfile::NamedTempFile {
        let tmp = tmp_png();
        let file = std::fs::File::create(tmp.path()).unwrap();
        let mut encoder = png::Encoder::new(file, width, height);
        encoder.set_color(color);
        encoder.set_depth(depth);
        if let Some((plte, trns)) = palette {
            encoder.set_palette(plte.to_vec());
            if !trns.is_empty() {
                encoder.set_trns(trns.to_vec());
            }
        }
        let mut writer = encoder.write_header().unwrap();
        writer.write_image_data(data).unwrap();
        writer.finish().unwrap();
        tmp
    }

    // A PNG assembled chunk by chunk, for layouts the encoder does not
    // write: interlaced rows and malformed palettes. `scanlines` are the
    // unfiltered rows; each gets filter type 0.
    fn assembled(ihdr: &[u8; 13], plte: &[u8], scanlines: &[Vec<u8>]) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut out = b"\x89PNG\r\n\x1a\n".to_vec();
        png_chunk(&mut out, b"IHDR", ihdr);
        if !plte.is_empty() {
            png_chunk(&mut out, b"PLTE", plte);
        }
        let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        for line in scanlines {
            z.write_all(&[0]).unwrap();
            z.write_all(line).unwrap();
        }
        png_chunk(&mut out, b"IDAT", &z.finish().unwrap());
        png_chunk(&mut out, b"IEND", &[]);
        let tmp = tmp_png();
        std::fs::write(tmp.path(), out).unwrap();
        tmp
    }

    fn ihdr((width, height): (u32, u32), depth: u8, color: u8, interlace: u8) -> [u8; 13] {
        let mut ihdr = [0u8; 13];
        ihdr[..4].copy_from_slice(&width.to_be_bytes());
        ihdr[4..8].copy_from_slice(&height.to_be_bytes());
        ihdr[8..].copy_from_slice(&[depth, color, 0, 0, interlace]);
        ihdr
    }

    // Rows of `pixel(x, y)` packed at `depth` bits, in Adam7 pass order.
    fn adam7_rows(
        (width, height): (u32, u32),
        depth: u8,
        pixel: impl Fn(u32, u32) -> u8,
    ) -> Vec<Vec<u8>> {
        let passes = [(0, 0, 8, 8), (4, 0, 8, 8), (0, 4, 4, 8), (2, 0, 4, 4)];
        let passes = passes
            .into_iter()
            .chain([(0, 2, 2, 4), (1, 0, 2, 2), (0, 1, 1, 2)]);
        let mut rows = Vec::new();
        for (x0, y0, dx, dy) in passes {
            for y in (y0..height).step_by(dy) {
                let xs: Vec<u32> = (x0..width).step_by(dx).collect();
                let mut row = vec![0u8; (xs.len() * depth as usize).div_ceil(8)];
                for (i, &x) in xs.iter().enumerate() {
                    let bit = i * depth as usize;
                    row[bit / 8] |= pixel(x, y) << (8 - depth as usize - bit % 8);
                }
                rows.push(row);
            }
        }
        rows
    }

    #[test]
    fn open_rejects_a_size_the_file_cannot_inflate_to() {
        let tmp = tmp_png();
        std::fs::write(tmp.path(), declared_rgba16(4_000_000, 4_000_000)).unwrap();
        let error = Image::open_png(tmp.path()).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some(
                "PNG error: 4000000x4000000 pixels inflate to 128000004000000 bytes, \
                 more than a 78-byte file can hold"
            )
        );
    }

    #[test]
    fn open_expands_one_and_two_bit_gray_to_eight_bits() {
        let rows = [0b1010_1010, 0b1100_0000, 0b0101_0101, 0b0000_0000];
        let tmp = encoded(
            (10, 2),
            png::ColorType::Grayscale,
            png::BitDepth::One,
            None,
            &rows,
        );
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.dimensions, vec![10usize, 2]);
        let (on, off) = (255u8, 0u8);
        let expected = vec![
            on, off, on, off, on, off, on, off, on, on, //
            off, on, off, on, off, on, off, on, off, off,
        ];
        assert_eq!(image.pixels.as_u8(), Some(&expected));

        let tmp = encoded(
            (4, 1),
            png::ColorType::Grayscale,
            png::BitDepth::Two,
            None,
            &[0b0001_1011],
        );
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.pixels.as_u8(), Some(&vec![0, 85, 170, 255]));
    }

    #[test]
    fn open_expands_palette_indices_to_rgb() {
        let plte = [10, 20, 30, 40, 50, 60, 70, 80, 90];
        let tmp = encoded(
            (3, 1),
            png::ColorType::Indexed,
            png::BitDepth::Four,
            Some((&plte, &[])),
            &[0x20, 0x10],
        );
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.dimensions, vec![3usize, 1, 3]);
        assert_eq!(
            image.pixels.as_u8(),
            Some(&vec![70, 80, 90, 10, 20, 30, 40, 50, 60])
        );
    }

    // tRNS makes the expanded image RGBA, and the alpha is dropped at load.
    #[test]
    fn open_expands_a_transparent_palette_to_rgb() {
        let plte = [1, 2, 3, 4, 5, 6];
        let tmp = encoded(
            (2, 1),
            png::ColorType::Indexed,
            png::BitDepth::Eight,
            Some((&plte, &[0, 128])),
            &[1, 0],
        );
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.dimensions, vec![2usize, 1, 3]);
        assert_eq!(image.pixels.as_u8(), Some(&vec![4, 5, 6, 1, 2, 3]));
    }

    #[test]
    fn open_rejects_an_index_past_the_palette() {
        let plte = [1, 2, 3, 4, 5, 6];
        let tmp = encoded(
            (2, 1),
            png::ColorType::Indexed,
            png::BitDepth::Eight,
            Some((&plte, &[])),
            &[1, 2],
        );
        let error = Image::open_png(tmp.path()).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("PNG error: palette index 2 is past the last of 2 entries")
        );
    }

    // 2000x2000 one-bit zeros compress about 500 to 1, below the DEFLATE
    // limit, so the size check accepts the file.
    #[test]
    fn open_accepts_a_highly_compressed_image() {
        let rows = vec![0u8; 250 * 2000];
        let tmp = encoded(
            (2000, 2000),
            png::ColorType::Grayscale,
            png::BitDepth::One,
            None,
            &rows,
        );
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.dimensions, vec![2000usize, 2000]);
        assert_eq!(image.pixels.as_u8(), Some(&vec![0u8; 2000 * 2000]));
    }

    #[test]
    fn open_scales_four_bit_gray_by_seventeen() {
        let tmp = encoded(
            (4, 1),
            png::ColorType::Grayscale,
            png::BitDepth::Four,
            None,
            &[0x01, 0xF8],
        );
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.pixels.as_u8(), Some(&vec![0, 17, 255, 136]));
    }

    // 5x5 leaves every one of the seven Adam7 passes non-empty.
    #[test]
    fn open_reads_an_interlaced_two_bit_image() {
        let pixel = |x: u32, y: u32| ((x + 2 * y) % 4) as u8;
        let rows = adam7_rows((5, 5), 2, pixel);
        let tmp = assembled(&ihdr((5, 5), 2, 0, 1), &[], &rows);
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.dimensions, vec![5usize, 5]);
        let expected: Vec<u8> = (0..5)
            .flat_map(|y| (0..5).map(move |x| pixel(x, y) * 85))
            .collect();
        assert_eq!(image.pixels.as_u8(), Some(&expected));
    }

    #[test]
    fn open_rejects_a_palette_of_partial_entries() {
        let tmp = assembled(&ihdr((1, 1), 8, 3, 0), &[1, 2, 3, 4], &[vec![0]]);
        let error = Image::open_png(tmp.path()).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("PNG error: the PLTE chunk is 4 bytes, not a whole number of 3-byte entries")
        );
    }

    // A 1-bit image can index only 2 entries. The extra entries are unused.
    #[test]
    fn open_accepts_more_palette_entries_than_the_depth_indexes() {
        let plte = [1, 2, 3, 4, 5, 6, 7, 8, 9];
        let tmp = assembled(&ihdr((2, 1), 1, 3, 0), &plte, &[vec![0b0100_0000]]);
        let image = Image::open_png(tmp.path()).unwrap();
        assert_eq!(image.pixels.as_u8(), Some(&vec![1, 2, 3, 4, 5, 6]));
    }

    // The expanded length is checked and reserved before any row is read.
    // An empty buffer with a forged size would otherwise expand to nothing.
    #[test]
    fn expansion_rejects_an_output_that_cannot_be_allocated() {
        let info = |color_type, bit_depth| png::OutputInfo {
            width: u32::MAX,
            height: u32::MAX,
            color_type,
            bit_depth,
            line_size: 1,
        };
        let gray = info(png::ColorType::Grayscale, png::BitDepth::One);
        let error = expand_gray(&[], &gray).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("PNG error: cannot allocate 18446744065119617025 samples for the image")
        );

        let indexed = info(png::ColorType::Indexed, png::BitDepth::One);
        let error = expand_palette(&[], &indexed, Some(&[0; 6]))
            .err()
            .map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("PNG error: 4294967295x4294967295 pixels of 3 samples do not fit in memory")
        );
    }

    #[test]
    fn decode_and_encode_errors_name_png() {
        let tmp = tmp_png();
        std::fs::write(tmp.path(), b"not a png").unwrap();
        let error = Image::open_png(tmp.path()).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("PNG error: decode: Invalid PNG signature.")
        );

        let image = Image::new(PixelData::U8(vec![0; 3]), vec![2usize, 2]);
        let error = image.save(tmp.path()).err().map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some("PNG error: encode: wrong data size, expected 4 got 3")
        );
    }

    #[test]
    fn save_rgba_grayscale_alpha_and_rgba_channel_types() {
        // 2-channel (grayscale-alpha)
        let tmp = tmp_png();
        let img = Image::new(PixelData::U8(vec![0; 32]), vec![4usize, 4, 2]);
        img.save(tmp.path()).unwrap();

        // 4-channel (rgba)
        let tmp = tmp_png();
        let img = Image::new(PixelData::U8(vec![0; 64]), vec![4usize, 4, 4]);
        img.save(tmp.path()).unwrap();
    }
}
