// Tile-compressed images written by celestial and read back through
// `Image::open`.
use super::Image;
use crate::core::errors::ImageError;
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::DataArray;
use crate::fits::errors::FitsError;
use crate::fits::header::keywords::Keyword;
use crate::fits::header::parser::Header;
use crate::fits::io::reader::FitsFile;
use crate::fits::io::test_support::FitsWriter;
use crate::formats::astro_image::AstroImage;
use crate::formats::pixel_data::PixelData;
use std::path::Path;

const WRITTEN: [CompressionAlgorithm; 3] = [
    CompressionAlgorithm::Rice,
    CompressionAlgorithm::Gzip,
    CompressionAlgorithm::Gzip2,
];

// Image dimensions and the requested tile size.
type Shape = (&'static [usize], Option<(usize, usize)>);

fn tmp_fits() -> tempfile::NamedTempFile {
    tempfile::Builder::new().suffix(".fits").tempfile().unwrap()
}

fn write<T: DataArray>(
    path: &Path,
    data: &[T],
    dims: &[usize],
    tile: Option<(usize, usize)>,
    algorithm: CompressionAlgorithm,
) -> crate::fits::errors::Result<()> {
    FitsWriter::create(path)
        .unwrap()
        .write_compressed_image(data, dims, tile, algorithm, &[])
}

fn round_trip<T: DataArray>(
    data: &[T],
    dims: &[usize],
    tile: Option<(usize, usize)>,
    algorithm: CompressionAlgorithm,
) -> Image {
    let tmp = tmp_fits();
    write(tmp.path(), data, dims, tile, algorithm)
        .unwrap_or_else(|e| panic!("{algorithm:?} {dims:?}: {e}"));
    Image::open(tmp.path()).unwrap_or_else(|e| panic!("{algorithm:?} {dims:?}: {e}"))
}

fn table_header(path: &Path) -> Header {
    FitsFile::open(path).unwrap().get_header(1).unwrap()
}

// A fixed pseudo-random sequence, so failures repeat.
fn noise(len: usize, seed: u64) -> impl Iterator<Item = u64> {
    let mut state = seed;
    (0..len).map(move |_| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state ^ (state >> 29)
    })
}

// The peak ends the first 32-pixel Rice block, so the block holds one large
// difference among zeros.
fn spike<T: Copy>(background: T, peak: T) -> Vec<T> {
    let mut data = vec![background; 64];
    data[31] = peak;
    data
}

#[test]
fn rice_writes_blocks_with_large_steps() {
    let i32s = spike(0i32, 1000);
    let img = round_trip(&i32s, &[8, 8], None, CompressionAlgorithm::Rice);
    assert_eq!(img.pixels.as_i32(), Some(&i32s));

    let i16s = spike(0i16, 1000);
    let img = round_trip(&i16s, &[8, 8], Some((8, 8)), CompressionAlgorithm::Rice);
    assert_eq!(img.pixels.as_i16(), Some(&i16s));

    let step: Vec<u8> = (0..64).map(|i| if i < 37 { 0 } else { 200 }).collect();
    let img = round_trip(&step, &[8, 8], Some((8, 8)), CompressionAlgorithm::Rice);
    assert_eq!(img.pixels.as_u8(), Some(&step));
}

#[test]
fn rice_writes_noise_with_a_cosmic_ray() {
    let mut data: Vec<i16> = noise(7000, 1)
        .map(|n| 1000 + (n % 41) as i16 - 20)
        .collect();
    let img = round_trip(
        &data,
        &[70, 100],
        Some((70, 100)),
        CompressionAlgorithm::Rice,
    );
    assert_eq!(img.pixels.as_i16(), Some(&data));

    data[3456] = 30000;
    let img = round_trip(
        &data,
        &[70, 100],
        Some((70, 100)),
        CompressionAlgorithm::Rice,
    );
    assert_eq!(img.pixels.as_i16(), Some(&data));
}

#[test]
fn rice_round_trips_u8_and_i16_with_bytepix() {
    let u8s: Vec<u8> = noise(23 * 17, 2).map(|n| n as u8).collect();
    let tmp = tmp_fits();
    write(
        tmp.path(),
        &u8s,
        &[23, 17],
        Some((7, 5)),
        CompressionAlgorithm::Rice,
    )
    .unwrap();
    let header = table_header(tmp.path());
    assert_eq!(header.get("ZNAME1").as_str(), Some("BLOCKSIZE"));
    assert_eq!(header.get("ZVAL1").as_i64(), Some(32));
    assert_eq!(header.get("ZNAME2").as_str(), Some("BYTEPIX"));
    assert_eq!(header.get("ZVAL2").as_i64(), Some(1));
    assert_eq!(Image::open(tmp.path()).unwrap().pixels.as_u8(), Some(&u8s));

    let i16s: Vec<i16> = noise(23 * 17, 3).map(|n| n as i16).collect();
    let tmp = tmp_fits();
    write(
        tmp.path(),
        &i16s,
        &[23, 17],
        Some((7, 5)),
        CompressionAlgorithm::Rice,
    )
    .unwrap();
    assert_eq!(table_header(tmp.path()).get("ZVAL2").as_i64(), Some(2));
    assert_eq!(
        Image::open(tmp.path()).unwrap().pixels.as_i16(),
        Some(&i16s)
    );
}

#[test]
fn cube_round_trips_every_plane() {
    let planes: Vec<i16> = (0..5 * 4 * 2).map(|i| i * 37 - 500).collect();
    let img = round_trip(
        &planes,
        &[5, 4, 2],
        Some((5, 4)),
        CompressionAlgorithm::Rice,
    );
    assert_eq!(img.dimensions, vec![5, 4, 2]);
    assert_eq!(img.pixels.as_i16(), Some(&planes));

    // Three planes open as an interleaved RGB image.
    let rgb: Vec<u8> = (0..4 * 3 * 3).map(|i| i as u8 * 7).collect();
    let tmp = tmp_fits();
    write(
        tmp.path(),
        &rgb,
        &[4, 3, 3],
        Some((4, 3)),
        CompressionAlgorithm::Rice,
    )
    .unwrap();
    assert_eq!(table_header(tmp.path()).get("ZTILE3").as_i64(), Some(1));
    let mut img = Image::open(tmp.path()).unwrap();
    img.interleaved_to_planar();
    assert_eq!(img.dimensions, vec![4, 3, 3]);
    assert_eq!(img.pixels.as_u8(), Some(&rgb));
}

// The first tile of the table, read through its heap descriptor.
fn first_tile(path: &Path) -> Vec<u8> {
    let mut fits = FitsFile::open(path).unwrap();
    let header = fits.get_header(1).unwrap();
    let data = fits.read_data_unit(1).unwrap();
    let rows = header.get("NAXIS2").as_i64().unwrap() as usize;
    let row_bytes = header.get("NAXIS1").as_i64().unwrap() as usize;
    let count = u32::from_be_bytes(data[0..4].try_into().unwrap()) as usize;
    let offset = u32::from_be_bytes(data[4..8].try_into().unwrap()) as usize;
    let heap = rows * row_bytes + offset;
    data[heap..heap + count].to_vec()
}

fn inflate(stream: &[u8]) -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    flate2::read::GzDecoder::new(stream)
        .read_to_end(&mut out)
        .unwrap();
    out
}

#[test]
fn gzip1_writes_a_gzip_stream() {
    let data: Vec<i16> = noise(23 * 17, 4).map(|n| (n % 3000) as i16).collect();
    let tmp = tmp_fits();
    write(
        tmp.path(),
        &data,
        &[23, 17],
        Some((7, 5)),
        CompressionAlgorithm::Gzip,
    )
    .unwrap();
    assert_eq!(
        table_header(tmp.path()).get("ZCMPTYPE").as_str(),
        Some("GZIP_1")
    );
    assert_eq!(first_tile(tmp.path())[..2], [0x1f, 0x8b]);
    assert_eq!(
        Image::open(tmp.path()).unwrap().pixels.as_i16(),
        Some(&data)
    );
}

// GZIP_2 stores the high bytes of a tile's values, then the low bytes.
#[test]
fn gzip2_writes_shuffled_bytes() {
    let data: Vec<i16> = vec![0x0102, 0x0304, 0x0506];
    let tmp = tmp_fits();
    write(tmp.path(), &data, &[3], None, CompressionAlgorithm::Gzip2).unwrap();
    assert_eq!(
        table_header(tmp.path()).get("ZCMPTYPE").as_str(),
        Some("GZIP_2")
    );
    assert_eq!(inflate(&first_tile(tmp.path())), vec![1, 3, 5, 2, 4, 6]);
    assert_eq!(
        Image::open(tmp.path()).unwrap().pixels.as_i16(),
        Some(&data)
    );
}

fn special_f32s() -> Vec<f32> {
    let mut data: Vec<f32> = noise(23 * 17, 5)
        .map(|n| (n % 100000) as f32 * 1.5e-3 - 7.0)
        .collect();
    data[..7].copy_from_slice(&[
        f32::NAN,
        f32::from_bits(0x7fc0_1234),
        -0.0,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
    ]);
    data
}

fn bits32(pixels: &PixelData) -> Vec<u32> {
    pixels
        .as_f32()
        .unwrap()
        .iter()
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn gzip_round_trips_floats_bit_exactly() {
    let data = special_f32s();
    let expected: Vec<u32> = data.iter().map(|v| v.to_bits()).collect();
    for algorithm in [CompressionAlgorithm::Gzip, CompressionAlgorithm::Gzip2] {
        let img = round_trip(&data, &[23, 17], Some((7, 5)), algorithm);
        assert_eq!(bits32(&img.pixels), expected, "{algorithm:?}");
    }

    let f64s: Vec<f64> = data.iter().map(|&v| f64::from(v) / 3.0).collect();
    let img = round_trip(&f64s, &[23, 17], Some((7, 5)), CompressionAlgorithm::Gzip);
    let bits: Vec<u64> = img
        .pixels
        .as_f64()
        .unwrap()
        .iter()
        .map(|v| v.to_bits())
        .collect();
    assert_eq!(bits, f64s.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
}

fn assert_unsupported(result: crate::fits::errors::Result<()>, text: &str) {
    match result {
        Err(FitsError::UnsupportedCompression(m)) => assert!(m.contains(text), "{m}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn rice_rejects_floats_and_64_bit_integers_and_names_gzip() {
    let tmp = tmp_fits();
    let rice = CompressionAlgorithm::Rice;
    assert_unsupported(
        write(tmp.path(), &[1.0f32; 16], &[4, 4], None, rice),
        "GZIP",
    );
    assert_unsupported(
        write(tmp.path(), &[1.0f64; 16], &[4, 4], None, rice),
        "GZIP",
    );
    assert_unsupported(write(tmp.path(), &[1i64; 16], &[4, 4], None, rice), "GZIP");
}

#[test]
fn hcompress_and_plio_are_not_written() {
    let tmp = tmp_fits();
    let cases = [
        (CompressionAlgorithm::HCompress, "HCOMPRESS_1"),
        (CompressionAlgorithm::Plio, "PLIO_1"),
    ];
    for (algorithm, name) in cases {
        let result = write(tmp.path(), &[1i16; 16], &[4, 4], None, algorithm);
        let expected = format!("{name} tiles cannot be written; use RICE_1, GZIP_1 or GZIP_2");
        match result {
            Err(FitsError::UnsupportedCompression(m)) => assert_eq!(m, expected),
            other => panic!("{other:?}"),
        }
    }
}

#[test]
fn failed_compressed_write_creates_no_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.fits");
    let rice = CompressionAlgorithm::Rice;
    assert!(write(&path, &[1.0f64; 16], &[4, 4], None, rice).is_err());
    assert!(write(&path, &[1i16; 16], &[4, 4], Some((0, 4)), rice).is_err());
    assert!(!path.exists());
}

#[test]
fn write_rejects_empty_images_and_tiles() {
    let tmp = tmp_fits();
    let gzip = CompressionAlgorithm::Gzip;
    let cases: [(Shape, &str); 4] = [
        ((&[], None), "a compressed image needs at least one axis"),
        (
            (&[0, 4], None),
            "a compressed image needs at least one pixel, not a 0x4 image",
        ),
        ((&[4, 4], Some((0, 1))), "tile size 0x1 has an empty axis"),
        ((&[4, 4], Some((1, 0))), "tile size 1x0 has an empty axis"),
    ];
    for ((dims, tile), expected) in cases {
        let len = if dims.is_empty() {
            0
        } else {
            dims.iter().product()
        };
        match write(tmp.path(), &vec![1u8; len], dims, tile, gzip) {
            Err(FitsError::InvalidFormat(m)) => assert_eq!(m, expected),
            other => panic!("{dims:?} {tile:?}: {other:?}"),
        }
    }
}

#[test]
fn tiles_are_rows_by_default_and_never_exceed_the_image() {
    let data: Vec<i32> = (0..5 * 4 * 2).collect();
    let tmp = tmp_fits();
    write(
        tmp.path(),
        &data,
        &[5, 4, 2],
        None,
        CompressionAlgorithm::Rice,
    )
    .unwrap();
    let header = table_header(tmp.path());
    let tiles: Vec<_> = (1..=3)
        .map(|i| header.get(&format!("ZTILE{i}")).as_i64())
        .collect();
    assert_eq!(tiles, vec![Some(5), Some(1), Some(1)]);
    assert_eq!(header.get("NAXIS2").as_i64(), Some(8));

    let rice = CompressionAlgorithm::Rice;
    write(tmp.path(), &data, &[5, 4, 2], Some((64, 64)), rice).unwrap();
    let header = table_header(tmp.path());
    assert_eq!(header.get("ZTILE1").as_i64(), Some(5));
    assert_eq!(header.get("ZTILE2").as_i64(), Some(4));
    assert_eq!(header.get("NAXIS2").as_i64(), Some(2));
    let img = Image::open(tmp.path()).unwrap();
    assert_eq!(img.pixels.as_i32(), Some(&data));
}

#[test]
fn header_keeps_image_keywords_and_drops_table_ones() {
    let keywords = [
        Keyword::string("OBJECT", "M31"),
        Keyword::integer("EXPTIME", 60),
        Keyword::string("ZCMPTYPE", "PLIO_1"),
        Keyword::string("TFORM1", "1E"),
        Keyword::integer("NAXIS1", 99),
        Keyword::integer("THEAP", 99),
    ];
    let tmp = tmp_fits();
    FitsWriter::create(tmp.path())
        .unwrap()
        .write_compressed_image(
            &[7i16; 12],
            &[4, 3],
            None,
            CompressionAlgorithm::Rice,
            &keywords,
        )
        .unwrap();
    let header = table_header(tmp.path());
    assert_eq!(header.get("ZCMPTYPE").as_str(), Some("RICE_1"));
    assert!(header.get("TFORM1").as_str().unwrap().starts_with("1PB("));
    assert_eq!(header.get("NAXIS1").as_i64(), Some(8));
    assert!(header.get_keyword("THEAP").is_none());
    let img = Image::open(tmp.path()).unwrap();
    assert_eq!(img.get("OBJECT").as_str(), Some("M31"));
    assert_eq!(img.get("EXPTIME").as_i64(), Some(60));
    assert_eq!(img.pixels.as_i16(), Some(&vec![7; 12]));
}

#[test]
fn astro_image_writes_uncompressed_by_default() {
    let data: Vec<i16> = (0..16).collect();
    let tmp = tmp_fits();
    AstroImage::new(&data, [4, 4])
        .write_fits(tmp.path())
        .unwrap();
    let mut fits = FitsFile::open(tmp.path()).unwrap();
    assert_eq!(fits.num_hdus(), 1);
    assert_eq!(fits.get_header(0).unwrap().get("NAXIS").as_i64(), Some(2));
}

// Values that cover each type's range: noise, the extremes, a flat run and
// a lone spike.
fn pixels_of(kind: &str, len: usize, seed: u64) -> PixelData {
    let mut n: Vec<u64> = noise(len, seed).collect();
    n[..4].copy_from_slice(&[0, u64::MAX, 1 << 63, (1 << 63) - 1]);
    for v in n.iter_mut().skip(40).take(40) {
        *v = 5;
    }
    n[len / 2] = u64::MAX / 3;
    match kind {
        "U8" => PixelData::U8(n.iter().map(|&v| v as u8).collect()),
        "I8" => PixelData::I8(n.iter().map(|&v| v as i8).collect()),
        "U16" => PixelData::U16(n.iter().map(|&v| v as u16).collect()),
        "I16" => PixelData::I16(n.iter().map(|&v| v as i16).collect()),
        "U32" => PixelData::U32(n.iter().map(|&v| v as u32).collect()),
        "I32" => PixelData::I32(n.iter().map(|&v| v as i32).collect()),
        "U64" => PixelData::U64(n),
        "I64" => PixelData::I64(n.iter().map(|&v| v as i64).collect()),
        "F32" => PixelData::F32(n.iter().map(|&v| f32::from_bits(v as u32)).collect()),
        _ => PixelData::F64(n.iter().map(|&v| f64::from_bits(v)).collect()),
    }
}

const KINDS: [&str; 10] = [
    "U8", "I8", "U16", "I16", "U32", "I32", "U64", "I64", "F32", "F64",
];

// U64 is refused for every algorithm, and RICE_1 also refuses I64 and
// floats.
fn written(kind: &str, algorithm: CompressionAlgorithm) -> bool {
    match kind {
        "U64" => false,
        "I64" | "F32" | "F64" => algorithm != CompressionAlgorithm::Rice,
        _ => true,
    }
}

fn bits(pixels: &PixelData) -> Vec<u64> {
    match pixels {
        PixelData::U8(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::I8(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::U16(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::I16(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::U32(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::I32(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::U64(v) => v.clone(),
        PixelData::I64(v) => v.iter().map(|&x| x as u64).collect(),
        PixelData::F32(v) => v.iter().map(|x| x.to_bits() as u64).collect(),
        PixelData::F64(v) => v.iter().map(|x| x.to_bits()).collect(),
    }
}

fn assert_refused(result: crate::core::errors::Result<()>, case: &str) {
    assert!(
        matches!(
            result,
            Err(ImageError::Fits(FitsError::UnsupportedCompression(_)))
        ),
        "{case}: {result:?}"
    );
}

// Every pixel type and algorithm, on uneven tiles, default rows, an RGB
// cube and a four-plane cube.
#[test]
fn image_round_trips_every_type_and_algorithm() {
    let shapes: [Shape; 4] = [
        (&[23, 17], Some((7, 5))),
        (&[23, 17], None),
        (&[11, 7, 3], Some((5, 3))),
        (&[11, 7, 4], Some((4, 7))),
    ];
    for (k, kind) in KINDS.iter().enumerate() {
        for algorithm in WRITTEN {
            for (s, (dims, tile)) in shapes.iter().enumerate() {
                let len = dims.iter().product();
                let image = Image::new(pixels_of(kind, len, (k * 10 + s) as u64), *dims);
                let dir = tempfile::tempdir().unwrap();
                let path = dir.path().join("image.fits");
                let result = image.save_fits_compressed(&path, algorithm, *tile);
                let case = format!("{kind} {algorithm:?} {dims:?} {tile:?}");
                if !written(kind, algorithm) {
                    assert_refused(result, &case);
                    assert!(!path.exists(), "{case}: file left behind");
                    continue;
                }
                result.unwrap_or_else(|e| panic!("{case}: {e}"));
                let back = Image::open(&path).unwrap_or_else(|e| panic!("{case}: {e}"));
                assert_eq!(back.dimensions, image.dimensions, "{case}");
                assert_eq!(
                    std::mem::discriminant(&back.pixels),
                    std::mem::discriminant(&image.pixels),
                    "{case}"
                );
                assert_eq!(bits(&back.pixels), bits(&image.pixels), "{case}");
            }
        }
    }
}

#[test]
fn offset_types_store_the_standard_bzero() {
    let cases = [
        ("I8", 8, -128.0),
        ("U16", 16, 32768.0),
        ("U32", 32, 2147483648.0),
    ];
    for (kind, zbitpix, bzero) in cases {
        let image = Image::new(pixels_of(kind, 12, 9), [4, 3]);
        let tmp = tmp_fits();
        image
            .save_fits_compressed(tmp.path(), CompressionAlgorithm::Rice, None)
            .unwrap();
        let header = table_header(tmp.path());
        assert_eq!(header.get("ZBITPIX").as_i64(), Some(zbitpix), "{kind}");
        assert_eq!(header.get("BZERO").as_f64(), Some(bzero), "{kind}");
        assert_eq!(header.get("BSCALE").as_f64(), Some(1.0), "{kind}");
    }
}

// Other FITS readers cannot decode 64-bit unsigned tiles, so the save is
// refused before any file is created.
#[test]
fn u64_compressed_save_is_refused_and_leaves_no_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("u64.fits");
    let image = Image::new(PixelData::U64(vec![0, 1, u64::MAX, 1 << 63]), [2, 2]);
    for algorithm in WRITTEN {
        match image.save_fits_compressed(&path, algorithm, None) {
            Err(ImageError::Fits(FitsError::UnsupportedCompression(m))) => {
                assert!(m.contains("uncompressed"), "{m}")
            }
            other => panic!("{algorithm:?}: {other:?}"),
        }
        assert!(!path.exists(), "{algorithm:?}");
    }
}
