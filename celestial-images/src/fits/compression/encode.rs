use super::restore::describes_table_or_compression;
use super::rice::DEFAULT_BLOCK_SIZE;
use super::tiles::TileGrid;
use super::{gzip, rice, CompressionAlgorithm};
use crate::core::types::BitPix;
use crate::fits::block::pad_to_block;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::keywords::Keyword;
use crate::fits::header::parser::Header;
use crate::fits::io::writer::is_checksum_keyword;
use crate::fits::util::image::shape_text;

// A compressed image as a BINTABLE extension (FITS 4.0 §10.1): the header
// without END, and the data unit padded to whole blocks.
pub(crate) struct EncodedImage {
    pub(crate) header: Header,
    pub(crate) data: Vec<u8>,
}

#[derive(Clone, Copy)]
enum Codec {
    Rice { bytepix: usize },
    Gzip1,
    Gzip2,
}

// The tile shape and codec are fixed before any tile is compressed, and every
// tile is compressed before the caller writes a byte, so an error leaves no
// partial output.
pub(crate) fn encode_image(
    pixels: &[u8],
    bitpix: BitPix,
    axes: &[usize],
    tile_size: Option<(usize, usize)>,
    algorithm: CompressionAlgorithm,
    keywords: &[Keyword],
) -> Result<EncodedImage> {
    let codec = codec(algorithm, bitpix)?;
    let tile = tile_shape(axes, tile_size)?;
    let grid = TileGrid::new(axes.to_vec(), tile.clone())?;
    let tiles = compress_tiles(&grid, codec, bitpix.bytes_per_pixel(), pixels)?;
    let table = Table::new(&tiles);
    let mut header = table.header();
    add_compression_keywords(&mut header, algorithm, codec, bitpix, axes, &tile);
    // FITS 4.0 §10.1.2 recommends copying the other keywords of the image
    // header, in order.
    for keyword in keywords.iter().filter(|k| is_image_keyword(k)) {
        header.add_keyword(keyword.clone());
    }
    Ok(EncodedImage {
        header,
        data: table.data(&tiles),
    })
}

fn add_compression_keywords(
    header: &mut Header,
    algorithm: CompressionAlgorithm,
    codec: Codec,
    bitpix: BitPix,
    axes: &[usize],
    tile: &[usize],
) {
    header.add_keyword(Keyword::logical("ZIMAGE", true));
    header.add_keyword(Keyword::string("ZCMPTYPE", algorithm.fits_name()));
    header.add_keyword(Keyword::integer("ZBITPIX", i64::from(bitpix.value())));
    header.add_keyword(Keyword::integer("ZNAXIS", count(axes.len())));
    for (i, &axis) in axes.iter().enumerate() {
        header.add_keyword(Keyword::integer(format!("ZNAXIS{}", i + 1), count(axis)));
    }
    for (i, &size) in tile.iter().enumerate() {
        header.add_keyword(Keyword::integer(format!("ZTILE{}", i + 1), count(size)));
    }
    if let Codec::Rice { bytepix } = codec {
        header.add_keyword(Keyword::string("ZNAME1", "BLOCKSIZE"));
        header.add_keyword(Keyword::integer("ZVAL1", count(DEFAULT_BLOCK_SIZE)));
        header.add_keyword(Keyword::string("ZNAME2", "BYTEPIX"));
        header.add_keyword(Keyword::integer("ZVAL2", count(bytepix)));
    }
}

// Structure keywords come from the data, and a checksum of the old HDU no
// longer holds. Both are dropped, as the uncompressed writer drops them.
fn is_image_keyword(keyword: &Keyword) -> bool {
    !keyword.is_mandatory()
        && !is_checksum_keyword(&keyword.name)
        && !describes_table_or_compression(&keyword.name)
}

// Sizes of buffers in memory are at most isize::MAX, which fits an i64.
fn count(n: usize) -> i64 {
    n as i64
}

fn codec(algorithm: CompressionAlgorithm, bitpix: BitPix) -> Result<Codec> {
    let unsupported = |reason: &str| {
        FitsError::UnsupportedCompression(format!(
            "{} tiles cannot hold {reason}; use GZIP_1 or GZIP_2",
            algorithm.fits_name()
        ))
    };
    match (algorithm, bitpix) {
        (CompressionAlgorithm::Rice, BitPix::U8 | BitPix::I16 | BitPix::I32) => Ok(Codec::Rice {
            bytepix: bitpix.bytes_per_pixel(),
        }),
        (CompressionAlgorithm::Rice, BitPix::I64) => Err(unsupported("64-bit integer pixels")),
        (CompressionAlgorithm::Rice, BitPix::F32 | BitPix::F64) => Err(unsupported(
            "floating-point pixels without lossy quantization, which is not written",
        )),
        (CompressionAlgorithm::Gzip, _) => Ok(Codec::Gzip1),
        (CompressionAlgorithm::Gzip2, _) => Ok(Codec::Gzip2),
        (CompressionAlgorithm::HCompress | CompressionAlgorithm::Plio, _) => {
            let choices = if matches!(bitpix, BitPix::U8 | BitPix::I16 | BitPix::I32) {
                "RICE_1, GZIP_1 or GZIP_2"
            } else {
                "GZIP_1 or GZIP_2"
            };
            Err(FitsError::UnsupportedCompression(format!(
                "{} tiles cannot be written; use {choices}",
                algorithm.fits_name()
            )))
        }
    }
}

// Rows of axis 1 by default (FITS 4.0 §10.1.2). A requested size applies to
// the first two axes, and a tile never extends past the image.
fn tile_shape(axes: &[usize], tile_size: Option<(usize, usize)>) -> Result<Vec<usize>> {
    check_image_axes(axes)?;
    let (first, second) = match tile_size {
        Some((first, second)) if first == 0 || second == 0 => {
            return Err(FitsError::InvalidFormat(format!(
                "tile size {first}x{second} has an empty axis"
            )))
        }
        Some(size) => size,
        None => (axes[0], 1),
    };
    Ok(axes
        .iter()
        .enumerate()
        .map(|(i, &axis)| match i {
            0 => first.min(axis),
            1 => second.min(axis),
            _ => 1,
        })
        .collect())
}

fn check_image_axes(axes: &[usize]) -> Result<()> {
    if axes.is_empty() {
        return Err(FitsError::InvalidFormat(
            "a compressed image needs at least one axis".to_string(),
        ));
    }
    if axes.contains(&0) {
        return Err(FitsError::InvalidFormat(format!(
            "a compressed image needs at least one pixel, not a {} image",
            shape_text(axes)
        )));
    }
    Ok(())
}

fn compress_tiles(
    grid: &TileGrid,
    codec: Codec,
    width: usize,
    pixels: &[u8],
) -> Result<Vec<Vec<u8>>> {
    (0..grid.tile_count())
        .map(|i| {
            let tile = grid.extract(&grid.region(i), width, pixels)?;
            compress(codec, &tile, width)
        })
        .collect()
}

fn compress(codec: Codec, tile: &[u8], width: usize) -> Result<Vec<u8>> {
    match codec {
        Codec::Rice { bytepix } => rice::compress(tile, bytepix),
        Codec::Gzip1 => gzip::compress(tile),
        Codec::Gzip2 => gzip::compress(&gzip::shuffle(tile, width)),
    }
}

// The one-column table of FITS 4.0 §10.1.3. Descriptors are 32-bit ('P')
// while the heap fits a signed 32-bit offset, and 64-bit ('Q') above that.
struct Table {
    rows: usize,
    wide: bool,
    heap: usize,
    longest: usize,
}

impl Table {
    fn new(tiles: &[Vec<u8>]) -> Self {
        let heap: usize = tiles.iter().map(Vec::len).sum();
        Self::with_descriptors(tiles, i32::try_from(heap).is_err())
    }

    fn with_descriptors(tiles: &[Vec<u8>], wide: bool) -> Self {
        Self {
            rows: tiles.len(),
            wide,
            heap: tiles.iter().map(Vec::len).sum(),
            longest: tiles.iter().map(Vec::len).max().unwrap_or(0),
        }
    }

    fn row_bytes(&self) -> usize {
        if self.wide {
            16
        } else {
            8
        }
    }

    fn header(&self) -> Header {
        let code = if self.wide { 'Q' } else { 'P' };
        let mut header = Header::new();
        header.add_keyword(Keyword::string("XTENSION", "BINTABLE"));
        header.add_keyword(Keyword::integer("BITPIX", 8));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::integer("NAXIS1", count(self.row_bytes())));
        header.add_keyword(Keyword::integer("NAXIS2", count(self.rows)));
        header.add_keyword(Keyword::integer("PCOUNT", count(self.heap)));
        header.add_keyword(Keyword::integer("GCOUNT", 1));
        header.add_keyword(Keyword::integer("TFIELDS", 1));
        header.add_keyword(Keyword::string("TTYPE1", "COMPRESSED_DATA"));
        header.add_keyword(Keyword::string(
            "TFORM1",
            &format!("1{code}B({})", self.longest),
        ));
        header
    }

    fn data(&self, tiles: &[Vec<u8>]) -> Vec<u8> {
        let mut data = Vec::with_capacity(self.rows * self.row_bytes() + self.heap);
        let mut offset = 0usize;
        for tile in tiles {
            if self.wide {
                data.extend_from_slice(&(tile.len() as u64).to_be_bytes());
                data.extend_from_slice(&(offset as u64).to_be_bytes());
            } else {
                data.extend_from_slice(&(tile.len() as u32).to_be_bytes());
                data.extend_from_slice(&(offset as u32).to_be_bytes());
            }
            offset += tile.len();
        }
        for tile in tiles {
            data.extend_from_slice(tile);
        }
        pad_to_block(&mut data, 0);
        data
    }
}

#[cfg(test)]
mod tests {
    use super::super::decode::decode_image;
    use super::*;

    #[test]
    fn tile_shape_defaults_to_rows_and_clamps() {
        assert_eq!(tile_shape(&[9, 4, 3], None).unwrap(), vec![9, 1, 1]);
        assert_eq!(tile_shape(&[9, 4, 3], Some((5, 8))).unwrap(), vec![5, 4, 1]);
        assert_eq!(tile_shape(&[9], Some((20, 3))).unwrap(), vec![9]);
        assert!(tile_shape(&[], None).is_err());
        assert!(tile_shape(&[9, 0], None).is_err());
        assert!(tile_shape(&[9, 4], Some((0, 1))).is_err());
    }

    #[test]
    fn codec_follows_the_pixel_type() {
        use CompressionAlgorithm::{Gzip, Gzip2, HCompress, Plio, Rice};
        for bitpix in [BitPix::U8, BitPix::I16, BitPix::I32] {
            assert!(codec(Rice, bitpix).is_ok());
        }
        for bitpix in [BitPix::I64, BitPix::F32, BitPix::F64] {
            assert!(matches!(
                codec(Rice, bitpix),
                Err(FitsError::UnsupportedCompression(_))
            ));
            assert!(codec(Gzip, bitpix).is_ok());
            assert!(codec(Gzip2, bitpix).is_ok());
        }
        for algorithm in [HCompress, Plio] {
            assert!(matches!(
                codec(algorithm, BitPix::I16),
                Err(FitsError::UnsupportedCompression(_))
            ));
        }
        match codec(HCompress, BitPix::F32) {
            Err(FitsError::UnsupportedCompression(m)) => {
                assert_eq!(
                    m,
                    "HCOMPRESS_1 tiles cannot be written; use GZIP_1 or GZIP_2"
                )
            }
            _ => panic!("HCOMPRESS_1 accepted F32 pixels"),
        }
    }

    // A heap past 2^31 bytes needs 64-bit descriptors. The test forces them
    // on a small image and reads it back.
    #[test]
    fn wide_descriptors_read_back() {
        let pixels: Vec<u8> = (0..60u8).flat_map(|v| [v, 255 - v]).collect();
        let axes = [6, 5, 2];
        let grid = TileGrid::new(axes.to_vec(), vec![4, 3, 1]).unwrap();
        let codec = Codec::Rice { bytepix: 2 };
        let tiles = compress_tiles(&grid, codec, 2, &pixels).unwrap();
        let table = Table::with_descriptors(&tiles, true);
        let mut header = table.header();
        add_compression_keywords(
            &mut header,
            CompressionAlgorithm::Rice,
            codec,
            BitPix::I16,
            &axes,
            &[4, 3, 1],
        );
        assert!(header.get("TFORM1").as_str().unwrap().starts_with("1QB("));
        assert_eq!(header.get("NAXIS1").as_i64(), Some(16));
        let image = decode_image(&header, &table.data(&tiles)).unwrap();
        assert_eq!(image.data, pixels);
    }
}
