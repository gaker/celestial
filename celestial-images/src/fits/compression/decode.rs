use super::columns::{Columns, Quantized};
use super::quantize::{self, TileScaling};
use super::restore::image_header;
use super::spec::{Codec, Spec};
use super::table::Table;
use super::tiles::TileGrid;
use super::{gzip, rice};
use crate::core::types::BitPix;
use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::binary_table::tform::Element;
use crate::fits::header::parser::Header;

// A decoded tile-compressed image: the restored header, and the pixels as the
// big-endian bytes an uncompressed data unit would hold.
pub(crate) struct DecodedImage {
    pub(crate) header: Header,
    pub(crate) dimensions: Vec<usize>,
    pub(crate) bitpix: BitPix,
    pub(crate) data: Vec<u8>,
}

// Where one tile's pixels are stored. fpack stores a tile that does not
// quantize or compress well in GZIP_COMPRESSED_DATA or UNCOMPRESSED_DATA.
enum Source<'a> {
    Compressed(&'a [u8]),
    Gzip(&'a [u8]),
    Raw(&'a [u8]),
}

struct Decoder<'a> {
    spec: Spec,
    table: Table<'a>,
    columns: Columns,
    width: usize,
}

// `header` and `data` are the header and data unit of a BINTABLE with
// ZIMAGE = T.
pub(crate) fn decode_image(header: &Header, data: &[u8]) -> Result<DecodedImage> {
    let (decoder, grid) = Decoder::new(header, data)?;
    let image = decoder.image(&grid)?;
    let blank = decoder
        .columns
        .integer_blank(header, &decoder.table, &decoder.spec)?;
    Ok(DecodedImage {
        header: image_header(header, &decoder.spec, blank)?,
        dimensions: decoder.spec.axes,
        bitpix: decoder.spec.bitpix,
        data: image,
    })
}

impl<'a> Decoder<'a> {
    fn new(header: &Header, data: &'a [u8]) -> Result<(Self, TileGrid)> {
        let spec = Spec::from_header(header)?;
        let grid = TileGrid::new(spec.axes.clone(), spec.tile.clone())?;
        let table = Table::new(header, data)?;
        if table.rows() != grid.tile_count() {
            return Err(invalid(format!(
                "the table has {} rows for {} tiles",
                table.rows(),
                grid.tile_count()
            )));
        }
        let columns = Columns::find(header, &table, &spec)?;
        let width = spec.bitpix.bytes_per_pixel();
        let decoder = Self {
            spec,
            table,
            columns,
            width,
        };
        Ok((decoder, grid))
    }

    // Every tile is checked against its stored bytes before the image is
    // allocated, so a bad table cannot request a large allocation.
    fn image(&self, grid: &TileGrid) -> Result<Vec<u8>> {
        let sources = (0..grid.tile_count())
            .map(|row| self.source(row, grid.region(row).pixel_count()))
            .collect::<Result<Vec<_>>>()?;
        let mut image = zeroed(grid.pixel_count(), self.width)?;
        for (row, source) in sources.into_iter().enumerate() {
            let region = grid.region(row);
            let tile = self.tile(row, source, region.pixel_count())?;
            grid.place(&region, self.width, &tile, &mut image)?;
        }
        Ok(image)
    }

    fn source(&self, row: usize, pixels: usize) -> Result<Source<'a>> {
        if let Some(mask) = self.columns.mask {
            if !self.table.array(row, mask)?.bytes.is_empty() {
                return Err(FitsError::UnsupportedCompression(
                    "tiles with a NULL_PIXEL_MASK".to_string(),
                ));
            }
        }
        let compressed = self.bytes(row, self.columns.compressed)?;
        let gzip = match self.columns.gzip {
            Some(c) => self.bytes(row, c)?,
            None => &[],
        };
        let source = if !compressed.is_empty() {
            Source::Compressed(compressed)
        } else if !gzip.is_empty() {
            Source::Gzip(gzip)
        } else if let Some(c) = self.columns.raw {
            Source::Raw(self.raw_bytes(row, c)?)
        } else {
            Source::Compressed(compressed)
        };
        self.check_size(row, &source, pixels)?;
        Ok(source)
    }

    fn bytes(&self, row: usize, column: usize) -> Result<&'a [u8]> {
        let array = self.table.array(row, column)?;
        match array.element {
            Element::U8 => Ok(array.bytes),
            _ => Err(invalid(format!(
                "row {} holds compressed data that are not bytes",
                row + 1
            ))),
        }
    }

    fn raw_bytes(&self, row: usize, column: usize) -> Result<&'a [u8]> {
        let array = self.table.array(row, column)?;
        let expected = match self.spec.bitpix {
            BitPix::U8 => Element::U8,
            BitPix::I16 => Element::I16,
            BitPix::I32 => Element::I32,
            BitPix::I64 => Element::I64,
            BitPix::F32 => Element::F32,
            BitPix::F64 => Element::F64,
        };
        if array.element != expected {
            return Err(invalid(format!(
                "UNCOMPRESSED_DATA of type {:?} in an image of BITPIX {}",
                array.element,
                self.spec.bitpix.value()
            )));
        }
        Ok(array.bytes)
    }

    // A tile cannot hold more pixels than its stored bytes can encode. The
    // check runs on every tile before the image is allocated.
    fn check_size(&self, row: usize, source: &Source, pixels: usize) -> Result<()> {
        let fits = match *source {
            Source::Raw(bytes) => pixels.checked_mul(self.width) == Some(bytes.len()),
            Source::Gzip(bytes) => fits_inflated(bytes, pixels, self.width),
            Source::Compressed(bytes) => match self.spec.codec {
                Codec::Rice(params) => pixels <= rice::max_pixels(bytes.len(), params),
                Codec::Gzip1 | Codec::Gzip2 => fits_inflated(bytes, pixels, self.stored_width()),
                Codec::NoCompress => false,
            },
        };
        if fits {
            return Ok(());
        }
        Err(invalid(format!(
            "tile {} of {pixels} pixels does not match its stored bytes",
            row + 1
        )))
    }

    // Quantized pixels are stored as 32-bit integers (FITS 4.0 §10.2).
    fn stored_width(&self) -> usize {
        if self.columns.quantized.is_some() {
            4
        } else {
            self.width
        }
    }

    fn tile(&self, row: usize, source: Source, pixels: usize) -> Result<Vec<u8>> {
        match (source, self.columns.quantized) {
            (Source::Raw(bytes), _) => Ok(bytes.to_vec()),
            (Source::Gzip(bytes), _) => gzip::inflate(bytes, pixels * self.width),
            (Source::Compressed(bytes), Some(q)) => {
                let values = self.quantized_values(bytes, pixels)?;
                let floats = quantize::unquantize(&values, &self.scaling(row, q)?, row)?;
                Ok(float_bytes(&floats, self.spec.bitpix))
            }
            (Source::Compressed(bytes), None) => match self.spec.codec {
                Codec::Rice(params) if !is_float(self.spec.bitpix) => {
                    integer_bytes(&rice::decode(bytes, pixels, params)?, self.spec.bitpix)
                }
                Codec::Gzip1 | Codec::Gzip2 => self.gzip_tile(bytes, pixels, self.width),
                _ => Err(invalid(format!(
                    "{:?} tiles cannot hold BITPIX {} pixels",
                    self.spec.codec,
                    self.spec.bitpix.value()
                ))),
            },
        }
    }

    fn quantized_values(&self, bytes: &[u8], pixels: usize) -> Result<Vec<i64>> {
        match self.spec.codec {
            Codec::Rice(params) => rice::decode(bytes, pixels, params),
            Codec::Gzip1 | Codec::Gzip2 => {
                let tile = self.gzip_tile(bytes, pixels, 4)?;
                Ok(tile
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|&c| i64::from(i32::from_be_bytes(c)))
                    .collect())
            }
            Codec::NoCompress => Err(invalid("NOCOMPRESS tile in COMPRESSED_DATA".to_string())),
        }
    }

    fn gzip_tile(&self, bytes: &[u8], pixels: usize, width: usize) -> Result<Vec<u8>> {
        let tile = gzip::inflate(bytes, pixels * width)?;
        Ok(match self.spec.codec {
            Codec::Gzip2 => gzip::unshuffle(&tile, width),
            _ => tile,
        })
    }

    fn scaling(&self, row: usize, q: Quantized) -> Result<TileScaling> {
        // FITS 4.0 §10.1.3: a ZBLANK column takes precedence over the keyword.
        let blank = match self.columns.blank {
            Some(c) => Some(self.table.integer(row, c)?),
            None => self.spec.blank,
        };
        Ok(TileScaling {
            scale: self.table.real(row, q.scale)?,
            zero: self.table.real(row, q.zero)?,
            blank,
            method: q.method,
            dither_seed: self.spec.dither_seed,
        })
    }
}

fn fits_inflated(bytes: &[u8], pixels: usize, width: usize) -> bool {
    pixels
        .checked_mul(width)
        .is_some_and(|n| n <= gzip::max_inflated(bytes.len()))
}

fn is_float(bitpix: BitPix) -> bool {
    matches!(bitpix, BitPix::F32 | BitPix::F64)
}

// NaN marks undefined pixels. `f64::NAN as f32` gives the quiet NaN
// 0x7FC00000, which is what astropy returns for them.
fn float_bytes(values: &[f64], bitpix: BitPix) -> Vec<u8> {
    match bitpix {
        BitPix::F32 => values
            .iter()
            .flat_map(|&v| (v as f32).to_be_bytes())
            .collect(),
        _ => values.iter().flat_map(|&v| v.to_be_bytes()).collect(),
    }
}

fn integer_bytes(values: &[i64], bitpix: BitPix) -> Result<Vec<u8>> {
    let out_of_range = |v: i64| {
        invalid(format!(
            "decoded value {v} does not fit BITPIX {}",
            bitpix.value()
        ))
    };
    let mut out = Vec::with_capacity(values.len() * bitpix.bytes_per_pixel());
    for &v in values {
        match bitpix {
            BitPix::U8 => out.push(u8::try_from(v).map_err(|_| out_of_range(v))?),
            BitPix::I16 => out.extend(i16::try_from(v).map_err(|_| out_of_range(v))?.to_be_bytes()),
            BitPix::I32 => out.extend(i32::try_from(v).map_err(|_| out_of_range(v))?.to_be_bytes()),
            BitPix::I64 => out.extend(v.to_be_bytes()),
            BitPix::F32 | BitPix::F64 => return Err(out_of_range(v)),
        }
    }
    Ok(out)
}

fn zeroed(pixels: usize, width: usize) -> Result<Vec<u8>> {
    let len = pixels
        .checked_mul(width)
        .ok_or_else(|| invalid(format!("{pixels} pixels are too many")))?;
    let mut image = Vec::new();
    image
        .try_reserve_exact(len)
        .map_err(|_| invalid(format!("cannot allocate a {len}-byte image")))?;
    image.resize(len, 0);
    Ok(image)
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_bytes_checks_the_range_of_each_type() {
        assert_eq!(integer_bytes(&[0, 255], BitPix::U8).unwrap(), vec![0, 255]);
        assert!(integer_bytes(&[256], BitPix::U8).is_err());
        assert!(integer_bytes(&[-1], BitPix::U8).is_err());
        assert_eq!(
            integer_bytes(&[-2, 32767], BitPix::I16).unwrap(),
            vec![0xFF, 0xFE, 0x7F, 0xFF]
        );
        assert!(integer_bytes(&[32768], BitPix::I16).is_err());
        assert!(integer_bytes(&[i64::from(i32::MIN) - 1], BitPix::I32).is_err());
        assert_eq!(integer_bytes(&[-1], BitPix::I64).unwrap(), vec![0xFF; 8]);
        assert!(integer_bytes(&[1], BitPix::F32).is_err());
    }

    #[test]
    fn float_bytes_writes_nan_as_the_quiet_nan() {
        assert_eq!(
            float_bytes(&[f64::NAN, 1.5], BitPix::F32),
            [0x7FC0_0000u32.to_be_bytes(), 1.5f32.to_bits().to_be_bytes()].concat()
        );
        assert_eq!(
            float_bytes(&[-2.0], BitPix::F64),
            (-2.0f64).to_be_bytes().to_vec()
        );
    }

    #[test]
    fn zeroed_rejects_sizes_that_overflow() {
        assert!(zeroed(usize::MAX, 2).is_err());
        assert_eq!(zeroed(3, 2).unwrap(), vec![0; 6]);
    }
}
