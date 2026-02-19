//! Memory-mapped catalog reader for HEALPix-indexed star catalogs.
//!
//! The catalog binary format has three contiguous sections:
//!
//! 1. **Header** (64 bytes) — magic, version, HEALPix parameters, star count, epoch
//! 2. **Pixel offset table** (`npix × 16` bytes) — byte offset and count per pixel
//! 3. **Star data** (`total_stars × 56` bytes) — [`StarRecord`] structs grouped by pixel
//!
//! Open a catalog with [`Catalog::open`], then query stars by pixel index
//! or use the higher-level cone search in [`super::cone`].

use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::fmt;
use std::fs::File;
use std::path::Path;

const CATALOG_MAGIC: &[u8; 4] = b"CCAT";
const CATALOG_VERSION: u32 = 1;
const HEADER_SIZE: usize = 64;
const PIXEL_ENTRY_SIZE: usize = 16;

/// Metadata parsed from the first 64 bytes of a catalog file.
#[derive(Debug, Clone)]
pub struct CatalogHeader {
    /// HEALPix order (nside = 2^order). Order 8 gives 786,432 pixels.
    pub order: u32,
    /// HEALPix nside parameter. Always `2.pow(order)`.
    pub nside: u32,
    /// Total number of HEALPix pixels (`12 * nside * nside`).
    pub npix: u64,
    /// Total number of star records in the catalog.
    pub total_stars: u64,
    /// Catalog epoch as a Julian year (e.g. 2016.0 for Gaia DR3).
    pub epoch: f64,
    /// Faintest magnitude included in the catalog.
    pub mag_limit: f32,
}

impl fmt::Display for CatalogHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let avg = self.total_stars as f64 / self.npix as f64;
        writeln!(f, "HEALPix order: {}", self.order)?;
        writeln!(f, "nside: {}", self.nside)?;
        writeln!(f, "npix: {}", self.npix)?;
        writeln!(f, "Total stars: {}", self.total_stars)?;
        writeln!(f, "Epoch: J{:.1}", self.epoch)?;
        writeln!(f, "Magnitude limit: {:.2}", self.mag_limit)?;
        write!(f, "Average stars per pixel: {:.1}", avg)
    }
}

/// A single star entry (56 bytes, `repr(C)`).
///
/// Laid out for zero-copy reads from the memory-mapped file. Fields are
/// stored in catalog-epoch coordinates; use [`super::cone::cone_search`]
/// with an observation epoch to get proper-motion-corrected positions.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct StarRecord {
    /// Gaia DR3 source ID, or negative Hipparcos ID for Hipparcos-only stars.
    pub source_id: i64,
    /// Right ascension at catalog epoch, in degrees.
    pub ra: f64,
    /// Declination at catalog epoch, in degrees.
    pub dec: f64,
    /// Proper motion in RA (μα*), including cos(δ) factor, in mas/yr.
    pub pmra: f64,
    /// Proper motion in declination, in mas/yr.
    pub pmdec: f64,
    /// Trigonometric parallax, in milliarcseconds.
    pub parallax: f64,
    /// Apparent magnitude (G-band for Gaia, V-band for Hipparcos).
    pub mag: f32,
    /// Bitfield of quality and source flags. See `FLAG_*` constants.
    pub flags: u16,
    pub(crate) _padding: u16,
}

/// Star has measured proper motion (pmra and pmdec are valid).
pub const FLAG_HAS_PROPER_MOTION: u16 = 1 << 0;
/// Star has measured parallax.
pub const FLAG_HAS_PARALLAX: u16 = 1 << 1;
/// Gaia RUWE > 1.4 — astrometric solution may be unreliable.
pub const FLAG_RUWE_SUSPECT: u16 = 1 << 2;
/// No 5-parameter astrometric solution in Gaia.
pub const FLAG_NO_5PARAM: u16 = 1 << 3;
/// BP/RP flux excess factor is suspect (possible blend or extended source).
pub const FLAG_BP_RP_EXCESS_SUSPECT: u16 = 1 << 4;
/// Star originates from Hipparcos, not Gaia DR3.
pub const FLAG_SOURCE_HIPPARCOS: u16 = 1 << 5;

const _: () = assert!(std::mem::size_of::<StarRecord>() == 56);

#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct PixelEntry {
    offset: u64,
    count: u32,
    _reserved: u32,
}

const _: () = assert!(std::mem::size_of::<PixelEntry>() == 16);

/// Memory-mapped handle to a HEALPix star catalog.
///
/// Created by [`Catalog::open`]. The underlying file stays mapped for the
/// lifetime of this value. Star slices returned by [`Catalog::stars_in_pixel`]
/// borrow directly from the map with no allocation or copying.
pub struct Catalog {
    mmap: Mmap,
    header: CatalogHeader,
}

impl Catalog {
    /// Open and memory-map a catalog file.
    ///
    /// Validates the header (magic bytes, version, HEALPix consistency) and
    /// returns immediately. No star data is read until you call
    /// [`Catalog::stars_in_pixel`] or run a cone search.
    ///
    /// # Errors
    /// Returns an error if the file cannot be opened, is too small, has an
    /// invalid header, or has inconsistent HEALPix parameters.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file =
            File::open(path).with_context(|| format!("Failed to open catalog file: {:?}", path))?;

        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to memory-map catalog file: {:?}", path))?;

        if mmap.len() < HEADER_SIZE {
            bail!("Catalog file too small: {} bytes", mmap.len());
        }

        let header = parse_header(&mmap)?;

        let expected_offset_table_size = header.npix as usize * PIXEL_ENTRY_SIZE;
        let min_size = HEADER_SIZE + expected_offset_table_size;
        if mmap.len() < min_size {
            bail!(
                "Catalog file too small for offset table: {} bytes, expected at least {}",
                mmap.len(),
                min_size
            );
        }

        Ok(Self { mmap, header })
    }

    /// Returns the catalog header (order, nside, star count, epoch, etc.).
    pub fn header(&self) -> &CatalogHeader {
        &self.header
    }

    /// Returns a zero-copy slice of all stars in the given HEALPix pixel.
    ///
    /// The returned slice borrows directly from the memory map. Returns an
    /// empty slice if `pixel_index` is out of range, the pixel contains no
    /// stars, or the underlying data is misaligned.
    pub fn stars_in_pixel(&self, pixel_index: u64) -> &[StarRecord] {
        if pixel_index >= self.header.npix {
            return &[];
        }

        let offset_table_start = HEADER_SIZE;
        let entry_offset = offset_table_start + (pixel_index as usize * PIXEL_ENTRY_SIZE);

        if entry_offset + PIXEL_ENTRY_SIZE > self.mmap.len() {
            return &[];
        }

        let entry_bytes = &self.mmap[entry_offset..entry_offset + PIXEL_ENTRY_SIZE];
        let offset = u64::from_le_bytes(entry_bytes[0..8].try_into().unwrap());
        let count = u32::from_le_bytes(entry_bytes[8..12].try_into().unwrap());

        if count == 0 {
            return &[];
        }

        let star_data_start = HEADER_SIZE + (self.header.npix as usize * PIXEL_ENTRY_SIZE);
        let star_offset = star_data_start + offset as usize;
        let star_size = count as usize * std::mem::size_of::<StarRecord>();

        if star_offset + star_size > self.mmap.len() {
            return &[];
        }

        let star_bytes = &self.mmap[star_offset..star_offset + star_size];
        let ptr = star_bytes.as_ptr();
        if !(ptr as usize).is_multiple_of(std::mem::align_of::<StarRecord>()) {
            return &[];
        }

        unsafe { std::slice::from_raw_parts(ptr as *const StarRecord, count as usize) }
    }

    /// Returns the total size of the memory-mapped file in bytes.
    pub fn file_size(&self) -> usize {
        self.mmap.len()
    }
}

fn parse_header(mmap: &Mmap) -> Result<CatalogHeader> {
    let header_bytes = &mmap[0..HEADER_SIZE];

    let magic = &header_bytes[0..4];
    if magic != CATALOG_MAGIC {
        bail!(
            "Invalid catalog magic: expected {:?}, got {:?}",
            CATALOG_MAGIC,
            magic
        );
    }

    let version = u32::from_le_bytes(header_bytes[4..8].try_into().unwrap());
    if version != CATALOG_VERSION {
        bail!(
            "Unsupported catalog version: expected {}, got {}",
            CATALOG_VERSION,
            version
        );
    }

    let order = u32::from_le_bytes(header_bytes[8..12].try_into().unwrap());
    let nside = u32::from_le_bytes(header_bytes[12..16].try_into().unwrap());
    let npix = u64::from_le_bytes(header_bytes[16..24].try_into().unwrap());
    let total_stars = u64::from_le_bytes(header_bytes[24..32].try_into().unwrap());
    let epoch = f64::from_le_bytes(header_bytes[32..40].try_into().unwrap());
    let mag_limit = f32::from_le_bytes(header_bytes[40..44].try_into().unwrap());

    let expected_nside = 1u32 << order;
    if nside != expected_nside {
        bail!(
            "Inconsistent nside: order {} implies nside {}, got {}",
            order,
            expected_nside,
            nside
        );
    }

    let expected_npix = 12u64 * (nside as u64) * (nside as u64);
    if npix != expected_npix {
        bail!(
            "Inconsistent npix: nside {} implies npix {}, got {}",
            nside,
            expected_npix,
            npix
        );
    }

    Ok(CatalogHeader {
        order,
        nside,
        npix,
        total_stars,
        epoch,
        mag_limit,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_star_record_size() {
        assert_eq!(std::mem::size_of::<StarRecord>(), 56);
    }

    #[test]
    fn test_pixel_entry_size() {
        assert_eq!(std::mem::size_of::<PixelEntry>(), 16);
    }

    #[test]
    fn test_star_record_alignment() {
        assert_eq!(std::mem::align_of::<StarRecord>(), 8);
    }

    fn make_star(source_id: i64, ra: f64, dec: f64, mag: f32, flags: u16) -> StarRecord {
        StarRecord {
            source_id,
            ra,
            dec,
            pmra: 0.0,
            pmdec: 0.0,
            parallax: 0.0,
            mag,
            flags,
            _padding: 0,
        }
    }

    fn star_to_bytes(star: &StarRecord) -> &[u8] {
        unsafe {
            std::slice::from_raw_parts(
                star as *const StarRecord as *const u8,
                std::mem::size_of::<StarRecord>(),
            )
        }
    }

    fn build_test_catalog(
        order: u32,
        epoch: f64,
        mag_limit: f32,
        pixel_stars: &[(u64, Vec<StarRecord>)],
    ) -> NamedTempFile {
        let nside = 1u32 << order;
        let npix = 12u64 * (nside as u64) * (nside as u64);
        let total_stars: u64 = pixel_stars.iter().map(|(_, s)| s.len() as u64).sum();

        let mut buf: Vec<u8> = Vec::new();

        // Header (64 bytes)
        buf.extend_from_slice(b"CCAT");
        buf.extend_from_slice(&1u32.to_le_bytes());
        buf.extend_from_slice(&order.to_le_bytes());
        buf.extend_from_slice(&nside.to_le_bytes());
        buf.extend_from_slice(&npix.to_le_bytes());
        buf.extend_from_slice(&total_stars.to_le_bytes());
        buf.extend_from_slice(&epoch.to_le_bytes());
        buf.extend_from_slice(&mag_limit.to_le_bytes());
        buf.extend_from_slice(&[0u8; 20]); // padding to 64 bytes

        assert_eq!(buf.len(), HEADER_SIZE);

        // Build sorted star data and compute offsets per pixel
        let mut offsets: Vec<(u64, u32)> = vec![(0, 0); npix as usize];
        let mut star_data: Vec<u8> = Vec::new();

        for &(pixel_idx, ref stars) in pixel_stars {
            let byte_offset = star_data.len() as u64;
            offsets[pixel_idx as usize] = (byte_offset, stars.len() as u32);
            for star in stars {
                star_data.extend_from_slice(star_to_bytes(star));
            }
        }

        // Pixel offset table (npix * 16 bytes)
        for &(offset, count) in &offsets {
            buf.extend_from_slice(&offset.to_le_bytes());
            buf.extend_from_slice(&count.to_le_bytes());
            buf.extend_from_slice(&0u32.to_le_bytes()); // reserved
        }

        // Star data section
        buf.extend_from_slice(&star_data);

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_open_valid_catalog() {
        let star = make_star(1001, 180.0, -45.0, 5.5, FLAG_HAS_PROPER_MOTION);
        let file = build_test_catalog(1, 2016.0, 21.0, &[(0, vec![star])]);

        let catalog = Catalog::open(file.path()).unwrap();
        let hdr = catalog.header();
        assert_eq!(hdr.order, 1);
        assert_eq!(hdr.nside, 2);
        assert_eq!(hdr.npix, 48);
        assert_eq!(hdr.total_stars, 1);
        assert_eq!(hdr.epoch, 2016.0);
        assert_eq!(hdr.mag_limit, 21.0);
    }

    #[test]
    fn test_open_truncated_file() {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&[0u8; 32]).unwrap();
        file.flush().unwrap();

        let result = Catalog::open(file.path());
        let msg = result.err().expect("expected error").to_string();
        assert!(msg.contains("too small"), "unexpected error: {}", msg);
    }

    #[test]
    fn test_open_bad_magic() {
        let mut buf = vec![0u8; HEADER_SIZE + 48 * PIXEL_ENTRY_SIZE];
        buf[0..4].copy_from_slice(b"XXXX");
        buf[4..8].copy_from_slice(&1u32.to_le_bytes());
        buf[8..12].copy_from_slice(&1u32.to_le_bytes()); // order=1
        buf[12..16].copy_from_slice(&2u32.to_le_bytes()); // nside=2
        buf[16..24].copy_from_slice(&48u64.to_le_bytes()); // npix=48

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();

        let result = Catalog::open(file.path());
        let msg = result.err().expect("expected error").to_string();
        assert!(
            msg.contains("Invalid catalog magic"),
            "unexpected error: {}",
            msg
        );
    }

    #[test]
    fn test_open_bad_version() {
        let mut buf = vec![0u8; HEADER_SIZE + 48 * PIXEL_ENTRY_SIZE];
        buf[0..4].copy_from_slice(b"CCAT");
        buf[4..8].copy_from_slice(&99u32.to_le_bytes()); // bad version
        buf[8..12].copy_from_slice(&1u32.to_le_bytes());
        buf[12..16].copy_from_slice(&2u32.to_le_bytes());
        buf[16..24].copy_from_slice(&48u64.to_le_bytes());

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();

        let result = Catalog::open(file.path());
        let msg = result.err().expect("expected error").to_string();
        assert!(
            msg.contains("Unsupported catalog version"),
            "unexpected error: {}",
            msg
        );
    }

    #[test]
    fn test_open_inconsistent_nside() {
        let mut buf = vec![0u8; HEADER_SIZE + 48 * PIXEL_ENTRY_SIZE];
        buf[0..4].copy_from_slice(b"CCAT");
        buf[4..8].copy_from_slice(&1u32.to_le_bytes());
        buf[8..12].copy_from_slice(&1u32.to_le_bytes()); // order=1
        buf[12..16].copy_from_slice(&7u32.to_le_bytes()); // nside=7 (should be 2)
        buf[16..24].copy_from_slice(&48u64.to_le_bytes());

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();

        let result = Catalog::open(file.path());
        let msg = result.err().expect("expected error").to_string();
        assert!(
            msg.contains("Inconsistent nside"),
            "unexpected error: {}",
            msg
        );
    }

    #[test]
    fn test_open_inconsistent_npix() {
        let mut buf = vec![0u8; HEADER_SIZE + 48 * PIXEL_ENTRY_SIZE];
        buf[0..4].copy_from_slice(b"CCAT");
        buf[4..8].copy_from_slice(&1u32.to_le_bytes());
        buf[8..12].copy_from_slice(&1u32.to_le_bytes()); // order=1
        buf[12..16].copy_from_slice(&2u32.to_le_bytes()); // nside=2 (correct)
        buf[16..24].copy_from_slice(&999u64.to_le_bytes()); // npix=999 (should be 48)

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();

        let result = Catalog::open(file.path());
        let msg = result.err().expect("expected error").to_string();
        assert!(
            msg.contains("Inconsistent npix"),
            "unexpected error: {}",
            msg
        );
    }

    #[test]
    fn test_stars_in_pixel_populated() {
        let star = make_star(42, 83.633, -5.375, 0.42, FLAG_SOURCE_HIPPARCOS);
        let file = build_test_catalog(1, 2016.0, 21.0, &[(5, vec![star])]);
        let catalog = Catalog::open(file.path()).unwrap();

        let stars = catalog.stars_in_pixel(5);
        assert_eq!(stars.len(), 1);
        assert_eq!(stars[0].source_id, 42);
    }

    #[test]
    fn test_stars_in_pixel_empty() {
        let star = make_star(1, 10.0, 20.0, 8.0, 0);
        let file = build_test_catalog(1, 2016.0, 21.0, &[(0, vec![star])]);
        let catalog = Catalog::open(file.path()).unwrap();

        let stars = catalog.stars_in_pixel(1);
        assert_eq!(stars.len(), 0);
    }

    #[test]
    fn test_stars_in_pixel_out_of_bounds() {
        let file = build_test_catalog(1, 2016.0, 21.0, &[]);
        let catalog = Catalog::open(file.path()).unwrap();

        assert_eq!(catalog.stars_in_pixel(48).len(), 0);
        assert_eq!(catalog.stars_in_pixel(999).len(), 0);
        assert_eq!(catalog.stars_in_pixel(u64::MAX).len(), 0);
    }

    #[test]
    fn test_stars_in_pixel_multiple_stars() {
        let stars = vec![
            make_star(100, 10.0, 20.0, 5.0, 0),
            make_star(101, 10.1, 20.1, 6.0, FLAG_HAS_PROPER_MOTION),
            make_star(102, 10.2, 20.2, 7.0, FLAG_HAS_PARALLAX),
        ];
        let file = build_test_catalog(1, 2016.0, 21.0, &[(12, stars)]);
        let catalog = Catalog::open(file.path()).unwrap();

        let result = catalog.stars_in_pixel(12);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].source_id, 100);
        assert_eq!(result[1].source_id, 101);
        assert_eq!(result[2].source_id, 102);
    }

    #[test]
    fn test_file_size_matches() {
        let star = make_star(1, 0.0, 0.0, 10.0, 0);
        let file = build_test_catalog(1, 2016.0, 21.0, &[(0, vec![star])]);
        let expected = HEADER_SIZE + 48 * PIXEL_ENTRY_SIZE + std::mem::size_of::<StarRecord>();

        let catalog = Catalog::open(file.path()).unwrap();
        assert_eq!(catalog.file_size(), expected);
    }

    #[test]
    fn test_star_fields_round_trip() {
        let star = StarRecord {
            source_id: -9_999_999,
            ra: 359.99999999,
            dec: -89.99999999,
            pmra: 1234.5678,
            pmdec: -8765.4321,
            parallax: 0.001,
            mag: 21.49,
            flags: FLAG_HAS_PROPER_MOTION | FLAG_HAS_PARALLAX | FLAG_SOURCE_HIPPARCOS,
            _padding: 0,
        };
        let file = build_test_catalog(1, 2016.0, 21.5, &[(47, vec![star])]);
        let catalog = Catalog::open(file.path()).unwrap();

        let result = catalog.stars_in_pixel(47);
        assert_eq!(result.len(), 1);
        let s = &result[0];
        assert_eq!(s.source_id, -9_999_999);
        assert_eq!(s.ra, 359.99999999);
        assert_eq!(s.dec, -89.99999999);
        assert_eq!(s.pmra, 1234.5678);
        assert_eq!(s.pmdec, -8765.4321);
        assert_eq!(s.parallax, 0.001);
        assert_eq!(s.mag, 21.49);
        assert_eq!(
            s.flags,
            FLAG_HAS_PROPER_MOTION | FLAG_HAS_PARALLAX | FLAG_SOURCE_HIPPARCOS
        );
    }

    #[test]
    fn test_catalog_header_display() {
        let header = CatalogHeader {
            order: 8,
            nside: 256,
            npix: 786432,
            total_stars: 37_000_000,
            epoch: 2016.0,
            mag_limit: 21.0,
        };
        let output = format!("{}", header);

        assert!(output.contains("HEALPix order: 8"), "missing order");
        assert!(output.contains("nside: 256"), "missing nside");
        assert!(output.contains("npix: 786432"), "missing npix");
        assert!(
            output.contains("Total stars: 37000000"),
            "missing total_stars"
        );
        assert!(output.contains("Epoch: J2016.0"), "missing epoch");
        assert!(
            output.contains("Magnitude limit: 21.00"),
            "missing mag_limit"
        );
        assert!(
            output.contains("Average stars per pixel: 47.0"),
            "missing avg"
        );
    }
}
