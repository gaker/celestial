//! Memory-mapped catalog reader for HEALPix-indexed star catalogs.
//!
//! The catalog binary format has three contiguous sections:
//!
//! 1. **Header** (64 bytes): magic, version, HEALPix parameters, star count, epoch
//! 2. **Pixel offset table** (`npix × 16` bytes): byte offset and count per pixel
//! 3. **Star data** (`total_stars × 56` bytes): [`StarRecord`] structs grouped by pixel
//!
//! Open a catalog with [`Catalog::open`], then query stars by pixel index
//! or use the higher-level cone search in [`super::cone`].

use memmap2::Mmap;
use std::fmt;
use std::fs::File;
use std::path::Path;

use super::format::{
    concat_le, le_field, PixelEntry, CATALOG_MAGIC, CATALOG_VERSION, HEADER_SIZE, PIXEL_ENTRY_SIZE,
    STAR_RECORD_SIZE,
};
use super::validate::{check_header, check_magic_and_version, check_pixel_table, Layout};
use crate::error::{CatalogError, CatalogResult};

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

impl CatalogHeader {
    // Writes the magic, the current format version and the fields. Bytes 44 to 63 are zero.
    pub fn to_bytes(&self) -> [u8; HEADER_SIZE] {
        concat_le(&[
            CATALOG_MAGIC,
            &CATALOG_VERSION.to_le_bytes(),
            &self.order.to_le_bytes(),
            &self.nside.to_le_bytes(),
            &self.npix.to_le_bytes(),
            &self.total_stars.to_le_bytes(),
            &self.epoch.to_le_bytes(),
            &self.mag_limit.to_le_bytes(),
        ])
    }

    fn parse(file: &[u8]) -> CatalogResult<Self> {
        let Some(bytes) = file.first_chunk::<HEADER_SIZE>() else {
            return Err(CatalogError::invalid_header(format!(
                "Catalog file too small: {} bytes, the header alone is {HEADER_SIZE}",
                file.len()
            )));
        };
        check_magic_and_version(le_field(bytes, 0), u32::from_le_bytes(le_field(bytes, 4)))?;
        let header = Self {
            order: u32::from_le_bytes(le_field(bytes, 8)),
            nside: u32::from_le_bytes(le_field(bytes, 12)),
            npix: u64::from_le_bytes(le_field(bytes, 16)),
            total_stars: u64::from_le_bytes(le_field(bytes, 24)),
            epoch: f64::from_le_bytes(le_field(bytes, 32)),
            mag_limit: f32::from_le_bytes(le_field(bytes, 40)),
        };
        check_header(&header)?;
        Ok(header)
    }
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
    /// Gaia DR3 source ID. A Hipparcos record with a Gaia DR3 cross-match carries that
    /// source ID; one without a cross-match carries its HIP number negated.
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
    /// Apparent magnitude: Gaia G, or Hipparcos Hp for records with [`FLAG_SOURCE_HIPPARCOS`].
    pub mag: f32,
    /// Bitfield of quality and source flags. See `FLAG_*` constants.
    pub flags: u16,
    pub(crate) _padding: u16,
}

impl StarRecord {
    // The byte layout is the in-memory repr(C) layout on a little-endian target, so these
    // agree with the records Catalog::stars_in_pixel reads in place.
    pub fn to_le_bytes(&self) -> [u8; STAR_RECORD_SIZE] {
        concat_le(&[
            &self.source_id.to_le_bytes(),
            &self.ra.to_le_bytes(),
            &self.dec.to_le_bytes(),
            &self.pmra.to_le_bytes(),
            &self.pmdec.to_le_bytes(),
            &self.parallax.to_le_bytes(),
            &self.mag.to_le_bytes(),
            &self.flags.to_le_bytes(),
            &self._padding.to_le_bytes(),
        ])
    }

    pub fn from_le_bytes(bytes: &[u8; STAR_RECORD_SIZE]) -> Self {
        Self {
            source_id: i64::from_le_bytes(le_field(bytes, 0)),
            ra: f64::from_le_bytes(le_field(bytes, 8)),
            dec: f64::from_le_bytes(le_field(bytes, 16)),
            pmra: f64::from_le_bytes(le_field(bytes, 24)),
            pmdec: f64::from_le_bytes(le_field(bytes, 32)),
            parallax: f64::from_le_bytes(le_field(bytes, 40)),
            mag: f32::from_le_bytes(le_field(bytes, 48)),
            flags: u16::from_le_bytes(le_field(bytes, 52)),
            _padding: u16::from_le_bytes(le_field(bytes, 54)),
        }
    }
}

/// Star has measured proper motion (pmra and pmdec are valid).
pub const FLAG_HAS_PROPER_MOTION: u16 = 1 << 0;
/// Star has measured parallax.
pub const FLAG_HAS_PARALLAX: u16 = 1 << 1;
/// Gaia RUWE > 1.4: the astrometric solution may be unreliable.
pub const FLAG_RUWE_SUSPECT: u16 = 1 << 2;
/// No 5-parameter astrometric solution in Gaia.
pub const FLAG_NO_5PARAM: u16 = 1 << 3;
/// BP/RP flux excess factor is suspect (possible blend or extended source).
pub const FLAG_BP_RP_EXCESS_SUSPECT: u16 = 1 << 4;
/// Star originates from Hipparcos, not Gaia DR3.
pub const FLAG_SOURCE_HIPPARCOS: u16 = 1 << 5;
// Gaia DR3 astrometric_params_solved is 31 (5-parameter solution) or 95 (6-parameter
// solution), so the record has a measured parallax and proper motion.
pub const FLAG_5_OR_6_PARAM: u16 = 1 << 6;
// Gaia DR3 duplicated_source is true: the source was duplicated during processing and this
// is the identifier that was kept.
pub const FLAG_DUPLICATED_SOURCE: u16 = 1 << 7;

// Writes one pixel table entry. The reserved bytes 12 to 15 are zero.
pub fn pixel_entry_bytes(offset: u64, count: u32) -> [u8; PIXEL_ENTRY_SIZE] {
    concat_le(&[&offset.to_le_bytes(), &count.to_le_bytes()])
}

/// Memory-mapped handle to a HEALPix star catalog.
///
/// Created by [`Catalog::open`]. The underlying file stays mapped for the
/// lifetime of this value. Star slices returned by [`Catalog::stars_in_pixel`]
/// borrow directly from the map with no allocation or copying.
pub struct Catalog {
    mmap: Mmap,
    header: CatalogHeader,
    layout: Layout,
}

impl Catalog {
    /// Open and memory-map a catalog file.
    ///
    /// Validates the header (magic bytes, version, HEALPix consistency, epoch
    /// and magnitude limit), the file length and every pixel table entry. No
    /// star data is read until you call [`Catalog::stars_in_pixel`] or run a
    /// cone search.
    ///
    /// # Errors
    /// Returns [`CatalogError::Io`] if the file cannot be opened or mapped,
    /// [`CatalogError::InvalidHeader`] if the header is invalid or the file
    /// length does not match it, and [`CatalogError::CorruptTable`] if a pixel
    /// table entry does not point at its own run of contiguous records.
    pub fn open(path: impl AsRef<Path>) -> CatalogResult<Self> {
        let path = path.as_ref();
        let mmap = map_file(path)?;
        let header = CatalogHeader::parse(&mmap)?;
        let layout = Layout::new(&header, mmap.len())?;
        let catalog = Self {
            mmap,
            header,
            layout,
        };
        check_pixel_table(catalog.pixel_table(), catalog.header.total_stars)?;
        Ok(catalog)
    }

    /// Returns the catalog header (order, nside, star count, epoch, etc.).
    pub fn header(&self) -> &CatalogHeader {
        &self.header
    }

    /// Returns a zero-copy slice of all stars in the given HEALPix pixel.
    ///
    /// The returned slice borrows directly from the memory map. Returns an
    /// empty slice if `pixel_index` is not below `npix` or the pixel contains
    /// no stars.
    pub fn stars_in_pixel(&self, pixel_index: u64) -> &[StarRecord] {
        let entry = usize::try_from(pixel_index)
            .ok()
            .and_then(|pixel| self.pixel_table().get(pixel));
        let Some(entry) = entry else {
            return &[];
        };
        // open checked that the entry's records lie inside the star data, so the offset fits
        // in usize and both slices are in bounds.
        let first = (entry.offset / STAR_RECORD_SIZE as u64) as usize;
        &self.records()[first..][..entry.count as usize]
    }

    /// Returns the total size of the memory-mapped file in bytes.
    pub fn file_size(&self) -> usize {
        self.mmap.len()
    }

    fn pixel_table(&self) -> &[PixelEntry] {
        // SAFETY: Layout::new checked that the file holds the header and `pixels` entries, and
        // map_file checked that the map starts 8-byte aligned, so the entries are in bounds and
        // aligned.
        // PixelEntry is repr(C) integers with no padding, so every bit pattern is a value.
        unsafe {
            let first = self.mmap.as_ptr().add(HEADER_SIZE).cast::<PixelEntry>();
            std::slice::from_raw_parts(first, self.layout.pixels)
        }
    }

    fn records(&self) -> &[StarRecord] {
        // SAFETY: Layout::new checked that the file ends with `stars` records starting at
        // `records_start`, a multiple of 8 bytes into a map that starts 8-byte aligned.
        // StarRecord is repr(C) integers and floats with no implicit padding, so every bit
        // pattern is a value.
        unsafe {
            let first = self.mmap.as_ptr().add(self.layout.records_start);
            std::slice::from_raw_parts(first.cast::<StarRecord>(), self.layout.stars)
        }
    }
}

fn map_file(path: &Path) -> CatalogResult<Mmap> {
    let io_error = |action: &str, source| CatalogError::Io {
        context: format!("Failed to {action} catalog file {path:?}"),
        source,
    };
    let file = File::open(path).map_err(|e| io_error("open", e))?;
    // SAFETY: the map is read-only. Another process that truncates or rewrites the file while
    // it is mapped changes the bytes behind the slices this module hands out; catalog files
    // are written once by forge and then only read.
    let mmap = unsafe { Mmap::map(&file) }.map_err(|e| io_error("memory-map", e))?;
    if mmap
        .as_ptr()
        .align_offset(std::mem::align_of::<StarRecord>())
        != 0
    {
        let source = std::io::Error::other("the map does not start on an 8-byte boundary");
        return Err(io_error("align", source));
    }
    Ok(mmap)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::query::test_support::{catalog_bytes, header, make_star};
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

    // Writes a catalog with every pixel's records in pixel order. A pixel listed twice keeps
    // the stars of both lists.
    pub(crate) fn build_test_catalog(
        order: u32,
        epoch: f64,
        mag_limit: f32,
        pixel_stars: &[(u64, Vec<StarRecord>)],
    ) -> NamedTempFile {
        write_file(&catalog_bytes(order, epoch, mag_limit, pixel_stars))
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
        let mut buf = header_only_bytes(CATALOG_VERSION, 1, 2, 48);
        buf[0..4].copy_from_slice(b"XXXX");
        let file = write_file(&buf);

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
        let file = header_only_file(99, 1, 2, 48);

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
        let file = header_only_file(CATALOG_VERSION, 1, 7, 48); // nside 7, should be 2

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
        let file = header_only_file(CATALOG_VERSION, 1, 2, 999); // npix 999, should be 48

        let result = Catalog::open(file.path());
        let msg = result.err().expect("expected error").to_string();
        assert!(
            msg.contains("Inconsistent npix"),
            "unexpected error: {}",
            msg
        );
    }

    // A header followed by a 48-entry table of zeros, with any version, order, nside and npix.
    fn header_only_bytes(version: u32, order: u32, nside: u32, npix: u64) -> Vec<u8> {
        let header = CatalogHeader {
            order,
            nside,
            npix,
            ..header(0, 0, 2016.0, 21.0)
        };
        let mut buf = header.to_bytes().to_vec();
        buf[4..8].copy_from_slice(&version.to_le_bytes());
        buf.resize(HEADER_SIZE + 48 * PIXEL_ENTRY_SIZE, 0);
        buf
    }

    fn header_only_file(version: u32, order: u32, nside: u32, npix: u64) -> NamedTempFile {
        write_file(&header_only_bytes(version, order, nside, npix))
    }

    #[test]
    fn test_open_version_1_asks_for_rebuild() {
        let file = header_only_file(1, 1, 2, 48);
        let msg = Catalog::open(file.path())
            .err()
            .expect("expected error")
            .to_string();
        assert!(
            msg.contains("old pixel numbering"),
            "unexpected error: {msg}"
        );
        assert!(msg.contains("forge build-index"), "unexpected error: {msg}");
    }

    #[test]
    fn test_open_rejects_order_above_29() {
        let cases = [
            (30, 1u32 << 30, 12u64 << 60),
            (31, 1u32 << 31, 12u64 << 62),
            (32, 1, 12),
        ];
        for (order, nside, npix) in cases {
            let file = header_only_file(CATALOG_VERSION, order, nside, npix);
            let msg = Catalog::open(file.path())
                .err()
                .expect("expected error")
                .to_string();
            assert!(
                msg.contains(&format!("HEALPix order {order} is above the maximum of 29")),
                "unexpected error: {msg}"
            );
        }
    }

    #[test]
    fn test_open_rejects_epoch_outside_model_range() {
        for epoch in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -0.5,
            4000.5,
            1e300,
        ] {
            let file = build_test_catalog(1, epoch, 21.0, &[]);
            let msg = Catalog::open(file.path())
                .err()
                .expect("expected error")
                .to_string();
            assert!(msg.contains("Catalog epoch"), "epoch {epoch}: {msg}");
        }
        for epoch in [0.0, 1991.25, 4000.0] {
            let file = build_test_catalog(1, epoch, 21.0, &[]);
            assert_eq!(Catalog::open(file.path()).unwrap().header().epoch, epoch);
        }
    }

    #[test]
    fn test_epoch_error_prints_large_epochs_in_exponent_form() {
        let file = build_test_catalog(1, 1e300, 21.0, &[]);
        let msg = Catalog::open(file.path())
            .err()
            .expect("expected error")
            .to_string();
        assert!(msg.contains("Catalog epoch J1e300 is not within"), "{msg}");
        let file = build_test_catalog(1, 5000.0, 21.0, &[]);
        let msg = Catalog::open(file.path())
            .err()
            .expect("expected error")
            .to_string();
        assert!(msg.contains("Catalog epoch J5000 is not within"), "{msg}");
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

    fn hostile_header(order: u32, total_stars: u64, mag_limit: f32) -> Vec<u8> {
        header(order, total_stars, 2016.0, mag_limit)
            .to_bytes()
            .to_vec()
    }

    // An order 0 file: header, 12 pixel entries given as (offset, count), then `records`
    // zeroed star records.
    fn hostile_file(total_stars: u64, entries: &[(u64, u32)], records: usize) -> NamedTempFile {
        let mut buf = hostile_header(0, total_stars, 21.0);
        for &(offset, count) in entries {
            buf.extend_from_slice(&pixel_entry_bytes(offset, count));
        }
        buf.resize(buf.len() + records * STAR_RECORD_SIZE, 0);
        write_file(&buf)
    }

    pub(crate) fn write_file(bytes: &[u8]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(bytes).unwrap();
        file.flush().unwrap();
        file
    }

    fn open_error(file: &NamedTempFile) -> CatalogError {
        match Catalog::open(file.path()) {
            Ok(_) => panic!("hostile catalog opened"),
            Err(e) => e,
        }
    }

    fn assert_header_error(err: CatalogError, text: &str) {
        let msg = err.to_string();
        assert!(matches!(err, CatalogError::InvalidHeader { .. }), "{msg}");
        assert!(msg.contains(text), "{msg}");
    }

    fn assert_table_error(err: CatalogError, text: &str) {
        let msg = err.to_string();
        assert!(matches!(err, CatalogError::CorruptTable { .. }), "{msg}");
        assert!(msg.contains(text), "{msg}");
    }

    // Every entry after the last non-empty one points at the end of the records.
    fn tail(entries: &[(u64, u32)], end: u64) -> Vec<(u64, u32)> {
        let mut all = entries.to_vec();
        all.resize(12, (end, 0));
        all
    }

    #[test]
    fn test_open_rejects_order_29_header_without_table() {
        let file = write_file(&hostile_header(29, 0, 21.0));
        assert_header_error(open_error(&file), "file length");
    }

    #[test]
    fn test_order_29_header_fails_at_open_before_cone_search() {
        use crate::query::cone::{cone_search, ConeSearchParams};
        use std::sync::mpsc;
        use std::time::Duration;

        let file = write_file(&hostile_header(29, 0, 21.0));
        let path = file.path().to_path_buf();
        let (sender, receiver) = mpsc::channel();
        std::thread::spawn(move || {
            let params = ConeSearchParams {
                ra_deg: 10.0,
                dec_deg: 10.0,
                radius_deg: 1.0,
                max_mag: None,
                max_results: None,
                epoch: None,
            };
            let outcome =
                Catalog::open(&path).map(|catalog| cone_search(&catalog, &params).is_ok());
            sender.send(outcome.is_err()).unwrap();
        });
        let open_failed = receiver
            .recv_timeout(Duration::from_secs(10))
            .expect("open and cone_search finished within 10 s without a panic");
        assert!(open_failed, "the order 29 header opened");
    }

    #[test]
    fn test_open_rejects_pixel_offset_near_usize_max() {
        let start = (HEADER_SIZE + 12 * PIXEL_ENTRY_SIZE) as u64;
        let offset = u64::MAX - start - 7;
        let file = hostile_file(1, &tail(&[(offset, 1)], 56), 1);
        assert_table_error(open_error(&file), "pixel 0");
    }

    #[test]
    fn test_open_rejects_count_of_u32_max() {
        let file = hostile_file(1, &tail(&[(0, u32::MAX)], 56), 1);
        assert_table_error(open_error(&file), "pixel 0: 4294967295 records");
    }

    #[test]
    fn test_open_rejects_offset_not_a_multiple_of_56() {
        let file = hostile_file(2, &tail(&[(0, 1), (52, 1)], 112), 2);
        assert_table_error(
            open_error(&file),
            "pixel 1: offset 52 is not a multiple of 56",
        );
    }

    #[test]
    fn test_open_rejects_entries_out_of_pixel_order() {
        let file = hostile_file(3, &tail(&[(56, 1), (112, 1), (0, 1)], 168), 3);
        assert_table_error(open_error(&file), "pixel 0: offset 56 should be 0");
    }

    #[test]
    fn test_open_rejects_overlapping_entries() {
        let file = hostile_file(3, &tail(&[(0, 2), (0, 1)], 168), 3);
        assert_table_error(open_error(&file), "pixel 1: offset 0 should be 112");
    }

    #[test]
    fn test_open_rejects_counts_that_do_not_sum_to_total_stars() {
        let file = hostile_file(2, &tail(&[(0, 1)], 56), 2);
        assert_table_error(open_error(&file), "sum to 1");
    }

    #[test]
    fn test_open_rejects_trailing_bytes() {
        let mut buf = hostile_header(0, 0, 21.0);
        buf.resize(HEADER_SIZE + 12 * PIXEL_ENTRY_SIZE + 8, 0);
        assert_header_error(open_error(&write_file(&buf)), "file length 264 bytes");
    }

    #[test]
    fn test_open_rejects_total_stars_that_overflow_the_file_size() {
        let mut buf = hostile_header(0, u64::MAX, 21.0);
        buf.resize(HEADER_SIZE + 12 * PIXEL_ENTRY_SIZE, 0);
        assert_header_error(open_error(&write_file(&buf)), "total_stars");
    }

    #[test]
    fn test_open_rejects_non_finite_mag_limit() {
        for mag_limit in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            let mut buf = hostile_header(0, 0, mag_limit);
            buf.resize(HEADER_SIZE + 12 * PIXEL_ENTRY_SIZE, 0);
            assert_header_error(open_error(&write_file(&buf)), "mag_limit");
        }
    }

    #[test]
    fn test_open_empty_file_is_too_small() {
        let msg = open_error(&write_file(&[])).to_string();
        assert!(msg.contains("too small: 0 bytes"), "{msg}");
    }

    // The OS message is the error's source, which error reporters print on their own line.
    #[test]
    fn test_io_error_message_leaves_the_os_message_to_source() {
        use std::error::Error;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("absent.bin");
        let err = Catalog::open(&path).err().unwrap();
        assert_eq!(
            err.to_string(),
            format!("Failed to open catalog file {path:?}")
        );
        let source = err.source().expect("the I/O error is the source");
        let io = source.downcast_ref::<std::io::Error>().unwrap();
        assert_eq!(io.kind(), std::io::ErrorKind::NotFound);
    }

    #[test]
    fn test_open_missing_file_is_io_error() {
        let dir = tempfile::tempdir().unwrap();
        let err = Catalog::open(dir.path().join("absent.bin")).err().unwrap();
        assert!(matches!(err, CatalogError::Io { .. }), "{err}");
        assert!(err.to_string().contains("absent.bin"), "{err}");
    }

    fn sample_star() -> StarRecord {
        StarRecord {
            source_id: -9_999_999,
            ra: 359.99999999,
            dec: -89.99999999,
            pmra: 1234.5678,
            pmdec: -8765.4321,
            parallax: 0.001,
            mag: 21.49,
            flags: FLAG_HAS_PROPER_MOTION | FLAG_SOURCE_HIPPARCOS,
            _padding: 0,
        }
    }

    #[test]
    fn test_star_record_le_bytes_layout() {
        let bytes = sample_star().to_le_bytes();
        assert_eq!(bytes[0..8], (-9_999_999i64).to_le_bytes());
        assert_eq!(bytes[8..16], 359.99999999f64.to_le_bytes());
        assert_eq!(bytes[40..48], 0.001f64.to_le_bytes());
        assert_eq!(bytes[48..52], 21.49f32.to_le_bytes());
        assert_eq!(bytes[52..54], 33u16.to_le_bytes());
        assert_eq!(bytes[54..56], [0, 0]);
        assert_eq!(StarRecord::from_le_bytes(&bytes).to_le_bytes(), bytes);
    }

    // stars_in_pixel reads records in place, so the file layout must be the memory layout.
    #[test]
    fn test_star_record_le_bytes_match_memory_layout() {
        let star = sample_star();
        let memory = unsafe {
            std::slice::from_raw_parts((&star as *const StarRecord).cast::<u8>(), STAR_RECORD_SIZE)
        };
        assert_eq!(star.to_le_bytes().as_slice(), memory);
    }

    #[test]
    fn test_catalog_header_bytes_round_trip() {
        let original = header(3, 17, 2016.0, 20.5);
        let bytes = original.to_bytes();
        assert_eq!(&bytes[0..8], b"CCAT\x02\0\0\0");
        assert_eq!(bytes[44..], [0u8; 20]);
        let parsed = CatalogHeader::parse(&bytes).unwrap();
        assert_eq!(
            (parsed.order, parsed.nside, parsed.npix, parsed.total_stars),
            (3, 8, 768, 17)
        );
        assert_eq!((parsed.epoch, parsed.mag_limit), (2016.0, 20.5));
    }

    #[test]
    fn test_pixel_entry_bytes_layout() {
        let bytes = pixel_entry_bytes(0x0102_0304_0506_0708, 0x0a0b_0c0d);
        assert_eq!(
            bytes,
            [8, 7, 6, 5, 4, 3, 2, 1, 0x0d, 0x0c, 0x0b, 0x0a, 0, 0, 0, 0]
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
