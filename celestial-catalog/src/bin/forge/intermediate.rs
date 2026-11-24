// Intermediate files passed between forge stages: hipparcos_ingest.bin, gaia_ingest.bin,
// merged.bin and hipparcos_shared_matches.bin. Each is a 32-byte header followed by `records`
// little-endian star records. In hipparcos_shared_matches.bin only each record's source_id is
// set: a Gaia DR3 source that more than one ingested Hipparcos star cross-matches.
//
// | Bytes  | Field     | Type            |
// |--------|-----------|-----------------|
// | 0..4   | magic     | `HIPP`, `GAIA`, `MERG` or `HSHM` |
// | 4..8   | version   | u32             |
// | 8..16  | records   | u64             |
// | 16..24 | epoch     | f64 Julian year |
// | 24..28 | mag_limit | f32             |
// | 28..32 | reserved  | zero            |

use anyhow::Context;
use celestial_catalog::query::catalog::StarRecord;
use celestial_catalog::query::validate::check_epoch;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

pub(crate) const RECORD_SIZE: usize = std::mem::size_of::<StarRecord>();
pub(crate) const HEADER_SIZE: usize = 32;
// Version 3 moved the epoch into its own f64 field. Versions 1 and 2 had a 24-byte header
// whose f32 slot held mag_limit in the ingest files and the epoch in merged.bin.
const VERSION: u32 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Hipparcos,
    Gaia,
    Merged,
    SharedMatches,
}

impl Kind {
    fn magic(self) -> &'static [u8; 4] {
        match self {
            Kind::Hipparcos => b"HIPP",
            Kind::Gaia => b"GAIA",
            Kind::Merged => b"MERG",
            Kind::SharedMatches => b"HSHM",
        }
    }

    pub(crate) fn file_name(self) -> &'static str {
        match self {
            Kind::Hipparcos => "hipparcos_ingest.bin",
            Kind::Gaia => "gaia_ingest.bin",
            Kind::Merged => "merged.bin",
            Kind::SharedMatches => "hipparcos_shared_matches.bin",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Header {
    pub(crate) kind: Kind,
    pub(crate) records: u64,
    pub(crate) epoch: f64,
    pub(crate) mag_limit: f32,
}

impl Header {
    pub(crate) fn to_bytes(self) -> [u8; HEADER_SIZE] {
        let fields: [&[u8]; 5] = [
            self.kind.magic(),
            &VERSION.to_le_bytes(),
            &self.records.to_le_bytes(),
            &self.epoch.to_le_bytes(),
            &self.mag_limit.to_le_bytes(),
        ];
        let mut out = [0u8; HEADER_SIZE];
        for (slot, byte) in out.iter_mut().zip(fields.into_iter().flatten()) {
            *slot = *byte;
        }
        out
    }

    // Parses the header of a `kind` file that is `file_len` bytes long. The file must hold
    // exactly `records` records after the header.
    pub(crate) fn parse(
        bytes: &[u8; HEADER_SIZE],
        kind: Kind,
        file_len: u64,
    ) -> anyhow::Result<Self> {
        check_magic_and_version(bytes, kind)?;
        let header = Self {
            kind,
            records: u64::from_le_bytes(field(bytes, 8)),
            epoch: f64::from_le_bytes(field(bytes, 16)),
            mag_limit: f32::from_le_bytes(field(bytes, 24)),
        };
        // Every later stage copies the epoch into the catalog, so an epoch Catalog::open
        // refuses fails here, before anything is built from the file.
        check_epoch(header.epoch)?;
        if !header.mag_limit.is_finite() {
            anyhow::bail!("mag_limit {} is not finite", header.mag_limit);
        }
        header.check_length(file_len)?;
        Ok(header)
    }

    fn check_length(&self, file_len: u64) -> anyhow::Result<()> {
        let expected = self
            .records
            .checked_mul(RECORD_SIZE as u64)
            .and_then(|len| len.checked_add(HEADER_SIZE as u64));
        if expected != Some(file_len) {
            anyhow::bail!(
                "the header claims {} records, but the file is {file_len} bytes; it should be \
                 {HEADER_SIZE} header bytes plus {RECORD_SIZE} bytes per record",
                self.records
            );
        }
        Ok(())
    }
}

fn check_magic_and_version(bytes: &[u8; HEADER_SIZE], kind: Kind) -> anyhow::Result<()> {
    let magic: [u8; 4] = field(bytes, 0);
    if &magic != kind.magic() {
        anyhow::bail!(
            "magic {:?} is not {:?}",
            String::from_utf8_lossy(&magic),
            String::from_utf8_lossy(kind.magic())
        );
    }
    let version = u32::from_le_bytes(field(bytes, 4));
    if version != VERSION {
        anyhow::bail!(
            "intermediate format version {version} is not {VERSION}; rerun the earlier forge \
             stages to regenerate it"
        );
    }
    Ok(())
}

fn field<const N: usize>(bytes: &[u8; HEADER_SIZE], at: usize) -> [u8; N] {
    std::array::from_fn(|i| bytes[at + i])
}

// Opens a `kind` file in `dir`, checks its header against the file length and returns the
// reader positioned at the first record.
pub(crate) fn open(dir: &Path, kind: Kind) -> anyhow::Result<(Header, BufReader<File>)> {
    open_path(&dir.join(kind.file_name()), kind)
}

// As `open`, for a `kind` file at `path` under any name.
pub(crate) fn open_path(path: &Path, kind: Kind) -> anyhow::Result<(Header, BufReader<File>)> {
    let context = || format!("Failed to read intermediate file {path:?}");
    let file = File::open(path).with_context(context)?;
    let file_len = file.metadata().with_context(context)?.len();
    let mut reader = BufReader::new(file);
    let mut bytes = [0u8; HEADER_SIZE];
    reader
        .read_exact(&mut bytes)
        .with_context(|| format!("{path:?} is shorter than the {HEADER_SIZE}-byte header"))?;
    let header = Header::parse(&bytes, kind, file_len).with_context(context)?;
    Ok((header, reader))
}

pub(crate) fn read_record<R: Read>(reader: &mut R) -> std::io::Result<StarRecord> {
    let mut bytes = [0u8; RECORD_SIZE];
    reader.read_exact(&mut bytes)?;
    Ok(StarRecord::from_le_bytes(&bytes))
}

pub(crate) fn write_record<W: Write>(writer: &mut W, record: &StarRecord) -> std::io::Result<()> {
    writer.write_all(&record.to_le_bytes())
}

// Writes `path` through `<path>.tmp` and a rename, so a stage that fails or is killed leaves
// no partial file under the final name. The temporary file is removed on error.
pub(crate) fn write_via_temp<T>(
    path: &Path,
    write: impl FnOnce(&mut BufWriter<File>) -> anyhow::Result<T>,
) -> anyhow::Result<T> {
    let temp_path = path.with_extension("bin.tmp");
    let result = File::create(&temp_path)
        .map_err(anyhow::Error::from)
        .and_then(|file| {
            let mut writer = BufWriter::new(file);
            let value = write(&mut writer)?;
            writer.into_inner()?.sync_all()?;
            std::fs::rename(&temp_path, path)?;
            Ok(value)
        });
    if result.is_err() {
        remove_temp_file(&temp_path);
    }
    result
}

// Writes `header` over the placeholder at the start of the file and syncs it. A writer
// fills in the record count this way once the records are written.
pub(crate) fn finish(writer: BufWriter<File>, header: &Header) -> anyhow::Result<()> {
    let mut file = writer.into_inner()?;
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&header.to_bytes())?;
    file.sync_all()?;
    Ok(())
}

pub(crate) fn zeroed_record() -> StarRecord {
    StarRecord::from_le_bytes(&[0; RECORD_SIZE])
}

// Removes a temporary file after a failed stage. A file that was never created is not an
// error; any other failure is reported and does not replace the stage's own error.
pub(crate) fn remove_temp_file(path: &Path) {
    match std::fs::remove_file(path) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => eprintln!("Warning: failed to remove temporary file {path:?}: {e}"),
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) fn star(source_id: i64, ra: f64, dec: f64, mag: f32, flags: u16) -> StarRecord {
        let mut record = zeroed_record();
        record.source_id = source_id;
        record.ra = ra;
        record.dec = dec;
        record.mag = mag;
        record.flags = flags;
        record
    }

    pub(crate) fn write_file(dir: &Path, header: Header, stars: &[StarRecord]) {
        let mut bytes = header.to_bytes().to_vec();
        for star in stars {
            write_record(&mut bytes, star).unwrap();
        }
        std::fs::write(dir.join(header.kind.file_name()), bytes).unwrap();
    }

    // Writes hipparcos_ingest.bin and the empty hipparcos_shared_matches.bin that
    // ingest-hipparcos writes beside it when no Gaia source is matched by more than one
    // Hipparcos star. merge requires both.
    pub(crate) fn write_hipparcos(dir: &Path, hipparcos: Header, stars: &[StarRecord]) {
        write_file(dir, hipparcos, stars);
        let shared = Header {
            kind: Kind::SharedMatches,
            records: 0,
            ..hipparcos
        };
        write_file(dir, shared, &[]);
    }

    pub(crate) fn header(kind: Kind, records: u64, epoch: f64, mag_limit: f32) -> Header {
        Header {
            kind,
            records,
            epoch,
            mag_limit,
        }
    }

    fn open_error(dir: &Path, kind: Kind) -> String {
        format!("{:#}", open(dir, kind).unwrap_err())
    }

    #[test]
    fn test_header_byte_layout() {
        let bytes = header(Kind::Merged, 0x0102, 2016.0, 18.5).to_bytes();
        assert_eq!(&bytes[0..4], b"MERG");
        assert_eq!(bytes[4..8], 3u32.to_le_bytes());
        assert_eq!(bytes[8..16], 0x0102u64.to_le_bytes());
        assert_eq!(bytes[16..24], 2016.0f64.to_le_bytes());
        assert_eq!(bytes[24..28], 18.5f32.to_le_bytes());
        assert_eq!(bytes[28..32], [0; 4]);
    }

    #[test]
    fn test_round_trip_through_a_file() {
        let dir = tempfile::tempdir().unwrap();
        let stars = [
            star(1, 10.0, 20.0, 5.5, 3),
            star(-2, 359.5, -89.0, 17.25, 32),
        ];
        let written = header(Kind::Gaia, 2, 2016.0, 18.5);
        write_file(dir.path(), written, &stars);
        let (read, mut reader) = open(dir.path(), Kind::Gaia).unwrap();
        assert_eq!(read, written);
        for expected in &stars {
            assert_eq!(
                read_record(&mut reader).unwrap().to_le_bytes(),
                expected.to_le_bytes()
            );
        }
    }

    #[test]
    fn test_open_rejects_another_kind_of_file() {
        let dir = tempfile::tempdir().unwrap();
        write_file(dir.path(), header(Kind::Hipparcos, 0, 2016.0, 7.0), &[]);
        std::fs::rename(
            dir.path().join("hipparcos_ingest.bin"),
            dir.path().join("gaia_ingest.bin"),
        )
        .unwrap();
        let err = open_error(dir.path(), Kind::Gaia);
        assert!(err.contains("magic \"HIPP\" is not \"GAIA\""), "{err}");
    }

    #[test]
    fn test_open_rejects_old_version() {
        let dir = tempfile::tempdir().unwrap();
        let mut bytes = header(Kind::Merged, 0, 2016.0, 18.5).to_bytes();
        bytes[4..8].copy_from_slice(&1u32.to_le_bytes());
        std::fs::write(dir.path().join("merged.bin"), bytes).unwrap();
        let err = open_error(dir.path(), Kind::Merged);
        assert!(err.contains("version 1 is not 3"), "{err}");
    }

    #[test]
    fn test_open_rejects_count_without_records() {
        let dir = tempfile::tempdir().unwrap();
        write_file(dir.path(), header(Kind::Merged, 10, 2016.0, 18.5), &[]);
        let err = open_error(dir.path(), Kind::Merged);
        assert!(
            err.contains("claims 10 records, but the file is 32 bytes"),
            "{err}"
        );
    }

    #[test]
    fn test_open_rejects_count_whose_size_overflows() {
        let dir = tempfile::tempdir().unwrap();
        write_file(
            dir.path(),
            header(Kind::Hipparcos, 1 << 60, 2016.0, 7.0),
            &[],
        );
        let err = open_error(dir.path(), Kind::Hipparcos);
        assert!(err.contains("claims 1152921504606846976 records"), "{err}");
    }

    #[test]
    fn test_open_rejects_trailing_partial_record() {
        let dir = tempfile::tempdir().unwrap();
        let mut bytes = header(Kind::Merged, 1, 2016.0, 18.5).to_bytes().to_vec();
        bytes.extend_from_slice(&[0; RECORD_SIZE + 1]);
        std::fs::write(dir.path().join("merged.bin"), bytes).unwrap();
        let err = open_error(dir.path(), Kind::Merged);
        assert!(
            err.contains("claims 1 records, but the file is 89 bytes"),
            "{err}"
        );
    }

    #[test]
    fn test_open_rejects_short_header() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("merged.bin"), b"MERG").unwrap();
        let err = open_error(dir.path(), Kind::Merged);
        assert!(err.contains("shorter than the 32-byte header"), "{err}");
    }

    #[test]
    fn test_open_rejects_non_finite_mag_limit() {
        let dir = tempfile::tempdir().unwrap();
        write_file(dir.path(), header(Kind::Gaia, 0, 2016.0, f32::NAN), &[]);
        let err = open_error(dir.path(), Kind::Gaia);
        assert!(err.contains("is not finite"), "{err}");
    }
}
