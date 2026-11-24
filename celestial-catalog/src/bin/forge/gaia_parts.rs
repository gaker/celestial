// Part files of ingest-gaia, one per input CSV, and their concatenation into gaia_ingest.bin.
// A part file has the gaia_ingest.bin layout: an intermediate header of kind Gaia, which
// carries the magic, the record count, the epoch and the magnitude limit, then the records.

use crate::gaia::{GaiaParser, RowCounts, GAIA_DR3_EPOCH};
use crate::intermediate::{self, Header, Kind, HEADER_SIZE, RECORD_SIZE};
use anyhow::Context;
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const PART_PREFIX: &str = "gaia_part_";

// Maps GaiaSource_000000-003111.csv.gz to gaia_part_000000-003111.bin in `output_dir`.
pub(crate) fn part_path(input: &Path, output_dir: &Path) -> Option<PathBuf> {
    let filename = input.file_name()?.to_str()?;
    let numeric = filename
        .strip_suffix(".csv.gz")?
        .strip_prefix("GaiaSource_")?;
    Some(output_dir.join(format!("{PART_PREFIX}{numeric}.bin")))
}

pub(crate) fn part_exists(input: &Path, output_dir: &Path) -> bool {
    part_path(input, output_dir)
        .and_then(|path| fs::metadata(path).ok())
        .is_some_and(|meta| meta.len() > 0)
}

// Writes the parser's rows to `output` through a temporary file, so an interrupted run leaves
// no part file that a later run would take as complete.
pub(crate) fn write_part<R: BufRead>(
    parser: GaiaParser<R>,
    output: &Path,
    mag_limit: f32,
) -> anyhow::Result<RowCounts> {
    let temp_path = output.with_extension("bin.tmp");
    let counts = write_records(parser, &temp_path, mag_limit).and_then(|counts| {
        fs::rename(&temp_path, output)?;
        Ok(counts)
    });
    if counts.is_err() {
        intermediate::remove_temp_file(&temp_path);
    }
    counts
}

fn write_records<R: BufRead>(
    mut parser: GaiaParser<R>,
    path: &Path,
    mag_limit: f32,
) -> anyhow::Result<RowCounts> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(&[0u8; HEADER_SIZE])?;
    while let Some(star) = parser.next() {
        let star = star.with_context(|| format!("Error at data row {}", parser.counts().read))?;
        intermediate::write_record(&mut writer, &star)?;
    }
    let counts = parser.counts();
    let header = gaia_header(counts.kept, mag_limit);
    intermediate::finish(writer, &header)?;
    Ok(counts)
}

fn gaia_header(records: u64, mag_limit: f32) -> Header {
    Header {
        kind: Kind::Gaia,
        records,
        epoch: GAIA_DR3_EPOCH,
        mag_limit,
    }
}

// Writes gaia_ingest.bin from every part file in `output_dir`, then deletes the part files.
// Any failure leaves the part files in place and removes the partial output.
pub(crate) fn concatenate(output_dir: &Path, mag_limit: f32) -> anyhow::Result<()> {
    let parts = collect_part_files(output_dir)?;
    if parts.is_empty() {
        anyhow::bail!("No part files found in {:?}", output_dir);
    }
    println!("\n=== Concatenating {} part files ===", parts.len());
    let final_path = output_dir.join(Kind::Gaia.file_name());
    let temp_path = final_path.with_extension("bin.tmp");
    let written = write_concatenated(&parts, &temp_path, mag_limit).and_then(|total| {
        fs::rename(&temp_path, &final_path)?;
        Ok(total)
    });
    if written.is_err() {
        intermediate::remove_temp_file(&temp_path);
    }
    println!("Written {} stars to {:?}", written?, final_path);
    delete_part_files(&parts)
}

fn collect_part_files(output_dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(output_dir)? {
        let path = entry?.path();
        if is_part_file(&path) {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

fn is_part_file(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.starts_with(PART_PREFIX) && n.ends_with(".bin"))
}

fn write_concatenated(parts: &[PathBuf], output: &Path, mag_limit: f32) -> anyhow::Result<u64> {
    let total = total_records(parts, mag_limit)?;
    let mut writer = BufWriter::new(File::create(output)?);
    writer.write_all(&gaia_header(total, mag_limit).to_bytes())?;
    let mut copied = 0u64;
    for path in parts {
        copied = copied.saturating_add(copy_part(&mut writer, path)?);
    }
    // A part file that changed between the two reads would leave the header count wrong.
    if copied != total {
        anyhow::bail!("Copied {copied} records, but the part headers claimed {total}");
    }
    writer.flush()?;
    writer.into_inner()?.sync_all()?;
    Ok(total)
}

// Reads only the headers, so the part files are not all open at once.
fn total_records(parts: &[PathBuf], mag_limit: f32) -> anyhow::Result<u64> {
    let mut total = 0u64;
    for path in parts {
        let records = read_part_header(path, mag_limit)?.records;
        total = total
            .checked_add(records)
            .with_context(|| format!("The part files hold more than {} records", u64::MAX))?;
    }
    Ok(total)
}

// open_path checks the magic, the version and the claimed count against the file length.
// Part files written before they had a header fail the magic check.
fn read_part_header(path: &Path, mag_limit: f32) -> anyhow::Result<Header> {
    let (header, _) = intermediate::open_path(path, Kind::Gaia).with_context(|| {
        format!(
            "{path:?} is not a part file of this version; delete the {PART_PREFIX}*.bin files \
             and rerun"
        )
    })?;
    if header.mag_limit != mag_limit || header.epoch != GAIA_DR3_EPOCH {
        anyhow::bail!(
            "Part file {path:?} was written with mag_limit {} at epoch J{}, not {mag_limit} at \
             J{GAIA_DR3_EPOCH}; delete the {PART_PREFIX}*.bin files and rerun",
            header.mag_limit,
            header.epoch
        );
    }
    Ok(header)
}

// Returns the number of records copied.
fn copy_part<W: Write>(writer: &mut W, path: &Path) -> anyhow::Result<u64> {
    let (header, reader) = intermediate::open_path(path, Kind::Gaia)?;
    // open_path checked that the file holds this many records, so the product fits.
    let bytes = header.records * RECORD_SIZE as u64;
    let copied = std::io::copy(&mut reader.take(bytes), writer)?;
    if copied != bytes {
        anyhow::bail!("Part file {path:?}: expected {bytes} bytes of records, read {copied}");
    }
    Ok(header.records)
}

fn delete_part_files(part_files: &[PathBuf]) -> anyhow::Result<()> {
    for path in part_files {
        fs::remove_file(path).with_context(|| format!("Failed to delete part file {path:?}"))?;
    }
    println!("Deleted {} part files", part_files.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gaia::tests::{csv, row};

    fn write_csv_part(path: &Path, rows: &[String]) -> RowCounts {
        let text = csv(rows);
        let parser = GaiaParser::new(std::io::Cursor::new(text), 18.5).unwrap();
        write_part(parser, path, 18.5).unwrap()
    }

    fn three_rows() -> [String; 3] {
        [
            row("1", "10.0", "20.0", "12.0", 31),
            row("2", "11.0", "21.0", "19.0", 31),
            row("3", "12.0", "22.0", "12.0", 31).replace("\"False\"", "\"True\""),
        ]
    }

    #[test]
    fn test_write_part_file_counts_rows_read_and_rows_kept() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("part.bin");
        let counts = write_csv_part(&path, &three_rows());
        assert_eq!((counts.read, counts.kept), (3, 2));
        let (header, _) = intermediate::open_path(&path, Kind::Gaia).unwrap();
        assert_eq!(header, gaia_header(2, 18.5));
        assert!(!path.with_extension("bin.tmp").exists());
    }

    #[test]
    fn test_concatenated_parts_open_as_gaia_ingest() {
        let dir = tempfile::tempdir().unwrap();
        for name in ["gaia_part_000.bin", "gaia_part_001.bin"] {
            write_csv_part(&dir.path().join(name), &three_rows());
        }
        concatenate(dir.path(), 18.5).unwrap();
        let (header, mut reader) = intermediate::open(dir.path(), Kind::Gaia).unwrap();
        assert_eq!(header, gaia_header(4, 18.5));
        let star = intermediate::read_record(&mut reader).unwrap();
        assert_eq!(
            (star.source_id, star.ra, star.dec, star.mag),
            (1, 10.0, 20.0, 12.0)
        );
        assert!(!dir.path().join("gaia_part_000.bin").exists());
    }

    fn assert_concatenate_fails_and_keeps_parts(dir: &Path, names: &[&str], expected: &str) {
        let err = format!("{:#}", concatenate(dir, 18.5).unwrap_err());
        assert!(err.contains(expected), "{err}");
        assert!(!dir.join("gaia_ingest.bin").exists());
        assert!(!dir.join("gaia_ingest.bin.tmp").exists());
        for name in names {
            assert!(dir.join(name).exists(), "{name} was deleted");
        }
    }

    // The format before part files had a header was an 8-byte record count. Two parts that
    // claim 2^63 records each overflowed the sum.
    #[test]
    fn test_part_without_header_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let names = ["gaia_part_000.bin", "gaia_part_001.bin"];
        for name in names {
            fs::write(dir.path().join(name), (1u64 << 63).to_le_bytes()).unwrap();
        }
        assert_concatenate_fails_and_keeps_parts(dir.path(), &names, "shorter than the 32-byte");
    }

    #[test]
    fn test_part_whose_count_disagrees_with_its_length_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let names = ["gaia_part_000.bin", "gaia_part_001.bin"];
        write_csv_part(&dir.path().join(names[0]), &three_rows());
        let mut bytes = gaia_header(u64::MAX / 2, 18.5).to_bytes().to_vec();
        bytes.extend_from_slice(&[0; RECORD_SIZE]);
        fs::write(dir.path().join(names[1]), bytes).unwrap();
        assert_concatenate_fails_and_keeps_parts(dir.path(), &names, "claims 9223372036854775807");
    }

    #[test]
    fn test_part_from_another_mag_limit_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let names = ["gaia_part_000.bin"];
        fs::write(dir.path().join(names[0]), gaia_header(0, 15.0).to_bytes()).unwrap();
        assert_concatenate_fails_and_keeps_parts(dir.path(), &names, "mag_limit 15");
    }

    // A part file written before part files had a header: an 8-byte count, then records.
    #[test]
    fn test_old_format_part_names_the_files_to_delete() {
        let dir = tempfile::tempdir().unwrap();
        let names = ["gaia_part_000.bin"];
        let mut bytes = 1u64.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0; RECORD_SIZE]);
        fs::write(dir.path().join(names[0]), bytes).unwrap();
        let expected = "delete the gaia_part_*.bin files";
        assert_concatenate_fails_and_keeps_parts(dir.path(), &names, expected);
    }
}
