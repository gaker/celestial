//! Catalog merging: Hipparcos + Gaia deduplication
//!
//! Hipparcos stars (cross-matched at ingest time) take priority over Gaia.
//! Streams through Gaia, skips duplicates, writes merged output.

use crate::cli::{Cli, MergeArgs};
use crate::gaia::FLAG_SOURCE_HIPPARCOS;
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

const RECORD_SIZE: usize = 56;
const HEADER_SIZE: usize = 24;
const MERGED_MAGIC: &[u8; 4] = b"MERG";
const MERGED_VERSION: u32 = 1;
const EPOCH_J2016: f64 = 2016.0;

#[repr(C)]
struct StarRecord {
    source_id: i64,
    ra: f64,
    dec: f64,
    pmra: f64,
    pmdec: f64,
    parallax: f64,
    mag: f32,
    flags: u16,
    _padding: u16,
}

struct MergeStats {
    hipparcos_count: u64,
    gaia_scanned: u64,
    gaia_skipped: u64,
    gaia_kept: u64,
}

pub fn run(args: &MergeArgs, cli: &Cli) -> anyhow::Result<()> {
    validate_paths(args)?;
    print_plan(args, cli);
    let hipparcos_ids = load_hipparcos_source_ids(args)?;
    let stats = merge_catalogs(args, &hipparcos_ids)?;
    print_stats(&stats);
    validate_output(args)?;
    Ok(())
}

fn validate_paths(args: &MergeArgs) -> anyhow::Result<()> {
    if !args.workdir.exists() {
        anyhow::bail!("Working directory does not exist: {:?}", args.workdir);
    }
    let hip_path = args.workdir.join("hipparcos_ingest.bin");
    if !hip_path.exists() {
        anyhow::bail!("Hipparcos ingest file not found: {:?}", hip_path);
    }
    let gaia_path = args.workdir.join("gaia_ingest.bin");
    if !gaia_path.exists() {
        anyhow::bail!("Gaia ingest file not found: {:?}", gaia_path);
    }
    Ok(())
}

fn print_plan(args: &MergeArgs, cli: &Cli) {
    println!("=== Catalog Merge ===");
    println!("Working directory: {:?}", args.workdir);
    println!("Verbose: {}", cli.verbose);
    println!();
}

fn load_hipparcos_source_ids(args: &MergeArgs) -> anyhow::Result<HashSet<i64>> {
    let path = args.workdir.join("hipparcos_ingest.bin");
    let (count, mut reader) = open_catalog_file(&path)?;
    println!("Loading {} Hipparcos source IDs...", count);
    let mut ids = HashSet::with_capacity(count as usize);
    let mut buf = [0u8; RECORD_SIZE];
    for _ in 0..count {
        reader.read_exact(&mut buf)?;
        let source_id = i64::from_le_bytes(buf[0..8].try_into().unwrap());
        ids.insert(source_id);
    }
    println!("Loaded {} Hipparcos source IDs", ids.len());
    Ok(ids)
}

fn open_catalog_file(path: &Path) -> anyhow::Result<(u64, BufReader<File>)> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = [0u8; HEADER_SIZE];
    reader.read_exact(&mut header)?;
    let count = u64::from_le_bytes(header[8..16].try_into().unwrap());
    Ok((count, reader))
}

fn merge_catalogs(args: &MergeArgs, hipparcos_ids: &HashSet<i64>) -> anyhow::Result<MergeStats> {
    let output_path = args.workdir.join("merged.bin");
    let temp_path = output_path.with_extension("bin.tmp");
    let file = File::create(&temp_path)?;
    let mut writer = BufWriter::new(file);
    write_placeholder_header(&mut writer)?;
    let hip_count = write_hipparcos_records(args, &mut writer)?;
    let (gaia_scanned, gaia_skipped, gaia_kept) =
        stream_gaia_records(args, &mut writer, hipparcos_ids)?;
    let total = hip_count + gaia_kept;
    finalize_output(&mut writer, total)?;
    drop(writer);
    fs::rename(&temp_path, &output_path)?;
    println!("Written {} stars to {:?}", total, output_path);
    Ok(MergeStats {
        hipparcos_count: hip_count,
        gaia_scanned,
        gaia_skipped,
        gaia_kept,
    })
}

fn write_placeholder_header<W: Write>(writer: &mut W) -> anyhow::Result<()> {
    writer.write_all(&[0u8; HEADER_SIZE])?;
    Ok(())
}

fn write_hipparcos_records(args: &MergeArgs, writer: &mut BufWriter<File>) -> anyhow::Result<u64> {
    let path = args.workdir.join("hipparcos_ingest.bin");
    let (count, mut reader) = open_catalog_file(&path)?;
    println!("Writing {} Hipparcos records...", count);
    copy_records(&mut reader, writer, count)?;
    Ok(count)
}

fn copy_records<R: Read, W: Write>(
    reader: &mut R,
    writer: &mut W,
    count: u64,
) -> anyhow::Result<()> {
    let mut buf = vec![0u8; 64 * 1024];
    let total_bytes = count * RECORD_SIZE as u64;
    let mut remaining = total_bytes;
    while remaining > 0 {
        let to_read = remaining.min(buf.len() as u64) as usize;
        reader.read_exact(&mut buf[..to_read])?;
        writer.write_all(&buf[..to_read])?;
        remaining -= to_read as u64;
    }
    Ok(())
}

fn stream_gaia_records(
    args: &MergeArgs,
    writer: &mut BufWriter<File>,
    hipparcos_ids: &HashSet<i64>,
) -> anyhow::Result<(u64, u64, u64)> {
    let path = args.workdir.join("gaia_ingest.bin");
    let (count, mut reader) = open_catalog_file(&path)?;
    println!("Streaming {} Gaia records...", count);
    let (skipped, kept) = filter_gaia_records(&mut reader, writer, count, hipparcos_ids)?;
    Ok((count, skipped, kept))
}

fn filter_gaia_records<R: Read, W: Write>(
    reader: &mut R,
    writer: &mut W,
    count: u64,
    hipparcos_ids: &HashSet<i64>,
) -> anyhow::Result<(u64, u64)> {
    let mut buf = [0u8; RECORD_SIZE];
    let mut skipped = 0u64;
    let mut kept = 0u64;
    for _ in 0..count {
        reader.read_exact(&mut buf)?;
        let source_id = i64::from_le_bytes(buf[0..8].try_into().unwrap());
        if hipparcos_ids.contains(&source_id) {
            skipped += 1;
        } else {
            writer.write_all(&buf)?;
            kept += 1;
        }
    }
    Ok((skipped, kept))
}

fn finalize_output(writer: &mut BufWriter<File>, total_stars: u64) -> anyhow::Result<()> {
    writer.flush()?;
    let file = writer.get_mut();
    file.seek(SeekFrom::Start(0))?;
    write_final_header(file, total_stars)?;
    Ok(())
}

fn write_final_header<W: Write + Seek>(writer: &mut W, total: u64) -> anyhow::Result<()> {
    writer.write_all(MERGED_MAGIC)?;
    writer.write_all(&MERGED_VERSION.to_le_bytes())?;
    writer.write_all(&total.to_le_bytes())?;
    writer.write_all(&(EPOCH_J2016 as f32).to_le_bytes())?;
    writer.write_all(&[0u8; 4])?;
    Ok(())
}

fn print_stats(stats: &MergeStats) {
    println!();
    println!("=== Merge Statistics ===");
    println!("Hipparcos stars loaded: {}", stats.hipparcos_count);
    println!("Gaia stars scanned: {}", stats.gaia_scanned);
    println!("Gaia stars skipped (duplicates): {}", stats.gaia_skipped);
    println!("Gaia stars kept: {}", stats.gaia_kept);
    println!("Total merged: {}", stats.hipparcos_count + stats.gaia_kept);
}

fn validate_output(args: &MergeArgs) -> anyhow::Result<()> {
    let path = args.workdir.join("merged.bin");
    let (count, mut reader) = open_catalog_file(&path)?;
    println!();
    println!("=== Validation (sample records) ===");
    println!("Total stars in merged catalog: {}", count);
    print_hipparcos_samples(&mut reader, count)?;
    let path = args.workdir.join("merged.bin");
    let (_, mut reader) = open_catalog_file(&path)?;
    print_gaia_samples(&mut reader, count)?;
    Ok(())
}

fn print_hipparcos_samples<R: Read>(reader: &mut R, count: u64) -> anyhow::Result<()> {
    println!();
    println!("Hipparcos-flagged records:");
    let mut found = 0;
    let mut buf = [0u8; RECORD_SIZE];
    for i in 0..count.min(10000) {
        reader.read_exact(&mut buf)?;
        let record = parse_record(&buf);
        if record.flags & FLAG_SOURCE_HIPPARCOS != 0 {
            print_record(i, &record);
            found += 1;
            if found >= 3 {
                break;
            }
        }
    }
    Ok(())
}

fn print_gaia_samples<R: Read>(reader: &mut R, count: u64) -> anyhow::Result<()> {
    println!();
    println!("Gaia records (no Hipparcos flag):");
    let mut found = 0;
    let mut buf = [0u8; RECORD_SIZE];
    for i in 0..count.min(100000) {
        reader.read_exact(&mut buf)?;
        let record = parse_record(&buf);
        if record.flags & FLAG_SOURCE_HIPPARCOS == 0 {
            print_record(i, &record);
            found += 1;
            if found >= 3 {
                break;
            }
        }
    }
    Ok(())
}

fn parse_record(buf: &[u8; RECORD_SIZE]) -> StarRecord {
    StarRecord {
        source_id: i64::from_le_bytes(buf[0..8].try_into().unwrap()),
        ra: f64::from_le_bytes(buf[8..16].try_into().unwrap()),
        dec: f64::from_le_bytes(buf[16..24].try_into().unwrap()),
        pmra: f64::from_le_bytes(buf[24..32].try_into().unwrap()),
        pmdec: f64::from_le_bytes(buf[32..40].try_into().unwrap()),
        parallax: f64::from_le_bytes(buf[40..48].try_into().unwrap()),
        mag: f32::from_le_bytes(buf[48..52].try_into().unwrap()),
        flags: u16::from_le_bytes(buf[52..54].try_into().unwrap()),
        _padding: 0,
    }
}

fn print_record(index: u64, record: &StarRecord) {
    println!(
        "  [{}] source_id={}, RA={:.6}, Dec={:.6}, mag={:.2}, flags=0x{:04x}",
        index, record.source_id, record.ra, record.dec, record.mag, record.flags
    );
}
