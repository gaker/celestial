//! Gaia DR3 catalog ingestion

use crate::cli::{Cli, IngestGaiaArgs};
use crate::gaia::{GaiaParser, GaiaStarRaw};
use anyhow::Context;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

const RECORD_SIZE: usize = 56;
const PART_HEADER_SIZE: usize = 8;
const FINAL_MAGIC: &[u8; 4] = b"GAIA";
const FINAL_VERSION: u32 = 1;

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

struct FileStats {
    kept: u64,
    scanned: u64,
}

pub fn run(args: &IngestGaiaArgs, cli: &Cli) -> anyhow::Result<()> {
    validate_args(args)?;
    let files = find_gzipped_csvs(&args.path)?;
    if files.is_empty() {
        anyhow::bail!("No .csv.gz files found in {:?}", args.path);
    }
    print_plan(args, files.len());
    configure_thread_pool(args.threads)?;
    let files_to_process = filter_already_processed(&files, &args.output)?;
    println!(
        "Files to process: {} (skipping {} already done)",
        files_to_process.len(),
        files.len() - files_to_process.len()
    );
    if !files_to_process.is_empty() {
        process_files(&files_to_process, args, cli)?;
    }
    if args.no_concat {
        println!("Skipping concatenation (--no-concat)");
        return Ok(());
    }
    concatenate_part_files(&args.output, args.mag_limit)
}

fn validate_args(args: &IngestGaiaArgs) -> anyhow::Result<()> {
    if !args.path.exists() {
        anyhow::bail!("Gaia directory does not exist: {:?}", args.path);
    }
    if !args.path.is_dir() {
        anyhow::bail!("Gaia path is not a directory: {:?}", args.path);
    }
    if !args.output.exists() {
        fs::create_dir_all(&args.output)?;
    }
    Ok(())
}

fn find_gzipped_csvs(dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        if is_gzipped_csv(&path) {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

fn is_gzipped_csv(path: &Path) -> bool {
    path.extension().is_some_and(|e| e == "gz")
        && path
            .file_stem()
            .and_then(|s| s.to_str())
            .is_some_and(|s| s.ends_with(".csv"))
}

fn filter_already_processed(files: &[PathBuf], output_dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    Ok(files
        .iter()
        .filter(|f| !part_file_exists(f, output_dir))
        .cloned()
        .collect())
}

fn part_file_exists(input: &Path, output_dir: &Path) -> bool {
    let Some(part_path) = compute_part_path(input, output_dir) else {
        return false;
    };
    part_path.exists()
        && fs::metadata(&part_path)
            .map(|m| m.len() > 0)
            .unwrap_or(false)
}

fn extract_numeric_portion(filename: &str) -> Option<&str> {
    let stem = filename.strip_suffix(".csv.gz")?;
    stem.strip_prefix("GaiaSource_")
}

fn compute_part_path(input: &Path, output_dir: &Path) -> Option<PathBuf> {
    let filename = input.file_name()?.to_str()?;
    let numeric = extract_numeric_portion(filename)?;
    Some(output_dir.join(format!("gaia_part_{}.bin", numeric)))
}

fn print_plan(args: &IngestGaiaArgs, file_count: usize) {
    println!("=== Gaia DR3 Ingestion ===");
    println!("Input directory: {:?}", args.path);
    println!("Files found: {}", file_count);
    println!("Magnitude limit: {:.1}", args.mag_limit);
    println!("Output directory: {:?}", args.output);
    println!("Threads: {}", resolve_threads(args.threads));
    println!();
}

fn resolve_threads(threads: usize) -> usize {
    if threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        threads
    }
}

fn configure_thread_pool(threads: usize) -> anyhow::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(resolve_threads(threads))
        .build_global()
        .ok();
    Ok(())
}

fn process_files(files: &[PathBuf], args: &IngestGaiaArgs, cli: &Cli) -> anyhow::Result<()> {
    let pb = create_progress_bar(files.len() as u64);
    let total_kept = AtomicU64::new(0);
    let total_scanned = AtomicU64::new(0);
    let errors = AtomicU64::new(0);
    let failed_files: Mutex<Vec<(PathBuf, String)>> = Mutex::new(Vec::new());

    files.par_iter().for_each(|file| {
        let result = process_single_file(file, args, cli, &total_kept, &total_scanned);
        if let Err(e) = result {
            eprintln!(
                "\nWarning: {:?}: {}",
                file.file_name().unwrap_or_default(),
                e
            );
            errors.fetch_add(1, Ordering::Relaxed);
            if let Ok(mut failed) = failed_files.lock() {
                failed.push((file.clone(), e.to_string()));
            }
        }
        pb.inc(1);
    });

    pb.finish_with_message("Done");
    write_failed_log(&args.output, &failed_files)?;
    print_summary(files.len(), &total_kept, &total_scanned, &errors)
}

fn write_failed_log(
    output_dir: &Path,
    failed: &Mutex<Vec<(PathBuf, String)>>,
) -> anyhow::Result<()> {
    let failed = failed
        .lock()
        .map_err(|e| anyhow::anyhow!("Lock error: {}", e))?;
    if failed.is_empty() {
        return Ok(());
    }
    let log_path = output_dir.join("failed_files.log");
    let mut file = File::create(&log_path)?;
    for (path, error) in failed.iter() {
        writeln!(file, "{}\t{}", path.display(), error)?;
    }
    println!("Failed files logged to: {:?}", log_path);
    Ok(())
}

fn create_progress_bar(total: u64) -> ProgressBar {
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    pb
}

fn process_single_file(
    path: &Path,
    args: &IngestGaiaArgs,
    cli: &Cli,
    total_kept: &AtomicU64,
    total_scanned: &AtomicU64,
) -> anyhow::Result<()> {
    let output_path =
        compute_part_path(path, &args.output).context("Failed to compute output path")?;
    let stats = stream_and_filter(path, &output_path, args.mag_limit)?;
    total_kept.fetch_add(stats.kept, Ordering::Relaxed);
    total_scanned.fetch_add(stats.scanned, Ordering::Relaxed);
    if cli.verbose {
        let name = path.file_name().unwrap_or_default();
        eprintln!(
            "\n{:?}: kept={}, scanned={}",
            name, stats.kept, stats.scanned
        );
    }
    Ok(())
}

fn stream_and_filter(input: &Path, output: &Path, mag_limit: f32) -> anyhow::Result<FileStats> {
    let file = File::open(input)?;
    let mut decoder = GzDecoder::new(BufReader::new(file));
    validate_gzip_header(&mut decoder, input)?;
    let buf_reader = BufReader::new(decoder);
    let parser = GaiaParser::new(buf_reader, mag_limit)?;
    let temp_path = output.with_extension("bin.tmp");
    let stats = write_part_file(parser, &temp_path)?;
    fs::rename(&temp_path, output)?;
    Ok(stats)
}

fn validate_gzip_header<R: Read>(decoder: &mut GzDecoder<R>, path: &Path) -> anyhow::Result<()> {
    let header = decoder.header();
    if header.is_none() {
        anyhow::bail!("Invalid or corrupt gzip file: {:?}", path);
    }
    Ok(())
}

fn write_part_file<R: std::io::BufRead>(
    parser: GaiaParser<R>,
    output: &Path,
) -> anyhow::Result<FileStats> {
    let file = File::create(output)?;
    let mut writer = BufWriter::new(file);
    writer.write_all(&[0u8; PART_HEADER_SIZE])?;
    let mut stats = FileStats {
        kept: 0,
        scanned: 0,
    };
    for result in parser {
        stats.scanned += 1;
        match result {
            Ok(star) => {
                write_star_record(&mut writer, &star)?;
                stats.kept += 1;
            }
            Err(e) => {
                // Decompression or parse error — bail immediately
                anyhow::bail!("Error at row {}: {}", stats.scanned, e);
            }
        }
    }
    finalize_part_file(writer, stats)
}

fn finalize_part_file(mut writer: BufWriter<File>, stats: FileStats) -> anyhow::Result<FileStats> {
    writer.flush()?;
    let mut file = writer.into_inner()?;
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&stats.kept.to_le_bytes())?;
    Ok(stats)
}

fn write_star_record<W: Write>(writer: &mut W, star: &GaiaStarRaw) -> anyhow::Result<()> {
    let record = StarRecord {
        source_id: star.source_id,
        ra: star.ra,
        dec: star.dec,
        pmra: star.pmra,
        pmdec: star.pmdec,
        parallax: star.parallax,
        mag: star.mag,
        flags: star.flags,
        _padding: 0,
    };
    let bytes = unsafe {
        std::slice::from_raw_parts(&record as *const StarRecord as *const u8, RECORD_SIZE)
    };
    writer.write_all(bytes)?;
    Ok(())
}

fn print_summary(
    file_count: usize,
    total_kept: &AtomicU64,
    total_scanned: &AtomicU64,
    errors: &AtomicU64,
) -> anyhow::Result<()> {
    let kept = total_kept.load(Ordering::Relaxed);
    let scanned = total_scanned.load(Ordering::Relaxed);
    let errs = errors.load(Ordering::Relaxed);
    println!();
    println!("=== Summary ===");
    println!("Files processed: {}", file_count);
    println!("Stars scanned: {}", scanned);
    println!("Stars kept: {}", kept);
    if errs > 0 {
        println!("Files with errors: {}", errs);
    }
    if errs == file_count as u64 {
        anyhow::bail!("All files failed to process");
    }
    Ok(())
}

fn concatenate_part_files(output_dir: &Path, mag_limit: f32) -> anyhow::Result<()> {
    let part_files = collect_part_files(output_dir)?;
    if part_files.is_empty() {
        anyhow::bail!("No part files found in {:?}", output_dir);
    }
    println!("\n=== Concatenating {} part files ===", part_files.len());
    let final_path = output_dir.join("gaia_ingest.bin");
    let temp_path = final_path.with_extension("bin.tmp");
    let total_stars = write_final_catalog(&part_files, &temp_path, mag_limit)?;
    fs::rename(&temp_path, &final_path)?;
    println!("Written {} stars to {:?}", total_stars, final_path);
    delete_part_files(&part_files)?;
    Ok(())
}

fn collect_part_files(output_dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = fs::read_dir(output_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| is_part_file(p))
        .collect();
    files.sort();
    Ok(files)
}

fn is_part_file(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .map(|n| n.starts_with("gaia_part_") && n.ends_with(".bin"))
        .unwrap_or(false)
}

fn write_final_catalog(
    part_files: &[PathBuf],
    output: &Path,
    mag_limit: f32,
) -> anyhow::Result<u64> {
    let total_stars = count_total_stars(part_files)?;
    let file = File::create(output)?;
    let mut writer = BufWriter::new(file);
    write_final_header(&mut writer, total_stars, mag_limit)?;
    copy_star_records(&mut writer, part_files)?;
    writer.flush()?;
    Ok(total_stars)
}

fn count_total_stars(part_files: &[PathBuf]) -> anyhow::Result<u64> {
    let mut total = 0u64;
    for path in part_files {
        total += read_part_star_count(path)?;
    }
    Ok(total)
}

fn read_part_star_count(path: &Path) -> anyhow::Result<u64> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 8];
    file.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn write_final_header<W: Write>(writer: &mut W, total: u64, mag_limit: f32) -> anyhow::Result<()> {
    writer.write_all(FINAL_MAGIC)?;
    writer.write_all(&FINAL_VERSION.to_le_bytes())?;
    writer.write_all(&total.to_le_bytes())?;
    writer.write_all(&mag_limit.to_le_bytes())?;
    writer.write_all(&[0u8; 4])?;
    Ok(())
}

fn copy_star_records<W: Write>(writer: &mut W, part_files: &[PathBuf]) -> anyhow::Result<()> {
    let mut buf = vec![0u8; 64 * 1024];
    for path in part_files {
        copy_single_part(writer, path, &mut buf)?;
    }
    Ok(())
}

fn copy_single_part<W: Write>(writer: &mut W, path: &Path, buf: &mut [u8]) -> anyhow::Result<()> {
    let mut file = File::open(path)?;
    file.seek(SeekFrom::Start(PART_HEADER_SIZE as u64))?;
    loop {
        let n = file.read(buf)?;
        if n == 0 {
            break;
        }
        writer.write_all(&buf[..n])?;
    }
    Ok(())
}

fn delete_part_files(part_files: &[PathBuf]) -> anyhow::Result<()> {
    for path in part_files {
        fs::remove_file(path)?;
    }
    println!("Deleted {} part files", part_files.len());
    Ok(())
}
