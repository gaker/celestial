//! Gaia DR3 catalog ingestion

use crate::cli::{Cli, IngestGaiaArgs};
use crate::gaia::{GaiaParser, RowCounts};
use crate::gaia_parts::{concatenate, part_exists, part_path, write_part};
use anyhow::Context;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};

const FAILED_LOG: &str = "failed_files.log";

pub(crate) fn run(args: &IngestGaiaArgs, cli: &Cli) -> anyhow::Result<()> {
    validate_args(args)?;
    let files = find_gzipped_csvs(&args.path)?;
    if files.is_empty() {
        anyhow::bail!("No .csv.gz files found in {:?}", args.path);
    }
    print_plan(args, files.len());
    let pending: Vec<PathBuf> = files
        .iter()
        .filter(|f| !part_exists(f, &args.output))
        .cloned()
        .collect();
    println!(
        "Files to process: {} (skipping {} already done)",
        pending.len(),
        files.len() - pending.len()
    );
    let failed = process_files(&pending, args, cli)?;
    finish(args, failed)
}

// Concatenates the part files unless --no-concat was given. A failed file has no part file,
// so concatenating after a failure would build an incomplete catalog.
fn finish(args: &IngestGaiaArgs, failed: usize) -> anyhow::Result<()> {
    let log_path = args.output.join(FAILED_LOG);
    if failed == 0 {
        remove_stale_log(&log_path)?;
    }
    if args.no_concat {
        println!("Skipping concatenation (--no-concat)");
        return Ok(());
    }
    if failed > 0 {
        anyhow::bail!(
            "{failed} input files failed, listed in {log_path:?}; gaia_ingest.bin was not \
             built. To build it, fix or remove those files and rerun; files that already have \
             a part file are skipped"
        );
    }
    concatenate(&args.output, args.mag_limit)
}

// A failed_files.log left by an earlier run lists files that have since succeeded.
fn remove_stale_log(path: &Path) -> anyhow::Result<()> {
    match fs::remove_file(path) {
        Ok(()) => println!("Removed {path:?} from an earlier run"),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).with_context(|| format!("Failed to remove {path:?}")),
    }
    Ok(())
}

fn validate_args(args: &IngestGaiaArgs) -> anyhow::Result<()> {
    if !args.path.is_dir() {
        anyhow::bail!("Gaia path is not a directory: {:?}", args.path);
    }
    if !args.mag_limit.is_finite() {
        anyhow::bail!("Magnitude limit {} is not finite", args.mag_limit);
    }
    fs::create_dir_all(&args.output)?;
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
        std::thread::available_parallelism().map_or(1, |n| n.get())
    } else {
        threads
    }
}

// Returns the number of files that failed. Every file failing is an error.
fn process_files(files: &[PathBuf], args: &IngestGaiaArgs, cli: &Cli) -> anyhow::Result<usize> {
    if files.is_empty() {
        return Ok(0);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(resolve_threads(args.threads))
        .build()
        .context("Failed to build the thread pool")?;
    let pb = create_progress_bar(files.len() as u64)?;
    let results: Vec<_> = pool.install(|| {
        files
            .par_iter()
            .map(|file| {
                let result = process_single_file(file, args, cli);
                pb.inc(1);
                result
            })
            .collect()
    });
    pb.finish_with_message("Done");
    let (counts, failed) = tally(files, results);
    write_failed_log(&args.output, &failed)?;
    print_summary(files.len(), &counts, failed.len());
    if failed.len() == files.len() {
        anyhow::bail!("All files failed to process");
    }
    Ok(failed.len())
}

type FileResult = anyhow::Result<RowCounts>;

fn tally(files: &[PathBuf], results: Vec<FileResult>) -> (RowCounts, Vec<(&Path, String)>) {
    let mut counts = RowCounts::default();
    let mut failed = Vec::new();
    for (file, result) in files.iter().zip(results) {
        match result {
            Ok(file_counts) => counts.add(file_counts),
            Err(e) => failed.push((file.as_path(), format!("{e:#}"))),
        }
    }
    (counts, failed)
}

fn write_failed_log(output_dir: &Path, failed: &[(&Path, String)]) -> anyhow::Result<()> {
    if failed.is_empty() {
        return Ok(());
    }
    let log_path = output_dir.join(FAILED_LOG);
    let mut file = File::create(&log_path)?;
    for (path, error) in failed {
        writeln!(file, "{}\t{}", path.display(), error)?;
    }
    println!("Failed files logged to: {:?}", log_path);
    Ok(())
}

fn create_progress_bar(total: u64) -> anyhow::Result<ProgressBar> {
    let style = ProgressStyle::default_bar()
        .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")?
        .progress_chars("#>-");
    Ok(ProgressBar::new(total).with_style(style))
}

fn process_single_file(path: &Path, args: &IngestGaiaArgs, cli: &Cli) -> FileResult {
    let result = ingest_file(path, &args.output, args.mag_limit);
    match &result {
        Err(e) => eprintln!("\nWarning: {}: {e:#}", path.display()),
        Ok(counts) if cli.verbose => eprintln!("\n{}: {counts:?}", path.display()),
        Ok(_) => {}
    }
    result
}

fn ingest_file(input: &Path, output_dir: &Path, mag_limit: f32) -> FileResult {
    let output = part_path(input, output_dir).with_context(|| {
        format!("{input:?} is not named GaiaSource_<range>.csv.gz, so it has no part file name")
    })?;
    let decoder = GzDecoder::new(BufReader::new(File::open(input)?));
    if decoder.header().is_none() {
        anyhow::bail!("Invalid or corrupt gzip file: {:?}", input);
    }
    let parser = GaiaParser::new(BufReader::new(decoder), mag_limit)?;
    write_part(parser, &output, mag_limit)
}

fn print_summary(file_count: usize, counts: &RowCounts, failed: usize) {
    let fainter = counts.read - counts.kept - counts.bad_magnitude - counts.bad_position;
    println!();
    println!("=== Summary ===");
    println!("Files processed: {}", file_count);
    println!("Rows read: {}", counts.read);
    println!("Stars kept: {}", counts.kept);
    println!("Rows fainter than the limit: {}", fainter);
    println!(
        "Rows skipped, missing or invalid magnitude: {}",
        counts.bad_magnitude
    );
    println!(
        "Rows skipped, missing or invalid source_id, ra or dec: {}",
        counts.bad_position
    );
    if failed > 0 {
        println!("Files with errors: {}", failed);
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::cli::Commands;
    use crate::gaia::tests::{csv, row};
    use crate::intermediate::{self, Kind};
    use flate2::write::GzEncoder;

    pub(crate) fn write_gz(path: &Path, text: &str) {
        let mut encoder = GzEncoder::new(File::create(path).unwrap(), Default::default());
        encoder.write_all(text.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    pub(crate) fn ingest(input: &Path, output: &Path, no_concat: bool) -> anyhow::Result<()> {
        let args = || IngestGaiaArgs {
            path: input.to_path_buf(),
            mag_limit: 18.5,
            output: output.to_path_buf(),
            no_concat,
            threads: 2,
        };
        let cli = Cli {
            verbose: false,
            command: Commands::IngestGaia(args()),
        };
        run(&args(), &cli)
    }

    fn setup() -> (tempfile::TempDir, PathBuf, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let (input, output) = (dir.path().join("in"), dir.path().join("out"));
        fs::create_dir(&input).unwrap();
        let good = csv(&[row("1", "10.0", "20.0", "12.0", 31)]);
        write_gz(&input.join("GaiaSource_000000-000001.csv.gz"), &good);
        fs::write(input.join("GaiaSource_000002-000003.csv.gz"), b"not gzip").unwrap();
        (dir, input, output)
    }

    #[test]
    fn test_failed_file_error_says_how_to_unblock_the_run() {
        let (_dir, input, output) = setup();
        let err = format!("{:#}", ingest(&input, &output, false).unwrap_err());
        assert!(err.contains("fix or remove those files and rerun"), "{err}");
        assert!(!err.contains("--no-concat"), "{err}");
    }

    #[test]
    fn test_successful_run_deletes_a_stale_failed_files_log() {
        let (_dir, input, output) = setup();
        assert!(ingest(&input, &output, false).is_err());
        fs::remove_file(input.join("GaiaSource_000002-000003.csv.gz")).unwrap();
        ingest(&input, &output, false).unwrap();
        assert!(!output.join("failed_files.log").exists());
    }

    #[test]
    fn test_failed_file_stops_concatenation_and_keeps_parts() {
        let (_dir, input, output) = setup();
        let err = format!("{:#}", ingest(&input, &output, false).unwrap_err());
        assert!(err.contains("1 input files failed"), "{err}");
        assert!(!output.join("gaia_ingest.bin").exists());
        assert!(output.join("gaia_part_000000-000001.bin").exists());
        assert!(output.join("failed_files.log").exists());
    }

    #[test]
    fn test_no_concat_keeps_parts_after_a_failed_file() {
        let (_dir, input, output) = setup();
        ingest(&input, &output, true).unwrap();
        assert!(!output.join("gaia_ingest.bin").exists());
        assert!(output.join("gaia_part_000000-000001.bin").exists());
    }

    #[test]
    fn test_rerun_after_fixing_the_file_concatenates_every_part() {
        let (_dir, input, output) = setup();
        assert!(ingest(&input, &output, false).is_err());
        let second = csv(&[row("2", "11.0", "21.0", "13.0", 31)]);
        write_gz(&input.join("GaiaSource_000002-000003.csv.gz"), &second);
        ingest(&input, &output, false).unwrap();
        let (header, mut reader) = intermediate::open(&output, Kind::Gaia).unwrap();
        assert_eq!(header.records, 2);
        let ids = [0, 1].map(|_| intermediate::read_record(&mut reader).unwrap().source_id);
        assert_eq!(ids, [1, 2]);
        assert!(!output.join("gaia_part_000000-000001.bin").exists());
    }
}
