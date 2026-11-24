//! Download Gaia DR3 source files from ESA CDN
//!
//! Lists the S3 bucket, builds a manifest with ETags,
//! and downloads files with resume support and parallel fetching.

use crate::cli::{Cli, DownloadGaiaArgs};
use crate::gaia_listing::{list_remote_files, RemoteFile};
use crate::http::Downloader;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

const CDN_BASE: &str = "https://cdn.gea.esac.esa.int/";
const MANIFEST_FILENAME: &str = "gaia_manifest.json";
const BYTES_PER_GB: f64 = 1_073_741_824.0;

#[derive(Debug, Default, Serialize, Deserialize)]
struct Manifest {
    files: HashMap<String, FileEntry>,
}

#[derive(Debug, Serialize, Deserialize)]
struct FileEntry {
    etag: String,
    size: u64,
    downloaded: bool,
}

struct DownloadJob<'a> {
    http: &'a Downloader,
    base_url: &'a str,
    output: &'a Path,
    concurrency: usize,
}

struct Progress {
    total: usize,
    completed: AtomicUsize,
    failed: AtomicUsize,
    bytes: AtomicU64,
}

pub(crate) fn run(args: &DownloadGaiaArgs, cli: &Cli) -> anyhow::Result<()> {
    fs::create_dir_all(&args.output)?;
    print_plan(args);
    let http = Downloader::new(args.retries)?;
    let remote_files = list_files(&http, args.limit, cli.verbose)?;
    let manifest_path = args.output.join(MANIFEST_FILENAME);
    let manifest = load_manifest(&manifest_path)?;
    let to_download = plan_downloads(&remote_files, &manifest, &args.output);
    if to_download.is_empty() {
        println!("All files already downloaded and verified.");
        return Ok(());
    }
    let skipped = remote_files.len() - to_download.len();
    print_download_plan(skipped, &to_download);
    let job = DownloadJob {
        http: &http,
        base_url: CDN_BASE,
        output: &args.output,
        concurrency: args.concurrency,
    };
    let progress = download_and_record(&job, &to_download, manifest, &manifest_path)?;
    print_summary(&progress, skipped)
}

fn print_plan(args: &DownloadGaiaArgs) {
    println!("=== Gaia DR3 Download ===");
    println!("Output: {:?}", args.output);
    println!("Concurrency: {}", args.concurrency);
    if let Some(limit) = args.limit {
        println!("Limit: {} files", limit);
    }
    println!();
}

fn list_files(
    http: &Downloader,
    limit: Option<usize>,
    verbose: bool,
) -> anyhow::Result<Vec<RemoteFile>> {
    println!("Listing files from ESA CDN...");
    let mut remote_files = list_remote_files(http, verbose)?;
    remote_files.sort_by(|a, b| a.filename.cmp(&b.filename));
    if let Some(limit) = limit {
        remote_files.truncate(limit);
    }
    let total_size: u64 = remote_files.iter().map(|f| f.size).sum();
    println!(
        "Found {} files ({:.1} GB total)",
        remote_files.len(),
        total_size as f64 / BYTES_PER_GB
    );
    Ok(remote_files)
}

fn print_download_plan(skipped: usize, to_download: &[&RemoteFile]) {
    let size: u64 = to_download.iter().map(|f| f.size).sum();
    println!(
        "Skipping {} already downloaded, {} to download ({:.1} GB)",
        skipped,
        to_download.len(),
        size as f64 / BYTES_PER_GB
    );
    println!();
}

// Records each file in the manifest as soon as it is downloaded and verified, so a run that
// stops early keeps the files it finished. A failed file keeps whatever entry it had.
fn download_and_record(
    job: &DownloadJob,
    files: &[&RemoteFile],
    manifest: Manifest,
    manifest_path: &Path,
) -> anyhow::Result<Progress> {
    let progress = Progress::new(files.len());
    let manifest = Mutex::new(manifest);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(job.concurrency)
        .build()
        .context("Failed to build thread pool")?;
    pool.scope(|s| {
        for &file in files {
            let (progress, manifest) = (&progress, &manifest);
            s.spawn(
                move |_| match download_one(job, file, manifest, manifest_path) {
                    Ok(()) => progress.record_success(file),
                    Err(e) => progress.record_failure(file, &e),
                },
            );
        }
    });
    Ok(progress)
}

fn download_one(
    job: &DownloadJob,
    file: &RemoteFile,
    manifest: &Mutex<Manifest>,
    manifest_path: &Path,
) -> anyhow::Result<()> {
    let url = format!("{}{}", job.base_url, file.key);
    job.http
        .download(&url, &job.output.join(&file.filename), Some(file.size))?;
    let mut manifest = manifest
        .lock()
        .map_err(|_| anyhow::anyhow!("Another download panicked while saving the manifest"))?;
    let entry = FileEntry {
        etag: file.etag.clone(),
        size: file.size,
        downloaded: true,
    };
    manifest.files.insert(file.filename.clone(), entry);
    save_manifest(&manifest, manifest_path)
}

impl Progress {
    fn new(total: usize) -> Self {
        Progress {
            total,
            completed: AtomicUsize::new(0),
            failed: AtomicUsize::new(0),
            bytes: AtomicU64::new(0),
        }
    }

    fn record_success(&self, file: &RemoteFile) {
        let done = self.completed.fetch_add(1, Ordering::Relaxed) + 1;
        let bytes = self.bytes.fetch_add(file.size, Ordering::Relaxed) + file.size;
        println!(
            "[{}/{}] {} ({:.1} MB) - {:.1} GB done",
            done,
            self.total,
            file.filename,
            file.size as f64 / 1_048_576.0,
            bytes as f64 / BYTES_PER_GB,
        );
    }

    fn record_failure(&self, file: &RemoteFile, error: &anyhow::Error) {
        self.failed.fetch_add(1, Ordering::Relaxed);
        eprintln!("FAILED {}: {:#}", file.filename, error);
    }
}

fn print_summary(progress: &Progress, skipped: usize) -> anyhow::Result<()> {
    let failed = progress.failed.load(Ordering::Relaxed);
    println!("\n=== Summary ===");
    println!("Downloaded: {}", progress.completed.load(Ordering::Relaxed));
    println!("Skipped: {}", skipped);
    println!("Failed: {}", failed);
    if failed > 0 {
        anyhow::bail!("{} downloads failed. Re-run to retry.", failed);
    }
    Ok(())
}

// A missing manifest means nothing has been downloaded yet.
fn load_manifest(path: &Path) -> anyhow::Result<Manifest> {
    let data = match fs::read_to_string(path) {
        Ok(data) => data,
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(Manifest::default()),
        Err(e) => return Err(e).with_context(|| format!("Failed to read {:?}", path)),
    };
    serde_json::from_str(&data).with_context(|| {
        format!(
            "{:?} is corrupt. Delete it and re-run; every file is then downloaded again",
            path
        )
    })
}

// The manifest is written to a temporary file and renamed, so an interrupted save leaves the
// previous manifest in place.
fn save_manifest(manifest: &Manifest, path: &Path) -> anyhow::Result<()> {
    let json = serde_json::to_string_pretty(manifest)?;
    // Synced before the rename: a manifest cut short by a power loss reads as corrupt, and
    // the only way out of that is downloading every file again.
    crate::intermediate::write_via_temp(path, |writer| {
        io::Write::write_all(writer, json.as_bytes())?;
        Ok(())
    })
    .with_context(|| format!("Failed to save {:?}", path))
}

fn plan_downloads<'a>(
    remote: &'a [RemoteFile],
    manifest: &Manifest,
    output_dir: &Path,
) -> Vec<&'a RemoteFile> {
    remote
        .iter()
        .filter(|f| !is_already_good(f, manifest, output_dir))
        .collect()
}

fn is_already_good(file: &RemoteFile, manifest: &Manifest, output_dir: &Path) -> bool {
    if !has_size(&output_dir.join(&file.filename), file.size) {
        return false;
    }
    manifest
        .files
        .get(&file.filename)
        .is_some_and(|entry| entry.etag == file.etag && entry.size == file.size && entry.downloaded)
}

fn has_size(path: &Path, size: u64) -> bool {
    fs::metadata(path).is_ok_and(|m| m.len() == size)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::http::tests::downloader;
    use crate::test_server::{Reply, Server};
    use std::time::Duration;

    fn remote(name: &str, etag: &str) -> RemoteFile {
        RemoteFile {
            key: format!("Gaia/{name}"),
            filename: name.to_string(),
            size: 3,
            etag: etag.to_string(),
        }
    }

    fn job<'a>(http: &'a Downloader, base_url: &'a str, output: &'a Path) -> DownloadJob<'a> {
        DownloadJob {
            http,
            base_url,
            output,
            concurrency: 1,
        }
    }

    fn planned(remote: &[RemoteFile], dir: &Path) -> Vec<String> {
        let manifest = load_manifest(&dir.join(MANIFEST_FILENAME)).unwrap();
        let plan = plan_downloads(remote, &manifest, dir);
        plan.iter().map(|f| f.filename.clone()).collect()
    }

    // The third request stalls and then fails, as a run stopped by Ctrl-C during its third
    // file would. A run that starts at that moment must skip the two finished files.
    #[test]
    fn test_manifest_records_each_file_as_it_finishes() {
        let dir = tempfile::tempdir().unwrap();
        let files: Vec<RemoteFile> = ["a", "b", "c"]
            .iter()
            .map(|n| remote(&format!("{n}.csv.gz"), n))
            .collect();
        let replies = vec![
            Reply::ok(b"abc"),
            Reply::ok(b"abc"),
            Reply::stalled(3, 1, Duration::from_secs(1)),
        ];
        let server = Server::start(replies);
        let http = downloader(5000, 0);
        let manifest_path = dir.path().join(MANIFEST_FILENAME);
        let progress = std::thread::scope(|s| {
            let run = s.spawn(|| {
                let to_download: Vec<&RemoteFile> = files.iter().collect();
                let job = job(&http, &server.url, dir.path());
                download_and_record(&job, &to_download, Manifest::default(), &manifest_path)
            });
            server.wait_for_requests(3);
            let left = planned(&files, dir.path());
            assert_eq!(left.len(), 1, "{left:?}");
            assert!(!dir.path().join(&left[0]).exists());
            run.join().unwrap().unwrap()
        });
        server.join();
        assert_eq!(progress.failed.load(Ordering::Relaxed), 1);
        assert_eq!(planned(&files, dir.path()).len(), 1);
    }

    // The local copy matches the new size, but downloading the new version failed, so the
    // manifest must keep the old ETag and the next run must download the file again.
    #[test]
    fn test_failed_update_is_not_recorded() {
        let dir = tempfile::tempdir().unwrap();
        let manifest_path = dir.path().join(MANIFEST_FILENAME);
        fs::write(dir.path().join("a.csv.gz"), b"old").unwrap();
        let mut manifest = Manifest::default();
        let entry = FileEntry {
            etag: "old".to_string(),
            size: 3,
            downloaded: true,
        };
        manifest.files.insert("a.csv.gz".to_string(), entry);
        save_manifest(&manifest, &manifest_path).unwrap();
        let files = [remote("a.csv.gz", "new")];
        let server = Server::start(vec![Reply::status("503 Service Unavailable")]);
        let http = downloader(1000, 0);
        let to_download: Vec<&RemoteFile> = files.iter().collect();
        let job = job(&http, &server.url, dir.path());
        let progress = download_and_record(&job, &to_download, manifest, &manifest_path).unwrap();
        server.join();
        assert_eq!(progress.failed.load(Ordering::Relaxed), 1);
        let saved = load_manifest(&manifest_path).unwrap();
        assert_eq!(saved.files["a.csv.gz"].etag, "old");
        assert_eq!(planned(&files, dir.path()), ["a.csv.gz"]);
    }

    #[test]
    fn test_missing_manifest_is_empty() {
        let dir = tempfile::tempdir().unwrap();
        let manifest = load_manifest(&dir.path().join(MANIFEST_FILENAME)).unwrap();
        assert!(manifest.files.is_empty());
    }

    #[test]
    fn test_corrupt_manifest_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(MANIFEST_FILENAME);
        fs::write(&path, "{\"files\": {").unwrap();
        let err = load_manifest(&path).unwrap_err();
        assert!(
            format!("{:#}", err).contains("is corrupt. Delete it"),
            "{err:#}"
        );
    }

    #[test]
    fn test_saved_manifest_loads() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(MANIFEST_FILENAME);
        let mut manifest = Manifest::default();
        let entry = FileEntry {
            etag: "e".to_string(),
            size: 3,
            downloaded: true,
        };
        manifest.files.insert("a.csv.gz".to_string(), entry);
        save_manifest(&manifest, &path).unwrap();
        let loaded = load_manifest(&path).unwrap();
        let entry = &loaded.files["a.csv.gz"];
        assert_eq!(
            (entry.etag.as_str(), entry.size, entry.downloaded),
            ("e", 3, true)
        );
        let names: Vec<_> = fs::read_dir(dir.path()).unwrap().collect();
        assert_eq!(names.len(), 1);
    }
}
