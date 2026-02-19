//! Download Gaia DR3 source files from ESA CDN
//!
//! Lists the S3 bucket, builds a manifest with ETags,
//! and downloads files with resume support and parallel fetching.

use crate::cli::{Cli, DownloadGaiaArgs};
use anyhow::Context;
use quick_xml::events::Event;
use quick_xml::Reader;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;

const LISTING_URL: &str =
    "https://gaia.eu-1.cdn77-storage.com/?prefix=Gaia/gdr3/gaia_source/&delimiter=/";
const CDN_BASE: &str = "https://cdn.gea.esac.esa.int/";
const MANIFEST_FILENAME: &str = "gaia_manifest.json";

#[derive(Debug, Serialize, Deserialize)]
struct Manifest {
    files: HashMap<String, FileEntry>,
}

#[derive(Debug, Serialize, Deserialize)]
struct FileEntry {
    etag: String,
    size: u64,
    downloaded: bool,
}

struct RemoteFile {
    key: String,
    filename: String,
    size: u64,
    etag: String,
}

pub fn run(args: &DownloadGaiaArgs, cli: &Cli) -> anyhow::Result<()> {
    fs::create_dir_all(&args.output)?;
    println!("=== Gaia DR3 Download ===");
    println!("Output: {:?}", args.output);
    println!("Concurrency: {}", args.concurrency);
    if let Some(limit) = args.limit {
        println!("Limit: {} files", limit);
    }
    println!();

    println!("Listing files from ESA CDN...");
    let mut remote_files = list_remote_files(cli.verbose)?;
    remote_files.sort_by(|a, b| a.filename.cmp(&b.filename));

    if let Some(limit) = args.limit {
        remote_files.truncate(limit);
    }

    let total_size: u64 = remote_files.iter().map(|f| f.size).sum();
    println!(
        "Found {} files ({:.1} GB total)",
        remote_files.len(),
        total_size as f64 / 1_073_741_824.0
    );

    let manifest_path = args.output.join(MANIFEST_FILENAME);
    let mut manifest = load_manifest(&manifest_path);
    let to_download = plan_downloads(&remote_files, &manifest, &args.output);

    if to_download.is_empty() {
        println!("All files already downloaded and verified.");
        return Ok(());
    }

    let skip_count = remote_files.len() - to_download.len();
    let dl_size: u64 = to_download.iter().map(|f| f.size).sum();
    println!(
        "Skipping {} already downloaded, {} to download ({:.1} GB)",
        skip_count,
        to_download.len(),
        dl_size as f64 / 1_073_741_824.0
    );
    println!();

    let completed = Arc::new(AtomicUsize::new(0));
    let bytes_done = Arc::new(AtomicU64::new(0));
    let failed = Arc::new(AtomicUsize::new(0));
    let total_count = to_download.len();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.concurrency)
        .build()
        .context("Failed to build thread pool")?;

    pool.scope(|s| {
        for file in &to_download {
            let output = &args.output;
            let retries = args.retries;
            let completed = Arc::clone(&completed);
            let bytes_done = Arc::clone(&bytes_done);
            let failed = Arc::clone(&failed);

            s.spawn(move |_| {
                let dest = output.join(&file.filename);
                match download_with_retry(&file.key, &dest, file.size, retries) {
                    Ok(etag) => {
                        let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                        let bytes = bytes_done.fetch_add(file.size, Ordering::Relaxed) + file.size;
                        println!(
                            "[{}/{}] {} ({:.1} MB) - {:.1} GB done",
                            done,
                            total_count,
                            file.filename,
                            file.size as f64 / 1_048_576.0,
                            bytes as f64 / 1_073_741_824.0,
                        );
                        let _ = etag; // used when saving manifest below
                    }
                    Err(e) => {
                        failed.fetch_add(1, Ordering::Relaxed);
                        eprintln!("FAILED {}: {}", file.filename, e);
                    }
                }
            });
        }
    });

    update_manifest(&mut manifest, &remote_files, &args.output);
    save_manifest(&manifest, &manifest_path)?;

    let fail_count = failed.load(Ordering::Relaxed);
    let ok_count = completed.load(Ordering::Relaxed);
    println!("\n=== Summary ===");
    println!("Downloaded: {}", ok_count);
    println!("Skipped: {}", skip_count);
    println!("Failed: {}", fail_count);

    if fail_count > 0 {
        anyhow::bail!("{} downloads failed. Re-run to retry.", fail_count);
    }
    Ok(())
}

fn list_remote_files(verbose: bool) -> anyhow::Result<Vec<RemoteFile>> {
    let mut files = Vec::new();
    let mut marker: Option<String> = None;

    loop {
        let url = match &marker {
            Some(m) => format!("{}&marker={}", LISTING_URL, m),
            None => LISTING_URL.to_string(),
        };

        if verbose {
            eprintln!("Listing: {}", url);
        }

        let body = reqwest::blocking::get(&url)
            .context("Failed to fetch bucket listing")?
            .text()
            .context("Failed to read listing response")?;

        let (batch, next_marker) = parse_listing(&body)?;
        let batch_len = batch.len();
        files.extend(batch);

        if verbose {
            eprintln!("  Got {} keys (total: {})", batch_len, files.len());
        }

        match next_marker {
            Some(m) => marker = Some(m),
            None => break,
        }
    }

    Ok(files)
}

fn parse_listing(xml: &str) -> anyhow::Result<(Vec<RemoteFile>, Option<String>)> {
    let mut reader = Reader::from_str(xml);
    let mut files = Vec::new();
    let mut next_marker: Option<String> = None;
    let mut buf = String::new();

    let mut in_contents = false;
    let mut in_key = false;
    let mut in_size = false;
    let mut in_etag = false;
    let mut in_next_marker = false;

    let mut cur_key = String::new();
    let mut cur_size = 0u64;
    let mut cur_etag = String::new();

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                match name.as_str() {
                    "Contents" => {
                        in_contents = true;
                        cur_key.clear();
                        cur_size = 0;
                        cur_etag.clear();
                    }
                    "Key" if in_contents => in_key = true,
                    "Size" if in_contents => in_size = true,
                    "ETag" if in_contents => in_etag = true,
                    "NextMarker" => in_next_marker = true,
                    _ => {}
                }
            }
            Ok(Event::Text(e)) => {
                buf.clear();
                buf.push_str(&e.unescape().unwrap_or_default());
                if in_key {
                    cur_key.push_str(&buf);
                } else if in_size {
                    cur_size = buf.trim().parse().unwrap_or(0);
                } else if in_etag {
                    cur_etag.push_str(buf.trim().trim_matches('"'));
                } else if in_next_marker {
                    next_marker = Some(buf.trim().to_string());
                }
            }
            Ok(Event::End(e)) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                match name.as_str() {
                    "Contents" => {
                        in_contents = false;
                        if cur_key.ends_with(".csv.gz") {
                            let filename = extract_filename(&cur_key);
                            files.push(RemoteFile {
                                key: cur_key.clone(),
                                filename,
                                size: cur_size,
                                etag: cur_etag.clone(),
                            });
                        }
                    }
                    "Key" => in_key = false,
                    "Size" => in_size = false,
                    "ETag" => in_etag = false,
                    "NextMarker" => in_next_marker = false,
                    _ => {}
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => anyhow::bail!("XML parse error: {}", e),
            _ => {}
        }
    }

    Ok((files, next_marker))
}

fn extract_filename(key: &str) -> String {
    key.rsplit('/').next().unwrap_or(key).to_string()
}

fn load_manifest(path: &Path) -> Manifest {
    match fs::read_to_string(path) {
        Ok(data) => serde_json::from_str(&data).unwrap_or(Manifest {
            files: HashMap::new(),
        }),
        Err(_) => Manifest {
            files: HashMap::new(),
        },
    }
}

fn save_manifest(manifest: &Manifest, path: &Path) -> anyhow::Result<()> {
    let json = serde_json::to_string_pretty(manifest)?;
    fs::write(path, json)?;
    Ok(())
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
    let local_path = output_dir.join(&file.filename);
    let meta = match fs::metadata(&local_path) {
        Ok(m) => m,
        Err(_) => return false,
    };
    if meta.len() != file.size {
        return false;
    }
    if let Some(entry) = manifest.files.get(&file.filename) {
        return entry.etag == file.etag && entry.size == file.size && entry.downloaded;
    }
    false
}

fn download_with_retry(
    key: &str,
    dest: &Path,
    expected_size: u64,
    max_retries: u32,
) -> anyhow::Result<String> {
    let url = format!("{}{}", CDN_BASE, key);
    let mut last_err = None;

    for attempt in 0..=max_retries {
        if attempt > 0 {
            eprintln!("  Retry {}/{} for {}", attempt, max_retries, key);
        }
        match download_file(&url, dest, expected_size) {
            Ok(etag) => return Ok(etag),
            Err(e) => {
                last_err = Some(e);
                if dest.exists() {
                    let _ = fs::remove_file(dest);
                }
            }
        }
    }

    Err(last_err.unwrap())
}

fn download_file(url: &str, dest: &Path, expected_size: u64) -> anyhow::Result<String> {
    let response =
        reqwest::blocking::get(url).with_context(|| format!("Failed to connect: {}", url))?;

    if !response.status().is_success() {
        anyhow::bail!("HTTP {}", response.status());
    }

    let etag = response
        .headers()
        .get("etag")
        .and_then(|v| v.to_str().ok())
        .map(|s| s.trim_matches('"').to_string())
        .unwrap_or_default();

    let tmp_path = dest.with_extension("csv.gz.tmp");
    let mut file =
        File::create(&tmp_path).with_context(|| format!("Failed to create {:?}", tmp_path))?;

    let bytes = response.bytes().context("Failed to read body")?;
    file.write_all(&bytes)?;
    file.flush()?;
    drop(file);

    let written = fs::metadata(&tmp_path)?.len();
    if written != expected_size {
        let _ = fs::remove_file(&tmp_path);
        anyhow::bail!(
            "Size mismatch: expected {} got {} for {}",
            expected_size,
            written,
            url
        );
    }

    fs::rename(&tmp_path, dest)
        .with_context(|| format!("Failed to rename {:?} -> {:?}", tmp_path, dest))?;

    Ok(etag)
}

fn update_manifest(manifest: &mut Manifest, remote: &[RemoteFile], output_dir: &Path) {
    for file in remote {
        let local_path = output_dir.join(&file.filename);
        let ok = match fs::metadata(&local_path) {
            Ok(m) => m.len() == file.size,
            Err(_) => false,
        };
        if ok {
            manifest.files.insert(
                file.filename.clone(),
                FileEntry {
                    etag: file.etag.clone(),
                    size: file.size,
                    downloaded: true,
                },
            );
        }
    }
}
