//! HTTP downloads shared by download-gaia and ingest-hipparcos.
//!
//! Integrity rests on the byte count and on gzip: every file forge fetches is gzip, whose
//! trailer holds a CRC32 of the data. ingest-gaia and ingest-hipparcos check it when they
//! decompress, and ingest-gaia stops before concatenation when any file fails.

use anyhow::Context;
use reqwest::blocking::{Client, ClientBuilder};
use reqwest::StatusCode;
use std::fs::{self, File};
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::Duration;

const CONNECT_TIMEOUT: Duration = Duration::from_secs(30);
// A Gaia file is about 200 MB, so a limit on the whole transfer would fail slow links. A
// request instead fails when one read waits this long for data.
const READ_TIMEOUT: Duration = Duration::from_secs(120);
const FIRST_RETRY_DELAY: Duration = Duration::from_secs(2);
const MAX_RETRY_DELAY: Duration = Duration::from_secs(60);

#[derive(Clone, Copy)]
struct Timeouts {
    connect: Duration,
    read: Duration,
}

#[derive(Clone, Copy)]
struct RetryPolicy {
    retries: u32,
    first_delay: Duration,
    max_delay: Duration,
}

pub(crate) struct Downloader {
    client: Client,
    retry: RetryPolicy,
}

impl Downloader {
    pub(crate) fn new(retries: u32) -> anyhow::Result<Self> {
        let timeouts = Timeouts {
            connect: CONNECT_TIMEOUT,
            read: READ_TIMEOUT,
        };
        let retry = RetryPolicy {
            retries,
            first_delay: FIRST_RETRY_DELAY,
            max_delay: MAX_RETRY_DELAY,
        };
        Self::with_config(timeouts, retry)
    }

    fn with_config(timeouts: Timeouts, retry: RetryPolicy) -> anyhow::Result<Self> {
        // The blocking client applies `timeout` to the wait for the response headers and to
        // each Read::read of the body, so download() streams the body through Read. A body
        // read with Response::bytes waits at most `timeout` in total. The async builder's
        // read_timeout cannot be used: reading such a body through blocking Read panics
        // outside a Tokio runtime in reqwest 0.13.3.
        let client = ClientBuilder::new()
            .connect_timeout(timeouts.connect)
            .timeout(timeouts.read)
            .build()
            .context("Failed to build the HTTP client")?;
        Ok(Downloader { client, retry })
    }

    // Returns the body of `url`, which must be UTF-8. The whole body must arrive within the
    // read timeout, which suits listing pages and not data files.
    pub(crate) fn get_text(&self, url: &str) -> anyhow::Result<String> {
        self.with_retries(url, || {
            let body = checked_get(&self.client, url)?
                .bytes()
                .context("Failed to read the response body")?;
            String::from_utf8(body.to_vec()).context("The response body is not UTF-8")
        })
    }

    // Writes the body of `url` to `dest` and returns its size. The body goes to a `.part` file
    // first, which becomes `dest` only when its size equals `expected_size`, or the
    // Content-Length header when `expected_size` is None.
    pub(crate) fn download(
        &self,
        url: &str,
        dest: &Path,
        expected_size: Option<u64>,
    ) -> anyhow::Result<u64> {
        self.download_and_check(url, dest, expected_size, |_| Ok(()))
    }

    // As download(), and `check` must also accept the `.part` file before it is renamed. A
    // rejected file counts as a failed attempt.
    pub(crate) fn download_and_check(
        &self,
        url: &str,
        dest: &Path,
        expected_size: Option<u64>,
        check: impl Fn(&Path) -> anyhow::Result<()>,
    ) -> anyhow::Result<u64> {
        let part = part_path(dest)?;
        self.with_retries(url, || {
            let result = download_once(&self.client, url, &part, expected_size);
            let result = result.and_then(|size| {
                check(&part)?;
                fs::rename(&part, dest)
                    .with_context(|| format!("Failed to rename {:?} to {:?}", part, dest))?;
                Ok(size)
            });
            if result.is_err() {
                remove_part(&part);
            }
            result
        })
    }

    fn with_retries<T>(
        &self,
        url: &str,
        mut attempt_once: impl FnMut() -> anyhow::Result<T>,
    ) -> anyhow::Result<T> {
        let mut retry = 0;
        loop {
            match attempt_once() {
                Ok(value) => return Ok(value),
                Err(e) if is_permanent(&e) => return Err(e),
                Err(e) if retry == self.retry.retries => {
                    return Err(e.context(format!("{} failed after {} tries", url, retry + 1)));
                }
                Err(e) => {
                    retry += 1;
                    let delay = backoff_delay(&self.retry, retry);
                    eprintln!(
                        "  {}: {:#}. Retry {}/{} in {:?}",
                        url, e, retry, self.retry.retries, delay
                    );
                    thread::sleep(delay);
                }
            }
        }
    }
}

fn checked_get(client: &Client, url: &str) -> anyhow::Result<reqwest::blocking::Response> {
    let response = client
        .get(url)
        .send()
        .with_context(|| format!("Failed to request {}", url))?;
    if !response.status().is_success() {
        let status = response.status();
        let url = url.to_string();
        return Err(StatusError { status, url }.into());
    }
    Ok(response)
}

#[derive(Debug, thiserror::Error)]
#[error("HTTP {status} for {url}")]
struct StatusError {
    status: StatusCode,
    url: String,
}

// A client error other than 408 Request Timeout or 429 Too Many Requests repeats on every
// attempt, so it is not retried.
fn is_permanent(error: &anyhow::Error) -> bool {
    error.downcast_ref::<StatusError>().is_some_and(|e| {
        e.status.is_client_error()
            && e.status != StatusCode::REQUEST_TIMEOUT
            && e.status != StatusCode::TOO_MANY_REQUESTS
    })
}

fn download_once(
    client: &Client,
    url: &str,
    part: &Path,
    expected_size: Option<u64>,
) -> anyhow::Result<u64> {
    let mut response = checked_get(client, url)?;
    let expected = expected_size.or(response.content_length());
    let file = File::create(part).with_context(|| format!("Failed to create {:?}", part))?;
    let mut writer = BufWriter::new(file);
    let written = io::copy(&mut response, &mut writer)
        .with_context(|| format!("Failed to download {}", url))?;
    let file = writer.into_inner().map_err(|e| e.into_error())?;
    file.sync_all()?;
    if let Some(expected) = expected.filter(|&n| n != written) {
        anyhow::bail!(
            "Size mismatch for {}: expected {} bytes, got {}",
            url,
            expected,
            written
        );
    }
    Ok(written)
}

fn part_path(dest: &Path) -> anyhow::Result<PathBuf> {
    let name = dest
        .file_name()
        .with_context(|| format!("{:?} has no file name", dest))?;
    let mut part = name.to_os_string();
    part.push(".part");
    Ok(dest.with_file_name(part))
}

fn remove_part(part: &Path) {
    if let Err(e) = fs::remove_file(part) {
        if e.kind() != io::ErrorKind::NotFound {
            eprintln!("Warning: could not remove {:?}: {}", part, e);
        }
    }
}

// The delay before retry number `retry` (1-based): the first delay, doubled for each earlier
// retry, at most the maximum.
fn backoff_delay(policy: &RetryPolicy, retry: u32) -> Duration {
    let factor = retry
        .checked_sub(1)
        .and_then(|n| 1u32.checked_shl(n))
        .unwrap_or(u32::MAX);
    policy
        .first_delay
        .saturating_mul(factor)
        .min(policy.max_delay)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::test_server::{Reply, Server};

    pub(crate) fn downloader(read_ms: u64, retries: u32) -> Downloader {
        let timeouts = Timeouts {
            connect: Duration::from_secs(5),
            read: Duration::from_millis(read_ms),
        };
        let retry = RetryPolicy {
            retries,
            first_delay: Duration::from_millis(10),
            max_delay: Duration::from_millis(40),
        };
        Downloader::with_config(timeouts, retry).unwrap()
    }

    fn file_names(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    // The body takes 1.2 s, twice the read timeout, so only a timeout on each read, with none
    // on the whole transfer, lets it complete.
    #[test]
    fn test_body_slower_than_read_timeout_downloads() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let body = vec![b'g'; 40];
        let server = Server::start(vec![Reply::slow(&body, 4, Duration::from_millis(300))]);
        let size = downloader(600, 0).download(&server.url, &dest, Some(40));
        server.join();
        assert_eq!(size.unwrap(), 40);
        assert_eq!(fs::read(&dest).unwrap(), body);
        assert_eq!(file_names(dir.path()), ["a.csv.gz"]);
    }

    #[test]
    fn test_stalled_body_times_out_and_leaves_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let server = Server::start(vec![Reply::stalled(20, 10, Duration::from_secs(1))]);
        let result = downloader(200, 0).download(&server.url, &dest, Some(20));
        server.join();
        assert!(result.is_err());
        assert!(file_names(dir.path()).is_empty());
    }

    #[test]
    fn test_size_mismatch_fails_and_leaves_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let server = Server::start(vec![Reply::ok(b"0123456789")]);
        let err = downloader(1000, 0)
            .download(&server.url, &dest, Some(12))
            .unwrap_err();
        server.join();
        assert!(
            format!("{:#}", err).contains("expected 12 bytes, got 10"),
            "{err:#}"
        );
        assert!(file_names(dir.path()).is_empty());
    }

    #[test]
    fn test_size_defaults_to_content_length() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("hip2.dat.gz");
        let server = Server::start(vec![Reply::ok(b"0123456789")]);
        let size = downloader(1000, 0).download(&server.url, &dest, None);
        server.join();
        assert_eq!(size.unwrap(), 10);
        assert_eq!(fs::read(&dest).unwrap(), b"0123456789");
    }

    #[test]
    fn test_retries_after_failure() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let replies = vec![
            Reply::status("500 Internal Server Error"),
            Reply::ok(b"abc"),
        ];
        let server = Server::start(replies);
        let size = downloader(1000, 1).download(&server.url, &dest, Some(3));
        assert_eq!(server.join(), 2);
        assert_eq!(size.unwrap(), 3);
        assert_eq!(file_names(dir.path()), ["a.csv.gz"]);
    }

    #[test]
    fn test_gives_up_after_the_last_retry() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let replies = vec![Reply::status("503 Service Unavailable"); 3];
        let server = Server::start(replies);
        let err = downloader(1000, 2)
            .download(&server.url, &dest, Some(3))
            .unwrap_err();
        assert_eq!(server.join(), 3);
        let message = format!("{:#}", err);
        assert!(message.contains("failed after 3 tries"), "{message}");
        assert!(message.contains("HTTP 503"), "{message}");
        assert!(file_names(dir.path()).is_empty());
    }

    #[test]
    fn test_get_text_rejects_invalid_utf8() {
        let server = Server::start(vec![Reply::ok(b"<Key>a\xff</Key>")]);
        let err = downloader(1000, 0).get_text(&server.url).unwrap_err();
        server.join();
        assert!(format!("{:#}", err).contains("not UTF-8"), "{err:#}");
    }

    #[test]
    fn test_backoff_doubles_up_to_the_maximum() {
        let policy = RetryPolicy {
            retries: 8,
            first_delay: Duration::from_secs(2),
            max_delay: Duration::from_secs(60),
        };
        let delays: Vec<u64> = (1..=7)
            .map(|r| backoff_delay(&policy, r).as_secs())
            .collect();
        assert_eq!(delays, [2, 4, 8, 16, 32, 60, 60]);
        assert_eq!(backoff_delay(&policy, 40), Duration::from_secs(60));
    }

    #[test]
    fn test_not_found_is_not_retried() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let server = Server::start(vec![Reply::status("404 Not Found"), Reply::ok(b"abc")]);
        let err = downloader(1000, 1)
            .download(&server.url, &dest, Some(3))
            .unwrap_err();
        assert_eq!(server.finish(), 1);
        assert!(
            format!("{:#}", err).contains("HTTP 404 Not Found"),
            "{err:#}"
        );
        assert!(file_names(dir.path()).is_empty());
    }

    #[test]
    fn test_too_many_requests_is_retried() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("a.csv.gz");
        let replies = vec![Reply::status("429 Too Many Requests"), Reply::ok(b"abc")];
        let server = Server::start(replies);
        let size = downloader(1000, 1).download(&server.url, &dest, Some(3));
        assert_eq!(server.join(), 2);
        assert_eq!(size.unwrap(), 3);
    }
}
