//! Hipparcos New Reduction (van Leeuwen 2007) catalog ingestion
//!
//! Parses hip2.dat fixed-width format from I/311.
//! Auto-downloads from CDS if not present locally.
//! Epoch-propagates positions from J1991.25 to J2016.0.

use crate::cli::{Cli, IngestHipparcosArgs};
use crate::gaia::GAIA_DR3_EPOCH;
use crate::hip2::Hip2Star;
use crate::http::Downloader;
use crate::intermediate::{self, Header, Kind};
use anyhow::Context;
use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

const HIP2_URL: &str = "https://cdsarc.cds.unistra.fr/ftp/I/311/hip2.dat.gz";
const CROSSMATCH_URL: &str =
    "https://cdn.gea.esac.esa.int/Gaia/gedr3/cross_match/hipparcos2_best_neighbour/Hipparcos2BestNeighbour.csv.gz";
const HIP2_FILENAME: &str = "hip2.dat";
const DOWNLOAD_RETRIES: u32 = 3;
const CROSSMATCH_FILENAME: &str = "Hipparcos2BestNeighbour.csv";

struct ParsedHip2 {
    stars: Vec<Hip2Star>,
    // Lines that parsed to a star, before the magnitude filter.
    rows_parsed: u64,
    // Non-blank lines with a field that does not parse, or a place that does not propagate.
    rows_skipped: u64,
}

struct IngestStats {
    with_gaia_match: u64,
    without_match: u64,
}

pub(crate) fn run(args: &IngestHipparcosArgs, cli: &Cli) -> anyhow::Result<()> {
    fs::create_dir_all(&args.workdir)?;
    fs::create_dir_all(&args.output)?;
    let http = Downloader::new(DOWNLOAD_RETRIES)?;
    let hip_path = ensure_file(
        &http,
        &args.workdir.join(HIP2_FILENAME),
        HIP2_URL,
        "hip2.dat (van Leeuwen 2007)",
    )?;
    let crossmatch_path = ensure_file(
        &http,
        &args.workdir.join(CROSSMATCH_FILENAME),
        CROSSMATCH_URL,
        "Hipparcos2BestNeighbour.csv (Gaia eDR3)",
    )?;
    print_plan(args, cli);
    let crossmatch = load_crossmatch(&crossmatch_path)?;
    println!("Loaded {} cross-match entries", crossmatch.len());
    let parsed = parse_hip2(&hip_path, args.mag_limit)?;
    let (unique, shared) = split_matches(&parsed.stars, &crossmatch);
    let stats = write_output(&parsed.stars, &unique, &args.output, args.mag_limit)?;
    write_shared_matches(&shared, &args.output, args.mag_limit)?;
    print_validation(&parsed.stars, &unique);
    print_summary(&parsed, &stats, shared.len());
    Ok(())
}

// Uses `path` if it exists. Otherwise decompresses `path` + ".gz", downloading it from `url`
// first when it is missing too.
fn ensure_file(http: &Downloader, path: &Path, url: &str, label: &str) -> anyhow::Result<PathBuf> {
    if path.exists() {
        println!("Found {}: {:?}", label, path);
        return Ok(path.to_path_buf());
    }
    let gz_path = path.with_extension(gz_extension(path));
    if gz_path.exists() {
        println!("Decompressing {:?}...", gz_path);
        decompress_gz(&gz_path, path)
            .with_context(|| format!("{:?} is corrupt. Delete it to download it again", gz_path))?;
        return Ok(path.to_path_buf());
    }
    println!("{} not found at {:?}", label, path);
    println!("Downloading and decompressing: {}", url);
    // A body without a Content-Length that ends early passes the size check, so the
    // download is kept only when the whole gzip stream decodes.
    let size = http.download_and_check(url, &gz_path, None, |part| decompress_gz(part, path))?;
    println!("Downloaded {:.1} MB", size as f64 / 1_048_576.0);
    Ok(path.to_path_buf())
}

fn gz_extension(path: &Path) -> String {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    format!("{}.gz", ext)
}

// The decoder checks the gzip CRC32 at the end of the stream, so a corrupt file fails here
// and leaves no output.
fn decompress_gz(gz_path: &Path, dest: &Path) -> anyhow::Result<()> {
    let gz_file = File::open(gz_path).with_context(|| format!("Failed to open {:?}", gz_path))?;
    let mut decoder = flate2::read::GzDecoder::new(BufReader::new(gz_file));
    intermediate::write_via_temp(dest, |writer| {
        std::io::copy(&mut decoder, writer)
            .with_context(|| format!("Failed to decompress {:?}", gz_path))?;
        Ok(())
    })
}

fn print_plan(args: &IngestHipparcosArgs, cli: &Cli) {
    println!("\n=== Hipparcos New Reduction Ingestion ===");
    println!("Workdir: {:?}", args.workdir);
    println!("Magnitude limit: {:.1}", args.mag_limit);
    println!("Output directory: {:?}", args.output);
    println!("Verbose: {}", cli.verbose);
    println!();
}

fn load_crossmatch(path: &Path) -> anyhow::Result<HashMap<u32, i64>> {
    let file = File::open(path).context("Failed to open cross-match file")?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        if line_num == 0 && line.contains("source_id") {
            continue;
        }
        if let Some((hip, gaia_id)) = parse_crossmatch_line(&line) {
            map.insert(hip, gaia_id);
        }
    }
    Ok(map)
}

fn parse_crossmatch_line(line: &str) -> Option<(u32, i64)> {
    let mut fields = line.split(',');
    let gaia_id: i64 = fields.next()?.trim().parse().ok()?;
    let hip: u32 = fields.next()?.trim().parse().ok()?;
    Some((hip, gaia_id))
}

fn parse_hip2(path: &Path, mag_limit: f32) -> anyhow::Result<ParsedHip2> {
    let file = File::open(path).context("Failed to open hip2.dat")?;
    let mut parsed = ParsedHip2 {
        stars: Vec::with_capacity(120_000),
        rows_parsed: 0,
        rows_skipped: 0,
    };
    for (line_num, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let Some(star) = Hip2Star::parse(&line) else {
            eprintln!("Warning: skipped hip2.dat line {}", line_num + 1);
            parsed.rows_skipped += 1;
            continue;
        };
        parsed.rows_parsed += 1;
        if star.hpmag <= mag_limit {
            parsed.stars.push(star);
        }
    }
    Ok(parsed)
}

// Splits the cross-matches of the kept stars into those whose Gaia DR3 source no other kept
// star matches, keyed by HIP number, and the sources that more than one kept star matches.
fn split_matches(
    stars: &[Hip2Star],
    crossmatch: &HashMap<u32, i64>,
) -> (HashMap<u32, i64>, BTreeSet<i64>) {
    let mut stars_per_source: HashMap<i64, u32> = HashMap::new();
    for gaia_id in stars.iter().filter_map(|s| crossmatch.get(&s.hip)) {
        *stars_per_source.entry(*gaia_id).or_default() += 1;
    }
    let shared: BTreeSet<i64> = stars_per_source
        .into_iter()
        .filter_map(|(gaia_id, count)| (count > 1).then_some(gaia_id))
        .collect();
    let unique = stars
        .iter()
        .filter_map(|s| Some((s.hip, *crossmatch.get(&s.hip)?)))
        .filter(|(_, gaia_id)| !shared.contains(gaia_id))
        .collect();
    (unique, shared)
}

// A star without a unique Gaia DR3 cross-match gets its HIP number negated, which no Gaia DR3
// source_id uses. Two stars that match one source would otherwise share its ID.
fn source_id(hip: u32, unique_matches: &HashMap<u32, i64>) -> i64 {
    unique_matches.get(&hip).copied().unwrap_or(-i64::from(hip))
}

// Merge drops the Gaia record of each of these sources and keeps every Hipparcos star that
// matched it.
fn write_shared_matches(
    shared: &BTreeSet<i64>,
    output_dir: &Path,
    mag_limit: f32,
) -> anyhow::Result<()> {
    let header = Header {
        kind: Kind::SharedMatches,
        records: shared.len() as u64,
        epoch: GAIA_DR3_EPOCH,
        mag_limit,
    };
    let path = output_dir.join(Kind::SharedMatches.file_name());
    intermediate::write_via_temp(&path, |writer| {
        writer.write_all(&header.to_bytes())?;
        for &gaia_id in shared {
            let mut record = intermediate::zeroed_record();
            record.source_id = gaia_id;
            intermediate::write_record(writer, &record)?;
        }
        Ok(())
    })
}

fn write_output(
    stars: &[Hip2Star],
    unique_matches: &HashMap<u32, i64>,
    output_dir: &Path,
    mag_limit: f32,
) -> anyhow::Result<IngestStats> {
    let output_path = output_dir.join(Kind::Hipparcos.file_name());
    let header = Header {
        kind: Kind::Hipparcos,
        records: stars.len() as u64,
        epoch: GAIA_DR3_EPOCH,
        mag_limit,
    };
    let stats = intermediate::write_via_temp(&output_path, |writer| {
        writer.write_all(&header.to_bytes())?;
        write_records(writer, stars, unique_matches)
    })?;
    println!("Written {} stars to {:?}", header.records, output_path);
    Ok(stats)
}

fn write_records<W: Write>(
    writer: &mut W,
    stars: &[Hip2Star],
    unique_matches: &HashMap<u32, i64>,
) -> anyhow::Result<IngestStats> {
    let mut stats = IngestStats {
        with_gaia_match: 0,
        without_match: 0,
    };
    for star in stars {
        if unique_matches.contains_key(&star.hip) {
            stats.with_gaia_match += 1;
        } else {
            stats.without_match += 1;
        }
        intermediate::write_record(writer, &star.record(source_id(star.hip, unique_matches)))?;
    }
    Ok(stats)
}

fn print_validation(stars: &[Hip2Star], unique_matches: &HashMap<u32, i64>) {
    println!("\n=== Validation (known stars) ===");
    validate_star(stars, unique_matches, 32349, "Sirius");
    validate_star(stars, unique_matches, 91262, "Vega");
}

fn validate_star(stars: &[Hip2Star], unique_matches: &HashMap<u32, i64>, hip: u32, name: &str) {
    let Some(star) = stars.iter().find(|s| s.hip == hip) else {
        println!("HIP {} ({}): not found in filtered catalog", hip, name);
        return;
    };
    let place = star.at_gaia_epoch;
    let match_status = if unique_matches.contains_key(&hip) {
        "unique Gaia match"
    } else {
        "no unique Gaia match"
    };
    println!(
        "HIP {} ({}): RA={:.6} Dec={:.6} deg (J2016.0), Hp={:.3}, {}",
        hip, name, place.ra, place.dec, star.hpmag, match_status
    );
    print_source_details(star, source_id(hip, unique_matches));
}

fn print_source_details(star: &Hip2Star, source_id: i64) {
    let measured = &star.at_hipparcos_epoch;
    println!(
        "  Original (J1991.25): RA={:.6} Dec={:.6} deg",
        measured.ra.degrees(),
        measured.dec.degrees()
    );
    println!(
        "  PM: pmRA={:.2} mas/yr, pmDec={:.2} mas/yr",
        measured.pm_ra_cos_dec_mas_per_year, measured.pm_dec_mas_per_year
    );
    if let Some(bv) = star.b_v {
        println!("  B-V={:.3}", bv);
    }
    println!("  Source ID: {}", source_id);
}

fn print_summary(parsed: &ParsedHip2, stats: &IngestStats, shared_sources: usize) {
    println!("\n=== Summary ===");
    println!("Total parsed: {}", parsed.rows_parsed);
    println!(
        "Skipped (unparseable or not propagatable): {}",
        parsed.rows_skipped
    );
    println!("Kept after mag filter: {}", parsed.stars.len());
    println!("With a unique Gaia match: {}", stats.with_gaia_match);
    println!(
        "Without a unique match (negated HIP number as ID): {}",
        stats.without_match
    );
    println!(
        "Gaia sources matched by more than one kept star: {}",
        shared_sources
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hip2::tests::hip2_line;
    use crate::test_server::{Reply, Server};
    use celestial_catalog::query::catalog::{
        FLAG_HAS_PARALLAX, FLAG_HAS_PROPER_MOTION, FLAG_SOURCE_HIPPARCOS,
    };

    const ASTROMETRY: [&str; 5] = ["1.0000000000", "0.5000000000", "10.00", "10.00", "-10.00"];

    fn parse_lines(dir: &Path, lines: &[String]) -> ParsedHip2 {
        let path = dir.join("hip2.dat");
        std::fs::write(&path, lines.join("\n") + "\n").unwrap();
        parse_hip2(&path, 7.0).unwrap()
    }

    fn parse_two_lines(dir: &Path) -> ParsedHip2 {
        let lines = [
            hip2_line(1, ASTROMETRY, "5.0000"),
            hip2_line(2, ASTROMETRY, "9.0000"),
        ];
        parse_lines(dir, &lines)
    }

    #[test]
    fn test_stats_count_rows_parsed_and_rows_kept() {
        let dir = tempfile::tempdir().unwrap();
        let parsed = parse_two_lines(dir.path());
        assert_eq!((parsed.rows_parsed, parsed.stars.len()), (2, 1));
    }

    #[test]
    fn test_write_output_opens_as_hipparcos_ingest() {
        let dir = tempfile::tempdir().unwrap();
        let parsed = parse_two_lines(dir.path());
        let crossmatch = HashMap::from([(1, 77)]);
        write_output(&parsed.stars, &crossmatch, dir.path(), 7.0).unwrap();
        let (header, mut reader) = intermediate::open(dir.path(), Kind::Hipparcos).unwrap();
        let expected = Header {
            kind: Kind::Hipparcos,
            records: 1,
            epoch: 2016.0,
            mag_limit: 7.0,
        };
        assert_eq!(header, expected);
        let star = intermediate::read_record(&mut reader).unwrap();
        assert_eq!((star.source_id, star.mag), (77, 5.0));
        assert_eq!(
            star.flags,
            FLAG_SOURCE_HIPPARCOS | FLAG_HAS_PROPER_MOTION | FLAG_HAS_PARALLAX
        );
    }

    #[test]
    fn test_star_without_gaia_match_gets_negated_hip_number() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [hip2_line(50176, ASTROMETRY, "5.0000")];
        let parsed = parse_lines(dir.path(), &lines);
        let mut bytes = Vec::new();
        let stats = write_records(&mut bytes, &parsed.stars, &HashMap::new()).unwrap();
        let star = intermediate::read_record(&mut bytes.as_slice()).unwrap();
        assert_eq!(star.source_id, -50176);
        assert_eq!((stats.with_gaia_match, stats.without_match), (0, 1));
    }

    #[test]
    fn test_unparseable_line_is_skipped_and_counted() {
        let dir = tempfile::tempdir().unwrap();
        let bad_parallax = ["1.0000000000", "0.5000000000", "x", "10.00", "-10.00"];
        let lines = [
            hip2_line(1, ASTROMETRY, "5.0000"),
            hip2_line(2, bad_parallax, "5.0000"),
            String::new(),
        ];
        let parsed = parse_lines(dir.path(), &lines);
        let hips: Vec<u32> = parsed.stars.iter().map(|s| s.hip).collect();
        assert_eq!(hips, [1]);
        assert_eq!((parsed.rows_parsed, parsed.rows_skipped), (1, 1));
    }

    #[test]
    fn test_failed_write_leaves_no_temp_file() {
        let dir = tempfile::tempdir().unwrap();
        // A non-empty directory at the final name makes the rename fail after the write.
        let target = dir.path().join(Kind::SharedMatches.file_name());
        std::fs::create_dir_all(target.join("x")).unwrap();
        assert!(write_shared_matches(&BTreeSet::from([5]), dir.path(), 7.0).is_err());
        assert!(!target.with_extension("bin.tmp").exists());
    }

    fn corrupt_gzip() -> Vec<u8> {
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(b"hip2 line\n").unwrap();
        let mut gz = gz.finish().unwrap();
        let crc_at = gz.len() - 8;
        gz[crc_at] ^= 0xff;
        gz
    }

    fn file_names(dir: &Path) -> Vec<std::ffi::OsString> {
        let mut names: Vec<_> = std::fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn test_corrupt_download_leaves_no_files() {
        let dir = tempfile::tempdir().unwrap();
        let server = Server::start(vec![Reply::ok(&corrupt_gzip())]);
        let http = Downloader::new(0).unwrap();
        let dest = dir.path().join(HIP2_FILENAME);
        let err = ensure_file(&http, &dest, &server.url, "hip2.dat").unwrap_err();
        server.join();
        assert!(
            format!("{:#}", err).contains("Failed to decompress"),
            "{err:#}"
        );
        assert!(file_names(dir.path()).is_empty());
    }

    #[test]
    fn test_corrupt_local_gzip_is_kept_and_reported() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("hip2.dat.gz"), corrupt_gzip()).unwrap();
        let http = Downloader::new(0).unwrap();
        let dest = dir.path().join(HIP2_FILENAME);
        let err = ensure_file(&http, &dest, "http://127.0.0.1:9/", "hip2.dat").unwrap_err();
        assert!(format!("{:#}", err).contains("Delete it"), "{err:#}");
        assert_eq!(file_names(dir.path()), ["hip2.dat.gz"]);
    }

    #[test]
    fn test_missing_file_is_downloaded_and_decompressed() {
        let dir = tempfile::tempdir().unwrap();
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(b"hip2 line\n").unwrap();
        let server = Server::start(vec![Reply::ok(&gz.finish().unwrap())]);
        let http = Downloader::new(0).unwrap();
        let dest = dir.path().join(HIP2_FILENAME);
        let path = ensure_file(&http, &dest, &server.url, "hip2.dat").unwrap();
        assert_eq!(server.join(), 1);
        assert_eq!(std::fs::read(path).unwrap(), b"hip2 line\n");
    }

    // A server that sends no Content-Length and closes early delivers a short body that no
    // size check catches. The gzip decode must reject it before either file gets its name.
    #[test]
    fn test_truncated_download_without_length_leaves_no_files() {
        let dir = tempfile::tempdir().unwrap();
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(&[b'h'; 4000]).unwrap();
        let gz = gz.finish().unwrap();
        let server = Server::start(vec![Reply::without_length(&gz[..gz.len() - 6])]);
        let http = Downloader::new(0).unwrap();
        let dest = dir.path().join(HIP2_FILENAME);
        let result = ensure_file(&http, &dest, &server.url, "hip2.dat");
        server.join();
        assert!(result.is_err(), "{:?}", result);
        let names = file_names(dir.path());
        assert!(names.is_empty(), "{names:?}");
    }
}
