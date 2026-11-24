//! Catalog merging: Hipparcos + Gaia deduplication
//!
//! A star in both catalogs keeps its Gaia DR3 record when that record has a 5- or 6-parameter
//! astrometric solution, and its Hipparcos record otherwise. Gaia DR3 measures the place at
//! J2016.0; the Hipparcos place is a 24.75-year extrapolation from J1991.25.
//!
//! A Gaia source that more than one Hipparcos star cross-matches loses its Gaia record, and
//! every one of those Hipparcos stars stays, under its negated HIP number.

use crate::cli::{Cli, MergeArgs};
use crate::intermediate::{self, Header, Kind, HEADER_SIZE};
use celestial_catalog::query::catalog::{StarRecord, FLAG_5_OR_6_PARAM, FLAG_SOURCE_HIPPARCOS};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

#[derive(Debug, Default, PartialEq, Eq)]
struct MergeStats {
    hipparcos_count: u64,
    hipparcos_replaced: u64,
    gaia_scanned: u64,
    gaia_skipped: u64,
    gaia_kept: u64,
}

// Source IDs of the Hipparcos records, the subset whose Gaia record replaces them, and the
// Gaia sources that more than one Hipparcos star matches.
struct CrossMatch {
    hipparcos_ids: HashSet<i64>,
    gaia_wins: HashSet<i64>,
    shared: HashSet<i64>,
}

impl CrossMatch {
    fn keeps_hipparcos(&self, star: &StarRecord) -> bool {
        !self.gaia_wins.contains(&star.source_id)
    }

    fn keeps_gaia(&self, star: &StarRecord) -> bool {
        let id = &star.source_id;
        !self.shared.contains(id)
            && (!self.hipparcos_ids.contains(id) || self.gaia_wins.contains(id))
    }
}

pub(crate) fn run(args: &MergeArgs, cli: &Cli) -> anyhow::Result<()> {
    validate_paths(args)?;
    print_plan(args, cli);
    let (hip_header, hipparcos_ids) = load_hipparcos_source_ids(&args.workdir)?;
    let shared = load_shared_matches(&args.workdir, &hip_header, &hipparcos_ids)?;
    let gaia_wins = find_gaia_wins(&args.workdir, &hipparcos_ids)?;
    let crossmatch = CrossMatch {
        hipparcos_ids,
        gaia_wins,
        shared,
    };
    let stats = merge_catalogs(&args.workdir, &crossmatch)?;
    print_stats(&stats, crossmatch.shared.len());
    print_samples(&args.workdir)?;
    Ok(())
}

fn validate_paths(args: &MergeArgs) -> anyhow::Result<()> {
    if !args.workdir.exists() {
        anyhow::bail!("Working directory does not exist: {:?}", args.workdir);
    }
    for kind in [Kind::Hipparcos, Kind::Gaia] {
        let path = args.workdir.join(kind.file_name());
        if !path.exists() {
            anyhow::bail!("Ingest file not found: {:?}", path);
        }
    }
    // ingest-hipparcos writes this file beside hipparcos_ingest.bin. Without it the Hipparcos
    // file comes from an older forge, whose IDs merge cannot use.
    let shared = args.workdir.join(Kind::SharedMatches.file_name());
    if !shared.exists() {
        anyhow::bail!("{shared:?} not found; rerun ingest-hipparcos to write it");
    }
    Ok(())
}

fn print_plan(args: &MergeArgs, cli: &Cli) {
    println!("=== Catalog Merge ===");
    println!("Working directory: {:?}", args.workdir);
    println!("Verbose: {}", cli.verbose);
    println!();
}

fn load_hipparcos_source_ids(dir: &Path) -> anyhow::Result<(Header, HashSet<i64>)> {
    let (header, mut reader) = intermediate::open(dir, Kind::Hipparcos)?;
    println!("Loading {} Hipparcos source IDs...", header.records);
    // open checked the count against the file length, so the capacity is bounded by the file.
    let mut ids = HashSet::with_capacity(usize::try_from(header.records)?);
    for _ in 0..header.records {
        ids.insert(intermediate::read_record(&mut reader)?.source_id);
    }
    println!("Loaded {} Hipparcos source IDs", ids.len());
    Ok((header, ids))
}

// The shared-matches file must come from the same ingest-hipparcos run as the Hipparcos
// records: same mag_limit and epoch, and no shared source among the record IDs, because that
// run gave every star on a shared source its negated HIP number.
fn load_shared_matches(
    dir: &Path,
    hipparcos: &Header,
    hipparcos_ids: &HashSet<i64>,
) -> anyhow::Result<HashSet<i64>> {
    let (header, mut reader) = intermediate::open(dir, Kind::SharedMatches)?;
    let rerun = "the two files come from different runs; rerun ingest-hipparcos";
    if header.mag_limit != hipparcos.mag_limit || header.epoch != hipparcos.epoch {
        anyhow::bail!(
            "{} has mag_limit {} at J{} and {} has mag_limit {} at J{}; {rerun}",
            Kind::SharedMatches.file_name(),
            header.mag_limit,
            header.epoch,
            Kind::Hipparcos.file_name(),
            hipparcos.mag_limit,
            hipparcos.epoch
        );
    }
    let mut ids = HashSet::new();
    for _ in 0..header.records {
        let id = intermediate::read_record(&mut reader)?.source_id;
        if hipparcos_ids.contains(&id) {
            anyhow::bail!("Shared Gaia source {id} is also a Hipparcos record ID; {rerun}");
        }
        ids.insert(id);
    }
    Ok(ids)
}

// The Hipparcos records come first in merged.bin, so the Gaia file is read once here to find
// which of them a Gaia record replaces.
fn find_gaia_wins(dir: &Path, hipparcos_ids: &HashSet<i64>) -> anyhow::Result<HashSet<i64>> {
    let (header, mut reader) = intermediate::open(dir, Kind::Gaia)?;
    println!(
        "Scanning {} Gaia records for cross-matches...",
        header.records
    );
    let mut wins = HashSet::new();
    for _ in 0..header.records {
        let star = intermediate::read_record(&mut reader)?;
        if star.flags & FLAG_5_OR_6_PARAM != 0 && hipparcos_ids.contains(&star.source_id) {
            wins.insert(star.source_id);
        }
    }
    Ok(wins)
}

// The merged file holds stars down to the fainter of the two limits, so its mag_limit is the
// larger value. Both inputs must be at the same epoch, because merge copies positions as is.
// The record count is filled in once the Gaia duplicates are known.
fn merged_header(hipparcos: &Header, gaia: &Header) -> anyhow::Result<Header> {
    if hipparcos.epoch != gaia.epoch {
        anyhow::bail!(
            "Hipparcos positions are at epoch J{} and Gaia positions at J{}; merge needs one epoch",
            hipparcos.epoch,
            gaia.epoch
        );
    }
    Ok(Header {
        kind: Kind::Merged,
        records: 0,
        epoch: gaia.epoch,
        mag_limit: hipparcos.mag_limit.max(gaia.mag_limit),
    })
}

type Input = (Header, BufReader<File>);

fn merge_catalogs(dir: &Path, crossmatch: &CrossMatch) -> anyhow::Result<MergeStats> {
    let hipparcos = intermediate::open(dir, Kind::Hipparcos)?;
    let gaia = intermediate::open(dir, Kind::Gaia)?;
    let header = merged_header(&hipparcos.0, &gaia.0)?;
    let output_path = dir.join(Kind::Merged.file_name());
    let temp_path = output_path.with_extension("bin.tmp");
    let stats = write_merged(&temp_path, header, hipparcos, gaia, crossmatch).and_then(|stats| {
        fs::rename(&temp_path, &output_path)?;
        Ok(stats)
    });
    if stats.is_err() {
        intermediate::remove_temp_file(&temp_path);
    }
    let stats = stats?;
    let total = stats.hipparcos_count - stats.hipparcos_replaced + stats.gaia_kept;
    println!("Written {} stars to {:?}", total, output_path);
    Ok(stats)
}

fn write_merged(
    path: &Path,
    mut header: Header,
    (hip_header, mut hip_reader): Input,
    (gaia_header, mut gaia_reader): Input,
    crossmatch: &CrossMatch,
) -> anyhow::Result<MergeStats> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(&[0u8; HEADER_SIZE])?;
    println!("Writing {} Hipparcos records...", hip_header.records);
    let keep = |star: &StarRecord| crossmatch.keeps_hipparcos(star);
    let hip_kept = copy_kept(&mut hip_reader, &mut writer, hip_header.records, keep)?;
    println!("Streaming {} Gaia records...", gaia_header.records);
    let keep = |star: &StarRecord| crossmatch.keeps_gaia(star);
    let gaia_kept = copy_kept(&mut gaia_reader, &mut writer, gaia_header.records, keep)?;
    header.records = hip_kept + gaia_kept;
    intermediate::finish(writer, &header)?;
    Ok(MergeStats {
        hipparcos_count: hip_header.records,
        hipparcos_replaced: hip_header.records - hip_kept,
        gaia_scanned: gaia_header.records,
        gaia_skipped: gaia_header.records - gaia_kept,
        gaia_kept,
    })
}

// Copies the records among the next `count` for which `keep` is true. Returns how many.
fn copy_kept<R: Read, W: Write>(
    reader: &mut R,
    writer: &mut W,
    count: u64,
    keep: impl Fn(&StarRecord) -> bool,
) -> anyhow::Result<u64> {
    let mut kept = 0u64;
    for _ in 0..count {
        let record = intermediate::read_record(reader)?;
        if keep(&record) {
            intermediate::write_record(writer, &record)?;
            kept += 1;
        }
    }
    Ok(kept)
}

fn print_stats(stats: &MergeStats, shared_sources: usize) {
    let hipparcos_kept = stats.hipparcos_count - stats.hipparcos_replaced;
    println!();
    println!("=== Merge Statistics ===");
    println!("Hipparcos stars loaded: {}", stats.hipparcos_count);
    println!(
        "Hipparcos stars replaced by a 5- or 6-parameter Gaia record: {}",
        stats.hipparcos_replaced
    );
    println!(
        "Gaia sources matched by more than one Hipparcos star (Gaia record dropped): {}",
        shared_sources
    );
    println!("Gaia stars scanned: {}", stats.gaia_scanned);
    println!("Gaia stars skipped (duplicates): {}", stats.gaia_skipped);
    println!("Gaia stars kept: {}", stats.gaia_kept);
    println!("Total merged: {}", hipparcos_kept + stats.gaia_kept);
}

fn print_samples(dir: &Path) -> anyhow::Result<()> {
    let (header, mut reader) = intermediate::open(dir, Kind::Merged)?;
    println!();
    println!("=== Validation (sample records) ===");
    println!("Total stars in merged catalog: {}", header.records);
    println!(
        "Epoch: J{}, magnitude limit: {:.2}",
        header.epoch, header.mag_limit
    );
    println!();
    println!("Hipparcos-flagged records:");
    print_matching(&mut reader, header.records.min(10_000), true)?;
    let (_, mut reader) = intermediate::open(dir, Kind::Merged)?;
    println!();
    println!("Gaia records (no Hipparcos flag):");
    print_matching(&mut reader, header.records.min(100_000), false)
}

// Prints the first three of `limit` records whose Hipparcos flag equals `hipparcos`.
fn print_matching(reader: &mut BufReader<File>, limit: u64, hipparcos: bool) -> anyhow::Result<()> {
    let mut found = 0;
    for index in 0..limit {
        let record = intermediate::read_record(reader)?;
        if (record.flags & FLAG_SOURCE_HIPPARCOS != 0) == hipparcos {
            print_record(index, &record);
            found += 1;
            if found >= 3 {
                break;
            }
        }
    }
    Ok(())
}

fn print_record(index: u64, record: &StarRecord) {
    println!(
        "  [{}] source_id={}, RA={:.6}, Dec={:.6}, mag={:.2}, flags=0x{:04x}",
        index, record.source_id, record.ra, record.dec, record.mag, record.flags
    );
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::cli::Commands;
    use crate::intermediate::tests::{header, star, write_file, write_hipparcos};
    use celestial_catalog::query::catalog::FLAG_NO_5PARAM;

    fn cli() -> Cli {
        Cli {
            verbose: false,
            command: Commands::Merge(MergeArgs {
                workdir: "unused".into(),
            }),
        }
    }

    pub(crate) fn merge(dir: &Path) -> anyhow::Result<()> {
        let args = MergeArgs {
            workdir: dir.to_path_buf(),
        };
        run(&args, &cli())
    }

    fn write_inputs(dir: &Path, hip: &[StarRecord], gaia: &[StarRecord]) {
        write_hipparcos(
            dir,
            header(Kind::Hipparcos, hip.len() as u64, 2016.0, 7.0),
            hip,
        );
        write_file(
            dir,
            header(Kind::Gaia, gaia.len() as u64, 2016.0, 18.5),
            gaia,
        );
    }

    #[test]
    fn test_merge_writes_fainter_mag_limit_and_input_epoch() {
        let dir = tempfile::tempdir().unwrap();
        let hip = [star(100, 10.0, 20.0, 3.0, FLAG_SOURCE_HIPPARCOS)];
        let gaia = [
            star(100, 10.0, 20.0, 3.1, 0),
            star(200, 200.0, -40.0, 15.0, 0),
        ];
        write_inputs(dir.path(), &hip, &gaia);
        merge(dir.path()).unwrap();
        let (merged, mut reader) = intermediate::open(dir.path(), Kind::Merged).unwrap();
        assert_eq!(merged, header(Kind::Merged, 2, 2016.0, 18.5));
        let ids = [0, 1].map(|_| intermediate::read_record(&mut reader).unwrap().source_id);
        assert_eq!(ids, [100, 200]);
    }

    fn assert_rejects_shared_matches(shared: Header, records: &[StarRecord]) {
        let dir = tempfile::tempdir().unwrap();
        let hip = [star(100, 10.0, 20.0, 3.0, FLAG_SOURCE_HIPPARCOS)];
        write_inputs(
            dir.path(),
            &hip,
            &[star(100, 10.1, 20.1, 3.1, FLAG_5_OR_6_PARAM)],
        );
        write_file(dir.path(), shared, records);
        let err = format!("{:#}", merge(dir.path()).unwrap_err());
        assert!(err.contains("rerun ingest-hipparcos"), "{err}");
        assert!(!dir.path().join("merged.bin").exists());
    }

    #[test]
    fn test_merge_rejects_shared_matches_from_another_mag_limit() {
        let shared = header(Kind::SharedMatches, 0, 2016.0, 12.0);
        assert_rejects_shared_matches(shared, &[]);
    }

    #[test]
    fn test_merge_rejects_shared_matches_at_another_epoch() {
        let shared = header(Kind::SharedMatches, 0, 1991.25, 7.0);
        assert_rejects_shared_matches(shared, &[]);
    }

    // After ingest-hipparcos, every star on a shared source carries its negated HIP number.
    #[test]
    fn test_merge_rejects_shared_source_that_a_hipparcos_record_carries() {
        let shared = header(Kind::SharedMatches, 1, 2016.0, 7.0);
        assert_rejects_shared_matches(shared, &[star(100, 0.0, 0.0, 0.0, 0)]);
    }

    #[test]
    fn test_merge_requires_the_shared_matches_file() {
        let dir = tempfile::tempdir().unwrap();
        write_inputs(dir.path(), &[], &[]);
        fs::remove_file(dir.path().join("hipparcos_shared_matches.bin")).unwrap();
        let err = format!("{:#}", merge(dir.path()).unwrap_err());
        assert!(err.contains("rerun ingest-hipparcos"), "{err}");
        assert!(!dir.path().join("merged.bin").exists());
    }

    #[test]
    fn test_merge_rejects_header_count_of_2_pow_60() {
        let dir = tempfile::tempdir().unwrap();
        write_inputs(dir.path(), &[], &[]);
        write_hipparcos(
            dir.path(),
            header(Kind::Hipparcos, 1 << 60, 2016.0, 7.0),
            &[],
        );
        let err = format!("{:#}", merge(dir.path()).unwrap_err());
        assert!(err.contains("claims 1152921504606846976 records"), "{err}");
    }

    #[test]
    fn test_merge_rejects_wrong_magic() {
        let dir = tempfile::tempdir().unwrap();
        write_inputs(dir.path(), &[], &[]);
        write_file(dir.path(), header(Kind::Merged, 0, 2016.0, 18.5), &[]);
        fs::rename(
            dir.path().join("merged.bin"),
            dir.path().join("gaia_ingest.bin"),
        )
        .unwrap();
        let err = format!("{:#}", merge(dir.path()).unwrap_err());
        assert!(err.contains("magic \"MERG\" is not \"GAIA\""), "{err}");
    }

    #[test]
    fn test_merge_rejects_inputs_at_different_epochs() {
        let dir = tempfile::tempdir().unwrap();
        write_inputs(dir.path(), &[], &[]);
        write_hipparcos(dir.path(), header(Kind::Hipparcos, 0, 1991.25, 7.0), &[]);
        let err = format!("{:#}", merge(dir.path()).unwrap_err());
        assert!(err.contains("epoch J1991.25"), "{err}");
    }

    #[test]
    fn test_merge_checks_epochs_before_creating_output() {
        let dir = tempfile::tempdir().unwrap();
        write_inputs(dir.path(), &[], &[]);
        write_hipparcos(dir.path(), header(Kind::Hipparcos, 0, 1991.25, 7.0), &[]);
        // A directory where the temporary file goes makes any attempt to create it fail.
        fs::create_dir(dir.path().join("merged.bin.tmp")).unwrap();
        let err = format!("{:#}", merge(dir.path()).unwrap_err());
        assert!(err.contains("epoch J1991.25"), "{err}");
    }

    #[test]
    fn test_merge_failure_leaves_no_temp_file() {
        let dir = tempfile::tempdir().unwrap();
        let hip = [star(100, 10.0, 20.0, 3.0, FLAG_SOURCE_HIPPARCOS)];
        let gaia = [star(200, 200.0, -40.0, 15.0, 0)];
        write_hipparcos(dir.path(), header(Kind::Hipparcos, 1, 1991.25, 7.0), &hip);
        write_file(dir.path(), header(Kind::Gaia, 1, 2016.0, 18.5), &gaia);
        assert!(merge(dir.path()).is_err());
        assert!(!dir.path().join("merged.bin.tmp").exists());
        assert!(!dir.path().join("merged.bin").exists());
    }

    #[test]
    fn test_merge_failure_after_writing_removes_temp_file() {
        let dir = tempfile::tempdir().unwrap();
        write_inputs(
            dir.path(),
            &[star(1, 10.0, 20.0, 3.0, 0)],
            &[star(2, 5.0, 6.0, 9.0, 0)],
        );
        // A non-empty directory at merged.bin makes the final rename fail after the write.
        fs::create_dir_all(dir.path().join("merged.bin").join("x")).unwrap();
        assert!(merge(dir.path()).is_err());
        assert!(!dir.path().join("merged.bin.tmp").exists());
    }

    // Cross-matches with a 5-parameter, a 6-parameter and a 2-parameter Gaia solution, a
    // Hipparcos star without a match and a Gaia star without a match.
    #[test]
    fn test_gaia_record_replaces_hipparcos_only_with_5_or_6_parameters() {
        let dir = tempfile::tempdir().unwrap();
        let hip_flag = FLAG_SOURCE_HIPPARCOS;
        let hip = [
            star(100, 10.0, 20.0, 3.0, hip_flag),
            star(101, 11.0, 21.0, 4.0, hip_flag),
            star(102, 12.0, 22.0, 5.0, hip_flag),
            star(-5, 13.0, 23.0, 6.0, hip_flag),
        ];
        let gaia = [
            star(100, 10.1, 20.1, 3.1, FLAG_5_OR_6_PARAM),
            star(101, 11.1, 21.1, 4.1, FLAG_5_OR_6_PARAM | FLAG_NO_5PARAM),
            star(102, 12.1, 22.1, 5.1, FLAG_NO_5PARAM),
            star(300, 14.1, 24.1, 15.0, FLAG_5_OR_6_PARAM),
        ];
        write_inputs(dir.path(), &hip, &gaia);
        merge(dir.path()).unwrap();
        let (header, mut reader) = intermediate::open(dir.path(), Kind::Merged).unwrap();
        let merged: Vec<(i64, f32, u16)> = (0..header.records)
            .map(|_| intermediate::read_record(&mut reader).unwrap())
            .map(|s| (s.source_id, s.mag, s.flags))
            .collect();
        let expected = [
            (102, 5.0, hip_flag),
            (-5, 6.0, hip_flag),
            (100, 3.1, FLAG_5_OR_6_PARAM),
            (101, 4.1, FLAG_5_OR_6_PARAM | FLAG_NO_5PARAM),
            (300, 15.0, FLAG_5_OR_6_PARAM),
        ];
        assert_eq!(merged, expected);
    }
}
