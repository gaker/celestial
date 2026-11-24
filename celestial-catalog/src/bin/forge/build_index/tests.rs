use super::*;
use crate::cli::{Commands, MergeArgs};
use crate::intermediate::tests::{header, star, write_file, write_hipparcos};
use crate::intermediate::Kind;
use crate::merge::tests::merge;
use celestial_catalog::query::catalog::{Catalog, StarRecord, FLAG_SOURCE_HIPPARCOS};
use celestial_catalog::query::healpix::ang2pix_nest;
use clap::Parser;

fn cli() -> Cli {
    Cli {
        verbose: false,
        command: Commands::Merge(MergeArgs {
            workdir: "unused".into(),
        }),
    }
}

fn build(dir: &Path, order: u32, max_per_cell: Option<u32>) -> anyhow::Result<Catalog> {
    let args = BuildIndexArgs {
        workdir: dir.to_path_buf(),
        healpix_order: order,
        output: dir.join("out").join("catalog.bin"),
        max_per_cell,
    };
    run(&args, &cli())?;
    Ok(Catalog::open(&args.output)?)
}

fn build_error(dir: &Path, order: u32) -> String {
    match build(dir, order, None) {
        Ok(_) => panic!("order {order} built a catalog"),
        Err(e) => format!("{e:#}"),
    }
}

fn write_merged(dir: &Path, stars: &[StarRecord]) {
    let records = stars.len() as u64;
    write_file(dir, header(Kind::Merged, records, 2016.0, 18.5), stars);
}

// Every catalog star, as (pixel, source_id, mag), in file order.
fn contents(catalog: &Catalog) -> Vec<(u64, i64, f32)> {
    (0..catalog.header().npix)
        .flat_map(|pixel| {
            let stars = catalog.stars_in_pixel(pixel);
            stars.iter().map(move |s| (pixel, s.source_id, s.mag))
        })
        .collect()
}

fn pixel(order: u32, ra: f64, dec: f64) -> u64 {
    ang2pix_nest(order, ra, dec).unwrap()
}

#[test]
fn test_ingest_merge_build_and_open() {
    let dir = tempfile::tempdir().unwrap();
    let hip = [
        star(100, 10.0, 20.0, 3.0, FLAG_SOURCE_HIPPARCOS),
        star(101, 359.9, -89.9, 6.5, FLAG_SOURCE_HIPPARCOS),
    ];
    let gaia = [
        star(100, 10.0, 20.0, 3.1, 0),
        star(200, 200.0, -40.0, 15.0, 0),
        star(201, 10.0001, 20.0001, 12.0, 0),
        star(202, 90.0, 89.99, 18.25, 0),
    ];
    write_hipparcos(dir.path(), header(Kind::Hipparcos, 2, 2016.0, 7.0), &hip);
    write_file(dir.path(), header(Kind::Gaia, 4, 2016.0, 18.5), &gaia);
    merge(dir.path()).unwrap();
    let catalog = build(dir.path(), 3, None).unwrap();
    let mut expected = vec![
        (pixel(3, 10.0, 20.0), 100, 3.0),
        (pixel(3, 10.0001, 20.0001), 201, 12.0),
        (pixel(3, 359.9, -89.9), 101, 6.5),
        (pixel(3, 200.0, -40.0), 200, 15.0),
        (pixel(3, 90.0, 89.99), 202, 18.25),
    ];
    expected.sort_by_key(|&(pixel, ..)| pixel);
    assert_eq!(contents(&catalog), expected);
    assert_eq!(catalog.header().total_stars, 5);
    assert_eq!(catalog.header().mag_limit, 18.5);
    assert_eq!(catalog.header().epoch, 2016.0);
}

#[test]
fn test_catalog_epoch_and_mag_limit_come_from_merged_header() {
    let dir = tempfile::tempdir().unwrap();
    let stars = [star(1, 30.0, 40.0, 9.0, 0)];
    write_file(dir.path(), header(Kind::Merged, 1, 2000.0, 12.25), &stars);
    let catalog = build(dir.path(), 2, None).unwrap();
    assert_eq!(catalog.header().epoch, 2000.0);
    assert_eq!(catalog.header().mag_limit, 12.25);
}

#[test]
fn test_pixels_are_sorted_brightest_first() {
    let dir = tempfile::tempdir().unwrap();
    let mags = [14.0, 2.5, 9.0, -1.0, 9.0];
    let stars: Vec<StarRecord> = (0..5)
        .map(|i| star(i, 45.0, 45.0, mags[i as usize], 0))
        .collect();
    write_merged(dir.path(), &stars);
    let catalog = build(dir.path(), 1, None).unwrap();
    let p = pixel(1, 45.0, 45.0);
    let expected = vec![
        (p, 3, -1.0),
        (p, 1, 2.5),
        (p, 2, 9.0),
        (p, 4, 9.0),
        (p, 0, 14.0),
    ];
    assert_eq!(contents(&catalog), expected);
}

#[test]
fn test_max_per_cell_keeps_the_brightest_stars() {
    let dir = tempfile::tempdir().unwrap();
    let stars = [
        star(1, 45.0, 45.0, 8.0, 0),
        star(2, 45.0, 45.0, 4.0, 0),
        star(3, 45.0, 45.0, 6.0, 0),
        star(4, 250.0, -30.0, 11.0, 0),
    ];
    write_merged(dir.path(), &stars);
    let catalog = build(dir.path(), 1, Some(2)).unwrap();
    let (north, south) = (pixel(1, 45.0, 45.0), pixel(1, 250.0, -30.0));
    let expected = vec![(north, 2, 4.0), (north, 3, 6.0), (south, 4, 11.0)];
    assert_eq!(contents(&catalog), expected);
    assert_eq!(catalog.header().total_stars, 3);
}

#[test]
fn test_empty_merged_catalog_builds_an_empty_catalog() {
    let dir = tempfile::tempdir().unwrap();
    write_merged(dir.path(), &[]);
    let catalog = build(dir.path(), 0, None).unwrap();
    assert_eq!(
        (catalog.header().npix, catalog.header().total_stars),
        (12, 0)
    );
}

#[test]
fn test_build_index_rejects_merged_count_without_records() {
    let dir = tempfile::tempdir().unwrap();
    write_file(dir.path(), header(Kind::Merged, 10, 2016.0, 18.5), &[]);
    let err = build_error(dir.path(), 2);
    assert!(
        err.contains("claims 10 records, but the file is 32 bytes"),
        "{err}"
    );
}

#[test]
fn test_build_index_rejects_merged_file_shorter_than_header() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("merged.bin"), b"MERG").unwrap();
    let err = build_error(dir.path(), 2);
    assert!(err.contains("shorter than the 32-byte header"), "{err}");
}

#[test]
fn test_build_index_rejects_ingest_file_renamed_to_merged() {
    let dir = tempfile::tempdir().unwrap();
    write_file(dir.path(), header(Kind::Gaia, 0, 2016.0, 18.5), &[]);
    std::fs::rename(
        dir.path().join("gaia_ingest.bin"),
        dir.path().join("merged.bin"),
    )
    .unwrap();
    let err = build_error(dir.path(), 2);
    assert!(err.contains("magic \"GAIA\" is not \"MERG\""), "{err}");
}

#[test]
fn test_build_index_stores_ra_725_reduced() {
    let dir = tempfile::tempdir().unwrap();
    write_merged(
        dir.path(),
        &[star(7, 725.0, 10.0, 9.0, 0), star(8, -30.0, 5.0, 9.5, 0)],
    );
    let catalog = build(dir.path(), 2, None).unwrap();
    let first = |ra, dec| catalog.stars_in_pixel(pixel(2, ra, dec))[0];
    assert_eq!((first(5.0, 10.0).source_id, first(5.0, 10.0).ra), (7, 5.0));
    assert_eq!(
        (first(330.0, 5.0).source_id, first(330.0, 5.0).ra),
        (8, 330.0)
    );
}

#[test]
fn test_build_index_error_names_the_bad_record() {
    let dir = tempfile::tempdir().unwrap();
    write_merged(
        dir.path(),
        &[star(1, 10.0, 10.0, 9.0, 0), star(4242, 10.0, 95.0, 9.0, 0)],
    );
    let err = build_error(dir.path(), 2);
    assert!(
        err.contains("Merged record 1 (source_id 4242) has no HEALPix pixel"),
        "{err}"
    );
    assert!(err.contains("declination 95 is outside [-90, 90]"), "{err}");
}

#[test]
fn test_build_index_rejects_non_finite_position() {
    let dir = tempfile::tempdir().unwrap();
    write_merged(dir.path(), &[star(9, f64::NAN, 10.0, 9.0, 0)]);
    let err = build_error(dir.path(), 2);
    assert!(err.contains("source_id 9"), "{err}");
}

#[test]
fn test_build_index_run_rejects_orders_above_10() {
    let dir = tempfile::tempdir().unwrap();
    write_merged(dir.path(), &[]);
    for order in [11, 14, 29, 32] {
        let err = build_error(dir.path(), order);
        let expected = format!("HEALPix order {order} is above the build maximum of 10");
        assert!(err.contains(&expected), "{err}");
    }
}

fn parse_order(order: &str) -> Result<Cli, clap::Error> {
    let argv = ["forge", "build-index", "--workdir", "w", "--output", "o"];
    Cli::try_parse_from(argv.into_iter().chain(["--healpix-order", order]))
}

#[test]
fn test_build_order_above_maximum_is_rejected_at_parsing() {
    for order in ["11", "14", "32"] {
        let err = parse_order(order).err().unwrap().to_string();
        let expected = format!("HEALPix order {order} is above the build maximum of 10");
        assert!(err.contains(&expected), "{err}");
    }
    assert!(parse_order("-1").is_err());
}

#[test]
fn test_build_order_up_to_maximum_parses() {
    for order in ["0", "8", "10"] {
        let Commands::BuildIndex(args) = parse_order(order).unwrap().command else {
            panic!("build-index parsed as another command");
        };
        assert_eq!(args.healpix_order.to_string(), order);
    }
}

fn temp_files(dir: &Path) -> Vec<String> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.ends_with(".tmp"))
        .collect()
}

#[test]
fn test_build_index_rejects_non_finite_magnitude() {
    let dir = tempfile::tempdir().unwrap();
    let mags = [-f32::NAN, 1.0, 5.0, f32::NAN];
    let stars: Vec<StarRecord> = (0..4)
        .map(|i| star(i, 45.0, 45.0, mags[i as usize], 0))
        .collect();
    write_merged(dir.path(), &stars);
    let err = match build(dir.path(), 1, Some(2)) {
        Ok(_) => panic!("built a catalog with NaN magnitudes"),
        Err(e) => format!("{e:#}"),
    };
    assert!(
        err.contains("Merged record 0 (source_id 0) has magnitude NaN"),
        "{err}"
    );
    write_merged(
        dir.path(),
        &[
            star(5, 45.0, 45.0, 2.0, 0),
            star(6, 45.0, 45.0, f32::NEG_INFINITY, 0),
        ],
    );
    let err = build_error(dir.path(), 1);
    assert!(
        err.contains("Merged record 1 (source_id 6) has magnitude -inf"),
        "{err}"
    );
}

#[test]
fn test_build_index_rejects_refused_epoch_and_keeps_existing_output() {
    let dir = tempfile::tempdir().unwrap();
    write_merged(dir.path(), &[star(1, 30.0, 40.0, 9.0, 0)]);
    build(dir.path(), 2, None).unwrap();
    let stars = [star(1, 30.0, 40.0, 9.0, 0), star(2, 31.0, 41.0, 8.0, 0)];
    write_file(dir.path(), header(Kind::Merged, 2, 5000.0, 18.5), &stars);
    let err = build_error(dir.path(), 2);
    assert!(err.contains("Catalog epoch J5000 is not within"), "{err}");
    let output = dir.path().join("out");
    let kept = Catalog::open(output.join("catalog.bin")).unwrap();
    assert_eq!(
        (kept.header().total_stars, kept.header().epoch),
        (1, 2016.0)
    );
    assert_eq!(temp_files(&output), Vec::<String>::new());
}

#[test]
fn test_build_index_epoch_error_is_readable() {
    let dir = tempfile::tempdir().unwrap();
    write_file(dir.path(), header(Kind::Merged, 0, 1e300, 18.5), &[]);
    let err = build_error(dir.path(), 2);
    assert!(err.contains("Catalog epoch J1e300 is not within"), "{err}");
    assert!(err.len() < 300, "{err}");
}

#[test]
fn test_build_index_failure_after_writing_removes_temp_files() {
    for cap in [None, Some(1)] {
        let dir = tempfile::tempdir().unwrap();
        write_merged(
            dir.path(),
            &[star(1, 30.0, 40.0, 9.0, 0), star(2, 30.0, 40.0, 8.0, 0)],
        );
        // A directory at the output path makes the final rename fail after the build.
        let output = dir.path().join("out");
        std::fs::create_dir_all(output.join("catalog.bin")).unwrap();
        let err = format!("{:#}", build(dir.path(), 2, cap).err().unwrap());
        assert!(err.contains("Failed to move"), "{err}");
        assert_eq!(temp_files(&output), Vec::<String>::new(), "cap {cap:?}");
    }
}
