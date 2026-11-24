// hip2.dat lines and Gaia DR3 CSV rows through ingest-hipparcos, ingest-gaia, merge and
// build-index, then Catalog::open.

use crate::cli::{BuildIndexArgs, Cli, Commands, IngestHipparcosArgs};
use crate::gaia::tests::{csv, row};
use crate::hip2::tests::hip2_line;
use crate::ingest_gaia::tests::{ingest, write_gz};
use crate::intermediate::{self, Kind};
use crate::merge::tests::merge;
use crate::{build_index, ingest_hipparcos};
use celestial_catalog::query::catalog::{
    Catalog, FLAG_5_OR_6_PARAM, FLAG_DUPLICATED_SOURCE, FLAG_HAS_PARALLAX, FLAG_HAS_PROPER_MOTION,
    FLAG_NO_5PARAM, FLAG_SOURCE_HIPPARCOS,
};
use std::path::Path;

const CROSSMATCH: &str = "source_id,original_ext_source_id,angular_distance,\
number_of_neighbours,xm_flag
1000,32349,1.0,1,0
2000,91262,1.0,1,0
3000,11767,1.0,1,0
";

fn cli() -> Cli {
    Cli {
        verbose: false,
        command: Commands::Merge(crate::cli::MergeArgs {
            workdir: "unused".into(),
        }),
    }
}

// Sirius matches a 2-parameter Gaia record, Vega a 5-parameter one and Polaris a 6-parameter
// one. Alpha Centauri A has no match.
fn write_hipparcos_inputs(dir: &Path) {
    let lines = [
        hip2_line(
            32349,
            [
                "1.7678185359",
                "-0.2916993748",
                "379.21",
                "-546.01",
                "-1223.07",
            ],
            "-1.0876",
        ),
        hip2_line(
            91262,
            ["4.8735545728", "0.6768909262", "130.23", "200.94", "286.23"],
            "0.0868",
        ),
        hip2_line(
            71683,
            [
                "3.8383352142",
                "-1.0617735850",
                "754.81",
                "-3679.25",
                "473.67",
            ],
            "0.1373",
        ),
        hip2_line(
            11767,
            ["0.6622851337", "1.5579531082", "7.54", "44.48", "-11.85"],
            "2.1077",
        ),
    ];
    std::fs::write(dir.join("hip2.dat"), lines.join("\n") + "\n").unwrap();
    std::fs::write(dir.join("Hipparcos2BestNeighbour.csv"), CROSSMATCH).unwrap();
}

fn write_gaia_inputs(dir: &Path) {
    std::fs::create_dir(dir).unwrap();
    let rows = [
        row("1000", "101.287", "-16.716", "8.5", 3),
        row("2000", "279.2347", "38.7836", "0.5", 31),
        row("3000", "37.9546", "89.2641", "1.9", 95),
        row("4000", "200.0", "-40.0", "12.0", 31).replace("\"False\"", "\"True\""),
        row("5000", "null", "-40.0", "12.0", 31),
    ];
    write_gz(&dir.join("GaiaSource_000000-003111.csv.gz"), &csv(&rows));
}

fn run_pipeline(dir: &Path) -> Catalog {
    let hip_args = IngestHipparcosArgs {
        workdir: dir.join("hip"),
        mag_limit: 7.0,
        output: dir.join("work"),
    };
    std::fs::create_dir(&hip_args.workdir).unwrap();
    write_hipparcos_inputs(&hip_args.workdir);
    write_gaia_inputs(&dir.join("gaia"));
    ingest_hipparcos::run(&hip_args, &cli()).unwrap();
    ingest(&dir.join("gaia"), &dir.join("work"), false).unwrap();
    merge(&dir.join("work")).unwrap();
    let build_args = BuildIndexArgs {
        workdir: dir.join("work"),
        healpix_order: 3,
        output: dir.join("catalog.bin"),
        max_per_cell: None,
    };
    build_index::run(&build_args, &cli()).unwrap();
    Catalog::open(&build_args.output).unwrap()
}

#[test]
fn test_hip2_and_gaia_rows_build_the_expected_catalog() {
    let dir = tempfile::tempdir().unwrap();
    let catalog = run_pipeline(dir.path());
    let mut stars: Vec<(i64, f64, f64, f32, u16)> = (0..catalog.header().npix)
        .flat_map(|pixel| catalog.stars_in_pixel(pixel).to_vec())
        .map(|s| (s.source_id, s.ra, s.dec, s.mag, s.flags))
        .collect();
    stars.sort_by_key(|s| s.0);
    let hip = FLAG_SOURCE_HIPPARCOS | FLAG_HAS_PROPER_MOTION | FLAG_HAS_PARALLAX;
    let gaia = FLAG_HAS_PROPER_MOTION | FLAG_HAS_PARALLAX | FLAG_5_OR_6_PARAM;
    let expected = [
        (-71683, 219.86850784424806, -60.83187873405317, 0.1373, hip),
        (1000, 101.28462149217447, -16.721551631891643, -1.0876, hip),
        (2000, 279.2347, 38.7836, 0.5, gaia),
        (3000, 37.9546, 89.2641, 1.9, gaia | FLAG_NO_5PARAM),
        (4000, 200.0, -40.0, 12.0, gaia | FLAG_DUPLICATED_SOURCE),
    ];
    assert_eq!(stars, expected);
    assert_eq!(catalog.header().total_stars, 5);
    assert_eq!(catalog.header().mag_limit, 18.5);
}

// HIP 5 and HIP 6 both match Gaia source 1005, and HIP 7 alone matches 1007.
fn merge_shared_match(dir: &Path, params: u8) -> Vec<(i64, u16)> {
    let hipdir = dir.join("hip");
    std::fs::create_dir(&hipdir).unwrap();
    let astrometry = |ra| [ra, "0.2000000000", "10.00", "5.00", "-5.00"];
    let lines = [
        hip2_line(5, astrometry("0.5000000000"), "7.0000"),
        hip2_line(6, astrometry("0.6000000000"), "7.1000"),
        hip2_line(7, astrometry("0.7000000000"), "7.2000"),
    ];
    std::fs::write(hipdir.join("hip2.dat"), lines.join("\n") + "\n").unwrap();
    let crossmatch = "source_id,original_ext_source_id\n1005,5\n1005,6\n1007,7\n";
    std::fs::write(hipdir.join("Hipparcos2BestNeighbour.csv"), crossmatch).unwrap();
    let hip_args = IngestHipparcosArgs {
        workdir: hipdir,
        mag_limit: 9.0,
        output: dir.join("work"),
    };
    ingest_hipparcos::run(&hip_args, &cli()).unwrap();
    std::fs::create_dir(dir.join("gaia")).unwrap();
    let rows = [
        row("1005", "28.6", "11.4", "7.0", params),
        row("1007", "40.1", "11.4", "7.2", 3),
    ];
    write_gz(
        &dir.join("gaia/GaiaSource_000000-003111.csv.gz"),
        &csv(&rows),
    );
    ingest(&dir.join("gaia"), &dir.join("work"), false).unwrap();
    merge(&dir.join("work")).unwrap();
    let (header, mut reader) = intermediate::open(&dir.join("work"), Kind::Merged).unwrap();
    (0..header.records)
        .map(|_| intermediate::read_record(&mut reader).unwrap())
        .map(|s| (s.source_id, s.flags & FLAG_SOURCE_HIPPARCOS))
        .collect()
}

fn assert_shared_match_keeps_both_hipparcos_records(params: u8) {
    let dir = tempfile::tempdir().unwrap();
    let hip = FLAG_SOURCE_HIPPARCOS;
    let expected = [(-5, hip), (-6, hip), (1007, hip)];
    assert_eq!(
        merge_shared_match(dir.path(), params),
        expected,
        "params {params}"
    );
}

#[test]
fn test_two_hip_stars_on_a_5_or_6_parameter_source_keep_their_records() {
    assert_shared_match_keeps_both_hipparcos_records(31);
    assert_shared_match_keeps_both_hipparcos_records(95);
}

#[test]
fn test_two_hip_stars_on_a_2_parameter_source_get_unique_ids() {
    assert_shared_match_keeps_both_hipparcos_records(3);
}
