# celestial-catalog

HEALPix-indexed star catalog combining Gaia DR3 and Hipparcos. Memory-mapped for fast cone searches.

[![Crates.io](https://img.shields.io/crates/v/celestial-catalog)](https://crates.io/crates/celestial-catalog)
[![Documentation](https://docs.rs/celestial-catalog/badge.svg)](https://docs.rs/celestial-catalog)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-catalog)](https://github.com/gaker/celestial)

## Installation

As a library:

```toml
[dependencies]
celestial-catalog = "0.1"
```

For the CLI tools (`forge` and `query-catalog`):

```sh
cargo install celestial-catalog --features cli
```

## Example

```rust
use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::cone::{cone_search, ConeSearchParams};

let catalog = Catalog::open("/path/to/catalog.bin").unwrap();

let results = cone_search(&catalog, &ConeSearchParams {
    ra_deg: 83.633,
    dec_deg: 22.014,
    radius_deg: 0.5,
    max_mag: Some(12.0),
    max_results: Some(100),
    epoch: None,
}).unwrap();

for r in &results {
    println!("{} mag={:.2} dist={:.4}°", r.star.source_id, r.star.mag, r.distance_deg);
}
```

## Modules

| Module           | Purpose                                                    |
|------------------|------------------------------------------------------------|
| `query::catalog` | Memory-mapped catalog reader (mmap, zero-copy star access) |
| `query::cone`    | Cone search with optional proper-motion propagation        |
| `query::healpix` | HEALPix pixel math (`ang2pix_nest`)                        |

## Download a Pre-Built Catalog

You may [download the latest catalog](https://drive.google.com/drive/folders/1akV1qbERKQETLn6smW3-K0vGzVFgEqsB) from Google Drive. The pre-built catalog has a magnitude cutoff at 18.5.

## Data Sources & Attribution

The pre-built catalog is derived from:

- **Gaia DR3**: European Space Agency (ESA) mission Gaia ([gaia.esa.int](https://www.cosmos.esa.int/gaia)), processed by the Gaia Data Processing and Analysis Consortium ([DPAC](https://www.cosmos.esa.int/web/gaia/dpac/consortium)). Gaia DR3 is licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). See the [Gaia credit page](https://gaia.aip.de/cms/credit/) for full citation requirements.
- **Hipparcos**: ESA Hipparcos and Tycho Catalogues (ESA, 1997). Public domain.

If you use the pre-built catalog in published work, please cite Gaia DR3 per ESA's guidelines.

## Building the Catalog from Source

Building your own catalog from raw survey data requires significant disk space (~700 GB for raw Gaia CSVs) and time. The `forge` CLI handles the full pipeline: download, ingest, merge, and index. Every command takes `--verbose` (`-v`) and `--help`.

The ingest steps write to one working directory, and `merge` and `build-index` read from it. The examples use `/path/to/workdir`.

### Step 1: Download the Gaia Catalog

```sh
forge download-gaia --output /path/to/gaia/dir
```

Optional flags: `--concurrency` (parallel downloads, default 4), `--retries` (default 3) and `--limit N` (download only the first N files). Re-running skips the files that `gaia_manifest.json` in `--output` records as downloaded.

### Step 2: Ingest Raw Gaia Catalog

```sh
forge ingest-gaia \
    --path /path/to/gaia/dir \
    --output /path/to/workdir \
    --mag-limit 18.5
```

`--mag-limit` defaults to 15.0. The pre-built catalog uses 18.5. Optional flags: `--threads` (`-t`, default 0 for all cores) and `--no-concat`, which keeps the per-file `gaia_part_*.bin` files and skips writing `gaia_ingest.bin`.

Ingest skips any input file whose `gaia_part_*.bin` file already exists in `--output`. Before a rebuild, for example with a different `--mag-limit`, delete the old `gaia_part_*.bin` files.

Outputs `gaia_ingest.bin`.

#### Choosing a Magnitude Limit

Most applications don't need stars fainter than ~18.5. The pre-built catalog uses this cutoff with `--max-per-cell 150`, and its HEALPix binary is 4.5 GB. Going deeper balloons quickly.

### Step 3: Ingest Hipparcos Catalog

Hipparcos epochs are propagated to J2016.0 to match Gaia DR3.

```sh
forge ingest-hipparcos \
    --workdir /path/to/hipparcos/dir \
    --output /path/to/workdir \
    --mag-limit 7.0
```

`--mag-limit` defaults to 7.0. The command reads `hip2.dat` (Hipparcos new reduction, van Leeuwen 2007) and `Hipparcos2BestNeighbour.csv` (the Gaia cross-match) from `--workdir`. When a file is missing, it decompresses `hip2.dat.gz` or `Hipparcos2BestNeighbour.csv.gz` from `--workdir`, and downloads that `.gz` file only when it is missing too.

Outputs `hipparcos_ingest.bin` and `hipparcos_shared_matches.bin`.

### Step 4: Merge Catalogs

```sh
forge merge --verbose --workdir /path/to/workdir
```

Outputs `merged.bin`.

### Step 5: Build the HEALPix Index

```sh
forge build-index \
    --workdir /path/to/workdir \
    --output ./catalog.20260217.bin \
    --healpix-order 8 \
    --max-per-cell 150
```

`--healpix-order` sets nside = 2^order. It defaults to 8 and is at most 10. The `--max-per-cell` flag caps the number of stars per HEALPix pixel, dropping faint stars first. This trims dense regions along the galactic plane. Without it, every star is kept.

## Features

- **`cli`**: Enables the `forge` and `query-catalog` binaries (adds clap, rayon, reqwest, etc.)

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
