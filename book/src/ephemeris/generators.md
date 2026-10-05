# Regenerating the Tables

Two tools behind the `cli` feature generate the coefficient tables:

| Tool           | Writes                                                                            |
|----------------|-----------------------------------------------------------------------------------|
| `vsop2013-gen` | `celestial-ephemeris/src/planetary_coefficients/`, one file per body and `mod.rs` |
| `elpmpp02-gen` | `celestial-ephemeris/src/lunar_coefficients/moon.rs`                              |

Both tools have the same three subcommands:

- `download` fetches the authors' data files.
- `analyze` reports how many terms a threshold keeps.
- `generate` writes Rust source.

## Data Files

VSOP2013 comes as `VSOP2013p1.dat` to `VSOP2013p9.dat`, one file per body, from
IMCCE at <https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/>. Each
file's number is the body's number in the table below; 3 is the Earth-Moon
barycenter.
ELP/MPP02 comes as `ELP_MAIN.S1`–`S3` and `ELP_PERT.S1`–`S3` from
<http://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/>.

```sh
cargo run --release -p celestial-ephemeris --features cli --bin vsop2013-gen -- \
    download --output vsop2013
cargo run --release -p celestial-ephemeris --features cli --bin elpmpp02-gen -- \
    download --output elpmpp02
```

`vsop2013-gen download --planet 5` fetches a single body. The nine VSOP2013
files total about 305 MB.

## Thresholds

A threshold keeps only the terms whose amplitude passes it:

- **VSOP2013:** `sqrt(S² + C²) > threshold`. Units follow the variable: AU for
  the semi-major axis, radians for the mean longitude, and none for the four
  eccentricity and inclination variables.
- **ELP/MPP02:** `amplitude >= threshold`. For main-problem terms the amplitude
  is the absolute value of the first coefficient; for perturbation terms it is
  `sqrt(S² + C²)`. Units are arcseconds for longitude and latitude, and km for
  distance.

The shipped tables were cut at these thresholds, chosen as the
[Accuracy](./accuracy.md) page describes:

| Body                    | Threshold | Terms kept | Terms in the series |
|-------------------------|----------:|-----------:|--------------------:|
| 1 Mercury               |      1e-9 |      1,434 |             272,360 |
| 2 Venus                 |     1e-10 |      5,242 |             289,647 |
| 3 Earth-Moon barycenter |     1e-10 |      7,813 |             294,426 |
| 4 Mars                  |     1e-10 |     16,287 |             309,140 |
| 5 Jupiter               |      1e-9 |      7,592 |             324,608 |
| 6 Saturn                |     3e-10 |     27,689 |             350,525 |
| 7 Uranus                |      3e-8 |      5,124 |             330,581 |
| 8 Neptune               |      1e-8 |      4,870 |             322,572 |
| 9 Pluto                 |      1e-7 |      7,729 |             114,088 |
| Moon (ELP/MPP02)        |      1e-4 |     10,754 |              35,901 |

The generators use these thresholds unless `--threshold` overrides them.
Before changing one, see what it keeps:

```sh
cargo run --release -p celestial-ephemeris --features cli --bin vsop2013-gen -- \
    analyze --input vsop2013 --planet 4 --threshold 3e-11
cargo run --release -p celestial-ephemeris --features cli --bin elpmpp02-gen -- \
    analyze --input elpmpp02 --threshold 1e-5
```

## Generating

```sh
cargo run --release -p celestial-ephemeris --features cli --bin vsop2013-gen -- \
    generate --input vsop2013 --output celestial-ephemeris/src/planetary_coefficients
cargo run --release -p celestial-ephemeris --features cli --bin elpmpp02-gen -- \
    generate --input elpmpp02 --output celestial-ephemeris/src/lunar_coefficients
rustfmt --edition 2021 celestial-ephemeris/src/planetary_coefficients/*.rs \
    celestial-ephemeris/src/lunar_coefficients/moon.rs
```

`vsop2013-gen generate --planet N` regenerates a single body. The bodies
already in the output directory stay listed in the new `mod.rs`.

Each generated file starts with the command and threshold that produced it:

```text
//! VSOP2013 coefficients for Mercury
//!
//! Generated from VSOP2013p1.dat by
//! `vsop2013-gen generate --input <dir> --output <dir> --planet 1 --threshold 1e-9`
//! Terms retained: 1434 of 272360 (0.5%)
```

The output is deterministic. Running the commands above with the default
thresholds, then `rustfmt`, reproduces the shipped files byte for byte.

## Shipping New Thresholds

The tests pin the shipped tables, so a new threshold means updating the
following. Paths are relative to `celestial-ephemeris/`.

- `default_threshold` in `src/bin/vsop2013_gen/generate/mod.rs`, or
  `DEFAULT_THRESHOLD` in `src/bin/elpmpp02_gen/generate/mod.rs`. A test checks
  that these match the headers of the shipped files.
- The exact values in `src/planets/tests/state.rs`.
- The bounds in `src/planets/tests/ctl.rs` and `src/planets/tests/de432s.rs`.
  Each bound is the largest miss of the shipped table, rounded up.
- For the Moon, the vectors in `src/moon/tests/authors/vectors.rs`. They come
  from the authors' Fortran run on the same truncated series.
- The numbers on the [Accuracy](./accuracy.md) page.
