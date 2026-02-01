# astro-time

Astronomical time scales and sidereal time calculations.

[![Crates.io](https://img.shields.io/crates/v/astro-time)](https://crates.io/crates/astro-time)
[![Documentation](https://docs.rs/astro-time/badge.svg)](https://docs.rs/astro-time)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/astro-time)](https://github.com/gaker/celestial)

Pure Rust implementation of 8 astronomical time scales (UTC, TAI, TT, UT1, GPS, TDB, TCB, TCG) with nanosecond-precision Julian Date handling, leap second support, and IAU-standard sidereal time calculations. No runtime FFI.

## Installation

```toml
[dependencies]
astro-time = "0.1"
```

## Modules

| Module       | Purpose                                                    |
|--------------|------------------------------------------------------------|
| `julian`     | Split Julian Date for microsecond precision                |
| `scales`     | Time scale types (UTC, TAI, TT, UT1, GPS, TDB, TCB, TCG)   |
| `sidereal`   | GMST, GAST, LMST, LAST, Earth Rotation Angle               |
| `transforms` | IAU 2000A/2000B/2006A nutation and precession models       |
| `parsing`    | ISO 8601 and calendar date parsing                         |
| `constants`  | Leap second table, epoch offsets, time scale constants     |

## Example

```rust
use astro_time::{utc_from_calendar, ToTAI, ToTT, GMST};

// Create UTC from calendar date
let utc = utc_from_calendar(2024, 6, 15, 12, 0, 0.0).unwrap();

// Convert through time scales: UTC -> TAI -> TT
let tai = utc.to_tai();
let tt = tai.to_tt();

// Compute Greenwich Mean Sidereal Time
let gmst = GMST::from_ut1_tt(utc.jd1(), utc.jd2(), tt.jd1(), tt.jd2()).unwrap();
println!("GMST = {:.6} rad", gmst.radians());
```

## Features

- **`serde`** — Enables serialization for time scale types and `JulianDate`

## Design Notes

- **Split Julian Dates**: All time scales store `(jd1, jd2)` internally. Julian Dates are ~2.4 million, but f64 has only ~15 decimal digits. Split storage preserves microsecond accuracy by keeping high-magnitude integer parts separate from fractional parts.
- **Conversions chain through TAI**: UTC -> TAI -> TT -> TCG. This ensures consistent handling of leap seconds and fixed offsets.
- **Leap seconds handled internally**: The leap second table covers 1972-present. UTC/TAI conversions account for discontinuities automatically.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
