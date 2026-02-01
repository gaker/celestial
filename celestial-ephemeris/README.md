# celestial-ephemeris

Planetary and lunar ephemerides for astronomical calculations.

[![Crates.io](https://img.shields.io/crates/v/celestial-ephemeris)](https://crates.io/crates/celestial-ephemeris)
[![Documentation](https://docs.rs/celestial-ephemeris/badge.svg)](https://docs.rs/celestial-ephemeris)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-ephemeris)](https://github.com/gaker/celestial)

Pure Rust implementation of VSOP2013 planetary theory and ELP/MPP02
lunar theory, plus JPL SPK kernel support for high-precision
ephemerides. No runtime FFI.

## Installation

```toml
[dependencies]
celestial-ephemeris = "0.1"
```

## Modules

| Module    | Purpose                                                    |
|-----------|------------------------------------------------------------|
| `planets` | VSOP2013 planetary positions (Mercury through Pluto)       |
| `earth`   | VSOP2013 Earth heliocentric position                       |
| `sun`     | Sun position (geocentric, derived from Earth ephemeris)    |
| `moon`    | ELP/MPP02 lunar theory (geocentric Moon position)          |
| `jpl`     | JPL SPK kernel reader (DE440, DE432s, etc.)                |

## Example

```rust
use celestial_ephemeris::{Vsop2013Earth, Vsop2013Sun};
use celestial_ephemeris::planets::{Vsop2013Mars, Vsop2013Jupiter};
use celestial_ephemeris::moon::ElpMpp02Moon;
use celestial_time::{TDB, JulianDate};

// Create a TDB epoch
let tdb = TDB::from_julian_date(JulianDate::j2000());

// Planetary positions (heliocentric ICRS, AU)
let mars = Vsop2013Mars;
let (pos, vel) = mars.heliocentric_state(&tdb)?;
println!("Mars heliocentric: ({:.4}, {:.4}, {:.4}) AU", pos.x, pos.y, pos.z);

// Geocentric position
let geo_pos = mars.geocentric_position(&tdb)?;
println!("Mars geocentric: ({:.4}, {:.4}, {:.4}) AU", geo_pos.x, geo_pos.y, geo_pos.z);

// Moon position (geocentric ICRS, km)
let moon = ElpMpp02Moon::new();
let moon_pos = moon.geocentric_position_icrs(&tdb)?;
println!("Moon: ({:.1}, {:.1}, {:.1}) km", moon_pos[0], moon_pos[1], moon_pos[2]);
```

## JPL SPK Kernels

For higher precision, use JPL Development Ephemerides via SPK files:

```rust
use celestial_ephemeris::jpl::{SpkFile, bodies};

let spk = SpkFile::open("de440.bsp")?;

// Get Mars barycenter state relative to Solar System Barycenter
let jd = 2451545.0; // J2000.0
let (pos_km, vel_kms) = spk.compute_state(
    bodies::MARS_BARYCENTER,
    bodies::SOLAR_SYSTEM_BARYCENTER,
    jd
)?;

println!("Mars position: [{:.3}, {:.3}, {:.3}] km", pos_km[0], pos_km[1], pos_km[2]);
println!("Mars velocity: [{:.6}, {:.6}, {:.6}] km/s", vel_kms[0], vel_kms[1], vel_kms[2]);
```

### Supported Bodies

```rust
use celestial_ephemeris::jpl::bodies;

bodies::SOLAR_SYSTEM_BARYCENTER  // 0
bodies::MERCURY_BARYCENTER       // 1
bodies::VENUS_BARYCENTER         // 2
bodies::EARTH_MOON_BARYCENTER    // 3
bodies::MARS_BARYCENTER          // 4
bodies::JUPITER_BARYCENTER       // 5
bodies::SATURN_BARYCENTER        // 6
bodies::URANUS_BARYCENTER        // 7
bodies::NEPTUNE_BARYCENTER       // 8
bodies::PLUTO_BARYCENTER         // 9
bodies::SUN                      // 10
bodies::MOON                     // 301
bodies::EARTH                    // 399
```

## VSOP2013 Planets

All planets use the VSOP2013 theory with embedded coefficients:

| Struct            | Body                    |
|-------------------|-------------------------|
| `Vsop2013Mercury` | Mercury                 |
| `Vsop2013Venus`   | Venus                   |
| `Vsop2013Earth`   | Earth                   |
| `Vsop2013Mars`    | Mars                    |
| `Vsop2013Jupiter` | Jupiter                 |
| `Vsop2013Saturn`  | Saturn                  |
| `Vsop2013Uranus`  | Uranus                  |
| `Vsop2013Neptune` | Neptune                 |
| `Vsop2013Pluto`   | Pluto                   |
| `Vsop2013Emb`     | Earth-Moon Barycenter   |

Each provides:

- `heliocentric_position(&tdb)` - Position relative to Sun (AU)
- `heliocentric_state(&tdb)` - Position and velocity (AU, AU/day)
- `geocentric_position(&tdb)` - Position relative to Earth (AU)
- `geocentric_state(&tdb)` - Position and velocity (AU, AU/day)

## ELP/MPP02 Moon

```rust
use celestial_ephemeris::moon::ElpMpp02Moon;

// Standard theory
let moon = ElpMpp02Moon::new();

// With DE405 fit corrections
let moon_de405 = ElpMpp02Moon::with_de405_fit();

// Geocentric position in mean ecliptic frame (km)
let pos_ecl = moon.geocentric_position(&tdb)?;

// Geocentric position in ICRS frame (km)
let pos_icrs = moon.geocentric_position_icrs(&tdb)?;

// Full state with velocity (km, km/day)
let state = moon.geocentric_state_icrs(&tdb)?;
```

## Features

- **`serde`** - Serialization support
- **`cli`** - Build coefficient generation tools (vsop2013-gen, elpmpp02-gen)

## Coefficient Ablation

The embedded VSOP2013 coefficients are ablated (filtered) from the
full theory to reduce binary size while maintaining accuracy. The full
VSOP2013 dataset contains millions of terms; we retain only terms with
amplitude above a configurable threshold.

The `vsop2013-gen` tool (enabled with `--features cli`) handles this:

```bash
# Download original VSOP2013 data files
cargo run --features cli --bin vsop2013-gen -- download -o ./vsop2013

# Analyze term distribution at different thresholds
cargo run --features cli --bin vsop2013-gen -- analyze -i ./vsop2013 -t 1e-10

# Generate ablated Rust coefficients
cargo run --features cli --bin vsop2013-gen -- generate -i ./vsop2013 -o ./src/planetary_coefficients -t 1e-10
```

The analyze command shows reduction statistics:

```text
Planet                         Total Terms  Above Thresh    Reduction         %
--------------------------------------------------------------------------------
Mercury                            531456         12847       518609       2.4%
Venus                              420328          8234       412094       2.0%
Earth-Moon Barycenter              460874         15632       445242       3.4%
...
```

Terms are filtered by amplitude: `sqrt(S² + C²) > threshold`.
Lower thresholds retain more terms (higher accuracy, larger binary).
The embedded coefficients are tuned for ±50 years from 2026 (1976-2076),
balancing binary size against accuracy for typical observatory use.

## Accuracy

Tested against reference ephemerides:

| Theory    | Valid Range      | Position Error                              |
|-----------|------------------|---------------------------------------------|
| VSOP2013  | ±50 years        | < 5,000 km (inner), < 50,000 km (outer)     |
| ELP/MPP02 | ±50 years        | < 5 km vs JPL DE441                         |
| JPL DE440 | Kernel-dependent | Sub-meter                                   |

Note: VSOP2013 errors are relative to full-precision VSOP2013 reference data,
not JPL DE. For pointing applications, these translate to sub-arcsecond accuracy
for inner planets and ~1" for outer planets at typical observing distances.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
