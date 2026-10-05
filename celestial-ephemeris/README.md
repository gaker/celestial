# celestial-ephemeris

Planetary and lunar positions, and a reader for JPL SPK kernels.

[![Crates.io](https://img.shields.io/crates/v/celestial-ephemeris)](https://crates.io/crates/celestial-ephemeris)
[![Documentation](https://docs.rs/celestial-ephemeris/badge.svg)](https://docs.rs/celestial-ephemeris)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-ephemeris)](https://github.com/gaker/celestial)

Pure Rust implementation of the VSOP2013 planetary theory and the ELP/MPP02
lunar theory, plus a JPL SPK kernel reader. The theories' series are built in,
so they need no data files. No runtime FFI.

## Installation

```toml
[dependencies]
celestial-ephemeris = "0.1"
```

## Modules

| Module    | Purpose                                                             |
|-----------|---------------------------------------------------------------------|
| `planets` | VSOP2013 Mercury to Pluto and the Earth-Moon barycenter             |
| `earth`   | Heliocentric Earth: the Earth-Moon barycenter less the Moon's share |
| `sun`     | Geocentric Sun, given the Earth                                     |
| `moon`    | ELP/MPP02 geocentric Moon                                           |
| `jpl`     | JPL SPK kernel reader for DE440, DE432s and the other DE kernels    |

## Example

```rust
use celestial_ephemeris::earth::Vsop2013Earth;
use celestial_ephemeris::moon::ElpMpp02Moon;
use celestial_ephemeris::planets::Vsop2013Mars;
use celestial_ephemeris::sun::Vsop2013Sun;
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;

    // One Earth state serves every geocentric call.
    let earth = Vsop2013Earth::new().heliocentric_state(&tdb)?;
    let (mars, mars_velocity) = Vsop2013Mars.geocentric_state(&tdb, &earth)?;
    let sun = Vsop2013Sun.geocentric_position(&tdb, &earth.0)?;
    let moon = ElpMpp02Moon::new().geocentric_position(&tdb)?;

    println!("Mars: {:.6} AU, {:.6} AU/day", mars.magnitude(), mars_velocity.magnitude());
    println!("Sun:  {:.6} AU", sun.magnitude());
    println!("Moon: {:.8} AU", moon.magnitude());
    Ok(())
}
```

## JPL SPK Kernels

```rust
use celestial_ephemeris::jpl::bodies;
use celestial_ephemeris::jpl::spk::SpkFile;
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spk = SpkFile::open("de440.bsp")?;
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;

    // Kilometers and km/s. Pairs the kernel doesn't store directly are chained.
    let (earth, earth_velocity) =
        spk.compute_state(bodies::EARTH, bodies::SOLAR_SYSTEM_BARYCENTER, &tdb)?;
    let moon = spk.compute_position(bodies::MOON, bodies::EARTH, &tdb)?;

    println!("Earth: {:.3} km, {:.6} km/s", earth.magnitude(), earth_velocity.magnitude());
    println!("Moon:  {:.3} km", moon.magnitude());
    Ok(())
}
```

`open` reads the whole kernel into memory and validates it. Only type 2
segments in the J2000 frame can be evaluated, which covers the DE kernels.

## Features

- **`serde`** — Enables `serde` in celestial-core and celestial-time, whose types this crate takes and returns.
- **`cli`** — Builds `vsop2013-gen` and `elpmpp02-gen`, which download the theories' data files and regenerate the built-in tables.

## Design Notes

- **ICRS, AU and TDB**: the theories take a `&TDB` and return ICRS positions in AU and velocities in AU/day. SPK kernels return km and km/s.
- **Geometric positions**: no light time, aberration or light deflection is applied.
- **Heliocentric, not barycentric**: "heliocentric" means relative to the Sun's center. The Sun moves about the solar system barycenter at up to 16 m/s, so aberration computed from heliocentric velocities is off by up to 11 mas. Use a kernel for the barycentric Earth.
- **Two Moon fits**: `ElpMpp02Moon::new()` uses the constants fitted to lunar laser ranging, and `Vsop2013Earth` uses it too. `with_de405_fit()` uses the constants fitted to DE405.
- **Date ranges are enforced**: VSOP2013 accepts −4000 to +8000 (Pluto 0 to +4000), and ELP/MPP02 and the Earth accept −3000 to +3000. Epochs outside the range, and NaN or infinite epochs, return an error.
- **Truncated series**: the crate keeps 83,780 of VSOP2013's 2.6 million terms and 10,754 of ELP/MPP02's 35,901. Against DE432s over 1950–2050, directions seen from the Earth's center are within 7 to 30 mas for Mercury through Saturn, 0.2″ for Neptune, 0.8″ for Uranus, 3″ for Pluto, 38 mas for the Moon and 2.3 mas for the Sun.

The [Ephemeris chapter](https://github.com/gaker/celestial/tree/main/book/src/ephemeris)
of the book has the per-body accuracy, the kernel lookup rules, and how to
regenerate the tables at other thresholds.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
