# Celestial

[![CI](https://github.com/gaker/celestial/workflows/Rust%20CI/badge.svg)](https://github.com/gaker/celestial/actions)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/astro-core)](https://github.com/gaker/celestial)

Astronomical computation library for Rust.

## Status

Early development. Only `astro-core` for now

## What's Available

### astro-core [![Crates.io](https://img.shields.io/crates/v/astro-core)](https://crates.io/crates/astro-core)

Basic types for astronomical calculations:

- Earth locations using WGS84 ellipsoid
- Error types for astronomical calculations
- Physical constants (speed of light, Earth radius, etc.)

```toml
[dependencies]
astro-core = "0.1"
```

```rust
use astro_core::Location;

// Create a location for Mauna Kea Observatory
let mauna_kea = Location::from_degrees(19.8283, -155.4783, 4207.0)?;

// Get geocentric coordinates for Earth rotation calculations
let (u, v) = mauna_kea.to_geocentric_km()?;
println!("Distance from Earth's axis: {:.3} km", u);
```

### erfa-sys

Low-level FFI bindings to the ERFA astronomical library. Used internally for validation and reference implementations.

## Planned

- `astro-time`: High-precision time handling (nanosecond accuracy)
- `astro-coords`: Coordinate transformations
- `astro-earth`: Earth rotation models

## Development

This is early stage development. The API will change. Don't use this. If you want to help. Reach out.

Each crate has its own README with detailed documentation.

## Requirements

Rust 1.70 or later.

## License

Licensed under either Apache 2.0 or MIT at your option.