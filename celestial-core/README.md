# celestial-core

Low-level astronomical calculations for coordinate transformations.

[![Crates.io](https://img.shields.io/crates/v/celestial-core)](https://crates.io/crates/celestial-core)
[![Documentation](https://docs.rs/celestial-core/badge.svg)](https://docs.rs/celestial-core)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-core)](https://github.com/gaker/celestial)

Pure Rust implementation of IAU 2000/2006 standards for celestial mechanics: rotation matrices, nutation/precession models, angle handling, and geodetic conversions. No runtime FFI.

## Installation

```toml
[dependencies]
celestial-core = "0.1"
```

## Modules

| Module        | Purpose                                                   |
|---------------|-----------------------------------------------------------|
| `angle`       | Angle types, parsing (HMS/DMS), normalization, validation |
| `matrix`      | 3×3 rotation matrices and 3D vectors                      |
| `nutation`    | IAU 2000A/2000B/2006A nutation models                     |
| `precession`  | IAU 2000/2006 precession (Fukushima-Williams angles)      |
| `cio`         | CIO-based GCRS↔CIRS transformations                       |
| `obliquity`   | Mean obliquity of the ecliptic (IAU 1980, 2006)           |
| `location`    | Observer geodetic coordinates, geocentric conversion      |
| `constants`   | Astronomical constants (J2000, WGS84, unit conversions)   |

## Example

```rust
use celestial_core::constants::J2000_JD;
use celestial_core::nutation::NutationIAU2006A;
use celestial_core::{angle::Angle, errors::AstroError};

fn main() -> Result<(), AstroError> {
    // Nutation at J2000.0 TT
    let nutation = NutationIAU2006A::new().compute(J2000_JD, 0.0)?;
    let dpsi = Angle::from_radians(nutation.delta_psi).arcseconds();
    let deps = Angle::from_radians(nutation.delta_eps).arcseconds();
    println!("Δψ = {dpsi:.6}″, Δε = {deps:.6}″");
    Ok(())
}
```

## Features

- **`serde`** — `Serialize`/`Deserialize` for `Angle`, `Location`, `Vector3` and `RotationMatrix3`. `Location` and `RotationMatrix3` are validated on deserialize, as they are by their constructors.
- **`test-support`** — Exposes test helpers and nutation internals to the workspace's own test suites. Not a stable API.

## Design Notes

- **Two-part Julian Dates**: Functions accept `(jd1, jd2)` to preserve precision. Typically `jd1 = 2451545.0` (J2000.0) and `jd2` is days from epoch.
- **Radians internally**: All angular computations use radians. The `Angle` type provides conversion methods for degrees/HMS/DMS display.
- **Stateless models**: Nutation and precession calculators have no internal state. Call `compute(jd1, jd2)` with any epoch.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
