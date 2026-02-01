# astro-coords

Type-safe astronomical coordinate transformations between reference frames.

[![Crates.io](https://img.shields.io/crates/v/astro-coords)](https://crates.io/crates/astro-coords)
[![Documentation](https://docs.rs/astro-coords/badge.svg)](https://docs.rs/astro-coords)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/astro-coords)](https://github.com/gaker/celestial)

Pure Rust implementation of coordinate frame transformations with full aberration, light deflection, and Earth orientation support. Each frame is a distinct type to prevent accidental mixing. ICRS serves as the pivot for all transformations.

## Installation

```toml
[dependencies]
astro-coords = "0.1"
```

## Coordinate Frames

| Frame                    | Description                                                            |
|--------------------------|------------------------------------------------------------------------|
| `ICRSPosition`           | International Celestial Reference System (catalog positions, J2000)    |
| `CIRSPosition`           | Celestial Intermediate Reference System (precession + nutation + bias) |
| `GCRSPosition`           | Geocentric Celestial Reference System                                  |
| `TIRSPosition`           | Terrestrial Intermediate Reference System                              |
| `ITRSPosition`           | International Terrestrial Reference System (ECEF)                      |
| `GalacticPosition`       | Galactic coordinates (l, b) with IAU standard pole                     |
| `EclipticPosition`       | Ecliptic coordinates with IAU 2006 obliquity                           |
| `TopocentricPosition`    | Observer-specific azimuth/elevation                                    |
| `HourAnglePosition`      | Hour angle + declination for a given observer                          |
| `HeliographicCarrington` | Solar surface coordinates (Carrington rotation)                        |
| `HeliographicStonyhurst` | Solar surface coordinates (fixed grid)                                 |
| `SelenographicPosition`  | Lunar surface coordinates                                              |

## Modules

| Module       | Purpose                                                     |
|--------------|-------------------------------------------------------------|
| `frames`     | Coordinate frame types and conversions                      |
| `transforms` | `CoordinateFrame` trait, Cartesian utilities                |
| `distance`   | Distance type with parsec/AU/ly/km conversions              |
| `eop`        | Earth Orientation Parameters (polar motion, UT1-UTC)        |
| `aberration` | Stellar aberration and gravitational light deflection       |
| `lighttime`  | Light-time correction for proper motion and radial velocity |
| `solar`      | Solar orientation (B0, L0, P angle, Carrington rotation)    |
| `lunar`      | Lunar libration and orientation                             |

## Example

```rust
use astro_coords::{ICRSPosition, GalacticPosition, Distance};
use astro_coords::transforms::CoordinateFrame;
use astro_time::TT;

// Create a position in ICRS (catalog coordinates)
let sirius = ICRSPosition::from_hours_degrees(6.752, -16.716)?;

// Transform to Galactic coordinates
let epoch = TT::j2000();
let galactic = sirius.to_galactic(&epoch)?;
println!("l = {:.2}°, b = {:.2}°", galactic.longitude().degrees(), galactic.latitude().degrees());

// With distance (parallax-derived)
let distance = Distance::from_parallax_milliarcsec(379.21)?;
let proxima = ICRSPosition::from_degrees_with_distance(217.42, -62.68, distance)?;
println!("Distance: {:.2} pc", proxima.distance().unwrap().parsecs());
```

## Transformation Chain

```
ICRS (catalog)
  ↓ frame bias + precession + nutation (IAU 2006A)
  ↓ stellar aberration (~20.5")
  ↓ gravitational light deflection (~1.75" max)
CIRS (geocentric apparent)
  ↓ Earth rotation (GAST)
TIRS
  ↓ polar motion (EOP)
ITRS (terrestrial)
```

All transformations use ICRS as the pivot frame. The `CoordinateFrame` trait provides:

```rust
pub trait CoordinateFrame: Sized {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition>;
    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self>;
}
```

## Earth Orientation Parameters

Required for CIRS to ITRS transformations (polar motion, UT1-UTC):

```rust
use astro_coords::eop::{EopManager, EopBuilder};

// Load embedded baseline data (1962-present)
let mut eop = EopBuilder::default().build()?;

// Get parameters for a specific MJD
let params = eop.get(60000.0)?;
println!("UT1-UTC = {:.4} s", params.ut1_utc);
println!("Polar motion: x={:.3}\", y={:.3}\"",
    params.x_p * 206264.806,
    params.y_p * 206264.806);
```

## Features

- **`serde`** - Serialization for coordinate types and EOP records
- **`eop-download`** - Network download of EOP data (requires `serde`)

## Topocentric Observations

```rust
use astro_coords::{TopocentricPosition, Distance};
use astro_core::{Angle, Location};
use astro_time::TT;

let observer = Location::from_degrees(19.8283, -155.4783, 4145.0)?; // Keck
let epoch = TT::j2000();

let moon_distance = Distance::from_kilometers(384400.0)?;
let moon = TopocentricPosition::with_distance(
    Angle::from_degrees(180.0),
    Angle::from_degrees(45.0),
    observer,
    epoch,
    moon_distance,
)?;

// Airmass (Rozenberg formula)
println!("Airmass: {:.2}", moon.air_mass());

// Atmospheric refraction (standard conditions)
let refraction = moon.atmospheric_refraction(1013.25, 15.0, 0.5, 0.574);
println!("Refraction: {:.1}\"", refraction.arcseconds());

// Diurnal parallax
let parallax = moon.diurnal_parallax().unwrap();
println!("Parallax: {:.1}'", parallax.arcminutes());
```

## Solar and Lunar Coordinates

```rust
use astro_coords::solar::{compute_solar_orientation, carrington_rotation_number};
use astro_coords::lunar::compute_optical_libration;
use astro_time::TT;

let epoch = TT::j2000();

// Solar orientation
let solar = compute_solar_orientation(&epoch);
println!("B0 = {:.2}°", solar.b0.degrees());
println!("L0 = {:.2}°", solar.l0.degrees());
println!("Carrington rotation: {}", carrington_rotation_number(&epoch));

// Lunar libration
let (lib_lon, lib_lat) = compute_optical_libration(&epoch);
println!("Libration: lon={:.2}°, lat={:.2}°", lib_lon.degrees(), lib_lat.degrees());
```

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
