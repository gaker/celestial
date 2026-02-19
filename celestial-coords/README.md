# celestial-coords

Type-safe astronomical coordinate transformations between reference frames.

[![Crates.io](https://img.shields.io/crates/v/celestial-coords)](https://crates.io/crates/celestial-coords)
[![Documentation](https://docs.rs/celestial-coords/badge.svg)](https://docs.rs/celestial-coords)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-coords)](https://github.com/gaker/celestial)

Pure Rust implementation of coordinate frame transformations with full aberration, light deflection, and Earth orientation support. Each frame is a distinct type to prevent accidental mixing. ICRS serves as the pivot for all transformations.

## Installation

```toml
[dependencies]
celestial-coords = "0.1"
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
use celestial_coords::{ICRSPosition, GalacticPosition, Distance};
use celestial_coords::transforms::CoordinateFrame;
use celestial_time::TT;

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

The full IAU 2000/2006 transformation from catalog to telescope has two paths from ICRS.

**CIRS path** (full pipeline -- precession, nutation, aberration, light deflection):

```text
ICRS (catalog)
  | frame bias + precession + nutation (IAU 2006A)
  | stellar aberration (~20.5")
  | gravitational light deflection (~1.75" max)
  v
CIRS (geocentric apparent)
  | Earth Rotation Angle
  v
TIRS
  | polar motion (EOP)
  v
ITRS (terrestrial)
```

**GCRS path** (aberration only -- no light deflection, no precession/nutation):

```text
ICRS (catalog)
  | stellar aberration only
  v
GCRS (geocentric, no Earth rotation applied)
```

Use the CIRS path for telescope pointing and observational work. GCRS is used for intermediate calculations where you need aberration correction without the full pipeline.

All transformations route through ICRS as the pivot frame. The `CoordinateFrame` trait provides:

```rust
pub trait CoordinateFrame: Sized {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition>;
    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self>;
}
```

Eight frames implement `CoordinateFrame`: `ICRSPosition`, `GCRSPosition`, `CIRSPosition`, `GalacticPosition`, `EclipticPosition`, `HeliographicStonyhurst`, `HeliographicCarrington`, `SelenographicPosition`.

Four frames do **not** implement it: `TIRSPosition`, `ITRSPosition`, `TopocentricPosition`, `HourAnglePosition`. These require Earth Orientation Parameters or an observer location beyond just an epoch, so they use dedicated conversion methods instead.

## Earth Orientation Parameters

Required for CIRS to ITRS transformations (polar motion, UT1-UTC):

```rust
use celestial_coords::eop::EopProvider;

// Bundled IERS C04 + finals2000A data (1962-present + predictions)
let provider = EopProvider::bundled()?;

// Get parameters for a specific MJD
let params = provider.get(60000.0)?;
println!("UT1-UTC = {:.4} s", params.ut1_utc);
println!("Polar motion: x={:.6}\", y={:.6}\"", params.x_p, params.y_p);
```

Additional constructors:

```rust
// C04 data only (no finals2000A predictions)
let provider = EopProvider::bundled_c04()?;

// Parse a finals2000A file from disk
let provider = EopProvider::from_finals_file("/path/to/finals2000A.data")?;

// Bundled data merged with a newer finals file (extends prediction range)
let provider = EopProvider::bundled_with_update("/path/to/finals2000A.data")?;

// From raw finals2000A text
let provider = EopProvider::from_finals_str(&text_content)?;
```

## Topocentric Observations

```rust
use celestial_coords::{TopocentricPosition, Distance};
use celestial_core::{Angle, Location};
use celestial_time::TT;

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
use celestial_coords::solar::{compute_solar_orientation, carrington_rotation_number};
use celestial_coords::lunar::compute_optical_libration;
use celestial_time::TT;

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

## Features

- **`serde`** - Serialization for coordinate types and EOP records

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
