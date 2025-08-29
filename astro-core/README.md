# astro-core

Basic types for astronomical calculations in Rust.

[![Crates.io](https://img.shields.io/crates/v/astro-core)](https://crates.io/crates/astro-core)
[![Documentation](https://docs.rs/astro-core/badge.svg)](https://docs.rs/astro-core)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/astro-core)](https://github.com/gaker/celestial)

## Overview

Common types used in astronomical calculations. Includes Earth location handling, error types, and physical constants.

## What's included

- Earth locations using WGS84 ellipsoid
- Error types for astronomical calculations 
- Physical constants (speed of light, Earth radius, etc.)
- Geocentric coordinate conversions

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
astro-core = "0.1"
```

## Examples

### Working with Observatory Locations

```rust
use astro_core::Location;

// Create a location for Mauna Kea Observatory
let mauna_kea = Location::from_degrees(19.8283, -155.4783, 4207.0);

// Get geocentric coordinates for Earth rotation calculations
let (u, v) = mauna_kea.to_geocentric_km();
println!("Distance from Earth's axis: {:.3} km", u);
println!("Distance north of equator: {:.3} km", v);

// Validate coordinates
if mauna_kea.is_valid() {
    println!("Valid observatory location");
}
```

### Error Handling

```rust
use astro_core::{AstroError, AstroResult, MathErrorKind};

fn calculate_something() -> AstroResult<f64> {
    // Example calculation that might fail
    if some_condition {
        return Err(AstroError::math_error(
            "calculation_name",
            MathErrorKind::InvalidInput,
            "input out of range"
        ));
    }
    Ok(42.0)
}

match calculate_something() {
    Ok(result) => println!("Result: {}", result),
    Err(e) => eprintln!("Error: {}", e),
}
```

### Using Physical Constants

```rust
use astro_core::constants::*;

// Earth's equatorial radius in meters
println!("Earth radius: {} m", EARTH_RADIUS_EQUATORIAL_M);

// Speed of light
println!("Speed of light: {} m/s", SPEED_OF_LIGHT_M_PER_S);
```

## Types

### Location

Earth surface positions with WGS84 ellipsoid support:

```rust
pub struct Location {
    latitude_rad: f64,   // Geodetic latitude in radians
    longitude_rad: f64,  // Longitude in radians  
    altitude_m: f64,     // Height above WGS84 ellipsoid in meters
}
```

Key methods:
- `from_degrees(lat, lon, alt)` - Create from degrees and meters
- `to_geocentric_km()` - Get geocentric coordinates for Earth rotation
- `is_valid()` - Validate coordinate ranges

### Error Types

Standardized error handling with specific contexts:

- `AstroError::InvalidDate` - Calendar date validation failures
- `AstroError::MathError` - Numerical computation problems
- `AstroError::ExternalLibraryError` - C library function failures
- `AstroError::DataError` - External data file issues
- `AstroError::CalculationError` - General calculation failures

## Precision

- Geographic coordinates use f64 precision
- Geocentric conversions validated against ERFA/SOFA
- Physical constants match IAU recommendations

## Standards

- WGS84 ellipsoid (NIMA TR8350.2)
- IAU/CODATA physical constants
- Standard coordinate system definitions

## Related crates

- `astro-time` - Time calculations
- `astro-coords` - Coordinate transformations (planned)
- `astro-earth` - Earth rotation models (planned)

## Minimum Supported Rust Version

Rust 1.70 or later.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

at your option.

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.