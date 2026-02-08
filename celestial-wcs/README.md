# celestial-wcs

Pure Rust implementation of World Coordinate System (WCS) transformations.

[![Crates.io](https://img.shields.io/crates/v/celestial-wcs)](https://crates.io/crates/celestial-wcs)
[![Documentation](https://docs.rs/celestial-wcs/badge.svg)](https://docs.rs/celestial-wcs)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-wcs)](https://github.com/gaker/celestial)

Convert between pixel coordinates and celestial coordinates (RA/Dec) for astronomical images. Supports all standard FITS WCS projections, distortion models (SIP, TPV, TNX), and the complete spherical rotation pipeline. No runtime FFI.

## Installation

```toml
[dependencies]
celestial-wcs = "0.1"
```

## Modules

| Module       | Purpose                                                       |
|--------------|---------------------------------------------------------------|
| `builder`    | WcsBuilder for constructing WCS from FITS headers             |
| `coordinate` | PixelCoord, IntermediateCoord, NativeCoord, CelestialCoord    |
| `linear`     | CD matrix, CRPIX, PC+CDELT linear transformations             |
| `spherical`  | 25 projection types (TAN, SIN, ARC, ZEA, etc.) and rotations  |
| `distortion` | SIP, TPV, and TNX optical distortion corrections              |
| `header`     | KeywordProvider trait for FITS header integration             |

## Projections

| Family            | Codes                                       |
|-------------------|---------------------------------------------|
| Zenithal          | TAN, SIN, ARC, STG, ZEA, AZP, SZP, ZPN, AIR |
| Cylindrical       | CAR, MER, CEA, CYP                          |
| Pseudocylindrical | SFL, PAR, MOL, AIT                          |
| Conic             | COP, COE, COD, COO                          |
| Polyconic         | BON, PCO                                    |
| Quadcube          | TSC, CSC, QSC                               |

## Example

```rust
use celestial_wcs::{Wcs, WcsBuilder, PixelCoord};

// Build WCS from FITS keywords
let wcs = WcsBuilder::new()
    .crpix([512.0, 512.0])
    .crval([180.0, 45.0])
    .cdelt([-0.001, 0.001])
    .ctype(["RA---TAN", "DEC--TAN"])
    .build()?;

// Convert pixel to sky coordinates
let pixel = PixelCoord::new(256.0, 256.0);
let sky = wcs.pixel_to_world(pixel)?;
println!("RA: {:.4}°, Dec: {:.4}°", sky.alpha().degrees(), sky.delta().degrees());

// Round-trip back to pixels
let recovered = wcs.world_to_pixel(sky)?;
```

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
