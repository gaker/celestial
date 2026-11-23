# celestial-wcs

Pure Rust implementation of World Coordinate System (WCS) transformations.

[![Crates.io](https://img.shields.io/crates/v/celestial-wcs)](https://crates.io/crates/celestial-wcs)
[![Documentation](https://docs.rs/celestial-wcs/badge.svg)](https://docs.rs/celestial-wcs)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-wcs)](https://github.com/gaker/celestial)

Convert between pixel coordinates and celestial coordinates (RA/Dec) for astronomical images. Supports all standard FITS WCS projections, distortion models (SIP, TPV), and the complete spherical rotation pipeline. No runtime FFI.

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
| `spherical`  | 26 projection types (TAN, SIN, ARC, ZEA, etc.) and rotations  |
| `distortion` | SIP and TPV optical distortion corrections                    |
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
use celestial_wcs::builder::WcsBuilder;
use celestial_wcs::coordinate::PixelCoord;
use celestial_wcs::spherical::Projection;

let wcs = WcsBuilder::new()
    .crpix(512.0, 512.0)
    .crval(180.0, 45.0)
    .cd_matrix([[-0.001, 0.0], [0.0, 0.001]])
    .projection(Projection::tan())
    .build()?;

// Pixel to sky, in degrees
let sky = wcs.pixel_to_celestial(PixelCoord::new(256.0, 256.0))?;
println!("RA: {:.4}°, Dec: {:.4}°", sky.alpha_deg(), sky.delta_deg());

// Sky back to pixel
let pixel = wcs.celestial_to_pixel(sky)?;
```

From FITS keywords, through any `KeywordProvider`:

```rust
use celestial_wcs::builder::WcsBuilder;
use celestial_wcs::header::KeywordMap;

let mut header = KeywordMap::new();
header
    .set_string("CTYPE1", "RA---TAN")
    .set_string("CTYPE2", "DEC--TAN")
    .set_float("CRPIX1", 512.0)
    .set_float("CRPIX2", 512.0)
    .set_float("CRVAL1", 180.0)
    .set_float("CRVAL2", 45.0)
    .set_float("CDELT1", -0.001)
    .set_float("CDELT2", 0.001);

let wcs = WcsBuilder::from_header(&header)?.build()?;
let (ra, dec) = wcs.pix2world(256.0, 256.0)?;
```

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
