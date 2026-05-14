# celestial-images

Pure Rust astronomical image format library (FITS and XISF support).

[![Crates.io](https://img.shields.io/crates/v/celestial-images)](https://crates.io/crates/celestial-images)
[![Documentation](https://docs.rs/celestial-images/badge.svg)](https://docs.rs/celestial-images)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-images)](https://github.com/gaker/celestial)

Read, write, and process FITS, XISF, and SER scientific image formats with compression support (Gzip, Rice), binary/ASCII tables, and Bayer demosaicing. No runtime FFI.

## Installation

```toml
[dependencies]
celestial-images = "0.1"
```

## Modules

| Module     | Purpose                                                      |
|------------|--------------------------------------------------------------|
| `core`     | BitPix, ByteOrder, error types, Result alias                 |
| `fits`     | FITS reader/writer (Primary, Image, ASCII/Binary Table HDUs) |
| `xisf`     | XISF (Extensible Image Serialization Format) reader          |
| `ser`      | SER video format reader/writer with frame timestamps         |
| `formats`  | Unified AstroImage abstraction across formats                |
| `debayer`  | Bayer pattern demosaicing (bilinear interpolation)           |
| `ricecomp` | Rice compression/decompression codec                         |

## Example

Read, inspect, and write images through the unified [`formats`](src/formats/) API. Format
is detected from the file extension — FITS and XISF work out of the box, PNG/TIFF behind
the `standard-formats` feature.

```rust,ignore
use celestial_images::formats::{AstroImage, Image};

// Open any supported format — FITS, XISF, (PNG/TIFF with `standard-formats`)
let img = Image::open("m31.fits")?;

let object = img.get_keyword("OBJECT").and_then(|k| k.value.as_ref());
let exposure = img.get_keyword("EXPTIME").and_then(|k| k.value.as_ref());
println!("{:?}: {:?}s exposure", object, exposure);
println!("{}x{}, {} channel(s)", img.width(), img.height(), img.channels());

// Build and write a new image with fluent metadata setters
let pixels: Vec<f32> = vec![0.0; 1024 * 1024];
AstroImage::new(&pixels, [1024, 1024])
    .object("M31")
    .telescope("Planewave CDK14")
    .filter("Ha")
    .exposure(300.0)
    .temperature(-10.0)
    .gain(1.0)
    .binning(1, 1)
    .date_obs("2026-04-21T20:00:00")
    .write_to("m31_ha.fits")?;       // dispatches on extension
```

For lower-level access (HDU iteration, binary tables, raw compression), the `fits` and
`xisf` modules expose the format-specific types directly.

## Features

- **`parallel`** (default) — Enables parallel processing via rayon
- **`simd`** — SIMD acceleration for image operations via wide
- **`standard-formats`** — PNG/TIFF export support

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
