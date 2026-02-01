# astro-images

Pure Rust astronomical image format library (FITS and XISF support).

[![Crates.io](https://img.shields.io/crates/v/astro-images)](https://crates.io/crates/astro-images)
[![Documentation](https://docs.rs/astro-images/badge.svg)](https://docs.rs/astro-images)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/astro-images)](https://github.com/gaker/celestial)

Read, write, and process FITS, XISF, and SER scientific image formats with compression support (Gzip, Rice), binary/ASCII tables, and Bayer demosaicing. No runtime FFI.

## Installation

```toml
[dependencies]
astro-images = "0.1"
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

```rust
use astro_images::{FitsFile, BitPix};

// Open a FITS file and read the primary HDU
let mut fits = FitsFile::open("m31.fits")?;
let primary = fits.primary_hdu()?;

// Access header keywords
let object = primary.header().get_string("OBJECT")?;
let exposure = primary.header().get_f64("EXPTIME")?;
println!("{}: {:.1}s exposure", object, exposure);

// Read image data as f32
let (header, data) = fits.primary_hdu_with_data::<f32>()?;
let width = header.get_i64("NAXIS1")? as usize;
let height = header.get_i64("NAXIS2")? as usize;
println!("Image: {}x{} pixels", width, height);
```

## Features

- **`parallel`** (default) — Enables parallel processing via rayon
- **`simd`** — SIMD acceleration for image operations via wide
- **`standard-formats`** — PNG/TIFF export support

## Design Notes

- **Memory-mapped I/O**: Large files use memory mapping for efficient random access without loading entire files into RAM.
- **Strict FITS compliance**: 2880-byte block alignment is validated. Non-compliant files produce errors, not silent corruption.
- **Type-safe data access**: BitPix enum and DataArray trait prevent accidental type mismatches when reading image data.
- **Streaming headers**: HDU headers are parsed on demand and cached, avoiding upfront parsing of multi-extension files.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
