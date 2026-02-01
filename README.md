# Celestial

[![CI](https://github.com/gaker/celestial/workflows/Rust%20CI/badge.svg)](https://github.com/gaker/celestial/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Pure Rust astronomical computation library. No runtime FFI.

## Crates

### celestial-core

Low-level astronomical calculations: IAU 2000/2006 nutation/precession models, rotation matrices, angle handling, geodetic conversions.

```toml
[dependencies]
celestial-core = "0.1"
```

### celestial-time

8 astronomical time scales (UTC, TAI, TT, UT1, GPS, TDB, TCB, TCG) with nanosecond-precision Julian Dates, leap second support, and IAU-standard sidereal time.

```toml
[dependencies]
celestial-time = "0.1"
```

### celestial-coords

Type-safe coordinate frame transformations (ICRS, CIRS, GCRS, TIRS, ITRS, Galactic, Ecliptic, Topocentric) with aberration, light deflection, and Earth orientation support.

```toml
[dependencies]
celestial-coords = "0.1"
```

### celestial-ephemeris

Planetary and lunar ephemerides using VSOP2013 and ELP/MPP02 theories. JPL SPK kernel support.

```toml
[dependencies]
celestial-ephemeris = "0.1"
```

### celestial-images

FITS, XISF, and SER image format support with compression (Gzip, Rice), binary/ASCII tables, and Bayer demosaicing.

```toml
[dependencies]
celestial-images = "0.1"
```

### celestial-wcs

World Coordinate System (WCS) transformations for FITS images. Pixel to celestial coordinate mapping with distortion support.

```toml
[dependencies]
celestial-wcs = "0.1"
```

## Development

Early stage. API will change. Each crate has its own README with detailed documentation.

## Requirements

Rust 1.70 or later.

## License

Licensed under either Apache 2.0 or MIT.
