# Celestial

Astronomical computation in Rust. Eight crates covering time scales, coordinate
transforms, ephemerides, image I/O, WCS projections, telescope pointing models,
and star catalog queries.

Pure Rust at runtime. No FFI dependencies in production code. Algorithms follow
IAU 2000/2006 standards validated against ERFA/SOFA test vectors and JPL
Horizons.

## Crates

| Crate                 | What it does                                                                                 |
|-----------------------|----------------------------------------------------------------------------------------------|
| `celestial-catalog`   | HEALPix-indexed star catalog (Gaia DR3 + Hipparcos), memory-mapped, cone search              |
| `celestial-core`      | Angles, 3D vectors, rotation matrices, precession/nutation models (IAU 2000A/B, 2006A)       |
| `celestial-time`      | Eight time scales (UTC, TAI, TT, UT1, GPS, TDB, TCB, TCG), Julian dates, sidereal time       |
| `celestial-coords`    | Coordinate frames with full IAU transformation chain                                         |
| `celestial-ephemeris` | Planetary positions (VSOP2013), lunar positions (ELP/MPP02), JPL SPK kernel reader           |
| `celestial-images`    | FITS, XISF, SER reading and writing with compression, binary/ASCII tables, Bayer demosaicing |
| `celestial-wcs`       | WCS pixel-to-sky and sky-to-pixel transforms, 26 projections                                 |
| `celestial-pointing`  | TPOINT-compatible telescope pointing models, weighted least-squares fitting                  |

## Dependencies

```text
celestial-core              (no internal deps)
    |
celestial-time              (core)
    |
celestial-coords            (core, time)
    |
    +-- celestial-ephemeris  (core, time, coords)
    +-- celestial-wcs        (core, coords)
    +-- celestial-pointing   (core, time, coords)
    +-- celestial-catalog    (core, time, coords)
    |
celestial-images            (core, time, wcs)
```

No circular dependencies. `core` depends on nothing internal. `time` depends
only on `core`. `coords` depends on `core` and `time`. Everything else builds
on those three.

## Example

Convert a catalog position (ICRS) to where it appears in the local sky.

```rust
use celestial_core::{Angle, Location};
use celestial_time::{tt_from_calendar, TT};
use celestial_time::scales::conversions::ToUT1WithDeltaT;
use celestial_time::sidereal::GAST;
use celestial_coords::{ICRSPosition, CoordinateFrame, CIRSPosition};

// Sirius in ICRS (catalog coordinates)
let sirius = ICRSPosition::new(
    Angle::from_degrees(101.287),   // RA 6h 45m 8.9s
    Angle::from_degrees(-16.716),   // Dec -16d 42' 58"
).unwrap();

// Observation epoch in TT
let tt = tt_from_calendar(2024, 6, 15, 22, 30, 0.0);

// ICRS -> CIRS (applies precession, nutation, aberration, light deflection)
let cirs = CIRSPosition::from_icrs(&sirius, &tt).unwrap();

// Observer location
let observatory = Location::from_degrees(33.0, -117.0, 100.0).unwrap();

// CIRS -> hour angle -> topocentric (needs Delta-T for Earth rotation)
let delta_t = 69.2; // TT - UT1 in seconds (from IERS)
let ha = cirs.to_hour_angle(&observatory, delta_t).unwrap();
let topo = ha.to_topocentric().unwrap();

println!("Azimuth:   {}", topo.azimuth());
println!("Elevation: {}", topo.elevation());
```

The transformation chain is explicit. Each step requires its physical inputs:
TT epoch for precession/nutation, observer location for the terrestrial
conversion, Delta-T for Earth rotation angle. The type system enforces the
correct order.
