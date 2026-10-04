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
use celestial_core::{angle::Angle, location::Location};
use celestial_time::scales::tt::{TT, tt_from_calendar};
use celestial_coords::eop::record::EopRecord;
use celestial_coords::frames::cirs::CIRSPosition;
use celestial_coords::frames::icrs::ICRSPosition;
use celestial_coords::transforms::CoordinateFrame;

// Sirius in ICRS (catalog coordinates)
let sirius = ICRSPosition::new(
    Angle::from_degrees(101.287),   // RA 6h 45m 8.9s
    Angle::from_degrees(-16.716),   // Dec -16d 42' 58"
).unwrap();

// Observation epoch in TT
let tt = tt_from_calendar(2024, 6, 15, 22, 30, 0.0).unwrap();

// ICRS -> CIRS (applies precession, nutation, aberration, light deflection)
let cirs = CIRSPosition::from_icrs(&sirius, &tt).unwrap();

// Observer location
let observatory = Location::from_degrees(33.0, -117.0, 100.0).unwrap();

// CIRS -> hour angle -> topocentric (needs UT1-UTC for Earth rotation).
// Values are illustrative; real ones come from IERS data via EopProvider.
// Arguments: MJD, x_p and y_p (arcsec), UT1-UTC and LOD (seconds).
let eop = EopRecord::new(60477.0, 0.2, 0.4, 0.01).unwrap().to_parameters();
let ha = cirs.to_hour_angle(&observatory, &eop).unwrap();
let topo = ha.to_topocentric().unwrap();

println!("Azimuth:   {}", topo.azimuth());
println!("Elevation: {}", topo.elevation());
```

The transformation chain is explicit. Each step requires its physical inputs:
TT epoch for precession/nutation, observer location for the terrestrial
conversion, Earth orientation parameters (UT1-UTC) for the Earth rotation
angle. The type system enforces the correct order.
