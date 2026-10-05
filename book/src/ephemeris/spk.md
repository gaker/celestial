# JPL SPK Kernels

JPL's Development Ephemerides, such as DE440 and DE432s, are distributed as SPK
kernels, which store JPL's numerical integration as Chebyshev polynomials.
Kernels are more accurate than the built-in series and also give barycentric
positions. NAIF publishes them at
<https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/>.

## Opening a Kernel

```rust,ignore
use celestial_ephemeris::jpl::spk::SpkFile;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spk = SpkFile::open("de440.bsp")?;
    for segment in spk.segments() {
        println!(
            "{:>3} relative to {:>2}: type {}, frame {}, JD {} to {}",
            segment.body(),
            segment.center(),
            segment.data_type(),
            segment.frame(),
            segment.start().to_julian_date().to_f64(),
            segment.end().to_julian_date().to_f64(),
        );
    }
    Ok(())
}
```

`open` reads the whole file into memory and validates every segment it can
evaluate. A truncated or damaged kernel therefore fails at `open`, not at the
first lookup. Kernels in either byte order can be read.

## Positions and Velocities

`compute_state(body, center, &tdb)` returns the position of `body` relative to
`center` in km and its velocity in km/s. `compute_position` returns the
position alone. Both use the kernel's axes, which for the DE kernels are the
ICRF.

The units differ from the built-in theories, which use AU and AU/day. This
example converts the Earth's barycentric state, the one to use for aberration:

```rust,ignore
use celestial_core::constants::{AU_KM, SECONDS_PER_DAY_F64};
use celestial_ephemeris::jpl::bodies;
use celestial_ephemeris::jpl::spk::SpkFile;
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spk = SpkFile::open("de440.bsp")?;
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;

    let (position, velocity) =
        spk.compute_state(bodies::EARTH, bodies::SOLAR_SYSTEM_BARYCENTER, &tdb)?;
    let position = position / AU_KM;
    let velocity = velocity * (SECONDS_PER_DAY_F64 / AU_KM);

    println!("Earth: {:.9} AU from the barycenter", position.magnitude());
    println!("Speed: {:.9} AU/day", velocity.magnitude());
    Ok(())
}
```

## Chained Lookups

A kernel stores each body relative to one center. In the DE kernels:

- the planetary barycenters (1–9) and the Sun (10) are relative to the solar
  system barycenter (0);
- the Moon (301) and the Earth (399) are relative to the Earth-Moon barycenter (3);
- Mercury (199) and Venus (299) are relative to their own barycenters (1 and 2).

`compute_state` can relate any two bodies that the kernel connects, through
their shared centers. For example, the Moon relative to the Earth is the
Moon's offset from the Earth-Moon barycenter minus the Earth's. From each body
it follows the last segment in the file that covers the epoch, as SPICE does,
so later segments override earlier ones. A chain can be at most 32 links long.

```rust,ignore
use celestial_ephemeris::jpl::bodies;
use celestial_ephemeris::jpl::spk::SpkFile;
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spk = SpkFile::open("de440.bsp")?;
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;

    let moon = spk.compute_position(bodies::MOON, bodies::EARTH, &tdb)?;
    let mars = spk.compute_position(bodies::MARS_BARYCENTER, bodies::EARTH, &tdb)?;

    println!("Moon: {:.3} km", moon.magnitude());
    println!("Mars: {:.3} km", mars.magnitude());
    Ok(())
}
```

## Body Codes

`celestial_ephemeris::jpl::bodies` names the NAIF codes the DE kernels use:

| Constant                  | Code |
|---------------------------|------|
| `SOLAR_SYSTEM_BARYCENTER` | 0    |
| `MERCURY_BARYCENTER`      | 1    |
| `VENUS_BARYCENTER`        | 2    |
| `EARTH_MOON_BARYCENTER`   | 3    |
| `MARS_BARYCENTER`         | 4    |
| `JUPITER_BARYCENTER`      | 5    |
| `SATURN_BARYCENTER`       | 6    |
| `URANUS_BARYCENTER`       | 7    |
| `NEPTUNE_BARYCENTER`      | 8    |
| `PLUTO_BARYCENTER`        | 9    |
| `SUN`                     | 10   |
| `MERCURY`                 | 199  |
| `VENUS`                   | 299  |
| `MOON`                    | 301  |
| `EARTH`                   | 399  |

The methods take plain `i32` codes, so codes outside this list work too if the
kernel has them.

## What Can Be Evaluated

Only type 2 segments in frame 1 can be evaluated. Type 2 stores Chebyshev
position polynomials, and the velocity is their derivative. Frame 1 is the one
SPICE calls J2000, which the DE kernels use for the ICRF. Every DE planetary
kernel is made of these.

`segments()` also lists segments of other types and frames. Looking up a pair
whose chain passes through one returns `UnsupportedType` or
`UnsupportedFrame`.

## Errors

Every method returns `Result<_, SpkError>`:

| Variant                                | Cause                                                   |
|----------------------------------------|---------------------------------------------------------|
| `Io`                                   | the file can't be read                                  |
| `InvalidFormat`                        | not an SPK file, or damaged (found at `open`)           |
| `InvalidData`                          | the segments covering the epoch loop or exceed 32 links |
| `InvalidEpoch { jd }`                  | the epoch is NaN or infinite                            |
| `SegmentNotFound { body, center, jd }` | no segments connect the pair at that epoch              |
| `UnsupportedType`                      | the chain needs a segment of a type other than 2        |
| `UnsupportedFrame`                     | the chain needs a segment in a frame other than J2000   |

`SpkError` converts into `AstroError`, so `?` works in functions that return
`AstroResult`.
