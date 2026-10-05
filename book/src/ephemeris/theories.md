# Planets, Sun and Moon

`celestial-ephemeris` computes the planets from VSOP2013 (Simon et al. 2013) and the
Moon from ELP/MPP02 (Chapront & Francou 2003). Both theories are built into the
crate as truncated series, so they need no data files at run time. For JPL
kernels, see [JPL SPK Kernels](./spk.md).

## Conventions

Every body on this page follows the same conventions:

- **Input:** a `&TDB` epoch.
- **Output:** ICRS axes, positions in AU and velocities in AU/day.
- **Errors:** every method returns `AstroResult`.

The epoch's two parts are kept apart: the series see
`(jd1 − 2451545.0) + jd2` days from J2000, so a time of day in `jd2` keeps its
precision. `tdb_from_calendar` splits the date this way already.

Positions are **geometric**: where the body is at that instant. They include no
light time, aberration or light deflection.

**Heliocentric** means relative to the Sun's center, not the solar system
barycenter. Between 1950 and 2050 the Sun moves about the barycenter at 8.5 to
16.1 m/s. Aberration computed from heliocentric velocities is therefore off by
up to 11 mas. When that matters, take the Earth's barycentric velocity from a
JPL kernel.

## Planets

Each planet is a unit struct in `celestial_ephemeris::planets`:

| Type              | Body                  |
|-------------------|-----------------------|
| `Vsop2013Mercury` | Mercury               |
| `Vsop2013Venus`   | Venus                 |
| `Vsop2013Emb`     | Earth-Moon barycenter |
| `Vsop2013Mars`    | Mars                  |
| `Vsop2013Jupiter` | Jupiter               |
| `Vsop2013Saturn`  | Saturn                |
| `Vsop2013Uranus`  | Uranus                |
| `Vsop2013Neptune` | Neptune               |
| `Vsop2013Pluto`   | Pluto                 |

All of them have the same four methods:

| Method                                 | Returns                    |
|----------------------------------------|----------------------------|
| `heliocentric_position(&tdb)`          | position                   |
| `heliocentric_state(&tdb)`             | `(position, velocity)`     |
| `geocentric_position(&tdb, &earth)`    | position minus the Earth's |
| `geocentric_state(&tdb, &earth_state)` | state minus the Earth's    |

The geocentric methods take the Earth's heliocentric position or state rather
than computing it, so one Earth evaluation per epoch can serve every planet:

```rust,ignore
use celestial_ephemeris::earth::Vsop2013Earth;
use celestial_ephemeris::planets::{Vsop2013Jupiter, Vsop2013Mars};
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;
    let earth = Vsop2013Earth::new().heliocentric_state(&tdb)?;

    let (mars, mars_velocity) = Vsop2013Mars.geocentric_state(&tdb, &earth)?;
    let jupiter = Vsop2013Jupiter.geocentric_position(&tdb, &earth.0)?;

    println!("Mars:    {:.6} AU, {:.6} AU/day", mars.magnitude(), mars_velocity.magnitude());
    println!("Jupiter: {:.6} AU", jupiter.magnitude());
    Ok(())
}
```

## Earth and Sun

`Vsop2013Earth` has only `heliocentric_position` and `heliocentric_state`.
VSOP2013 has no Earth series. Instead the crate computes the Earth as the
Earth-Moon barycenter minus the Moon's share:

```text
earth = emb − moon × MOON_EMB_MASS_RATIO
```

Here `moon` is `ElpMpp02Moon::new()` and the ratio, the Moon's fraction of the
Earth-Moon mass, comes from `celestial_core::constants`. Because the Earth
depends on the Moon, it accepts only dates in the Moon's range (below).

`Vsop2013Sun` has all four methods, so code can treat the Sun like any other
body. Its heliocentric position and velocity are zero, and its geocentric
values are the negated Earth.

## Moon

`ElpMpp02Moon` returns the geocentric Moon:

```rust,ignore
use celestial_ephemeris::moon::ElpMpp02Moon;
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;

    let moon = ElpMpp02Moon::new();
    let position = moon.geocentric_position(&tdb)?;
    let (_, velocity) = moon.geocentric_state(&tdb)?;

    println!("Moon: {:.8} AU, {:.8} AU/day", position.magnitude(), velocity.magnitude());
    Ok(())
}
```

The authors fitted the theory's constants two ways, and each way has its own
constructor:

| Constructor        | Constants fitted to | Rotation to ICRS                                   |
|--------------------|---------------------|----------------------------------------------------|
| `new()`            | lunar laser ranging | ε = 84381.406″ (IAU 2006 obliquity)                |
| `with_de405_fit()` | JPL's DE405         | ε = 84381.4096″, φ = −0.05028″ (fitted with DE405) |

`Vsop2013Earth` uses `new()`. Over 1950–2050, the full series is within 60 m of
DE432s with the laser-ranging fit and within 21 m with the DE405 fit. The
shipped tables drop terms worth at most about 19 m, so they miss by 67 m and
36 m (see [Accuracy](./accuracy.md)). Given the same truncated series, the
crate reproduces the authors' Fortran, compiled without FMA, bit for bit.

## Date Ranges

Each theory rejects epochs outside its published range:

| Theory    | Bodies                                         | Years          | TDB Julian dates       |
|-----------|------------------------------------------------|----------------|------------------------|
| VSOP2013  | Mercury to Neptune, Earth-Moon barycenter, Sun | −4000 to +8000 | 260045.0 to 4643045.0  |
| VSOP2013  | Pluto                                          | 0 to +4000     | 1721045.0 to 3182045.0 |
| ELP/MPP02 | Moon, Earth                                    | −3000 to +3000 | 625295.0 to 2816795.0  |

Both ends are included. An epoch outside the range returns
`AstroError::MathError` with kind `OutOfRange`, and the message names the
theory and its years. A NaN or infinite epoch returns kind `NotFinite`:

```rust,ignore
use celestial_core::errors::{AstroError, MathErrorKind};
use celestial_ephemeris::planets::Vsop2013Pluto;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

fn main() {
    let year_5000 = TDB::from_julian_date(JulianDate::new(3547295.0, 0.0));
    match Vsop2013Pluto.heliocentric_position(&year_5000) {
        Err(AstroError::MathError { kind: MathErrorKind::OutOfRange, message, .. }) => {
            println!("{message}");
        }
        other => println!("unexpected: {other:?}"),
    }
}
```

The ranges are the theories' own. The shipped tables were measured only over
1950–2050; see [Accuracy](./accuracy.md#outside-19502050).

## Light Time

Because the output is geometric, light time is the caller's job. To see a planet
from the Earth at `t`, evaluate the planet at `t − τ`, where `τ` is the
light-travel time. The usual fixed-point loop converges in a few passes:

```rust,ignore
use celestial_core::constants::SPEED_OF_LIGHT_AU_PER_DAY;
use celestial_ephemeris::earth::Vsop2013Earth;
use celestial_ephemeris::planets::Vsop2013Mars;
use celestial_time::scales::tdb::tdb_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tdb = tdb_from_calendar(2026, 10, 5, 0, 0, 0.0)?;
    let earth = Vsop2013Earth::new().heliocentric_position(&tdb)?;

    let mut mars = Vsop2013Mars.geocentric_position(&tdb, &earth)?;
    for _ in 0..3 {
        let tau = mars.magnitude() / SPEED_OF_LIGHT_AU_PER_DAY;
        let emitted = tdb.add_days(-tau);
        mars = Vsop2013Mars.heliocentric_position(&emitted)? - earth;
    }
    println!("Mars, corrected for light time: {:.6} AU", mars.magnitude());
    Ok(())
}
```

Both positions here are heliocentric, so this ignores how far the Sun moves
during `τ`. That error is at most the Sun's speed over the speed of light,
11 mas.
