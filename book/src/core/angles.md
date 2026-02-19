# Angles

The `Angle` type is the most-used type in the Celestial workspace. It wraps an `f64`
storing radians and provides conversions to degrees, hours, arcminutes, and arcseconds.

All trigonometric methods use `libm` internally for cross-platform deterministic results.

## Creating Angles

Five constructors, one for each unit:

```rust,ignore
use celestial_core::Angle;
use celestial_core::constants::QUARTER_PI;

let a = Angle::from_degrees(45.0);
let b = Angle::from_radians(QUARTER_PI);
let c = Angle::from_hours(3.0);           // 3h = 45°
let d = Angle::from_arcseconds(162000.0); // 3600" per degree × 45
let e = Angle::from_arcminutes(2700.0);   // 60' per degree × 45
```

`from_radians` is the only `const fn` constructor, since radians are the internal
representation and no conversion is needed:

```rust,ignore
use celestial_core::Angle;
use celestial_core::constants::HALF_PI;

const RIGHT_ANGLE: Angle = Angle::from_radians(HALF_PI);
```

Use constants from `celestial_core::constants` rather than `std::f64::consts`.
The library defines its own `PI`, `HALF_PI`, `TWOPI`, etc. for the same reason
it uses `libm` — deterministic, cross-platform values that don't depend on the
std implementation.

### Constants

Three constants for common values:

```rust,ignore
use celestial_core::Angle;

Angle::ZERO;    // 0 radians
Angle::PI;      // π radians (180°)
Angle::HALF_PI; // π/2 radians (90°)
```

### Shorthand Functions

For terser code, free functions are re-exported at the crate root:

```rust,ignore
use celestial_core::{deg, rad, hours, arcsec, arcmin};
use celestial_core::constants::QUARTER_PI;

let a = deg(45.0);
let b = rad(QUARTER_PI);
let c = hours(3.0);
let d = arcsec(162000.0);
let e = arcmin(2700.0);
```

These are identical to calling `Angle::from_degrees`, `Angle::from_radians`, etc.

## Unit Conversions

Every constructor has a corresponding accessor:

```rust,ignore
use celestial_core::Angle;

let angle = Angle::from_degrees(45.0);

angle.radians();    // 0.7853981633974483
angle.degrees();    // 45.0
angle.hours();      // 3.0
angle.arcseconds(); // 162000.0
angle.arcminutes(); // 2700.0
```

### Hours and Degrees

Astronomy uses hours for right ascension. The relationship is:

- 24 hours = 360 degrees
- 1 hour = 15 degrees
- 1 minute of time = 15 arcminutes
- 1 second of time = 15 arcseconds

```rust,ignore
use celestial_core::Angle;

// Sirius RA: 6h 45m 8.9s
let ra = Angle::from_hours(6.0 + 45.0 / 60.0 + 8.9 / 3600.0);
ra.degrees(); // ~101.287
ra.hours();   // ~6.7525
```

## Trigonometry

Four methods, all using `libm` internally:

```rust,ignore
use celestial_core::Angle;

let angle = Angle::from_degrees(30.0);

angle.sin();     // 0.5
angle.cos();     // 0.8660254037844387
angle.tan();     // 0.5773502691896258

// When you need both sin and cos (common in rotation matrices):
let (sin, cos) = angle.sin_cos();
```

`sin_cos()` returns the tuple `(sin, cos)`. This is the method you'll reach for most
often in coordinate transforms where you need both values.

## Normalization

Two methods for wrapping angles into standard ranges:

```rust,ignore
use celestial_core::Angle;

let angle = Angle::from_degrees(270.0);

// Wraps to [0, 360°) — use for right ascension
let n = angle.normalized();
n.degrees(); // 270.0 (already in range)

// Wraps to [-180°, +180°) — use for hour angles, longitude differences
let w = angle.wrapped();
w.degrees(); // -90.0
```

The difference matters. 270° and -90° are the same direction, but RA is conventionally
positive (use `normalized`), while hour angles and longitude differences use the
shortest-arc representation (use `wrapped`).

```rust,ignore
use celestial_core::Angle;

let angle = Angle::from_degrees(-90.0);

angle.normalized().degrees(); // 270.0
angle.wrapped().degrees();    // -90.0
```

`abs()` returns the absolute value:

```rust,ignore
let neg = Angle::from_degrees(-45.0);
neg.abs().degrees(); // 45.0
```

### Free Functions on Raw Radians

If you're working with raw `f64` values in radians, three free functions are available:

```rust,ignore
use celestial_core::angle::{wrap_0_2pi, wrap_pm_pi, clamp_dec};

wrap_0_2pi(-1.0);  // wraps to [0, 2π)
wrap_pm_pi(5.0);   // wraps to [-π, +π)
clamp_dec(2.0);    // clamps to [-π/2, +π/2]
```

`clamp_dec` clamps rather than wraps — values beyond the poles are pinned to ±π/2.

## Validation

Four methods check that an angle falls within a physically meaningful range.
All return `Result<Angle, AstroError>` and reject NaN/Infinity inputs.

### Right Ascension

Cyclic — any finite angle is normalized to [0, 360°):

```rust,ignore
use celestial_core::Angle;

let ra = Angle::from_degrees(400.0);
let valid = ra.validate_right_ascension().unwrap();
valid.degrees(); // 40.0

let nan = Angle::from_radians(f64::NAN);
nan.validate_right_ascension().is_err(); // true
```

### Declination

Bounded. Standard range is [-90°, +90°]. The `beyond_pole` flag extends to [-180°, +180°]
for German equatorial mounts that can track past the pole:

```rust,ignore
use celestial_core::Angle;

let dec = Angle::from_degrees(45.0);
dec.validate_declination(false).unwrap(); // ok

let bad = Angle::from_degrees(100.0);
bad.validate_declination(false).is_err(); // true — outside [-90, +90]
bad.validate_declination(true).unwrap();  // ok — within [-180, +180]
```

### Latitude

Same as `validate_declination(false)` — range [-90°, +90°]:

```rust,ignore
use celestial_core::Angle;

let lat = Angle::from_degrees(33.0); // San Diego
lat.validate_latitude().unwrap();

let bad = Angle::from_degrees(91.0);
bad.validate_latitude().is_err(); // true
```

### Longitude

With `normalize: true`, wraps to [0, 360°) and always succeeds (for finite inputs).
With `normalize: false`, requires the angle to be within [-180°, +180°]:

```rust,ignore
use celestial_core::Angle;

let lon = Angle::from_degrees(200.0);
lon.validate_longitude(true).unwrap();   // ok — wraps
lon.validate_longitude(false).is_err();  // true — outside [-180, +180]
```

### Standalone Functions

The same validations are available as free functions in `celestial_core::angle`:

```rust,ignore
use celestial_core::Angle;
use celestial_core::angle::{
    validate_right_ascension,
    validate_declination,
    validate_latitude,
    validate_longitude,
};

let angle = Angle::from_degrees(45.0);

validate_right_ascension(angle).unwrap();
validate_declination(angle, false).unwrap();
validate_latitude(angle).unwrap();
validate_longitude(angle, false).unwrap();
```

## Parsing & Formatting

### Parsing from Strings

Two parsing systems are available. The `parse` module provides regex-based parsing
with verbose format support. The `format` module provides a lightweight splitter.

#### Explicit Unit Parsing (AngleUnits trait)

The `AngleUnits` trait is implemented on `str` and provides methods for each unit.
All return `Result<Angle, AstroError>`:

```rust,ignore
use celestial_core::angle::AngleUnits;

// Decimal values with explicit units
let a = "45.5".deg().unwrap();
let b = "0.785".rad().unwrap();
let c = "12.5".hours().unwrap();
let d = "60.0".arcmin().unwrap();
let e = "3600.0".arcsec().unwrap();

// Sexagesimal formats
let ra  = "12:34:56".hms().unwrap();    // colon-separated
let dec = "-45:30:15".dms().unwrap();   // negative DMS

// Verbose formats
let ra2  = "12h34m56s".hms().unwrap();
let ra3  = "12 hours 34 minutes 56 seconds".hms().unwrap();
let dec2 = "45d30m15s".dms().unwrap();
let dec3 = "45 degrees 30 arcmin 15 arcsec".dms().unwrap();
```

#### Auto-Detection (ParseAngle trait)

The `ParseAngle` trait provides `.to_angle()` which tries formats in order:
HMS, then DMS, then decimal degrees:

```rust,ignore
use celestial_core::angle::ParseAngle;

// Detected as HMS (colon format tries HMS first)
let a = "12:34:56".to_angle().unwrap();
a.hours(); // ~12.582

// Detected as DMS (degree marker gives it away)
let b = "45d30m15s".to_angle().unwrap();
b.degrees(); // ~45.504

// Falls back to decimal degrees
let c = "45.5".to_angle().unwrap();
c.degrees(); // 45.5
```

Colon-separated values like `12:34:56` are ambiguous — they could be HMS or DMS.
Auto-detection tries HMS first. If you know the format, use `.hms()` or `.dms()`.

#### Lightweight Parsing (parse_angle)

`parse_angle` in the `format` module is a simpler parser that returns `ParsedAngle`
(a wrapper containing the `Angle`). It tries HMS then DMS, splitting on delimiter
characters rather than using regex:

```rust,ignore
use celestial_core::angle::parse_angle;

let parsed = parse_angle("05h14m32.27s").unwrap();
let angle = parsed.angle; // extract the Angle

let dec = parse_angle("-08°12'05.9\"").unwrap();
dec.angle.degrees(); // ~-8.2016
```

### Formatting

Two formatters for astronomical notation, plus the default `Display` impl.

#### DMS (Degrees-Minutes-Seconds)

Used for declination, latitude, altitude:

```rust,ignore
use celestial_core::Angle;
use celestial_core::angle::DmsFmt;

let dec = Angle::from_degrees(-23.4392);

let fmt0 = DmsFmt { frac_digits: 0 };
fmt0.fmt(dec); // "-23° 26' 21\""

let fmt2 = DmsFmt { frac_digits: 2 };
fmt2.fmt(dec); // "-23° 26' 21.12\""
```

Sign is always shown (+ or -). Degrees and arcminutes are whole numbers.
`frac_digits` controls decimal places on the arcseconds component.

#### HMS (Hours-Minutes-Seconds)

Used for right ascension and hour angles:

```rust,ignore
use celestial_core::Angle;
use celestial_core::angle::HmsFmt;

let ra = Angle::from_hours(14.5); // 14h 30m 00s

let fmt = HmsFmt { frac_digits: 1 };
fmt.fmt(ra); // "14ʰ 30ᵐ 0.0ˢ"

// Negative angles wrap to [0, 24h)
let neg = Angle::from_hours(-1.5);
fmt.fmt(neg); // "22ʰ 30ᵐ 0.0ˢ"
```

Uses Unicode superscript markers: ʰ, ᵐ, ˢ.

#### Display Trait

The default `Display` formats as decimal degrees with 6 decimal places:

```rust,ignore
use celestial_core::Angle;

let a = Angle::from_degrees(45.123456789);
format!("{}", a); // "45.123457°"
```

## Arithmetic

Angles support addition, subtraction, scalar multiplication, scalar division,
and negation:

```rust,ignore
use celestial_core::Angle;

let a = Angle::from_degrees(30.0);
let b = Angle::from_degrees(15.0);

let sum  = a + b;   // 45°
let diff = a - b;   // 15°
let scaled = a * 2.0; // 60°
let half = a / 2.0;   // 15°
let neg  = -a;       // -30°
```

Results are not automatically normalized. If you need the result in a standard
range, call `.normalized()` or `.wrapped()` on the result.
