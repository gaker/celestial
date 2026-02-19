pub const J2000_JD: f64 = 2451545.0;

pub const DAYS_PER_JULIAN_CENTURY: f64 = 36525.0;

pub const DAYS_PER_JULIAN_MILLENNIUM: f64 = 365250.0;

pub const CIRCULAR_ARCSECONDS: f64 = 1296000.0;

/// WGS84 semi-major axis in kilometers.
pub const WGS84_SEMI_MAJOR_AXIS_KM: f64 = 6378.137;

/// WGS84 first eccentricity squared: e² = (a² - b²) / a².
pub const WGS84_ECCENTRICITY_SQUARED: f64 = 6.6943799901413165e-3;

pub const NANOSECONDS_PER_SECOND: u32 = 1_000_000_000;

pub const NANOSECONDS_PER_SECOND_F64: f64 = 1_000_000_000.0;

pub const SECONDS_PER_DAY: i64 = 86_400;

pub const SECONDS_PER_DAY_F64: f64 = 86_400.0;

pub const HOURS_PER_DAY: f64 = 24.0;

pub const MINUTES_PER_DAY: f64 = 1440.0;

#[allow(clippy::excessive_precision)]
pub const ARCMIN_TO_RAD: f64 = 2.908882086657215961539535e-4;

#[allow(clippy::excessive_precision)]
pub const ARCSEC_TO_RAD: f64 = 4.848136811095359935899141e-6;

#[allow(clippy::excessive_precision)]
pub const MILLIARCSEC_TO_RAD: f64 = 4.848136811095359935899141e-9;

#[allow(clippy::excessive_precision)]
pub const MICROARCSEC_TO_RAD: f64 = 4.848136811095359935899141e-13;

pub const MJD_ZERO_POINT: f64 = 2_400_000.5;

#[allow(clippy::excessive_precision)]
#[allow(clippy::approx_constant)]
pub const PI: f64 = 3.141592653589793238462643;

#[allow(clippy::excessive_precision)]
#[allow(clippy::approx_constant)]
pub const HALF_PI: f64 = 1.5707963267948966192313216;

#[allow(clippy::excessive_precision)]
#[allow(clippy::approx_constant)]
pub const QUARTER_PI: f64 = 0.7853981633974483096156608;

#[allow(clippy::excessive_precision)]
#[allow(clippy::approx_constant)]
pub const TWOPI: f64 = 6.283185307179586476925287;

#[allow(clippy::excessive_precision)]
pub const DEG_TO_RAD: f64 = 1.745329251994329576923691e-2;

#[allow(clippy::excessive_precision)]
pub const RAD_TO_DEG: f64 = 57.29577951308232087679815;

#[allow(clippy::excessive_precision)]
pub const ARCSEC_PER_RAD: f64 = 206264.8062470963551564734;

pub const WGS84_SEMI_MAJOR_AXIS: f64 = 6_378_137.0;

pub const WGS84_FLATTENING: f64 = 0.0033528106647474805;

#[allow(clippy::excessive_precision)]
#[allow(clippy::approx_constant)]
pub const SQRT2: f64 = 1.4142135623730950488;

/// Astronomical Unit in meters (IAU 2012 definition, exact)
pub const AU_M: f64 = 149_597_870_700.0;

/// Astronomical Unit in kilometers (derived from IAU 2012 definition)
pub const AU_KM: f64 = 149_597_870.7;

pub const DAYS_PER_JULIAN_YEAR: f64 = 365.25;

pub const SPEED_OF_LIGHT_AU_PER_DAY: f64 = 173.1446326846693;

#[allow(clippy::excessive_precision)]
pub const J2000_OBLIQUITY_RAD: f64 = (23.0 + 26.0 / 60.0 + 21.41136 / 3600.0) * PI / 180.0;

#[allow(clippy::excessive_precision)]
pub const FRAME_BIAS_PHI_RAD: f64 = -0.05188 / 3600.0 * PI / 180.0;

#[allow(clippy::excessive_precision)]
pub const GM_EARTH_KM3S2: f64 = 398600.435507;

#[allow(clippy::excessive_precision)]
pub const GM_MOON_KM3S2: f64 = 4902.800118;

pub const MOON_EARTH_MASS_RATIO: f64 = GM_MOON_KM3S2 / (GM_EARTH_KM3S2 + GM_MOON_KM3S2);
