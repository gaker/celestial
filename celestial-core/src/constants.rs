pub const J2000_JD: f64 = 2451545.0;

pub const DAYS_PER_JULIAN_CENTURY: f64 = 36525.0;

// Validity span of the IAU 2000/2006 series; the CIO locator and nutation reject epochs beyond it.
pub const MAX_CENTURIES_FROM_J2000: f64 = 20.0;

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

pub const ARCMIN_TO_RAD: f64 = 0.0002908882086657216;

pub const ARCSEC_TO_RAD: f64 = 4.84813681109536e-6;

pub const MILLIARCSEC_TO_RAD: f64 = ARCSEC_TO_RAD / 1e3;

pub const TENTH_MICROARCSEC_TO_RAD: f64 = ARCSEC_TO_RAD / 1e7;

pub const MJD_ZERO_POINT: f64 = 2_400_000.5;

// The workspace takes this constant from here rather than from std, so approx_constant
// is a false positive at this one definition.
#[allow(clippy::approx_constant)]
pub const PI: f64 = 3.141592653589793;

pub const HALF_PI: f64 = PI / 2.0;

pub const QUARTER_PI: f64 = PI / 4.0;

pub const TWOPI: f64 = PI * 2.0;

pub const DEG_TO_RAD: f64 = 0.017453292519943295;

pub const RAD_TO_DEG: f64 = 57.29577951308232;

pub const ARCSEC_PER_RAD: f64 = 206264.80624709636;

// Unit conversion factors as head + tail: the head is the correctly rounded factor and the
// tail the remainder rounded to a double, so Angle's conversions multiply with one rounding.
// Three tails are the remainder's other neighbouring double, because the nearest one
// misrounds every input with one mantissa: HOURS_TO_RAD_LO (1.06888805882525, e.g.
// 17.102208941204 h), RAD_TO_DEG_LO (1.8557083013070463 rad) and DEG_TO_RAD_LO
// (1.5888359943899804, e.g. 50.84275182047937°).
pub(crate) const RAD_TO_DEG_LO: f64 = -1.9878495670576287e-15;
pub(crate) const DEG_TO_RAD_LO: f64 = 2.948652270870168e-19;
pub(crate) const ARCSEC_PER_RAD_LO: f64 = -8.06575314318039e-12;
pub(crate) const ARCSEC_TO_RAD_LO: f64 = 9.320078015422868e-23;
pub(crate) const RAD_TO_HOURS: f64 = 3.819718634205488;
pub(crate) const RAD_TO_HOURS_LO: f64 = -1.4099515177158535e-17;
pub(crate) const HOURS_TO_RAD: f64 = 0.26179938779914946;
pub(crate) const HOURS_TO_RAD_LO: f64 = -2.6802044161277278e-17;
pub(crate) const RAD_TO_ARCMIN: f64 = 3437.746770784939;
pub(crate) const RAD_TO_ARCMIN_LO: f64 = 1.0810270141977435e-13;
pub(crate) const ARCMIN_TO_RAD_LO: f64 = 1.5756442176305324e-20;

pub const WGS84_SEMI_MAJOR_AXIS: f64 = 6_378_137.0;

pub const WGS84_FLATTENING: f64 = 0.0033528106647474805;

// The workspace takes this constant from here rather than from std, so approx_constant
// is a false positive at this one definition.
#[allow(clippy::approx_constant)]
pub const SQRT2: f64 = 1.4142135623730951;

/// Astronomical Unit in meters (IAU 2012 definition, exact)
pub const AU_M: f64 = 149_597_870_700.0;

/// Astronomical Unit in kilometers (derived from IAU 2012 definition)
pub const AU_KM: f64 = 149_597_870.7;

pub const DAYS_PER_JULIAN_YEAR: f64 = 365.25;

pub const SPEED_OF_LIGHT_M_PER_S: f64 = 299_792_458.0;

pub const SPEED_OF_LIGHT_AU_PER_DAY: f64 = SECONDS_PER_DAY_F64 / (AU_M / SPEED_OF_LIGHT_M_PER_S);

pub const VSOP2013_OBLIQUITY_RAD: f64 = 84381.41136 * ARCSEC_TO_RAD;

pub const VSOP2013_PHI_RAD: f64 = -0.05188 * ARCSEC_TO_RAD;

pub const GM_EARTH_KM3S2: f64 = 398600.435507;

pub const GM_MOON_KM3S2: f64 = 4902.800118;

pub const MOON_EMB_MASS_RATIO: f64 = GM_MOON_KM3S2 / (GM_EARTH_KM3S2 + GM_MOON_KM3S2);

#[cfg(test)]
mod tests {
    use super::*;

    // Speed of light in au/day from the IAU 2012 au and SI c; equals ERFA's DC.
    #[test]
    fn test_speed_of_light_au_per_day_matches_erfa_dc() {
        assert_eq!(SPEED_OF_LIGHT_AU_PER_DAY, 173.1446326742403);
    }

    // Nutation amplitudes are tabulated in units of 0.1 µas; equals ERFA's U2R factor.
    #[test]
    fn test_tenth_microarcsec_to_rad_matches_erfa_u2r() {
        assert_eq!(TENTH_MICROARCSEC_TO_RAD, 4.848136811095359e-13);
    }

    // Expected values are the IAU/IERS constants as ERFA defines them (erfam.h, and eform
    // for WGS84), compared bit for bit.
    #[test]
    fn test_time_constants_match_erfa() {
        assert_eq!(J2000_JD, 2451545.0);
        assert_eq!(MJD_ZERO_POINT, 2400000.5);
        assert_eq!(DAYS_PER_JULIAN_YEAR, 365.25);
        assert_eq!(DAYS_PER_JULIAN_CENTURY, 36525.0);
        assert_eq!(DAYS_PER_JULIAN_MILLENNIUM, 365250.0);
        assert_eq!(SECONDS_PER_DAY_F64, 86400.0);
    }

    #[test]
    fn test_angle_constants_match_erfa() {
        assert_eq!(PI.to_bits(), 0x400921fb54442d18);
        assert_eq!(TWOPI.to_bits(), 0x401921fb54442d18);
        assert_eq!(DEG_TO_RAD, 0.017453292519943295);
        assert_eq!(RAD_TO_DEG, 57.29577951308232);
        assert_eq!(ARCSEC_TO_RAD, 4.84813681109536e-6);
        assert_eq!(MILLIARCSEC_TO_RAD, ARCSEC_TO_RAD / 1e3);
        assert_eq!(ARCSEC_PER_RAD, 206264.80624709636);
        assert_eq!(CIRCULAR_ARCSECONDS, 1296000.0);
    }

    // Correctly rounded π/2, π/4, √2 and π/10800; the first three are compared as bits
    // because clippy's approx_constant rejects them as literals.
    #[test]
    fn test_other_angle_constants_are_correctly_rounded() {
        assert_eq!(HALF_PI.to_bits(), 0x3ff921fb54442d18);
        assert_eq!(QUARTER_PI.to_bits(), 0x3fe921fb54442d18);
        assert_eq!(SQRT2.to_bits(), 0x3ff6a09e667f3bcd);
        assert_eq!(ARCMIN_TO_RAD, 0.0002908882086657216);
    }

    // Correctly rounded 84381.41136″ and -0.05188″.
    #[test]
    fn test_vsop2013_angles_are_correctly_rounded() {
        assert_eq!(VSOP2013_OBLIQUITY_RAD, 0.4090926265865962);
        assert_eq!(VSOP2013_PHI_RAD, -2.515213377596273e-7);
    }

    #[test]
    fn test_physical_constants_match_erfa() {
        assert_eq!(AU_M, 149597870700.0);
        assert_eq!(SPEED_OF_LIGHT_M_PER_S, 299792458.0);
        assert_eq!(WGS84_SEMI_MAJOR_AXIS, 6378137.0);
        assert_eq!(WGS84_FLATTENING, 1.0 / 298.257223563);
        assert_eq!(
            WGS84_ECCENTRICITY_SQUARED,
            (2.0 - WGS84_FLATTENING) * WGS84_FLATTENING
        );
    }
}
