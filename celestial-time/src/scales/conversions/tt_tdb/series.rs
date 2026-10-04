use crate::constants::FAIRHD;
use celestial_core::constants::{DAYS_PER_JULIAN_MILLENNIUM, DEG_TO_RAD, J2000_JD, TWOPI};
use celestial_core::math::fmod;
use std::ops::Range;

/// Compute TDB-TT difference using Fairhead & Bretagnon (1990) model.
///
/// This implements the full 787-term series for the TDB-TT difference, including:
/// - Fundamental arguments (Sun, Moon, planets mean longitudes)
/// - Diurnal terms from observer's geocentric position
/// - JPL planetary-mass adjustment terms
///
/// # Arguments
///
/// - `date1`, `date2`: Two-part Julian Date
/// - `ut`: UT1 fraction of day from midnight (for local solar time calculation)
/// - `elong`: Observer's east longitude in radians
/// - `u`, `v`: Geocentric cylindrical coordinates in km (from Location::to_geocentric_km)
///
/// # Returns
///
/// TDB - TT in seconds. Range is approximately -0.00166 to +0.00166 seconds.
pub(super) fn calculate_tdb_tt_difference(
    date1: f64,
    date2: f64,
    ut: f64,
    elong: f64,
    u: f64,
    v: f64,
) -> f64 {
    let t = ((date1 - J2000_JD) + date2) / DAYS_PER_JULIAN_MILLENNIUM;
    let tsol = fmod(ut, 1.0) * TWOPI + elong;
    topocentric_terms(t, tsol, u, v) + fairhead_series(t) + jpl_mass_terms(t)
}

// `tsol` is local solar time as an angle. `u` and `v` are the observer's distances
// from Earth's spin axis and north of the equatorial plane, in km.
fn topocentric_terms(t: f64, tsol: f64, u: f64, v: f64) -> f64 {
    let w = t / 3600.0;
    let elsun = fmod(280.46645683 + 1296027711.03429 * w, 360.0) * DEG_TO_RAD;
    let emsun = fmod(357.52910918 + 1295965810.481 * w, 360.0) * DEG_TO_RAD;
    let d = fmod(297.85019547 + 16029616012.090 * w, 360.0) * DEG_TO_RAD;
    let elj = fmod(34.35151874 + 109306899.89453 * w, 360.0) * DEG_TO_RAD;
    let els = fmod(50.07744430 + 44046398.47038 * w, 360.0) * DEG_TO_RAD;

    0.00029e-10 * u * libm::sin(tsol + elsun - els)
        + 0.00100e-10 * u * libm::sin(tsol - 2.0 * emsun)
        + 0.00133e-10 * u * libm::sin(tsol - d)
        + 0.00133e-10 * u * libm::sin(tsol + elsun - elj)
        - 0.00229e-10 * u * libm::sin(tsol + 2.0 * elsun + emsun)
        - 0.02200e-10 * v * libm::cos(elsun + emsun)
        + 0.05312e-10 * u * libm::sin(tsol - emsun)
        - 0.13677e-10 * u * libm::sin(tsol + 2.0 * elsun)
        - 1.31840e-10 * v * libm::cos(elsun)
        + 3.17679e-10 * u * libm::sin(tsol)
}

// Each block of FAIRHD is the coefficient of one power of t, t^0 through t^4.
fn fairhead_series(t: f64) -> f64 {
    let [w0, w1, w2, w3, w4] =
        [0..474, 474..679, 679..764, 764..784, 784..787].map(|terms| sine_sum(t, terms));
    t * (t * (t * (t * w4 + w3) + w2) + w1) + w0
}

// The blocks run roughly largest term first, so summing in reverse adds the
// small terms before the large ones.
fn sine_sum(t: f64, terms: Range<usize>) -> f64 {
    let mut sum = 0.0;
    for [amplitude, frequency, phase] in FAIRHD[terms].iter().rev() {
        sum += amplitude * libm::sin(frequency * t + phase);
    }
    sum
}

fn jpl_mass_terms(t: f64) -> f64 {
    0.00065e-6 * libm::sin(6069.776754 * t + 4.021194)
        + 0.00033e-6 * libm::sin(213.299095 * t + 5.543132)
        - 0.00196e-6 * libm::sin(6208.294251 * t + 5.696701)
        - 0.00173e-6 * libm::sin(74.781599 * t + 2.435900)
        + 0.03638e-6 * t * t
}
