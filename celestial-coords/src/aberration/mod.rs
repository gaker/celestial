mod coefficients;
#[cfg(test)]
mod tests;

use crate::constants::AU_LIGHT_TIME_DAYS;
use crate::errors::{CoordError, CoordResult};
use celestial_core::{
    angle::{wrap_0_2pi, Angle},
    constants::{DAYS_PER_JULIAN_YEAR, J2000_JD},
    matrix::Vector3,
};
use celestial_time::scales::tt::TT;
use coefficients::Coefficients;

pub struct EarthState {
    pub barycentric_velocity: Vector3,
    pub heliocentric_position: Vector3,
}

// Each axis is a cosine series, plus t and t² times two more, with t in Julian years from
// J2000. One set gives the Earth's heliocentric motion, the other the Sun's barycentric motion.
type Series = [&'static [[f64; 3]]; 3];

const EARTH_SERIES: [Series; 3] = [
    [Coefficients::E0X, Coefficients::E1X, Coefficients::E2X],
    [Coefficients::E0Y, Coefficients::E1Y, Coefficients::E2Y],
    [Coefficients::E0Z, Coefficients::E1Z, Coefficients::E2Z],
];

const SUN_SERIES: [Series; 3] = [
    [Coefficients::S0X, Coefficients::S1X, Coefficients::S2X],
    [Coefficients::S0Y, Coefficients::S1Y, Coefficients::S2Y],
    [Coefficients::S0Z, Coefficients::S1Z, Coefficients::S2Z],
];

pub fn compute_earth_state(tt: &TT) -> CoordResult<EarthState> {
    let t = julian_years_since_j2000(tt)?;
    let mut heliocentric_position = [0.0; 3];
    let mut barycentric_velocity = [0.0; 3];
    for axis in 0..3 {
        let mut sums = (0.0, 0.0);
        accumulate_series(t, &EARTH_SERIES[axis], &mut sums);
        heliocentric_position[axis] = sums.0;
        // The Earth's barycentric motion is its heliocentric motion plus the Sun's.
        accumulate_series(t, &SUN_SERIES[axis], &mut sums);
        barycentric_velocity[axis] = sums.1 / DAYS_PER_JULIAN_YEAR;
    }
    Ok(EarthState {
        barycentric_velocity: rotate_to_bcrs(barycentric_velocity),
        heliocentric_position: rotate_to_bcrs(heliocentric_position),
    })
}

fn julian_years_since_j2000(tt: &TT) -> CoordResult<f64> {
    let jd = tt.to_julian_date();
    let t = ((jd.jd1() - J2000_JD) + jd.jd2()) / DAYS_PER_JULIAN_YEAR;
    if t.is_finite() {
        Ok(t)
    } else {
        Err(CoordError::invalid_coordinate("Epoch is not finite"))
    }
}

// Adds one axis's series to a running (value, rate per year).
fn accumulate_series(t: f64, series: &Series, sums: &mut (f64, f64)) {
    for (power, terms) in series.iter().enumerate() {
        for &[a, b, c] in *terms {
            let ct = c * t;
            let (sin_p, cos_p) = libm::sincos(b + ct);
            let (value, rate) = match power {
                0 => (a * cos_p, -(a * c * sin_p)),
                1 => (a * t * cos_p, a * (cos_p - ct * sin_p)),
                _ => (a * (t * t) * cos_p, a * t * (2.0 * cos_p - ct * sin_p)),
            };
            sums.0 += value;
            sums.1 += rate;
        }
    }
}

// From the ecliptic and equinox of J2000 to the BCRS.
fn rotate_to_bcrs(v: [f64; 3]) -> Vector3 {
    const AM12: f64 = 0.000000211284;
    const AM13: f64 = -0.000000091603;
    const AM21: f64 = -0.000000230286;
    const AM22: f64 = 0.917482137087;
    const AM23: f64 = -0.397776982902;
    const AM32: f64 = 0.397776982902;
    const AM33: f64 = 0.917482137087;

    let (x, y, z) = (v[0], v[1], v[2]);
    Vector3::new(
        x + AM12 * y + AM13 * z,
        AM21 * x + AM22 * y + AM23 * z,
        AM32 * y + AM33 * z,
    )
}

const SCHWARZSCHILD_RADIUS_SUN_AU: f64 = 1.97412574336e-8;

/// Apply gravitational light deflection by the Sun.
///
/// This implements the relativistic bending of starlight as it passes near the Sun,
/// for a star far enough away that its directions from the Sun and from the observer coincide.
///
/// # Arguments
/// * `star_direction` - Unit vector from observer to star (BCRS)
/// * `sun_to_observer` - Unit vector from Sun to observer (BCRS)
/// * `sun_observer_distance_au` - Distance from Sun to observer in AU
///
/// # Returns
/// Deflected star direction (unit vector)
pub fn apply_light_deflection(
    star_direction: Vector3,
    sun_to_observer: Vector3,
    sun_observer_distance_au: f64,
) -> Vector3 {
    // Within about 5' of the Sun's centre as seen from 1 au, well inside the limb, the
    // deflection is held back so that it falls to zero rather than diverging. The threshold
    // shrinks with the Sun's apparent size for observers farther out.
    let distance_squared = sun_observer_distance_au * sun_observer_distance_au;
    let dlim = 1e-6 / libm::fmax(distance_squared, 1.0);
    let (p, e) = (star_direction, sun_to_observer);
    let w =
        SCHWARZSCHILD_RADIUS_SUN_AU / sun_observer_distance_au / libm::fmax(p.dot(&(p + e)), dlim);
    p + p.cross(&e.cross(&p)) * w
}

/// Remove gravitational light deflection by the Sun (inverse operation).
pub fn remove_light_deflection(
    deflected_direction: Vector3,
    sun_to_observer: Vector3,
    sun_observer_distance_au: f64,
) -> CoordResult<Vector3> {
    invert(deflected_direction, 5, |p| {
        apply_light_deflection(p, sun_to_observer, sun_observer_distance_au)
    })
}

/// Remove stellar aberration (inverse operation).
pub fn remove_aberration(
    apparent_direction: Vector3,
    velocity_au_day: Vector3,
    sun_earth_distance_au: f64,
) -> CoordResult<Vector3> {
    let (v, bm1) = observer_motion(velocity_au_day);
    unaberrate(apparent_direction, v, sun_earth_distance_au, bm1)
}

pub fn apply_aberration(
    direction: Vector3,
    velocity_au_day: Vector3,
    sun_earth_distance_au: f64,
) -> Vector3 {
    let (v, bm1) = observer_motion(velocity_au_day);
    aberrate(direction, v, sun_earth_distance_au, bm1)
}

// Observer velocity in units of c, and sqrt(1 - v²).
pub(crate) fn observer_motion(velocity_au_day: Vector3) -> (Vector3, f64) {
    let v = velocity_au_day * AU_LIGHT_TIME_DAYS;
    (v, libm::sqrt(1.0 - v.magnitude_squared()))
}

pub(crate) fn aberrate(direction: Vector3, v: Vector3, sun_distance_au: f64, bm1: f64) -> Vector3 {
    let pdv = direction.dot(&v);
    let w1 = 1.0 + pdv / (1.0 + bm1);
    let w2 = SCHWARZSCHILD_RADIUS_SUN_AU / sun_distance_au;

    let p2 = direction * bm1 + v * w1 + (v - direction * pdv) * w2;
    p2 / p2.magnitude()
}

pub(crate) fn unaberrate(
    apparent: Vector3,
    v: Vector3,
    sun_distance_au: f64,
    bm1: f64,
) -> CoordResult<Vector3> {
    invert(apparent, 2, |p| aberrate(p, v, sun_distance_au, bm1))
}

// Fixed-point inversion of a small direction shift: subtract the shift that
// `forward` produces at the current estimate. Estimates are normalised by
// dividing by the modulus, which is what makes the result bit-exact with ERFA.
fn invert(
    observed: Vector3,
    iterations: usize,
    forward: impl Fn(Vector3) -> Vector3,
) -> CoordResult<Vector3> {
    let mut shift = Vector3::zeros();
    for _ in 0..iterations {
        let before = unit(observed - shift)?;
        shift = forward(before) - before;
    }
    unit(observed - shift)
}

fn unit(v: Vector3) -> CoordResult<Vector3> {
    let r = v.magnitude();
    if r == 0.0 || !r.is_finite() {
        return Err(CoordError::invalid_coordinate(
            "Direction has zero or non-finite length",
        ));
    }
    Ok(v / r)
}

pub fn apply_annual_aberration(ra: Angle, dec: Angle, tt: &TT) -> CoordResult<(Angle, Angle)> {
    annual_shift(ra, dec, tt, true)
}

pub fn remove_annual_aberration(ra: Angle, dec: Angle, tt: &TT) -> CoordResult<(Angle, Angle)> {
    annual_shift(ra, dec, tt, false)
}

fn annual_shift(ra: Angle, dec: Angle, tt: &TT, apply: bool) -> CoordResult<(Angle, Angle)> {
    let state = compute_earth_state(tt)?;
    let sun_distance_au = state.heliocentric_position.magnitude();
    let direction = Vector3::from_spherical(ra.radians(), dec.radians());
    let corrected = if apply {
        apply_aberration(direction, state.barycentric_velocity, sun_distance_au)
    } else {
        remove_aberration(direction, state.barycentric_velocity, sun_distance_au)?
    };
    let (ra, dec) = corrected.to_spherical();
    Ok((
        Angle::from_radians(wrap_0_2pi(ra)?),
        Angle::from_radians(dec),
    ))
}
