use crate::julian::JulianDate;
use crate::{TimeError, TimeResult};
use celestial_core::angle::wrap_0_2pi;
use celestial_core::constants::{J2000_JD, TWOPI};
use celestial_core::math::fmod;

pub fn earth_rotation_angle(ut1_jd: &JulianDate) -> TimeResult<f64> {
    let (d1, d2) = if ut1_jd.jd1() < ut1_jd.jd2() {
        (ut1_jd.jd1(), ut1_jd.jd2())
    } else {
        (ut1_jd.jd2(), ut1_jd.jd1())
    };
    let t = days_from_j2000(d1, d2)?;
    let f = fmod(d1, 1.0) + fmod(d2, 1.0);
    Ok(wrap_0_2pi(
        TWOPI * (f + 0.7790572732640 + 0.00273781191135448 * t),
    )?)
}

fn days_from_j2000(d1: f64, d2: f64) -> TimeResult<f64> {
    let t = d1 + (d2 - J2000_JD);
    if !t.is_finite() || libm::fabs(t) > 1e12 {
        return Err(TimeError::CalculationError(format!(
            "Time value out of valid range: {} days from J2000",
            t
        )));
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_era_matches_erfa() {
        // eraEra00
        let cases = [
            (J2000_JD, 0.0, 4.894961212823756),
            (2451545.5, 0.0, 1.7619696490215855),
            (2440587.5, 0.0, 1.7560450788883983),
            (J2000_JD, 0.5, 1.7619696490215855),
        ];
        for (jd1, jd2, expected) in cases {
            let era = earth_rotation_angle(&JulianDate::new(jd1, jd2)).unwrap();
            assert_eq!(era, expected, "({jd1}, {jd2})");
        }
    }
}
