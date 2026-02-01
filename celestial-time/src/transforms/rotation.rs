use crate::{JulianDate, TimeResult};
use celestial_core::angle::wrap_0_2pi;
use celestial_core::constants::{J2000_JD, TWOPI};
use celestial_core::math::fmod;

pub fn earth_rotation_angle(ut1_jd: &JulianDate) -> TimeResult<f64> {
    let ut1_jd1 = ut1_jd.jd1();
    let ut1_jd2 = ut1_jd.jd2();

    let (d1, d2) = if ut1_jd1 < ut1_jd2 {
        (ut1_jd1, ut1_jd2)
    } else {
        (ut1_jd2, ut1_jd1)
    };

    let t = d1 + (d2 - J2000_JD);

    let f = fmod(d1, 1.0) + fmod(d2, 1.0);

    let theta = wrap_0_2pi(TWOPI * (f + 0.7790572732640 + 0.00273781191135448 * t));

    Ok(theta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::UT1;
    use celestial_core::constants::J2000_JD;

    #[test]
    fn test_era_j2000() {
        let ut1 = UT1::j2000();
        let era = earth_rotation_angle(&ut1.to_julian_date()).unwrap();

        assert!((era - 4.894961212823757).abs() < 1e-12);
    }

    #[test]
    fn test_era_precision() {
        let test_cases = [
            (J2000_JD, 0.0),
            (2451545.5, 0.0),
            (2440587.5, 0.0),
            (J2000_JD, 0.5),
        ];

        for (jd1, jd2) in test_cases {
            let ut1 = UT1::from_julian_date(JulianDate::new(jd1, jd2));
            let era = earth_rotation_angle(&ut1.to_julian_date()).unwrap();

            assert!((0.0..2.0 * celestial_core::constants::PI).contains(&era));
        }
    }
}
