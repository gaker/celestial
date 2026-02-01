use crate::{CoordResult, ICRSPosition};
use astro_core::constants::{ARCSEC_TO_RAD, DEG_TO_RAD, J2000_JD, TWOPI};
use astro_core::utils::{normalize_angle_rad, normalize_angle_to_positive};
use astro_core::Angle;
use astro_time::TT;

const SOLAR_EQUATOR_INCLINATION_DEG: f64 = 7.25;
const SOLAR_EQUATOR_INCLINATION_RAD: f64 = SOLAR_EQUATOR_INCLINATION_DEG * DEG_TO_RAD;

const SOLAR_ASCENDING_NODE_J2000_DEG: f64 = 75.76;

pub const CARRINGTON_EPOCH_JD: f64 = 2398220.0;
pub const CARRINGTON_SYNODIC_PERIOD: f64 = 27.2753;

pub struct SolarOrientation {
    pub b0: Angle,
    pub l0: Angle,
    pub p: Angle,
}

pub fn compute_solar_orientation(epoch: &TT) -> SolarOrientation {
    let jd = epoch.to_julian_date();
    let d = (jd.jd1() - J2000_JD) + jd.jd2();
    let t = d / astro_core::constants::DAYS_PER_JULIAN_CENTURY;

    let (sun_lon, sun_lat, obliquity) = solar_ecliptic_coords(t);
    let (b0, l0, p) = heliographic_coords(t, sun_lon, sun_lat, obliquity);

    SolarOrientation {
        b0: Angle::from_radians(b0),
        l0: Angle::from_radians(l0),
        p: Angle::from_radians(p),
    }
}

pub fn compute_b0(epoch: &TT) -> Angle {
    compute_solar_orientation(epoch).b0
}

pub fn compute_l0(epoch: &TT) -> Angle {
    compute_solar_orientation(epoch).l0
}

pub fn compute_p(epoch: &TT) -> Angle {
    compute_solar_orientation(epoch).p
}

pub fn carrington_rotation_number(epoch: &TT) -> u32 {
    let jd = epoch.to_julian_date();
    let jd_days = (jd.jd1() - CARRINGTON_EPOCH_JD) + jd.jd2();
    (jd_days / CARRINGTON_SYNODIC_PERIOD).floor() as u32 + 1
}

pub fn sun_earth_distance(epoch: &TT) -> f64 {
    let jd = epoch.to_julian_date();
    let d = (jd.jd1() - J2000_JD) + jd.jd2();
    let t = d / astro_core::constants::DAYS_PER_JULIAN_CENTURY;

    let m = (357.52911 + 35999.05029 * t - 0.0001537 * t * t) * DEG_TO_RAD;
    let e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t;

    let c_rad = (1.914602 - 0.004817 * t - 0.000014 * t * t) * DEG_TO_RAD * m.sin()
        + (0.019993 - 0.000101 * t) * DEG_TO_RAD * (2.0 * m).sin()
        + 0.000289 * DEG_TO_RAD * (3.0 * m).sin();

    let true_anomaly = m + c_rad;
    let a = 1.000001018; // semi-major axis in AU

    a * (1.0 - e * e) / (1.0 + e * true_anomaly.cos())
}

fn heliographic_coords(t: f64, sun_lon: f64, _sun_lat: f64, obliquity: f64) -> (f64, f64, f64) {
    let i = SOLAR_EQUATOR_INCLINATION_RAD;
    let k = (SOLAR_ASCENDING_NODE_J2000_DEG + 1.3958333 * t) * DEG_TO_RAD;

    let lambda = sun_lon;
    let theta = lambda - k;
    let (sin_theta, cos_theta) = theta.sin_cos();
    let (sin_i, cos_i) = i.sin_cos();
    let (_sin_obl, cos_obl) = obliquity.sin_cos();

    let b0 = (sin_theta * sin_i).asin();

    let eta = (sin_i * cos_theta).atan2(cos_i);
    let jd_days =
        t * astro_core::constants::DAYS_PER_JULIAN_CENTURY + J2000_JD - CARRINGTON_EPOCH_JD;
    let l0_raw = 360.0 / CARRINGTON_SYNODIC_PERIOD * jd_days;
    let l0 = normalize_angle_to_positive((l0_raw * DEG_TO_RAD - eta) % TWOPI);

    let rho = (cos_theta * sin_i / cos_obl).atan();
    let sigma = (sin_theta * cos_i).atan();
    let p = normalize_angle_rad(rho + sigma);

    (b0, l0, p)
}

fn solar_ecliptic_coords(t: f64) -> (f64, f64, f64) {
    let l0 = 280.46646 + 36000.76983 * t + 0.0003032 * t * t;
    let m = 357.52911 + 35999.05029 * t - 0.0001537 * t * t;
    let m_rad = m * DEG_TO_RAD;

    let c = (1.914602 - 0.004817 * t - 0.000014 * t * t) * m_rad.sin()
        + (0.019993 - 0.000101 * t) * (2.0 * m_rad).sin()
        + 0.000289 * (3.0 * m_rad).sin();

    let sun_true_lon = l0 + c;

    let omega = 125.04 - 1934.136 * t;
    let omega_rad = omega * DEG_TO_RAD;
    let apparent_lon = sun_true_lon - 0.00569 - 0.00478 * omega_rad.sin();

    let obliquity = mean_obliquity(t);

    (
        normalize_angle_to_positive(apparent_lon * DEG_TO_RAD),
        0.0,
        obliquity,
    )
}

fn mean_obliquity(t: f64) -> f64 {
    let eps0_arcsec = 84381.448 - 46.8150 * t - 0.00059 * t * t + 0.001813 * t * t * t;
    eps0_arcsec * ARCSEC_TO_RAD
}

pub(crate) fn get_sun_icrs(epoch: &TT) -> CoordResult<ICRSPosition> {
    let jd = epoch.to_julian_date();
    let d = (jd.jd1() - J2000_JD) + jd.jd2();
    let t = d / astro_core::constants::DAYS_PER_JULIAN_CENTURY;

    let l0 = 280.46646 + 36000.76983 * t + 0.0003032 * t * t;
    let m = 357.52911 + 35999.05029 * t - 0.0001537 * t * t;
    let m_rad = m * DEG_TO_RAD;

    let c = (1.914602 - 0.004817 * t - 0.000014 * t * t) * m_rad.sin()
        + (0.019993 - 0.000101 * t) * (2.0 * m_rad).sin()
        + 0.000289 * (3.0 * m_rad).sin();

    let sun_true_lon = l0 + c;
    let omega = 125.04 - 1934.136 * t;
    let omega_rad = omega * DEG_TO_RAD;
    let apparent_lon = sun_true_lon - 0.00569 - 0.00478 * omega_rad.sin();

    let lambda = apparent_lon * DEG_TO_RAD;
    let eps = (23.439291 - 0.0130042 * t) * DEG_TO_RAD;

    let (sin_lambda, cos_lambda) = lambda.sin_cos();
    let (sin_eps, cos_eps) = eps.sin_cos();

    let ra = (sin_lambda * cos_eps).atan2(cos_lambda);
    let dec = (sin_lambda * sin_eps).asin();

    ICRSPosition::new(
        Angle::from_radians(normalize_angle_to_positive(ra)),
        Angle::from_radians(dec),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_time::julian::JulianDate;

    #[test]
    fn test_b0_range() {
        let epochs = [
            TT::j2000(),
            TT::from_julian_date(JulianDate::new(J2000_JD + 91.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 182.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 273.0, 0.0)),
        ];

        for epoch in &epochs {
            let b0 = compute_b0(epoch);
            assert!(
                b0.degrees().abs() <= 7.3,
                "B0 = {} degrees exceeds expected range ±7.25°",
                b0.degrees()
            );
        }
    }

    #[test]
    fn test_l0_range() {
        let epoch = TT::j2000();
        let l0 = compute_l0(&epoch);
        assert!(
            l0.degrees() >= 0.0 && l0.degrees() < 360.0,
            "L0 = {} degrees outside [0, 360) range",
            l0.degrees()
        );
    }

    #[test]
    fn test_p_range() {
        let epochs = [
            TT::j2000(),
            TT::from_julian_date(JulianDate::new(J2000_JD + 91.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 182.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 273.0, 0.0)),
        ];

        for epoch in &epochs {
            let p = compute_p(epoch);
            assert!(
                p.degrees().abs() <= 45.0,
                "P = {} degrees exceeds expected range ±45°",
                p.degrees()
            );
        }
    }

    #[test]
    fn test_carrington_rotation_period() {
        let epoch1 = TT::j2000();
        let l0_1 = compute_l0(&epoch1);

        let epoch2 =
            TT::from_julian_date(JulianDate::new(J2000_JD + CARRINGTON_SYNODIC_PERIOD, 0.0));
        let l0_2 = compute_l0(&epoch2);

        let diff = (l0_2.degrees() - l0_1.degrees()).abs();
        assert!(
            (diff - 360.0).abs() < 5.0 || diff < 5.0,
            "L0 should change by ~360° in one Carrington rotation, got {} degrees",
            diff
        );
    }

    #[test]
    fn test_solar_orientation_combined() {
        let epoch = TT::j2000();
        let orientation = compute_solar_orientation(&epoch);

        assert!(
            orientation.b0.degrees().abs() <= 7.3,
            "B0 = {} out of range",
            orientation.b0.degrees()
        );
        assert!(
            orientation.l0.degrees() >= 0.0 && orientation.l0.degrees() < 360.0,
            "L0 = {} out of range",
            orientation.l0.degrees()
        );
        assert!(
            orientation.p.degrees().abs() <= 30.0,
            "P = {} out of range",
            orientation.p.degrees()
        );
    }
}
