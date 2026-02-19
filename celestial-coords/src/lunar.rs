use crate::{CoordResult, ICRSPosition};
use celestial_core::constants::{DEG_TO_RAD, J2000_JD};
use celestial_core::utils::{normalize_angle_rad, normalize_angle_to_positive};
use celestial_core::Angle;
use celestial_time::TT;

const LUNAR_AXIAL_INCLINATION_DEG: f64 = 1.5424;
const LUNAR_AXIAL_INCLINATION_RAD: f64 = LUNAR_AXIAL_INCLINATION_DEG * DEG_TO_RAD;

pub struct LunarLibration {
    pub longitude: Angle,
    pub latitude: Angle,
}

pub struct LunarOrientation {
    pub optical_libration: LunarLibration,
    pub sub_earth_point: LunarLibration,
    pub position_angle: Angle,
}

pub fn compute_lunar_orientation(epoch: &TT) -> LunarOrientation {
    let (lib_lon, lib_lat) = compute_optical_libration_internal(epoch);
    let position_angle = compute_position_angle_internal(epoch);

    LunarOrientation {
        optical_libration: LunarLibration {
            longitude: Angle::from_radians(lib_lon),
            latitude: Angle::from_radians(lib_lat),
        },
        sub_earth_point: LunarLibration {
            longitude: Angle::from_radians(lib_lon),
            latitude: Angle::from_radians(lib_lat),
        },
        position_angle: Angle::from_radians(position_angle),
    }
}

pub fn compute_optical_libration(epoch: &TT) -> (Angle, Angle) {
    let (lon, lat) = compute_optical_libration_internal(epoch);
    (Angle::from_radians(lon), Angle::from_radians(lat))
}

pub fn compute_sub_earth_point(epoch: &TT) -> (Angle, Angle) {
    compute_optical_libration(epoch)
}

fn compute_optical_libration_internal(epoch: &TT) -> (f64, f64) {
    let jd = epoch.to_julian_date();
    let d = (jd.jd1() - J2000_JD) + jd.jd2();
    let t = d / celestial_core::constants::DAYS_PER_JULIAN_CENTURY;

    let _mean_anomaly = moon_mean_anomaly(t);
    let mean_argument_latitude = moon_argument_latitude(t);
    let mean_elongation = moon_mean_elongation(t);
    let ascending_node = moon_ascending_node(t);

    let lib_lon = -0.02752 * libm::cos(ascending_node)
        - 0.02245 * libm::sin(mean_argument_latitude)
        + 0.00684 * libm::cos(mean_argument_latitude - 2.0 * mean_elongation)
        - 0.00293 * libm::cos(2.0 * mean_argument_latitude)
        - 0.00085 * libm::cos(2.0 * mean_argument_latitude - 2.0 * mean_elongation)
        - 0.00054 * libm::sin(mean_argument_latitude - 2.0 * mean_elongation)
        - 0.00020 * libm::sin(mean_argument_latitude + ascending_node)
        - 0.00020 * libm::cos(2.0 * mean_argument_latitude - mean_elongation)
        - 0.00020 * libm::sin(mean_argument_latitude - ascending_node);

    let argument = mean_argument_latitude - ascending_node;

    let lib_lat = -0.02816 * libm::sin(argument) + 0.02244 * libm::cos(ascending_node)
        - 0.00682 * libm::sin(argument - 2.0 * mean_elongation)
        - 0.00279 * libm::sin(2.0 * mean_argument_latitude - argument)
        - 0.00083 * libm::sin(2.0 * mean_argument_latitude - argument - 2.0 * mean_elongation)
        + 0.00069 * libm::sin(argument + 2.0 * mean_elongation)
        + 0.00040 * libm::cos(2.0 * ascending_node);

    let lib_lon_rad = lib_lon * 10.0 * DEG_TO_RAD;
    let lib_lat_rad = lib_lat * 10.0 * DEG_TO_RAD;

    (
        normalize_angle_rad(lib_lon_rad),
        normalize_angle_rad(lib_lat_rad),
    )
}

fn compute_position_angle_internal(epoch: &TT) -> f64 {
    let jd = epoch.to_julian_date();
    let d = (jd.jd1() - J2000_JD) + jd.jd2();
    let t = d / celestial_core::constants::DAYS_PER_JULIAN_CENTURY;

    let ascending_node = moon_ascending_node(t);
    let obliquity = mean_obliquity(t);
    let i_prime = LUNAR_AXIAL_INCLINATION_RAD;

    let v = ascending_node;
    let (vs, vc) = libm::sincos(v);

    let (is, ic) = libm::sincos(i_prime);
    let (os, oc) = libm::sincos(obliquity);

    let x = is * vs;
    let y = is * vc * oc - ic * os;

    libm::atan2(y, x)
}

fn moon_mean_anomaly(t: f64) -> f64 {
    let m_prime = 134.9633964 + 477198.8675055 * t + 0.0087414 * t * t + t * t * t / 69699.0
        - t * t * t * t / 14712000.0;
    normalize_angle_to_positive(m_prime * DEG_TO_RAD)
}

fn moon_argument_latitude(t: f64) -> f64 {
    let f = 93.272095 + 483202.0175233 * t - 0.0036539 * t * t - t * t * t / 3526000.0
        + t * t * t * t / 863310000.0;
    normalize_angle_to_positive(f * DEG_TO_RAD)
}

fn moon_mean_elongation(t: f64) -> f64 {
    let d = 297.8501921 + 445267.1114034 * t - 0.0018819 * t * t + t * t * t / 545868.0
        - t * t * t * t / 113065000.0;
    normalize_angle_to_positive(d * DEG_TO_RAD)
}

fn moon_ascending_node(t: f64) -> f64 {
    let omega = 125.0445479 - 1934.1362891 * t + 0.0020754 * t * t + t * t * t / 467441.0
        - t * t * t * t / 60616000.0;
    normalize_angle_to_positive(omega * DEG_TO_RAD)
}

fn mean_obliquity(t: f64) -> f64 {
    let eps0 = 23.439291 - 0.0130042 * t - 1.64e-7 * t * t + 5.04e-7 * t * t * t;
    eps0 * DEG_TO_RAD
}

pub(crate) fn get_moon_icrs(epoch: &TT) -> CoordResult<ICRSPosition> {
    let jd = epoch.to_julian_date();
    let d = (jd.jd1() - J2000_JD) + jd.jd2();
    let t = d / celestial_core::constants::DAYS_PER_JULIAN_CENTURY;

    let l_prime = 218.3164477 + 481267.88123421 * t;
    let l_prime = normalize_angle_to_positive(l_prime * DEG_TO_RAD);

    let d_moon = 297.8501921 + 445267.1114034 * t;
    let d_moon = normalize_angle_to_positive(d_moon * DEG_TO_RAD);

    let m = 357.5291092 + 35999.0502909 * t;
    let m = normalize_angle_to_positive(m * DEG_TO_RAD);

    let m_prime = 134.9633964 + 477198.8675055 * t;
    let m_prime = normalize_angle_to_positive(m_prime * DEG_TO_RAD);

    let f = 93.272095 + 483202.0175233 * t;
    let f = normalize_angle_to_positive(f * DEG_TO_RAD);

    let moon_lon = l_prime
        + 6.289 * DEG_TO_RAD * libm::sin(m_prime)
        + 1.274 * DEG_TO_RAD * libm::sin(2.0 * d_moon - m_prime)
        + 0.658 * DEG_TO_RAD * libm::sin(2.0 * d_moon)
        + 0.214 * DEG_TO_RAD * libm::sin(2.0 * m_prime)
        - 0.186 * DEG_TO_RAD * libm::sin(m)
        - 0.114 * DEG_TO_RAD * libm::sin(2.0 * f);

    let moon_lat = 5.128 * DEG_TO_RAD * libm::sin(f)
        + 0.281 * DEG_TO_RAD * libm::sin(m_prime + f)
        + 0.278 * DEG_TO_RAD * libm::sin(m_prime - f);

    let eps = (23.439291 - 0.0130042 * t) * DEG_TO_RAD;

    let (sin_lon, cos_lon) = libm::sincos(moon_lon);
    let (sin_lat, cos_lat) = libm::sincos(moon_lat);
    let (sin_eps, cos_eps) = libm::sincos(eps);

    let part = sin_lon * cos_eps - libm::tan(moon_lat) * sin_eps;
    let ra = libm::atan2(part, cos_lon);
    let dec = libm::asin(sin_lat * cos_eps + cos_lat * sin_eps * sin_lon);

    ICRSPosition::new(
        Angle::from_radians(normalize_angle_to_positive(ra)),
        Angle::from_radians(dec),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_time::julian::JulianDate;

    #[test]
    fn test_optical_libration_longitude_range() {
        let epochs = [
            TT::j2000(),
            TT::from_julian_date(JulianDate::new(J2000_JD + 7.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 14.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 21.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 28.0, 0.0)),
        ];

        for epoch in &epochs {
            let (lon, _) = compute_optical_libration(epoch);
            assert!(
                lon.degrees().abs() <= 8.5,
                "Libration longitude = {} degrees exceeds expected range ±7.9°",
                lon.degrees()
            );
        }
    }

    #[test]
    fn test_optical_libration_latitude_range() {
        let epochs = [
            TT::j2000(),
            TT::from_julian_date(JulianDate::new(J2000_JD + 7.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 14.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 21.0, 0.0)),
            TT::from_julian_date(JulianDate::new(J2000_JD + 28.0, 0.0)),
        ];

        for epoch in &epochs {
            let (_, lat) = compute_optical_libration(epoch);
            assert!(
                lat.degrees().abs() <= 7.5,
                "Libration latitude = {} degrees exceeds expected range ±6.7°",
                lat.degrees()
            );
        }
    }

    #[test]
    fn test_sub_earth_point_equals_optical_libration() {
        let epoch = TT::j2000();
        let (lib_lon, lib_lat) = compute_optical_libration(&epoch);
        let (sub_lon, sub_lat) = compute_sub_earth_point(&epoch);

        assert_eq!(lib_lon.radians(), sub_lon.radians());
        assert_eq!(lib_lat.radians(), sub_lat.radians());
    }

    #[test]
    fn test_lunar_orientation_combined() {
        let epoch = TT::j2000();
        let orientation = compute_lunar_orientation(&epoch);

        assert!(
            orientation.optical_libration.longitude.degrees().abs() <= 8.5,
            "lib_lon = {}",
            orientation.optical_libration.longitude.degrees()
        );
        assert!(
            orientation.optical_libration.latitude.degrees().abs() <= 7.5,
            "lib_lat = {}",
            orientation.optical_libration.latitude.degrees()
        );
        assert!(
            orientation.sub_earth_point.longitude.degrees().abs() <= 8.5,
            "sub_lon = {}",
            orientation.sub_earth_point.longitude.degrees()
        );
        assert!(
            orientation.sub_earth_point.latitude.degrees().abs() <= 7.5,
            "sub_lat = {}",
            orientation.sub_earth_point.latitude.degrees()
        );
    }

    #[test]
    fn test_libration_variation_over_month() {
        let start = TT::j2000();
        let end = TT::from_julian_date(JulianDate::new(J2000_JD + 27.3, 0.0));

        let (lon_start, lat_start) = compute_optical_libration(&start);
        let (lon_end, lat_end) = compute_optical_libration(&end);

        let lon_diff = (lon_end.degrees() - lon_start.degrees()).abs();
        let lat_diff = (lat_end.degrees() - lat_start.degrees()).abs();

        assert!(
            lon_diff < 16.0,
            "Longitude libration change over one month should be < 16°, got {}",
            lon_diff
        );
        assert!(
            lat_diff < 14.0,
            "Latitude libration change over one month should be < 14°, got {}",
            lat_diff
        );
    }
}
