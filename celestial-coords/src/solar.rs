use crate::aberration::compute_earth_state;
use crate::errors::CoordResult;
use crate::frames::ecliptic::ecm06_matrix;
use celestial_core::angle::{wrap_0_2pi, Angle};
use celestial_core::constants::{ARCSEC_TO_RAD, DAYS_PER_JULIAN_CENTURY, DEG_TO_RAD, TWOPI};
use celestial_core::matrix::{RotationMatrix3, Vector3};
use celestial_core::obliquity::iau_2006_mean_obliquity;
use celestial_time::scales::tt::TT;
use celestial_time::transforms::nutation::NutationCalculator;

// Carrington's elements as Meeus gives them (Astronomical Algorithms, 2nd ed., ch. 29): the
// solar equator is inclined I to the mean ecliptic of date, its ascending node K advances with
// general precession, and the prime meridian turns with the 25.38-day sidereal period.
const INCLINATION_RAD: f64 = 7.25 * DEG_TO_RAD;
const NODE_EPOCH_JD: f64 = 2396758.0;
const NODE_AT_EPOCH_DEG: f64 = 73.6667;
const NODE_RATE_DEG_PER_CENTURY: f64 = 1.3958333;
const MERIDIAN_EPOCH_JD: f64 = 2398220.0;
const SIDEREAL_PERIOD_DAYS: f64 = 25.38;

// Meeus eq. 29.1 without its periodic terms. It only has to pick the integer rotation number,
// because the fraction comes from L0.
const ROTATION_EPOCH_JD: f64 = 2398140.2270;
const SYNODIC_PERIOD_DAYS: f64 = 27.2752316;

// Meeus eq. 25.10: the Sun's apparent longitude is its geometric longitude less 20.4898"/R.
const ABERRATION_ARCSEC_AU: f64 = 20.4898;

pub struct SolarOrientation {
    pub b0: Angle,
    pub l0: Angle,
    pub p: Angle,
}

pub fn compute_solar_orientation(epoch: &TT) -> CoordResult<SolarOrientation> {
    let earth = apparent_earth_direction(epoch)?;
    let (b0, l0) = disk_center(epoch, earth)?;
    Ok(SolarOrientation {
        b0,
        l0,
        p: position_angle(epoch, earth)?,
    })
}

// Only P needs nutation, so the single-value functions don't go through the full orientation.
pub(crate) fn compute_b0(epoch: &TT) -> CoordResult<Angle> {
    Ok(disk_center(epoch, apparent_earth_direction(epoch)?)?.0)
}

pub(crate) fn compute_l0(epoch: &TT) -> CoordResult<Angle> {
    Ok(disk_center(epoch, apparent_earth_direction(epoch)?)?.1)
}

// A rotation starts when L0 passes 360°, so the fraction is how far L0 has fallen since.
pub fn carrington_rotation_number(epoch: &TT) -> CoordResult<f64> {
    let jd = epoch.to_julian_date();
    let mean = ((jd.jd1() - ROTATION_EPOCH_JD) + jd.jd2()) / SYNODIC_PERIOD_DAYS;
    let fraction = 1.0 - compute_l0(epoch)?.radians() / TWOPI;
    Ok(fraction + libm::round(mean - fraction))
}

pub fn sun_earth_distance(epoch: &TT) -> CoordResult<f64> {
    Ok(compute_earth_state(epoch)?
        .heliocentric_position
        .magnitude())
}

pub(crate) fn icrs_to_carrington(epoch: &TT) -> CoordResult<RotationMatrix3> {
    Ok(ecliptic_to_carrington(epoch).multiply(&ecm06_matrix(epoch)?))
}

fn ecliptic_to_carrington(epoch: &TT) -> RotationMatrix3 {
    let mut m = RotationMatrix3::identity();
    m.rotate_z(node_longitude(epoch));
    m.rotate_x(INCLINATION_RAD);
    m.rotate_z(prime_meridian(epoch));
    m
}

fn node_longitude(epoch: &TT) -> f64 {
    let jd = epoch.to_julian_date();
    let centuries = ((jd.jd1() - NODE_EPOCH_JD) + jd.jd2()) / DAYS_PER_JULIAN_CENTURY;
    (NODE_AT_EPOCH_DEG + NODE_RATE_DEG_PER_CENTURY * centuries) * DEG_TO_RAD
}

fn prime_meridian(epoch: &TT) -> f64 {
    let jd = epoch.to_julian_date();
    let days = (jd.jd1() - MERIDIAN_EPOCH_JD) + jd.jd2();
    libm::fmod(days * 360.0 / SIDEREAL_PERIOD_DAYS, 360.0) * DEG_TO_RAD
}

// Meeus displaces the Sun's longitude by aberration before finding B0 and L0, so the Earth's
// heliocentric direction (mean ecliptic of date) gets the same displacement.
fn apparent_earth_direction(epoch: &TT) -> CoordResult<Vector3> {
    let earth = compute_earth_state(epoch)?.heliocentric_position;
    let mut aberration = RotationMatrix3::identity();
    aberration.rotate_z(ABERRATION_ARCSEC_AU * ARCSEC_TO_RAD / earth.magnitude());
    Ok(aberration.multiply(&ecm06_matrix(epoch)?) * earth)
}

// B0 and L0 are the Earth's heliographic latitude and Carrington longitude.
fn disk_center(epoch: &TT, earth: Vector3) -> CoordResult<(Angle, Angle)> {
    let (l0, b0) = (ecliptic_to_carrington(epoch) * earth).to_spherical();
    Ok((
        Angle::from_radians(b0),
        Angle::from_radians(wrap_0_2pi(l0)?),
    ))
}

// Meeus ch. 29: P = x + y, where x tilts the ecliptic pole against the true celestial pole and
// y tilts the solar pole against the ecliptic pole.
fn position_angle(epoch: &TT, earth: Vector3) -> CoordResult<Angle> {
    let sun_longitude = libm::atan2(-earth.y, -earth.x);
    let jd = epoch.to_julian_date();
    let nutation = epoch.nutation_iau2006a()?;
    let obliquity = iau_2006_mean_obliquity(jd.jd1(), jd.jd2())? + nutation.nutation_obliquity();
    let apparent_longitude = sun_longitude + nutation.nutation_longitude();

    let x = libm::atan(-libm::cos(apparent_longitude) * libm::tan(obliquity));
    let node_distance = sun_longitude - node_longitude(epoch);
    let y = libm::atan(-libm::cos(node_distance) * libm::tan(INCLINATION_RAD));
    Ok(Angle::from_radians(x + y))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rounded;
    use celestial_time::julian::JulianDate;

    fn tt(jd: f64) -> TT {
        TT::from_julian_date(JulianDate::new(jd, 0.0))
    }

    #[test]
    fn test_meeus_example_29a() {
        let o = compute_solar_orientation(&tt(2448908.50068)).unwrap();
        let printed = [o.p, o.b0, o.l0].map(|a| rounded(a.degrees(), 2));
        assert_eq!(printed, [26.27, 5.99, 238.63]);
    }

    #[test]
    fn test_orientation_regression_at_meeus_example_29a() {
        let o = compute_solar_orientation(&tt(2448908.50068)).unwrap();
        assert_eq!(
            [o.p.radians(), o.b0.radians(), o.l0.radians()],
            [0.45855809341835996, 0.10451052877571221, 4.1649065093158235]
        );
    }

    #[test]
    fn test_carrington_rotation_1699_starts_at_meeus_jde() {
        let c = carrington_rotation_number(&tt(2444480.7230)).unwrap();
        assert_eq!(rounded(c, 4), 1699.0);
    }

    #[test]
    fn test_carrington_rotation_minus_10_starts_at_meeus_jde() {
        let c = carrington_rotation_number(&tt(2397867.4913)).unwrap();
        assert_eq!(rounded(c, 4), -10.0);
    }

    #[test]
    fn test_carrington_rotation_regression_at_j2000() {
        assert_eq!(
            carrington_rotation_number(&TT::j2000()).unwrap(),
            1957.9959467485958
        );
    }

    #[test]
    fn test_sun_earth_distance_is_epv00_heliocentric_distance() {
        // eraEpv00 output at its own test epoch.
        let epoch = TT::from_julian_date(JulianDate::new(2400000.5, 53411.52501161));
        let earth = Vector3::new(-0.7757238809297661, 0.5598052241363407, 0.24269984664817157);
        assert_eq!(sun_earth_distance(&epoch).unwrap(), earth.magnitude());
    }

    #[test]
    fn test_non_finite_epoch_is_err() {
        let epoch = tt(f64::NAN);
        assert!(compute_solar_orientation(&epoch).is_err());
        assert!(carrington_rotation_number(&epoch).is_err());
        assert!(sun_earth_distance(&epoch).is_err());
    }
}
