mod libration;
mod position;
mod tables;

use crate::errors::CoordResult;
use crate::frames::ecliptic::ecm06_matrix;
use celestial_core::angle::{wrap_0_2pi, Angle};
use celestial_core::constants::{DEG_TO_RAD, PI};
use celestial_core::matrix::{RotationMatrix3, Vector3};
use celestial_core::precession::PrecessionIAU2006;
use celestial_core::utils::jd_to_centuries;
use celestial_time::scales::tt::TT;
use celestial_time::transforms::nutation::NutationCalculator;
use libration::PhysicalLibration;
use position::{moon_gcrs, MeanArguments};

// Meeus ch. 53: the IAU inclination of the mean lunar equator to the ecliptic.
const INCLINATION_RAD: f64 = 1.54242 * DEG_TO_RAD;

pub struct LunarLibration {
    pub longitude: Angle,
    pub latitude: Angle,
}

impl LunarLibration {
    fn of(direction: Vector3) -> Self {
        let (longitude, latitude) = direction.to_spherical();
        Self {
            longitude: Angle::from_radians(longitude),
            latitude: Angle::from_radians(latitude),
        }
    }
}

pub struct LunarOrientation {
    pub optical_libration: LunarLibration,
    pub sub_earth_point: LunarLibration,
    pub position_angle: Angle,
}

pub fn compute_lunar_orientation(epoch: &TT) -> CoordResult<LunarOrientation> {
    let args = MeanArguments::at(epoch)?;
    let moon = moon_gcrs(&args);
    let to_ecliptic = ecm06_matrix(epoch)?;
    let frame = |libration: &PhysicalLibration| {
        selenographic_frame(&args, libration).multiply(&to_ecliptic)
    };
    let selenographic = frame(&PhysicalLibration::at(&args));

    Ok(LunarOrientation {
        optical_libration: LunarLibration::of(frame(&PhysicalLibration::default()) * -moon),
        sub_earth_point: LunarLibration::of(selenographic * -moon),
        position_angle: Angle::from_radians(position_angle(epoch, &selenographic, moon)?),
    })
}

// Only the position angle needs nutation, so the single-value functions don't go through the
// full orientation.
pub fn compute_optical_libration(epoch: &TT) -> CoordResult<(Angle, Angle)> {
    earth_coordinates(epoch, |_| PhysicalLibration::default())
}

pub fn compute_sub_earth_point(epoch: &TT) -> CoordResult<(Angle, Angle)> {
    earth_coordinates(epoch, PhysicalLibration::at)
}

pub(crate) fn icrs_to_selenographic(epoch: &TT) -> CoordResult<RotationMatrix3> {
    let args = MeanArguments::at(epoch)?;
    let frame = selenographic_frame(&args, &PhysicalLibration::at(&args));
    Ok(frame.multiply(&ecm06_matrix(epoch)?))
}

fn earth_coordinates(
    epoch: &TT,
    libration: fn(&MeanArguments) -> PhysicalLibration,
) -> CoordResult<(Angle, Angle)> {
    let args = MeanArguments::at(epoch)?;
    let frame = selenographic_frame(&args, &libration(&args)).multiply(&ecm06_matrix(epoch)?);
    let earth = LunarLibration::of(frame * -moon_gcrs(&args));
    Ok((earth.longitude, earth.latitude))
}

// Meeus ch. 53 as a rotation from the mean ecliptic of date. The lunar equator crosses that
// ecliptic at Ω + σ/sin I, inclined at I + ρ, and the prime meridian (towards the mean Earth)
// lies F + 180° + τ − σ cot I along the equator from there. With ρ = σ = τ = 0 the Earth's
// coordinates in this frame are Meeus' optical libration (eq. 53.1); with them, they reproduce
// eq. 53.2 to first order.
fn selenographic_frame(args: &MeanArguments, libration: &PhysicalLibration) -> RotationMatrix3 {
    let node = args.node + libration.sigma / libm::sin(INCLINATION_RAD);
    let meridian = args.argument_of_latitude + PI + libration.tau
        - libration.sigma / libm::tan(INCLINATION_RAD);
    let mut m = RotationMatrix3::identity();
    m.rotate_z(node);
    m.rotate_x(-(INCLINATION_RAD + libration.rho));
    m.rotate_z(meridian);
    m
}

// Meeus measures P from the north point of the disk, towards the true celestial pole of date.
fn position_angle(epoch: &TT, selenographic: &RotationMatrix3, moon: Vector3) -> CoordResult<f64> {
    let true_equator = gcrs_to_true_equator(epoch)?;
    let pole = true_equator * (selenographic.transpose() * Vector3::new(0.0, 0.0, 1.0));
    let u = (true_equator * moon).normalize()?;
    let east = Vector3::new(-u.y, u.x, 0.0);
    let north = Vector3::new(-u.z * u.x, -u.z * u.y, 1.0 - u.z * u.z);
    Ok(wrap_0_2pi(libm::atan2(pole.dot(&east), pole.dot(&north)))?)
}

fn gcrs_to_true_equator(epoch: &TT) -> CoordResult<RotationMatrix3> {
    let jd = epoch.to_julian_date();
    let nutation = epoch.nutation_iau2006a()?;
    Ok(PrecessionIAU2006::new().npb_matrix_iau2006a(
        jd_to_centuries(jd.jd1(), jd.jd2()),
        nutation.nutation_longitude(),
        nutation.nutation_obliquity(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rounded;
    use celestial_time::julian::JulianDate;

    fn tt(jd: f64) -> TT {
        TT::from_julian_date(JulianDate::new(jd, 0.0))
    }

    fn degrees(libration: &LunarLibration) -> [f64; 2] {
        [libration.longitude.degrees(), libration.latitude.degrees()]
    }

    #[test]
    fn test_meeus_example_53a_optical_libration() {
        let (l, b) = compute_optical_libration(&tt(2448724.5)).unwrap();
        assert_eq!(
            [rounded(l.degrees(), 3), rounded(b.degrees(), 3)],
            [-1.206, 4.194]
        );
    }

    #[test]
    fn test_meeus_example_53a_physical_libration() {
        let o = compute_lunar_orientation(&tt(2448724.5)).unwrap();
        let [l, b] = degrees(&o.sub_earth_point);
        let [l_optical, b_optical] = degrees(&o.optical_libration);
        assert_eq!(
            [rounded(l - l_optical, 3), rounded(b - b_optical, 3)],
            [-0.025, 0.006]
        );
    }

    #[test]
    fn test_meeus_example_53a_sub_earth_point_is_total_libration() {
        let (l, b) = compute_sub_earth_point(&tt(2448724.5)).unwrap();
        assert_eq!(
            [rounded(l.degrees(), 2), rounded(b.degrees(), 2)],
            [-1.23, 4.20]
        );
    }

    #[test]
    fn test_meeus_example_53a_position_angle() {
        let o = compute_lunar_orientation(&tt(2448724.5)).unwrap();
        assert_eq!(rounded(o.position_angle.degrees(), 2), 15.08);
    }

    #[test]
    fn test_orientation_regression_at_meeus_example_53a() {
        let o = compute_lunar_orientation(&tt(2448724.5)).unwrap();
        let radians = |l: &LunarLibration| [l.longitude.radians(), l.latitude.radians()];
        assert_eq!(
            radians(&o.optical_libration),
            [-0.02104138290164475, 0.07319973069234371]
        );
        assert_eq!(
            radians(&o.sub_earth_point),
            [-0.02148494418846901, 0.07329815269437237]
        );
        assert_eq!(o.position_angle.radians(), 0.2632686895950156);
    }

    #[test]
    fn test_non_finite_epoch_is_err() {
        let epoch = tt(f64::NAN);
        assert!(compute_lunar_orientation(&epoch).is_err());
        assert!(compute_optical_libration(&epoch).is_err());
        assert!(compute_sub_earth_point(&epoch).is_err());
        assert!(icrs_to_selenographic(&epoch).is_err());
    }
}
