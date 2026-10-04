use super::tables::{LATITUDE, LONGITUDE_AND_DISTANCE};
use crate::errors::{CoordError, CoordResult};
use celestial_core::constants::{AU_M, DEG_TO_RAD};
use celestial_core::matrix::{RotationMatrix3, Vector3};
use celestial_core::precession::PrecessionIAU2006;
use celestial_core::utils::jd_to_centuries;
use celestial_time::scales::tt::TT;

// Meeus, Astronomical Algorithms (2nd ed.), eq. 47.1-47.5 and 47.7: degrees, then degrees per power of T.
// L′ starts from Simon et al. (1994) rather than Meeus' 218.3164477, which folds in the Moon's
// 0.7″ of light-time, so the position is geometric.
const MEAN_LONGITUDE: [f64; 5] = [
    218.31665436,
    481267.88123421,
    -0.0015786,
    1.0 / 538841.0,
    -1.0 / 65194000.0,
];
// The T⁴ term of D and the T³ term of F have the opposite sign to Meeus' printed equations, as
// in SOFA's moon98, so that this position agrees with it exactly. The difference stays below
// 1e-8° within five centuries of J2000.
const MEAN_ELONGATION: [f64; 5] = [
    297.8501921,
    445267.1114034,
    -0.0018819,
    1.0 / 545868.0,
    1.0 / 113065000.0,
];
const SUN_MEAN_ANOMALY: [f64; 5] = [
    357.5291092,
    35999.0502909,
    -0.0001536,
    1.0 / 24490000.0,
    0.0,
];
const MOON_MEAN_ANOMALY: [f64; 5] = [
    134.9633964,
    477198.8675055,
    0.0087414,
    1.0 / 69699.0,
    -1.0 / 14712000.0,
];
const ARGUMENT_OF_LATITUDE: [f64; 5] = [
    93.2720950,
    483202.0175233,
    -0.0036539,
    1.0 / 3526000.0,
    1.0 / 863310000.0,
];
const NODE: [f64; 5] = [
    125.0445479,
    -1934.1362891,
    0.0020754,
    1.0 / 467441.0,
    -1.0 / 60616000.0,
];

// Meeus' further arguments A1 (Venus), A2 (Jupiter) and A3, and eq. 47.6 for the eccentricity
// factor E applied to terms in M.
const A1: [f64; 2] = [119.75, 131.849];
const A2: [f64; 2] = [53.09, 479264.29];
const A3: [f64; 2] = [313.45, 481266.484];
const ECCENTRICITY: [f64; 2] = [-0.002516, -0.0000074];

// Meeus' Δ = 385000.56 km + Σr, with Σr in metres.
const MEAN_DISTANCE_M: f64 = 385000560.0;

pub(super) struct MeanArguments {
    pub(super) t: f64,
    pub(super) mean_longitude: f64,
    pub(super) elongation: f64,
    pub(super) sun_anomaly: f64,
    pub(super) moon_anomaly: f64,
    pub(super) argument_of_latitude: f64,
    pub(super) node: f64,
    pub(super) eccentricity: f64,
}

impl MeanArguments {
    pub(super) fn at(epoch: &TT) -> CoordResult<Self> {
        let jd = epoch.to_julian_date();
        let t = jd_to_centuries(jd.jd1(), jd.jd2());
        if !t.is_finite() {
            return Err(CoordError::invalid_coordinate("Epoch is not finite"));
        }
        Ok(Self {
            t,
            mean_longitude: polynomial(MEAN_LONGITUDE, t),
            elongation: polynomial(MEAN_ELONGATION, t),
            sun_anomaly: polynomial(SUN_MEAN_ANOMALY, t),
            moon_anomaly: polynomial(MOON_MEAN_ANOMALY, t),
            argument_of_latitude: polynomial(ARGUMENT_OF_LATITUDE, t),
            node: polynomial(NODE, t),
            eccentricity: 1.0 + (ECCENTRICITY[0] + ECCENTRICITY[1] * t) * t,
        })
    }

    pub(super) fn a1(&self) -> f64 {
        linear(A1, self.t)
    }

    fn term(&self, [d, m, m_prime, f]: [i8; 4]) -> (f64, f64) {
        let argument = f64::from(d) * self.elongation
            + f64::from(m) * self.sun_anomaly
            + f64::from(m_prime) * self.moon_anomaly
            + f64::from(f) * self.argument_of_latitude;
        let factor = match m.unsigned_abs() {
            1 => self.eccentricity,
            2 => self.eccentricity * self.eccentricity,
            _ => 1.0,
        };
        (argument, factor)
    }
}

pub(super) fn moon_gcrs(args: &MeanArguments) -> Vector3 {
    let (longitude, distance) = ecliptic_longitude_and_distance(args);
    let latitude = ecliptic_latitude(args);
    let rcp = distance * libm::cos(latitude);
    let ecliptic = Vector3::new(
        rcp * libm::cos(longitude),
        rcp * libm::sin(longitude),
        distance * libm::sin(latitude),
    );
    ecliptic_to_gcrs(args.t) * ecliptic
}

// The sums run from the last table row up, as moon98's do, so the rounding matches it.
fn ecliptic_longitude_and_distance(args: &MeanArguments) -> (f64, f64) {
    let mut longitude = 0.003958 * libm::sin(args.a1())
        + 0.001962 * libm::sin(args.mean_longitude - args.argument_of_latitude)
        + 0.000318 * libm::sin(linear(A2, args.t));
    let mut distance = 0.0;
    for &(multiples, l, r) in LONGITUDE_AND_DISTANCE.iter().rev() {
        let (argument, factor) = args.term(multiples);
        longitude += l * (libm::sin(argument) * factor);
        distance += r * (libm::cos(argument) * factor);
    }
    (
        args.mean_longitude + DEG_TO_RAD * longitude,
        (distance + MEAN_DISTANCE_M) / AU_M,
    )
}

fn ecliptic_latitude(args: &MeanArguments) -> f64 {
    let (l, m_prime, f, a1) = (
        args.mean_longitude,
        args.moon_anomaly,
        args.argument_of_latitude,
        args.a1(),
    );
    let mut latitude = -0.002235 * libm::sin(l)
        + 0.000382 * libm::sin(linear(A3, args.t))
        + 0.000175 * libm::sin(a1 - f)
        + 0.000175 * libm::sin(a1 + f)
        + 0.000127 * libm::sin(l - m_prime)
        - 0.000115 * libm::sin(l + m_prime);
    for &(multiples, b) in LATITUDE.iter().rev() {
        let (argument, factor) = args.term(multiples);
        latitude += b * (libm::sin(argument) * factor);
    }
    latitude * DEG_TO_RAD
}

// Built from the Fukushima-Williams angles, as moon98 builds it, rather than as the transpose of
// the ICRS-to-ecliptic matrix: the two agree only to rounding.
fn ecliptic_to_gcrs(t: f64) -> RotationMatrix3 {
    let fw = PrecessionIAU2006::new().fukushima_williams_angles(t);
    let mut m = RotationMatrix3::identity();
    m.rotate_z(fw.psi_bar);
    m.rotate_x(-fw.phi_bar);
    m.rotate_z(-fw.gamma_bar);
    m
}

fn polynomial(c: [f64; 5], t: f64) -> f64 {
    let degrees = c[0] + (c[1] + (c[2] + (c[3] + c[4] * t) * t) * t) * t;
    libm::fmod(degrees, 360.0) * DEG_TO_RAD
}

fn linear(c: [f64; 2], t: f64) -> f64 {
    (c[0] + c[1] * t) * DEG_TO_RAD
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::J2000_JD;
    use celestial_time::julian::JulianDate;

    fn moon(jd1: f64, jd2: f64) -> [f64; 3] {
        let epoch = TT::from_julian_date(JulianDate::new(jd1, jd2));
        let v = moon_gcrs(&MeanArguments::at(&epoch).unwrap());
        [v.x, v.y, v.z]
    }

    // eraMoon98 outputs, au.
    #[test]
    fn test_matches_moon98() {
        assert_eq!(
            moon(2400000.5, 43999.9),
            [
                -0.0026012959599712084,
                0.0006139750944294706,
                0.0002640794528226106
            ]
        );
        assert_eq!(
            moon(2448724.5, 0.0),
            [
                -0.0016853456834658809,
                0.0016977178353056598,
                0.0005848854873104547
            ]
        );
        assert_eq!(
            moon(J2000_JD, 0.0),
            [
                -0.0019492621453406477,
                -0.001782881213717378,
                -0.0005086906382512421
            ]
        );
        assert_eq!(
            moon(2400000.5, 60000.0),
            [
                0.0020030934215161016,
                0.0014478563698763357,
                0.0006396537649144986
            ]
        );
    }

    #[test]
    fn test_non_finite_epoch_is_err() {
        for jd in [f64::NAN, f64::INFINITY] {
            let epoch = TT::from_julian_date(JulianDate::new(jd, 0.0));
            assert!(MeanArguments::at(&epoch).is_err());
        }
    }
}
