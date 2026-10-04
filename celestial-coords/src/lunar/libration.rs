use super::position::MeanArguments;
use celestial_core::constants::DEG_TO_RAD;
use libm::{cos, sin};

// Meeus' K2 (ch. 53), in degrees and degrees per century. His K1 is A1 of ch. 47.
const K2: [f64; 2] = [72.56, 20.186];

// Meeus ch. 53: the physical libration in inclination (ρ), node (σ) and longitude (τ), radians.
#[derive(Default)]
pub(super) struct PhysicalLibration {
    pub(super) rho: f64,
    pub(super) sigma: f64,
    pub(super) tau: f64,
}

impl PhysicalLibration {
    pub(super) fn at(args: &MeanArguments) -> Self {
        Self {
            rho: rho(args) * DEG_TO_RAD,
            sigma: sigma(args) * DEG_TO_RAD,
            tau: tau(args) * DEG_TO_RAD,
        }
    }
}

fn rho(args: &MeanArguments) -> f64 {
    let (d, mp, f) = (
        args.elongation,
        args.moon_anomaly,
        args.argument_of_latitude,
    );
    -0.02752 * cos(mp) - 0.02245 * sin(f) + 0.00684 * cos(mp - 2.0 * f)
        - 0.00293 * cos(2.0 * f)
        - 0.00085 * cos(2.0 * (f - d))
        - 0.00054 * cos(mp - 2.0 * d)
        - 0.00020 * sin(mp + f)
        - 0.00020 * cos(mp + 2.0 * f)
        - 0.00020 * cos(mp - f)
        + 0.00014 * cos(mp + 2.0 * (f - d))
}

fn sigma(args: &MeanArguments) -> f64 {
    let (d, mp, f) = (
        args.elongation,
        args.moon_anomaly,
        args.argument_of_latitude,
    );
    -0.02816 * sin(mp) + 0.02244 * cos(f)
        - 0.00682 * sin(mp - 2.0 * f)
        - 0.00279 * sin(2.0 * f)
        - 0.00083 * sin(2.0 * (f - d))
        + 0.00069 * sin(mp - 2.0 * d)
        + 0.00040 * cos(mp + f)
        - 0.00025 * sin(2.0 * mp)
        - 0.00023 * sin(mp + 2.0 * f)
        + 0.00020 * cos(mp - f)
        + 0.00019 * sin(mp - f)
        + 0.00013 * sin(mp + 2.0 * (f - d))
        - 0.00010 * cos(mp - 3.0 * f)
}

fn tau(args: &MeanArguments) -> f64 {
    let (d, m, mp, f) = (
        args.elongation,
        args.sun_anomaly,
        args.moon_anomaly,
        args.argument_of_latitude,
    );
    let k2 = (K2[0] + K2[1] * args.t) * DEG_TO_RAD;
    0.02520 * args.eccentricity * sin(m) + 0.00473 * sin(2.0 * (mp - f)) - 0.00467 * sin(mp)
        + 0.00396 * sin(args.a1())
        + 0.00276 * sin(2.0 * (mp - d))
        + 0.00196 * sin(args.node)
        - 0.00183 * cos(mp - f)
        + 0.00115 * sin(mp - 2.0 * d)
        - 0.00096 * sin(mp - d)
        + 0.00046 * sin(2.0 * (f - d))
        - 0.00039 * sin(mp - f)
        - 0.00032 * sin(mp - m - d)
        + 0.00027 * sin(2.0 * (mp - d) - m)
        + 0.00023 * sin(k2)
        - 0.00014 * sin(2.0 * d)
        + 0.00014 * cos(2.0 * (mp - f))
        - 0.00012 * sin(mp - 2.0 * f)
        - 0.00012 * sin(2.0 * mp)
        + 0.00011 * sin(2.0 * (mp - m - d))
}
