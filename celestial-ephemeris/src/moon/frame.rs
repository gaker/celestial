use celestial_core::constants::{ARCSEC_TO_RAD, DAYS_PER_JULIAN_CENTURY};

use super::fit::Poly5;

type Matrix = [[f64; 3]; 3];

// (ε, φ) of the rotation from the ELP frame to ICRS. The LLR fit keeps the
// IAU 2006 obliquity alone; of the published angle sets, that leaves it
// closest to DE432s. The DE405 fit takes the angles fitted with it (Chapront
// et al. 2002), eq. 5 of the VSOP2013 paper.
pub(super) const LLR: (f64, f64) = (84381.406 * ARCSEC_TO_RAD, 0.0);
pub(super) const DE405: (f64, f64) = (84381.4096 * ARCSEC_TO_RAD, -0.05028 * ARCSEC_TO_RAD);

// Laskar's P and Q series for the precession of the ecliptic.
const P: Poly5 = [
    0.10180391e-04,
    0.47020439e-06,
    -0.5417367e-09,
    -0.2507948e-11,
    0.463486e-14,
];
const Q: Poly5 = [
    -0.113469002e-03,
    0.12372674e-06,
    0.1265417e-08,
    -0.1371808e-11,
    -0.320334e-14,
];

pub(super) struct Frame {
    sin_eps: f64,
    cos_eps: f64,
    sin_phi: f64,
    cos_phi: f64,
}

impl Frame {
    pub(super) fn new((eps, phi): (f64, f64)) -> Self {
        let (sin_eps, cos_eps) = libm::sincos(eps);
        let (sin_phi, cos_phi) = libm::sincos(phi);
        Self {
            sin_eps,
            cos_eps,
            sin_phi,
            cos_phi,
        }
    }

    pub(super) fn to_icrs(&self, [x, y, z]: [f64; 3]) -> [f64; 3] {
        let y1 = y * self.cos_eps - z * self.sin_eps;
        let z1 = y * self.sin_eps + z * self.cos_eps;
        [
            x * self.cos_phi - y1 * self.sin_phi,
            x * self.sin_phi + y1 * self.cos_phi,
            z1,
        ]
    }
}

// Longitude, latitude and distance with their rates to rectangular
// coordinates in the mean ecliptic of date.
pub(super) fn rectangular(v: &[f64; 6]) -> ([f64; 3], [f64; 3]) {
    let (slamb, clamb) = libm::sincos(v[0]);
    let (sbeta, cbeta) = libm::sincos(v[1]);
    let cw = v[2] * cbeta;
    let sw = v[2] * sbeta;
    let x = [cw * clamb, cw * slamb, sw];
    let xp = [
        (v[5] * cbeta - v[4] * sw) * clamb - v[3] * x[1],
        (v[5] * cbeta - v[4] * sw) * slamb + v[3] * x[0],
        v[5] * sbeta + v[4] * cw,
    ];
    (x, xp)
}

// From the mean ecliptic of date to the inertial mean ecliptic of J2000.
// Velocities come in per century and leave per day.
pub(super) fn to_j2000(x: &[f64; 3], xp: &[f64; 3], t: &Poly5) -> [f64; 6] {
    let (m, dm) = precession(t);
    let pos = |i: usize| m[i][0] * x[0] + m[i][1] * x[1] + m[i][2] * x[2];
    let vel = |i: usize| {
        (m[i][0] * xp[0]
            + m[i][1] * xp[1]
            + m[i][2] * xp[2]
            + dm[i][0] * x[0]
            + dm[i][1] * x[1]
            + dm[i][2] * x[2])
            / DAYS_PER_JULIAN_CENTURY
    };
    [pos(0), pos(1), pos(2), vel(0), vel(1), vel(2)]
}

fn precession(t: &Poly5) -> (Matrix, Matrix) {
    let (pw, qw) = (pole(&P, t), pole(&Q, t));
    let ra = 2.0 * libm::sqrt(1.0 - pw * pw - qw * qw);
    let (pw2, qw2) = (1.0 - 2.0 * pw * pw, 1.0 - 2.0 * qw * qw);
    let (pwqw, pwra, qwra) = (2.0 * pw * qw, pw * ra, qw * ra);
    let m = [
        [pw2, pwqw, pwra],
        [pwqw, qw2, -qwra],
        [-pwra, qwra, pw2 + qw2 - 1.0],
    ];
    (m, precession_rate(pw, qw, ra, t))
}

fn precession_rate(pw: f64, qw: f64, ra: f64, t: &Poly5) -> Matrix {
    let (ppw, qpw) = (pole_rate(&P, t), pole_rate(&Q, t));
    let (ppw2, qpw2) = (-4.0 * pw * ppw, -4.0 * qw * qpw);
    let ppwqpw = 2.0 * (ppw * qw + pw * qpw);
    let rap = (ppw2 + qpw2) / ra;
    let (ppwra, qpwra) = (ppw * ra + pw * rap, qpw * ra + qw * rap);
    [
        [ppw2, ppwqpw, ppwra],
        [ppwqpw, qpw2, -qpwra],
        [-ppwra, qpwra, ppw2 + qpw2],
    ]
}

fn pole(c: &Poly5, t: &Poly5) -> f64 {
    (c[0] + c[1] * t[1] + c[2] * t[2] + c[3] * t[3] + c[4] * t[4]) * t[1]
}

fn pole_rate(c: &Poly5, t: &Poly5) -> f64 {
    c[0] + (2.0 * c[1] + 3.0 * c[2] * t[1] + 4.0 * c[3] * t[2] + 5.0 * c[4] * t[3]) * t[1]
}
