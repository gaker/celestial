//! ELP/MPP02 Lunar Ephemeris
//!
//! Computes geocentric rectangular coordinates of the Moon using the
//! ELP/MPP02 semi-analytical lunar theory by Chapront & Francou (2003).
//!
//! Output is in ICRS, with positions in AU and velocities in AU/day.

mod fit;
mod frame;
mod series;
#[cfg(test)]
mod tests;

use celestial_core::{
    constants::{ARCSEC_PER_RAD, AU_KM, DAYS_PER_JULIAN_CENTURY},
    errors::AstroResult,
    matrix::Vector3,
};
use celestial_time::scales::tdb::TDB;

use crate::validity::ELPMPP02;
use fit::{Fit, Poly5};

const A405: f64 = 384747.9613701725;
const AELP: f64 = 384747.980674318;

pub struct ElpMpp02Moon {
    fit: &'static Fit,
}

impl Default for ElpMpp02Moon {
    fn default() -> Self {
        Self::new()
    }
}

impl ElpMpp02Moon {
    pub fn new() -> Self {
        Self { fit: Fit::llr() }
    }

    pub fn with_de405_fit() -> Self {
        Self { fit: Fit::de405() }
    }

    pub fn geocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
        Ok(self.icrs(self.ecliptic_position(tdb)?))
    }

    pub fn geocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
        let [x, y, z, vx, vy, vz] = self.ecliptic_state(tdb)?;
        Ok((self.icrs([x, y, z]), self.icrs([vx, vy, vz])))
    }

    // From km in the theory's frame to AU in ICRS.
    fn icrs(&self, v: [f64; 3]) -> Vector3 {
        Vector3::from_array(self.fit.frame.to_icrs(v)) / AU_KM
    }

    // In km, in the mean ecliptic and equinox of J2000.
    fn ecliptic_position(&self, tdb: &TDB) -> AstroResult<[f64; 3]> {
        let [x, y, z, ..] = self.evaluate::<false>(ELPMPP02.days_from_j2000(tdb)?);
        Ok([x, y, z])
    }

    // In km and km/day, in the mean ecliptic and equinox of J2000.
    fn ecliptic_state(&self, tdb: &TDB) -> AstroResult<[f64; 6]> {
        Ok(self.evaluate::<true>(ELPMPP02.days_from_j2000(tdb)?))
    }

    // Without RATES the series rates are skipped, and only the position is
    // meaningful.
    fn evaluate<const RATES: bool>(&self, tj: f64) -> [f64; 6] {
        let t1 = tj / DAYS_PER_JULIAN_CENTURY;
        let t2 = t1 * t1;
        let t3 = t2 * t1;
        let t = [1.0, t1, t2, t3, t3 * t1];
        let (x, xp) = frame::rectangular(&self.spherical::<RATES>(&t));
        frame::to_j2000(&x, &xp, &t)
    }

    // Longitude and latitude in radians, distance in km, and their rates per
    // century, in the mean ecliptic of date.
    fn spherical<const RATES: bool>(&self, t: &Poly5) -> [f64; 6] {
        let w1 = &self.fit.w1;
        let mut v = self.fit.series.sums::<RATES>(t);
        v[0] = v[0] / ARCSEC_PER_RAD
            + w1[0]
            + w1[1] * t[1]
            + w1[2] * t[2]
            + w1[3] * t[3]
            + w1[4] * t[4];
        v[1] /= ARCSEC_PER_RAD;
        v[2] = v[2] * A405 / AELP;
        v[3] = v[3] / ARCSEC_PER_RAD
            + w1[1]
            + 2.0 * w1[2] * t[1]
            + 3.0 * w1[3] * t[2]
            + 4.0 * w1[4] * t[3];
        v[4] /= ARCSEC_PER_RAD;
        v
    }
}
