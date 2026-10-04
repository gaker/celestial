use crate::errors::{CoordError, CoordResult};
use celestial_core::constants::SPEED_OF_LIGHT_AU_PER_DAY as C;
use celestial_core::matrix::Vector3;

// Space motion faster than this fraction of c is rejected rather than used.
const MAX_SPEED: f64 = 0.5;
const MAX_ITERATIONS: usize = 100;

// Spherical coordinates and their rates: radians and au, per day.
pub(super) struct SphericalMotion {
    pub(super) ra: f64,
    pub(super) dec: f64,
    pub(super) distance: f64,
    pub(super) ra_rate: f64,
    pub(super) dec_rate: f64,
    pub(super) radial_rate: f64,
}

// Barycentric position (au) and inertial velocity (au/day).
pub(super) struct SpaceMotion {
    position: Vector3,
    velocity: Vector3,
}

impl SpaceMotion {
    // A catalog's proper motion and radial velocity are observed rates; the straight-line
    // motion needs the inertial ones, which differ by the special-relativistic terms below.
    pub(super) fn from_observed(observed: &SphericalMotion) -> CoordResult<Self> {
        let observed = Self::from_spherical(observed);
        if observed.velocity.magnitude() / C > MAX_SPEED {
            return Err(CoordError::invalid_coordinate(
                "Space velocity exceeds half the speed of light",
            ));
        }
        let (unit, radial, transverse) = observed.split_velocity()?;
        let (beta_radial, beta_transverse) = (radial / C, transverse.magnitude() / C);
        let (d, del) = observed_to_inertial(beta_radial, beta_transverse);
        Ok(Self {
            position: observed.position,
            velocity: unit * (C * (d * beta_radial + del)) + transverse * d,
        })
    }

    pub(super) fn to_observed(&self) -> CoordResult<SphericalMotion> {
        let (unit, radial, transverse) = self.split_velocity()?;
        let (beta_radial, beta_transverse) = (radial / C, transverse.magnitude() / C);
        let d = 1.0 + beta_radial;
        let w = beta_radial * beta_radial + beta_transverse * beta_transverse;
        if d == 0.0 || w > 1.0 {
            return Err(CoordError::invalid_coordinate(
                "Space velocity reaches the speed of light",
            ));
        }
        let del = -w / (libm::sqrt(1.0 - w) + 1.0);
        let velocity = unit * (C * (beta_radial - del) / d) + transverse * (1.0 / d);
        Ok(spherical(self.position, velocity))
    }

    // The catalog place is where the star was when the light seen at the first epoch left it.
    // Move it to the second epoch's geometric place, then back by that place's light time.
    pub(super) fn propagate(&self, days: f64) -> CoordResult<Self> {
        let light_time = self.position.magnitude() / C;
        let geometric = self.after(days + light_time);
        Ok(self.after(days + (light_time - geometric.arrival_light_time()?)))
    }

    fn after(&self, days: f64) -> Self {
        Self {
            position: self.position + self.velocity * days,
            velocity: self.velocity,
        }
    }

    fn arrival_light_time(&self) -> CoordResult<f64> {
        let r2 = self.position.dot(&self.position);
        let rdv = self.position.dot(&self.velocity);
        let c2mv2 = C * C - self.velocity.dot(&self.velocity);
        if c2mv2 <= 0.0 {
            return Err(CoordError::invalid_coordinate(
                "Space velocity reaches the speed of light",
            ));
        }
        Ok((-rdv + libm::sqrt(rdv * rdv + c2mv2 * r2)) / c2mv2)
    }

    // The unit position vector, the radial speed along it, and the transverse velocity.
    fn split_velocity(&self) -> CoordResult<(Vector3, f64, Vector3)> {
        let unit = self.position.normalize()?;
        let radial = unit.dot(&self.velocity);
        Ok((unit, radial, self.velocity - unit * radial))
    }

    fn from_spherical(s: &SphericalMotion) -> Self {
        let (st, ct) = (libm::sin(s.ra), libm::cos(s.ra));
        let (sp, cp) = (libm::sin(s.dec), libm::cos(s.dec));
        let rcp = s.distance * cp;
        let (x, y) = (rcp * ct, rcp * st);
        let rpd = s.distance * s.dec_rate;
        let w = rpd * sp - cp * s.radial_rate;
        Self {
            position: Vector3::new(x, y, s.distance * sp),
            velocity: Vector3::new(
                -y * s.ra_rate - w * ct,
                x * s.ra_rate - w * st,
                rpd * cp + sp * s.radial_rate,
            ),
        }
    }
}

// The position is never zero here: split_velocity has already normalized it.
fn spherical(p: Vector3, v: Vector3) -> SphericalMotion {
    let rxy2 = p.x * p.x + p.y * p.y;
    let r2 = rxy2 + p.z * p.z;
    let (distance, rxy) = (libm::sqrt(r2), libm::sqrt(rxy2));
    let xyp = p.x * v.x + p.y * v.y;
    let on_axis = rxy2 == 0.0;
    SphericalMotion {
        ra: if on_axis { 0.0 } else { libm::atan2(p.y, p.x) },
        dec: libm::atan2(p.z, rxy),
        distance,
        ra_rate: if on_axis {
            0.0
        } else {
            (p.x * v.y - p.y * v.x) / rxy2
        },
        dec_rate: if on_axis {
            0.0
        } else {
            (v.z * rxy2 - p.z * xyp) / (r2 * rxy)
        },
        radial_rate: (xyp + p.z * v.z) / distance,
    }
}

// Solves for the inertial radial and transverse speeds (as fractions of c) that the observed
// ones imply, iterating until the corrections stop shrinking, as SOFA's starpv does. Below
// half of c the iteration always converges, but the stopping test can cycle in the last bits
// without passing (about 1 star in 75,000 in a random sweep), so the final iterate is kept,
// as SOFA keeps it, rather than refused.
fn observed_to_inertial(beta_radial: f64, beta_transverse: f64) -> (f64, f64) {
    let (mut betr, mut bett, mut d, mut del) = (beta_radial, beta_transverse, 0.0, 0.0);
    let (mut od, mut odel, mut odd, mut oddel) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..MAX_ITERATIONS {
        d = 1.0 + betr;
        let w = betr * betr + bett * bett;
        del = -w / (libm::sqrt(1.0 - w) + 1.0);
        betr = d * beta_radial + del;
        bett = d * beta_transverse;
        if i > 0 {
            let (dd, ddel) = (libm::fabs(d - od), libm::fabs(del - odel));
            if i > 1 && dd >= odd && ddel >= oddel {
                break;
            }
            (odd, oddel) = (dd, ddel);
        }
        (od, odel) = (d, del);
    }
    (d, del)
}
