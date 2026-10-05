mod orbit;
mod series;
#[cfg(test)]
mod tests;

use celestial_core::constants::DAYS_PER_JULIAN_MILLENNIUM;
use celestial_core::ecliptic::vsop2013_to_icrs;
use celestial_core::errors::AstroResult;
use celestial_core::matrix::Vector3;
use celestial_time::scales::tdb::TDB;

use crate::planetary_coefficients::{
    emb, jupiter, mars, mercury, neptune, pluto, saturn, uranus, venus,
};
use crate::validity::{Window, VSOP2013, VSOP2013_PLUTO};
use orbit::Orbit;
use series::Series;

struct Theory {
    series: Series,
    window: &'static Window,
}

impl Theory {
    fn position(&self, tdb: &TDB) -> AstroResult<Vector3> {
        let (elements, _) = series::elements(&self.series, self.millennia(tdb)?);
        Ok(to_icrs(Orbit::new(elements)?.position()))
    }

    fn state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
        let (elements, rates) = series::elements(&self.series, self.millennia(tdb)?);
        let orbit = Orbit::new(elements)?;
        let velocity = to_icrs(orbit.velocity(&rates)) / DAYS_PER_JULIAN_MILLENNIUM;
        Ok((to_icrs(orbit.position()), velocity))
    }

    fn millennia(&self, tdb: &TDB) -> AstroResult<f64> {
        Ok(self.window.days_from_j2000(tdb)? / DAYS_PER_JULIAN_MILLENNIUM)
    }
}

fn to_icrs(v: [f64; 3]) -> Vector3 {
    vsop2013_to_icrs(&Vector3::from_array(v))
}

macro_rules! impl_vsop2013_planet {
    ($name:ident, $coeffs:ident, $window:expr) => {
        pub struct $name;

        impl $name {
            const THEORY: Theory = Theory {
                series: [
                    $coeffs::A,
                    $coeffs::LAMBDA,
                    $coeffs::K,
                    $coeffs::H,
                    $coeffs::Q,
                    $coeffs::P,
                ],
                window: &$window,
            };

            pub fn heliocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
                Self::THEORY.position(tdb)
            }

            pub fn heliocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
                Self::THEORY.state(tdb)
            }

            pub fn geocentric_position(&self, tdb: &TDB, earth: &Vector3) -> AstroResult<Vector3> {
                Ok(self.heliocentric_position(tdb)? - *earth)
            }

            pub fn geocentric_state(
                &self,
                tdb: &TDB,
                (earth_position, earth_velocity): &(Vector3, Vector3),
            ) -> AstroResult<(Vector3, Vector3)> {
                let (position, velocity) = self.heliocentric_state(tdb)?;
                Ok((position - *earth_position, velocity - *earth_velocity))
            }
        }
    };
}

impl_vsop2013_planet!(Vsop2013Mercury, mercury, VSOP2013);
impl_vsop2013_planet!(Vsop2013Venus, venus, VSOP2013);
impl_vsop2013_planet!(Vsop2013Mars, mars, VSOP2013);
impl_vsop2013_planet!(Vsop2013Jupiter, jupiter, VSOP2013);
impl_vsop2013_planet!(Vsop2013Saturn, saturn, VSOP2013);
impl_vsop2013_planet!(Vsop2013Uranus, uranus, VSOP2013);
impl_vsop2013_planet!(Vsop2013Neptune, neptune, VSOP2013);
impl_vsop2013_planet!(Vsop2013Pluto, pluto, VSOP2013_PLUTO);

impl_vsop2013_planet!(Vsop2013Emb, emb, VSOP2013);
