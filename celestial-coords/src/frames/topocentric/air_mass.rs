use super::TopocentricPosition;
use celestial_core::constants::DEG_TO_RAD;

impl TopocentricPosition {
    pub fn air_mass(&self) -> f64 {
        self.air_mass_rozenberg()
    }

    pub fn air_mass_rozenberg(&self) -> f64 {
        let zenith = self.zenith_angle();
        if zenith.degrees() >= 90.0 {
            return 40.0;
        }
        let cos_z = libm::cos(zenith.radians());
        let term = cos_z + 0.025 * libm::exp(-11.0 * cos_z);
        1.0 / term
    }

    /// Computes airmass using Pickering's (2002) empirical formula.
    ///
    /// # Valid Range
    /// - Returns `f64::INFINITY` at and below about -1.115° elevation, where
    ///   the formula's corrected altitude reaches zero
    /// - Values grow without bound as elevation approaches that limit
    ///
    /// # Numerical Stability
    /// Near the horizon (0° to 5°), results can be very large but remain finite.
    /// Use this method only if observations extend below the horizon; otherwise
    /// prefer `air_mass_rozenberg()` or `air_mass_kasten_young()`.
    ///
    /// Reference: Pickering, K. A. (2002). "The Southern Limits of the Ancient Star
    /// Catalog". DIO 12, 3-27.
    pub fn air_mass_pickering(&self) -> f64 {
        let h = self.elevation.degrees();
        let altitude = h + 244.0 / (165.0 + 47.0 * libm::pow(libm::fabs(h), 1.1));
        if altitude <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / libm::sin(altitude * DEG_TO_RAD)
    }

    pub fn air_mass_kasten_young(&self) -> f64 {
        let zenith_deg = self.zenith_angle().degrees();
        if zenith_deg >= 90.0 {
            return 38.0;
        }
        let cos_z = libm::cos(self.zenith_angle().radians());
        let term = libm::pow(96.07995 - zenith_deg, -1.6364);
        1.0 / (cos_z + 0.50572 * term)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::topocentric::test_observer;
    use crate::test_support::rounded;
    use celestial_time::scales::tt::TT;

    fn at(elevation_deg: f64) -> TopocentricPosition {
        TopocentricPosition::from_degrees(0.0, elevation_deg, test_observer(), TT::j2000()).unwrap()
    }

    // Each published formula evaluated independently, to 10 decimals:
    // [elevation, Rozenberg, Pickering, Kasten-Young].
    const AIR_MASS: [[f64; 4]; 5] = [
        [90.0, 0.9999995825, 1.0000001962, 0.9997119919],
        [60.0, 1.1546981080, 1.1540579206, 1.1539922334],
        [30.0, 1.9995914063, 1.9931538464, 1.9942928525],
        [15.0, 3.8421715212, 3.8081682505, 3.8129118692],
        [5.0, 10.3369439795, 10.3337055994, 10.3057913279],
    ];

    #[test]
    fn test_zenith_angle() {
        assert_eq!(at(90.0).zenith_angle().degrees(), 0.0);
        // pi/2 and pi/3 each round to the nearest double, which leaves 90 - 60 a
        // little over an ulp below pi/6.
        assert_eq!(rounded(at(60.0).zenith_angle().degrees(), 12), 30.0);
    }

    #[test]
    fn test_air_mass_formulas() {
        for [elevation, rozenberg, pickering, kasten_young] in AIR_MASS {
            let p = at(elevation);
            let computed = [
                p.air_mass_rozenberg(),
                p.air_mass_pickering(),
                p.air_mass_kasten_young(),
            ]
            .map(|x| rounded(x, 10));
            assert_eq!(
                computed,
                [rozenberg, pickering, kasten_young],
                "elevation {elevation}"
            );
        }
    }

    #[test]
    fn test_air_mass_is_rozenberg() {
        let p = at(30.0);
        assert_eq!(p.air_mass(), p.air_mass_rozenberg());
    }

    #[test]
    fn test_air_mass_at_and_below_horizon() {
        for elevation in [0.0, -1.0, -10.0] {
            let p = at(elevation);
            assert_eq!(p.air_mass_rozenberg(), 40.0, "elevation {elevation}");
            assert_eq!(p.air_mass_kasten_young(), 38.0, "elevation {elevation}");
        }
        assert_eq!(rounded(at(0.0).air_mass_pickering(), 10), 38.7493987558);
        assert_eq!(rounded(at(-1.0).air_mass_pickering(), 10), 379.5849783511);
    }

    #[test]
    fn test_air_mass_formulas_agree_and_rise_toward_horizon() {
        let elevations = [
            90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 5.0, 2.0, 0.0,
        ];
        let mut previous = [0.0; 3];
        for elevation in elevations {
            let p = at(elevation);
            let values = [
                p.air_mass_rozenberg(),
                p.air_mass_pickering(),
                p.air_mass_kasten_young(),
            ];
            for (value, before) in values.iter().zip(previous) {
                assert!(*value > before, "elevation {elevation}");
            }
            if elevation > 30.0 {
                let mean = values.iter().sum::<f64>() / 3.0;
                for value in values {
                    assert!(
                        libm::fabs(value - mean) / mean < 0.05,
                        "elevation {elevation}"
                    );
                }
            }
            previous = values;
        }
    }

    #[test]
    fn test_air_mass_pickering_below_its_pole_is_infinite() {
        for elevation in [-1.12, -1.5, -1.99, -2.0, -30.0, -90.0] {
            assert_eq!(
                at(elevation).air_mass_pickering(),
                f64::INFINITY,
                "elevation {elevation}"
            );
        }
        assert_eq!(rounded(at(-1.11).air_mass_pickering(), 10), 5345.0513018512);
    }
}
