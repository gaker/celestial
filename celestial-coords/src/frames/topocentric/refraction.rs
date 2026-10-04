use super::TopocentricPosition;
use crate::errors::{CoordError, CoordResult};
use celestial_core::angle::Angle;
use celestial_core::constants::HALF_PI;

// Near the horizon the model is held at its value for these cos/sin of
// elevation (about 2.9 degrees), as ERFA does.
const MIN_COS_ELEVATION: f64 = 1e-6;
const MIN_SIN_ELEVATION: f64 = 0.05;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Refraction {
    a: f64,
    b: f64,
}

impl Refraction {
    pub fn new(
        pressure_hpa: f64,
        temperature_c: f64,
        relative_humidity: f64,
        wavelength_um: f64,
    ) -> CoordResult<Self> {
        let conditions = [
            pressure_hpa,
            temperature_c,
            relative_humidity,
            wavelength_um,
        ];
        if conditions.iter().any(|x| !x.is_finite()) {
            return Err(CoordError::invalid_coordinate(
                "refraction conditions must be finite",
            ));
        }
        let (gamma, beta) = gamma_beta(conditions);
        Ok(Self {
            a: gamma * (1.0 - beta),
            b: -gamma * (beta - gamma / 2.0),
        })
    }

    pub fn a(&self) -> f64 {
        self.a
    }

    pub fn b(&self) -> f64 {
        self.b
    }
}

// The refractivity gamma and the beta term of the model, with each condition clamped to the
// range the formulae hold over. Wavelengths above 100 µm take the radio formulae.
fn gamma_beta([pressure_hpa, temperature_c, humidity, wavelength_um]: [f64; 4]) -> (f64, f64) {
    let p = bounded(pressure_hpa, 0.0, 10000.0);
    let t = bounded(temperature_c, -150.0, 200.0);
    let pw = water_vapour_pressure(p, t, bounded(humidity, 0.0, 1.0));
    let tk = t + 273.15;
    if wavelength_um <= 100.0 {
        let w = bounded(wavelength_um, 0.1, 1e6);
        (optical_gamma(p, pw, tk, w), 4.4474e-6 * tk)
    } else {
        (radio_gamma(p, pw, tk), radio_beta(pw, tk))
    }
}

// Max then min, not f64::clamp: an input equal to a bound takes the bound
// itself, so -0.0 becomes +0.0 as it does in eraRefco.
fn bounded(x: f64, low: f64, high: f64) -> f64 {
    let x = if x > low { x } else { low };
    if x < high {
        x
    } else {
        high
    }
}

fn water_vapour_pressure(p: f64, t: f64, r: f64) -> f64 {
    if p > 0.0 {
        let ps = libm::pow(10.0, (0.7859 + 0.03477 * t) / (1.0 + 0.00412 * t))
            * (1.0 + p * (4.5e-6 + 6e-10 * t * t));
        r * ps / (1.0 - (1.0 - r) * ps / p)
    } else {
        0.0
    }
}

fn optical_gamma(p: f64, pw: f64, tk: f64, w: f64) -> f64 {
    let wlsq = w * w;
    ((77.53484e-6 + (4.39108e-7 + 3.666e-9 / wlsq) / wlsq) * p - 11.2684e-6 * pw) / tk
}

fn radio_gamma(p: f64, pw: f64, tk: f64) -> f64 {
    (77.6890e-6 * p - (6.3938e-6 - 0.375463 / tk) * pw) / tk
}

fn radio_beta(pw: f64, tk: f64) -> f64 {
    let beta = 4.4474e-6 * tk;
    beta - 0.0074 * pw * beta
}

impl TopocentricPosition {
    // self is the observed place.
    pub fn atmospheric_refraction(&self, refraction: &Refraction) -> Angle {
        let (_, _, dref) = self.observed_refraction(refraction);
        Angle::from_radians(dref)
    }

    // self is the true place; the result is where it is observed.
    pub fn with_refraction(&self, refraction: &Refraction) -> Self {
        let (sa, ca) = self.azimuth.sin_cos();
        let (se, ce) = self.elevation.sin_cos();
        let (x, y, z) = (ca * ce, sa * ce, se);
        let r = libm::fmax(libm::sqrt(x * x + y * y), MIN_COS_ELEVATION);
        let zc = libm::fmax(z, MIN_SIN_ELEVATION);
        let tz = r / zc;
        let w = refraction.b * tz * tz;
        let del = (refraction.a + w) * tz / (1.0 + (refraction.a + 3.0 * w) / (zc * zc));
        let cosdel = 1.0 - del * del / 2.0;
        let f = cosdel - del * zc / r;
        let (xo, yo, zo) = (x * f, y * f, cosdel * z + del * r);
        let zd = libm::atan2(libm::sqrt(xo * xo + yo * yo), zo);
        Self {
            elevation: Angle::from_radians(HALF_PI - zd),
            ..*self
        }
    }

    // self is the observed place; the result is the true place.
    pub fn without_refraction(&self, refraction: &Refraction) -> Self {
        let (az, zdo, dref) = self.observed_refraction(refraction);
        let (st, ct) = libm::sincos(zdo + dref);
        let (sa, ca) = libm::sincos(az);
        let (x, y, z) = (ca * st, sa * st, ct);
        let elevation = if z == 0.0 {
            0.0
        } else {
            libm::atan2(z, libm::sqrt(x * x + y * y))
        };
        Self {
            elevation: Angle::from_radians(elevation),
            ..*self
        }
    }

    // The azimuth (south = 0) and zenith distance of the observed place
    // rebuilt as a vector, and the refraction to add to that distance.
    fn observed_refraction(&self, refraction: &Refraction) -> (f64, f64, f64) {
        let (sa, ca) = self.azimuth.sin_cos();
        let (sz, cz) = libm::sincos(HALF_PI - self.elevation.radians());
        let (x, y, z) = (-ca * sz, sa * sz, cz);
        let az = if x != 0.0 || y != 0.0 {
            libm::atan2(y, x)
        } else {
            0.0
        };
        let s = libm::sqrt(x * x + y * y);
        let zdo = libm::atan2(s, z);
        let tz = s / if z > MIN_SIN_ELEVATION {
            z
        } else {
            MIN_SIN_ELEVATION
        };
        let dref = (refraction.a + refraction.b * tz * tz) * tz;
        (az, zdo, dref)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::topocentric::test_observer;
    use celestial_time::scales::tt::TT;

    #[test]
    fn test_atmospheric_refraction_standard_conditions() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Standard conditions: sea level, 15°C, 50% humidity, optical (0.574 μm)
        let standard = Refraction::new(1013.25, 15.0, 0.5, 0.574).unwrap();

        let zenith = TopocentricPosition::from_degrees(0.0, 90.0, observer, epoch).unwrap();
        assert_eq!(zenith.atmospheric_refraction(&standard).radians(), 0.0);

        let pos_45 = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        let ref_45 = pos_45.atmospheric_refraction(&standard);
        // Typical refraction at 45° elevation ~60 arcsec
        assert!(ref_45.arcseconds() > 50.0 && ref_45.arcseconds() < 70.0);

        let pos_10 = TopocentricPosition::from_degrees(0.0, 10.0, observer, epoch).unwrap();
        let ref_10 = pos_10.atmospheric_refraction(&standard);
        // Refraction increases dramatically near horizon, ~5-6 arcmin
        assert!(ref_10.arcminutes() > 4.0 && ref_10.arcminutes() < 7.0);
    }

    #[test]
    fn test_atmospheric_refraction_with_without() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let standard = Refraction::new(1013.25, 15.0, 0.5, 0.574).unwrap();

        let true_pos = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        let apparent = true_pos.with_refraction(&standard);
        assert!(apparent.elevation().degrees() > true_pos.elevation().degrees());

        // The forward model takes a single Newton step, so the two directions are not exact
        // inverses; at 45 degrees they disagree by about 8 µas.
        let back_to_true = apparent.without_refraction(&standard);
        let error = back_to_true.elevation().radians() - true_pos.elevation().radians();
        assert!(libm::fabs(error) < 1e-10, "{error:e} rad");
    }

    #[test]
    fn test_atmospheric_refraction_zero_pressure() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Zero pressure = no atmosphere = no refraction
        let pos = TopocentricPosition::from_degrees(0.0, 30.0, observer, epoch).unwrap();
        let refraction =
            pos.atmospheric_refraction(&Refraction::new(0.0, 15.0, 0.5, 0.574).unwrap());

        assert_eq!(refraction.radians(), 0.0);
    }

    #[test]
    fn test_atmospheric_refraction_radio_vs_optical() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let pos = TopocentricPosition::from_degrees(0.0, 30.0, observer, epoch).unwrap();

        // Optical wavelength (0.574 μm)
        let optical =
            pos.atmospheric_refraction(&Refraction::new(1013.25, 15.0, 0.5, 0.574).unwrap());

        // Radio wavelength (>100 μm)
        let radio =
            pos.atmospheric_refraction(&Refraction::new(1013.25, 15.0, 0.5, 200.0).unwrap());

        assert!(optical.arcseconds() > 0.0);
        // Water vapour refracts radio waves far more strongly than light.
        assert!(radio.arcseconds() > optical.arcseconds());
    }

    #[test]
    fn test_refraction_constants_match_erfa() {
        // eraRefco for the same conditions.
        let optical = Refraction::new(1013.25, 15.0, 0.5, 0.574).unwrap();
        assert_eq!(
            [optical.a(), optical.b()],
            [0.00027675592203715746, -3.167276123549832e-7]
        );
        let radio = Refraction::new(1013.25, 15.0, 0.5, 200.0).unwrap();
        assert_eq!(
            [radio.a(), radio.b()],
            [0.00031164709116880424, -3.256447759334908e-7]
        );
    }

    #[test]
    fn test_refraction_rejects_non_finite_conditions() {
        let conditions = [
            (f64::NAN, 15.0, 0.5, 0.574),
            (1013.25, f64::NAN, 0.5, 0.574),
            (1013.25, 15.0, f64::INFINITY, 0.574),
            (1013.25, 15.0, 0.5, f64::NEG_INFINITY),
        ];
        for (p, t, rh, wl) in conditions {
            assert!(Refraction::new(p, t, rh, wl).is_err(), "{p} {t} {rh} {wl}");
        }
    }
}
