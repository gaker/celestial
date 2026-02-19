use super::{MountTypeFlags, Term};

pub struct IA;
pub struct IE;
pub struct CA;
pub struct NPAE;
pub struct AN;
pub struct AW;

impl Term for IA {
    fn name(&self) -> &str {
        "IA"
    }
    fn description(&self) -> &str {
        "Azimuth index error"
    }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (-1.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags {
        MountTypeFlags::ALTAZ
    }
}

impl Term for IE {
    fn name(&self) -> &str {
        "IE"
    }
    fn description(&self) -> &str {
        "Elevation index error"
    }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 1.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags {
        MountTypeFlags::ALTAZ
    }
}

impl Term for CA {
    fn name(&self) -> &str {
        "CA"
    }
    fn description(&self) -> &str {
        "Left-right collimation error"
    }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, el: f64, _lat: f64) -> (f64, f64) {
        (-1.0 / libm::cos(el), 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags {
        MountTypeFlags::ALTAZ
    }
}

impl Term for NPAE {
    fn name(&self) -> &str {
        "NPAE"
    }
    fn description(&self) -> &str {
        "Non-perpendicularity of az/el axes"
    }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, el: f64, _lat: f64) -> (f64, f64) {
        (-libm::tan(el), 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags {
        MountTypeFlags::ALTAZ
    }
}

impl Term for AN {
    fn name(&self) -> &str {
        "AN"
    }
    fn description(&self) -> &str {
        "Azimuth axis tilt north-south"
    }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn jacobian_altaz(&self, az: f64, el: f64, _lat: f64) -> (f64, f64) {
        (-libm::sin(az) * libm::tan(el), libm::cos(az))
    }
    fn applicable_mounts(&self) -> MountTypeFlags {
        MountTypeFlags::ALTAZ
    }
}

impl Term for AW {
    fn name(&self) -> &str {
        "AW"
    }
    fn description(&self) -> &str {
        "Azimuth axis tilt east-west"
    }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn jacobian_altaz(&self, az: f64, el: f64, _lat: f64) -> (f64, f64) {
        (-libm::cos(az) * libm::tan(el), -libm::sin(az))
    }
    fn applicable_mounts(&self) -> MountTypeFlags {
        MountTypeFlags::ALTAZ
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    #[test]
    fn ia_jacobian() {
        let (da, de) = IA.jacobian_altaz(0.0, 0.0, 0.0);
        assert_eq!(da, -1.0);
        assert_eq!(de, 0.0);
    }

    #[test]
    fn ie_jacobian() {
        let (da, de) = IE.jacobian_altaz(0.0, 0.0, 0.0);
        assert_eq!(da, 0.0);
        assert_eq!(de, 1.0);
    }

    #[test]
    fn ca_at_zero_el() {
        let (da, de) = CA.jacobian_altaz(0.0, 0.0, 0.0);
        assert_eq!(da, -1.0);
        assert_eq!(de, 0.0);
    }

    #[test]
    fn npae_at_zero_el() {
        let (da, de) = NPAE.jacobian_altaz(0.0, 0.0, 0.0);
        assert_eq!(da, 0.0);
        assert_eq!(de, 0.0);
    }

    #[test]
    fn npae_at_45_el() {
        let el = FRAC_PI_4;
        let (da, de) = NPAE.jacobian_altaz(0.0, el, 0.0);
        assert_eq!(da, -el.tan());
        assert_eq!(de, 0.0);
    }

    #[test]
    fn an_at_pi_over_2_az() {
        let az = FRAC_PI_2;
        let el = FRAC_PI_4;
        let (da, de) = AN.jacobian_altaz(az, el, 0.0);
        assert_eq!(da, -libm::sin(az) * el.tan());
        assert_eq!(de, libm::cos(az));
    }

    #[test]
    fn aw_at_zero_az() {
        let el = FRAC_PI_4;
        let (da, de) = AW.jacobian_altaz(0.0, el, 0.0);
        assert_eq!(da, -el.tan());
        assert_eq!(de, 0.0);
    }

    #[test]
    fn altaz_terms_return_zero_for_equatorial() {
        let terms: Vec<Box<dyn Term>> = vec![
            Box::new(IA),
            Box::new(IE),
            Box::new(CA),
            Box::new(NPAE),
            Box::new(AN),
            Box::new(AW),
        ];
        for term in &terms {
            let (dh, dd) = term.jacobian_equatorial(1.0, 0.5, 0.7, 1.0);
            assert_eq!(
                (dh, dd),
                (0.0, 0.0),
                "term {} should return (0,0) for equatorial",
                term.name()
            );
        }
    }

    #[test]
    fn no_altaz_terms_are_pier_sensitive() {
        assert!(!IA.pier_sensitive());
        assert!(!IE.pier_sensitive());
        assert!(!CA.pier_sensitive());
        assert!(!NPAE.pier_sensitive());
        assert!(!AN.pier_sensitive());
        assert!(!AW.pier_sensitive());
    }
}
