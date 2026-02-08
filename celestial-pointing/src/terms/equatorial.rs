use super::{MountTypeFlags, Term};

pub struct IH;
pub struct ID;
pub struct CH;
pub struct NP;
pub struct MA;
pub struct ME;
pub struct TF;
pub struct TX;
pub struct DAF;
pub struct FO;
pub struct HCES;
pub struct HCEC;
pub struct DCES;
pub struct DCEC;

impl Term for IH {
    fn name(&self) -> &str { "IH" }
    fn description(&self) -> &str { "Hour angle index error" }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (-1.0, 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for ID {
    fn name(&self) -> &str { "ID" }
    fn description(&self) -> &str { "Declination index error" }
    fn pier_sensitive(&self) -> bool { true }
    fn jacobian_equatorial(&self, _h: f64, _dec: f64, _lat: f64, pier: f64) -> (f64, f64) {
        (0.0, -pier)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for CH {
    fn name(&self) -> &str { "CH" }
    fn description(&self) -> &str { "East-west collimation error" }
    fn pier_sensitive(&self) -> bool { true }
    fn jacobian_equatorial(&self, _h: f64, dec: f64, _lat: f64, pier: f64) -> (f64, f64) {
        (-pier / dec.cos(), 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for NP {
    fn name(&self) -> &str { "NP" }
    fn description(&self) -> &str { "Non-perpendicularity of axes" }
    fn pier_sensitive(&self) -> bool { true }
    fn jacobian_equatorial(&self, _h: f64, dec: f64, _lat: f64, pier: f64) -> (f64, f64) {
        (-pier * dec.tan(), 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for MA {
    fn name(&self) -> &str { "MA" }
    fn description(&self) -> &str { "Polar axis azimuth misalignment" }
    fn jacobian_equatorial(&self, h: f64, dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (h.cos() * dec.tan(), -h.sin())
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for ME {
    fn name(&self) -> &str { "ME" }
    fn description(&self) -> &str { "Polar axis elevation misalignment" }
    fn jacobian_equatorial(&self, h: f64, dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (-h.sin() * dec.tan(), -h.cos())
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for TF {
    fn name(&self) -> &str { "TF" }
    fn description(&self) -> &str { "Tube flexure (sin zeta)" }
    fn jacobian_equatorial(&self, h: f64, dec: f64, lat: f64, _pier: f64) -> (f64, f64) {
        let dh = lat.cos() * h.sin() / dec.cos();
        let dd = lat.cos() * h.cos() * dec.sin() - lat.sin() * dec.cos();
        (dh, dd)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for TX {
    fn name(&self) -> &str { "TX" }
    fn description(&self) -> &str { "Tube flexure (tan zeta)" }
    fn jacobian_equatorial(&self, h: f64, dec: f64, lat: f64, _pier: f64) -> (f64, f64) {
        let sin_alt = lat.sin() * dec.sin() + lat.cos() * dec.cos() * h.cos();
        let alt = sin_alt.asin();
        let sin_a = alt.sin();
        if sin_a.abs() < 1e-10 {
            return (0.0, 0.0);
        }
        let dh_tf = lat.cos() * h.sin() / dec.cos();
        let dd_tf = lat.cos() * h.cos() * dec.sin() - lat.sin() * dec.cos();
        (dh_tf / sin_a, dd_tf / sin_a)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for DAF {
    fn name(&self) -> &str { "DAF" }
    fn description(&self) -> &str { "Declination axis flexure" }
    fn jacobian_equatorial(&self, h: f64, dec: f64, lat: f64, _pier: f64) -> (f64, f64) {
        (-(lat.sin() * dec.tan() + lat.cos() * h.cos()), 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for FO {
    fn name(&self) -> &str { "FO" }
    fn description(&self) -> &str { "Fork flexure" }
    fn jacobian_equatorial(&self, h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, h.cos())
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for HCES {
    fn name(&self) -> &str { "HCES" }
    fn description(&self) -> &str { "HA centering error (sin)" }
    fn jacobian_equatorial(&self, h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (h.sin(), 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for HCEC {
    fn name(&self) -> &str { "HCEC" }
    fn description(&self) -> &str { "HA centering error (cos)" }
    fn jacobian_equatorial(&self, h: f64, _dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (h.cos(), 0.0)
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for DCES {
    fn name(&self) -> &str { "DCES" }
    fn description(&self) -> &str { "Dec centering error (sin)" }
    fn jacobian_equatorial(&self, _h: f64, dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, dec.sin())
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

impl Term for DCEC {
    fn name(&self) -> &str { "DCEC" }
    fn description(&self) -> &str { "Dec centering error (cos)" }
    fn jacobian_equatorial(&self, _h: f64, dec: f64, _lat: f64, _pier: f64) -> (f64, f64) {
        (0.0, dec.cos())
    }
    fn jacobian_altaz(&self, _az: f64, _el: f64, _lat: f64) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn applicable_mounts(&self) -> MountTypeFlags { MountTypeFlags::EQUATORIAL }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    #[test]
    fn ih_jacobian() {
        let (dh, dd) = IH.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, -1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn id_pier_east() {
        let (dh, dd) = ID.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, -1.0);
    }

    #[test]
    fn id_pier_west() {
        let (dh, dd) = ID.jacobian_equatorial(0.0, 0.0, 0.0, -1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 1.0);
    }

    #[test]
    fn ch_at_zero_dec_east() {
        let (dh, dd) = CH.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, -1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn ch_at_zero_dec_west() {
        let (dh, dd) = CH.jacobian_equatorial(0.0, 0.0, 0.0, -1.0);
        assert_eq!(dh, 1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn np_at_zero_dec() {
        let (dh, dd) = NP.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn np_at_45_dec_east() {
        let dec = FRAC_PI_4;
        let (dh, dd) = NP.jacobian_equatorial(0.0, dec, 0.0, 1.0);
        assert_eq!(dh, -dec.tan());
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn ma_at_zero_ha() {
        let dec = FRAC_PI_4;
        let (dh, dd) = MA.jacobian_equatorial(0.0, dec, 0.0, 1.0);
        assert_eq!(dh, dec.tan());
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn ma_at_6h_ha() {
        let h = FRAC_PI_2;
        let dec = FRAC_PI_4;
        let (dh, dd) = MA.jacobian_equatorial(h, dec, 0.0, 1.0);
        assert_eq!(dh, h.cos() * dec.tan());
        assert_eq!(dd, -h.sin());
    }

    #[test]
    fn me_at_zero_ha() {
        let dec = FRAC_PI_4;
        let (dh, dd) = ME.jacobian_equatorial(0.0, dec, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, -1.0);
    }

    #[test]
    fn hces_at_pi_over_2() {
        let (dh, dd) = HCES.jacobian_equatorial(FRAC_PI_2, 0.0, 0.0, 1.0);
        assert_eq!(dh, 1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn hcec_at_zero() {
        let (dh, dd) = HCEC.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn dces_at_pi_over_2() {
        let (dh, dd) = DCES.jacobian_equatorial(0.0, FRAC_PI_2, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 1.0);
    }

    #[test]
    fn dcec_at_zero() {
        let (dh, dd) = DCEC.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 1.0);
    }

    #[test]
    fn fo_at_zero_ha() {
        let (dh, dd) = FO.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 1.0);
    }

    #[test]
    fn fo_at_pi() {
        let (dh, dd) = FO.jacobian_equatorial(PI, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, -1.0);
    }

    #[test]
    fn equatorial_terms_return_zero_for_altaz() {
        let terms: Vec<Box<dyn Term>> = vec![
            Box::new(IH), Box::new(ID), Box::new(CH), Box::new(NP),
            Box::new(MA), Box::new(ME), Box::new(TF), Box::new(TX),
            Box::new(DAF), Box::new(FO), Box::new(HCES), Box::new(HCEC),
            Box::new(DCES), Box::new(DCEC),
        ];
        for term in &terms {
            let (da, de) = term.jacobian_altaz(1.0, 0.5, 0.7);
            assert_eq!((da, de), (0.0, 0.0), "term {} should return (0,0) for altaz", term.name());
        }
    }

    #[test]
    fn pier_sensitivity_flags() {
        assert!(!IH.pier_sensitive());
        assert!(ID.pier_sensitive());
        assert!(CH.pier_sensitive());
        assert!(NP.pier_sensitive());
        assert!(!MA.pier_sensitive());
        assert!(!ME.pier_sensitive());
        assert!(!TF.pier_sensitive());
    }
}
