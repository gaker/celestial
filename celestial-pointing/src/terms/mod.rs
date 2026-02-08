pub mod equatorial;
pub mod altaz;
pub mod harmonic;

use bitflags::bitflags;
use crate::error::{Error, Result};

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct MountTypeFlags: u8 {
        const EQUATORIAL = 0b001;
        const ALTAZ = 0b010;
        const ALL = 0b111;
    }
}

pub trait Term: Send + Sync {
    fn name(&self) -> &str;
    fn description(&self) -> &str;

    fn jacobian_equatorial(&self, h: f64, dec: f64, lat: f64, pier: f64) -> (f64, f64);
    fn jacobian_altaz(&self, az: f64, el: f64, lat: f64) -> (f64, f64);

    fn pier_sensitive(&self) -> bool { false }
    fn applicable_mounts(&self) -> MountTypeFlags;
}

pub fn create_term(name: &str) -> Result<Box<dyn Term>> {
    match name.to_uppercase().as_str() {
        "IH" => Ok(Box::new(equatorial::IH)),
        "ID" => Ok(Box::new(equatorial::ID)),
        "CH" => Ok(Box::new(equatorial::CH)),
        "NP" => Ok(Box::new(equatorial::NP)),
        "MA" => Ok(Box::new(equatorial::MA)),
        "ME" => Ok(Box::new(equatorial::ME)),
        "TF" => Ok(Box::new(equatorial::TF)),
        "TX" => Ok(Box::new(equatorial::TX)),
        "DAF" => Ok(Box::new(equatorial::DAF)),
        "FO" => Ok(Box::new(equatorial::FO)),
        "HCES" => Ok(Box::new(equatorial::HCES)),
        "HCEC" => Ok(Box::new(equatorial::HCEC)),
        "DCES" => Ok(Box::new(equatorial::DCES)),
        "DCEC" => Ok(Box::new(equatorial::DCEC)),
        "IA" => Ok(Box::new(altaz::IA)),
        "IE" => Ok(Box::new(altaz::IE)),
        "CA" => Ok(Box::new(altaz::CA)),
        "NPAE" => Ok(Box::new(altaz::NPAE)),
        "AN" => Ok(Box::new(altaz::AN)),
        "AW" => Ok(Box::new(altaz::AW)),
        name if name.starts_with('H') => {
            let spec = harmonic::parse_harmonic(name)?;
            Ok(Box::new(spec))
        }
        _ => Err(Error::UnknownTerm(name.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_known_equatorial_terms() {
        let names = ["IH", "ID", "CH", "NP", "MA", "ME", "TF", "TX",
                      "DAF", "FO", "HCES", "HCEC", "DCES", "DCEC"];
        for name in &names {
            let term = create_term(name).unwrap();
            assert_eq!(term.name(), *name);
            assert_eq!(term.applicable_mounts(), MountTypeFlags::EQUATORIAL);
        }
    }

    #[test]
    fn create_known_altaz_terms() {
        let names = ["IA", "IE", "CA", "NPAE", "AN", "AW"];
        for name in &names {
            let term = create_term(name).unwrap();
            assert_eq!(term.name(), *name);
            assert_eq!(term.applicable_mounts(), MountTypeFlags::ALTAZ);
        }
    }

    #[test]
    fn create_harmonic_term() {
        let term = create_term("HDSH").unwrap();
        assert_eq!(term.name(), "HDSH");
        assert_eq!(term.applicable_mounts(), MountTypeFlags::EQUATORIAL);
    }

    #[test]
    fn unknown_term_returns_error() {
        let result = create_term("ZZZZ");
        assert!(result.is_err());
    }

    #[test]
    fn case_insensitive_lookup() {
        let term = create_term("ih").unwrap();
        assert_eq!(term.name(), "IH");
    }
}
