use super::{MountTypeFlags, Term};
use crate::error::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrigFunc {
    Sin,
    Cos,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Coordinate {
    H,
    D,
    A,
    E,
    Z,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResultType {
    H,
    D,
    X,
    P,
    A,
    E,
    Z,
    S,
}

#[derive(Debug, Clone)]
pub struct HarmonicComponent {
    pub func: TrigFunc,
    pub coord: Coordinate,
    pub frequency: u8,
}

#[derive(Debug, Clone)]
pub struct HarmonicTerm {
    pub result: ResultType,
    pub components: Vec<HarmonicComponent>,
    pub original_name: String,
}

pub fn parse_harmonic(name: &str) -> Result<HarmonicTerm> {
    let chars: Vec<char> = name.chars().collect();
    if chars.len() < 4 || chars[0] != 'H' {
        return Err(Error::InvalidHarmonic(format!(
            "too short or missing H prefix: {}",
            name
        )));
    }
    let result = parse_result_type(chars[1], name)?;
    let components = parse_components(&chars[2..], name)?;
    if components.is_empty() {
        return Err(Error::InvalidHarmonic(format!(
            "no components found: {}",
            name
        )));
    }
    Ok(HarmonicTerm {
        result,
        components,
        original_name: name.to_string(),
    })
}

fn parse_result_type(ch: char, name: &str) -> Result<ResultType> {
    match ch {
        'H' => Ok(ResultType::H),
        'D' => Ok(ResultType::D),
        'X' => Ok(ResultType::X),
        'P' => Ok(ResultType::P),
        'A' => Ok(ResultType::A),
        'E' => Ok(ResultType::E),
        'Z' => Ok(ResultType::Z),
        'S' => Ok(ResultType::S),
        _ => Err(Error::InvalidHarmonic(format!(
            "unknown result type '{}' in {}",
            ch, name
        ))),
    }
}

fn parse_func(ch: char, name: &str) -> Result<TrigFunc> {
    match ch {
        'S' => Ok(TrigFunc::Sin),
        'C' => Ok(TrigFunc::Cos),
        _ => Err(Error::InvalidHarmonic(format!(
            "unknown function '{}' in {}",
            ch, name
        ))),
    }
}

fn parse_coord(ch: char, name: &str) -> Result<Coordinate> {
    match ch {
        'H' => Ok(Coordinate::H),
        'D' => Ok(Coordinate::D),
        'A' => Ok(Coordinate::A),
        'E' => Ok(Coordinate::E),
        'Z' => Ok(Coordinate::Z),
        _ => Err(Error::InvalidHarmonic(format!(
            "unknown coordinate '{}' in {}",
            ch, name
        ))),
    }
}

fn parse_components(chars: &[char], name: &str) -> Result<Vec<HarmonicComponent>> {
    let mut components = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        if i + 1 >= chars.len() {
            return Err(Error::InvalidHarmonic(format!(
                "incomplete component in {}",
                name
            )));
        }
        let func = parse_func(chars[i], name)?;
        let coord = parse_coord(chars[i + 1], name)?;
        i += 2;
        let frequency = if i < chars.len() && chars[i].is_ascii_digit() {
            let d = chars[i].to_digit(10).unwrap() as u8;
            i += 1;
            d
        } else {
            1
        };
        if frequency == 0 {
            return Err(Error::InvalidHarmonic(format!(
                "frequency 0 not allowed in {}",
                name
            )));
        }
        components.push(HarmonicComponent {
            func,
            coord,
            frequency,
        });
    }
    Ok(components)
}

impl HarmonicTerm {
    fn evaluate_equatorial(&self, h: f64, dec: f64) -> f64 {
        self.components.iter().fold(1.0, |acc, comp| {
            let coord_val = match comp.coord {
                Coordinate::H => h,
                Coordinate::D => dec,
                Coordinate::A | Coordinate::E | Coordinate::Z => 0.0,
            };
            acc * trig_value(comp.func, coord_val, comp.frequency)
        })
    }

    fn evaluate_altaz(&self, az: f64, el: f64) -> f64 {
        self.components.iter().fold(1.0, |acc, comp| {
            let coord_val = match comp.coord {
                Coordinate::A => az,
                Coordinate::E => el,
                Coordinate::Z => std::f64::consts::FRAC_PI_2 - el,
                Coordinate::H | Coordinate::D => 0.0,
            };
            acc * trig_value(comp.func, coord_val, comp.frequency)
        })
    }
}

fn trig_value(func: TrigFunc, coord_val: f64, frequency: u8) -> f64 {
    let arg = coord_val * frequency as f64;
    match func {
        TrigFunc::Sin => libm::sin(arg),
        TrigFunc::Cos => libm::cos(arg),
    }
}

impl Term for HarmonicTerm {
    fn name(&self) -> &str {
        &self.original_name
    }
    fn description(&self) -> &str {
        "Harmonic correction term"
    }

    fn pier_sensitive(&self) -> bool {
        matches!(self.result, ResultType::P)
    }

    fn applicable_mounts(&self) -> MountTypeFlags {
        match self.result {
            ResultType::H | ResultType::D | ResultType::X | ResultType::P => {
                MountTypeFlags::EQUATORIAL
            }
            ResultType::A | ResultType::E | ResultType::Z | ResultType::S => MountTypeFlags::ALTAZ,
        }
    }

    fn jacobian_equatorial(&self, h: f64, dec: f64, _lat: f64, pier: f64) -> (f64, f64) {
        let val = self.evaluate_equatorial(h, dec);
        match self.result {
            ResultType::H => (val, 0.0),
            ResultType::D => (0.0, val),
            ResultType::X => (val / libm::cos(dec), 0.0),
            ResultType::P => (-pier * libm::tan(dec) * val, 0.0),
            _ => (0.0, 0.0),
        }
    }

    fn jacobian_altaz(&self, az: f64, el: f64, _lat: f64) -> (f64, f64) {
        let val = self.evaluate_altaz(az, el);
        match self.result {
            ResultType::A => (val, 0.0),
            ResultType::E => (0.0, val),
            ResultType::Z => (0.0, -val),
            ResultType::S => (val / libm::cos(el), 0.0),
            _ => (0.0, 0.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    #[test]
    fn parse_hdsh() {
        let term = parse_harmonic("HDSH").unwrap();
        assert_eq!(term.result, ResultType::D);
        assert_eq!(term.components.len(), 1);
        assert_eq!(term.components[0].func, TrigFunc::Sin);
        assert_eq!(term.components[0].coord, Coordinate::H);
        assert_eq!(term.components[0].frequency, 1);
    }

    #[test]
    fn parse_hxch3() {
        let term = parse_harmonic("HXCH3").unwrap();
        assert_eq!(term.result, ResultType::X);
        assert_eq!(term.components.len(), 1);
        assert_eq!(term.components[0].func, TrigFunc::Cos);
        assert_eq!(term.components[0].coord, Coordinate::H);
        assert_eq!(term.components[0].frequency, 3);
    }

    #[test]
    fn parse_hdch2cd4() {
        let term = parse_harmonic("HDCH2CD4").unwrap();
        assert_eq!(term.result, ResultType::D);
        assert_eq!(term.components.len(), 2);
        assert_eq!(term.components[0].func, TrigFunc::Cos);
        assert_eq!(term.components[0].coord, Coordinate::H);
        assert_eq!(term.components[0].frequency, 2);
        assert_eq!(term.components[1].func, TrigFunc::Cos);
        assert_eq!(term.components[1].coord, Coordinate::D);
        assert_eq!(term.components[1].frequency, 4);
    }

    #[test]
    fn parse_too_short() {
        assert!(parse_harmonic("HD").is_err());
    }

    #[test]
    fn parse_bad_prefix() {
        assert!(parse_harmonic("XDSH").is_err());
    }

    #[test]
    fn parse_bad_result_type() {
        assert!(parse_harmonic("HQSH").is_err());
    }

    #[test]
    fn parse_bad_function() {
        assert!(parse_harmonic("HDXH").is_err());
    }

    #[test]
    fn parse_bad_coordinate() {
        assert!(parse_harmonic("HDSQ").is_err());
    }

    #[test]
    fn parse_zero_frequency() {
        assert!(parse_harmonic("HDSH0").is_err());
    }

    #[test]
    fn hdsh_jacobian_at_pi_over_2() {
        let term = parse_harmonic("HDSH").unwrap();
        let (dh, dd) = term.jacobian_equatorial(FRAC_PI_2, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 1.0);
    }

    #[test]
    fn hdsh_jacobian_at_zero() {
        let term = parse_harmonic("HDSH").unwrap();
        let (dh, dd) = term.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn hhch_jacobian_at_zero() {
        let term = parse_harmonic("HHCH").unwrap();
        let (dh, dd) = term.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn hxch_at_zero_ha_zero_dec() {
        let term = parse_harmonic("HXCH").unwrap();
        let (dh, dd) = term.jacobian_equatorial(0.0, 0.0, 0.0, 1.0);
        assert_eq!(dh, 1.0);
        assert_eq!(dd, 0.0);
    }

    #[test]
    fn hp_result_pier_sensitive() {
        let term = parse_harmonic("HPSH").unwrap();
        assert!(term.pier_sensitive());
    }

    #[test]
    fn hd_result_not_pier_sensitive() {
        let term = parse_harmonic("HDSH").unwrap();
        assert!(!term.pier_sensitive());
    }

    #[test]
    fn equatorial_result_mount_flags() {
        for prefix in &["HH", "HD", "HX", "HP"] {
            let name = format!("{}SH", prefix);
            let term = parse_harmonic(&name).unwrap();
            assert_eq!(term.applicable_mounts(), MountTypeFlags::EQUATORIAL);
        }
    }

    #[test]
    fn altaz_result_mount_flags() {
        for prefix in &["HA", "HE", "HZ", "HS"] {
            let name = format!("{}SA", prefix);
            let term = parse_harmonic(&name).unwrap();
            assert_eq!(term.applicable_mounts(), MountTypeFlags::ALTAZ);
        }
    }

    #[test]
    fn altaz_harmonic_jacobian() {
        let term = parse_harmonic("HASA").unwrap();
        let az = FRAC_PI_2;
        let (da, de) = term.jacobian_altaz(az, 0.0, 0.0);
        assert_eq!(da, 1.0);
        assert_eq!(de, 0.0);
    }

    #[test]
    fn altaz_elevation_result() {
        let term = parse_harmonic("HESE").unwrap();
        let el = FRAC_PI_4;
        let (da, de) = term.jacobian_altaz(0.0, el, 0.0);
        assert_eq!(da, 0.0);
        assert_eq!(de, libm::sin(el));
    }

    #[test]
    fn altaz_zenith_result() {
        let term = parse_harmonic("HZSE").unwrap();
        let el = FRAC_PI_4;
        let (da, de) = term.jacobian_altaz(0.0, el, 0.0);
        assert_eq!(da, 0.0);
        assert_eq!(de, -libm::sin(el));
    }

    #[test]
    fn two_component_product() {
        let term = parse_harmonic("HDCH2CD4").unwrap();
        let h = 0.0;
        let dec = 0.0;
        let (dh, dd) = term.jacobian_equatorial(h, dec, 0.0, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(dd, 1.0);
    }

    #[test]
    fn equatorial_harmonic_returns_zero_for_altaz() {
        let term = parse_harmonic("HDSH").unwrap();
        let (da, de) = term.jacobian_altaz(1.0, 0.5, 0.7);
        assert_eq!((da, de), (0.0, 0.0));
    }

    #[test]
    fn altaz_harmonic_returns_zero_for_equatorial() {
        let term = parse_harmonic("HASA").unwrap();
        let (dh, dd) = term.jacobian_equatorial(1.0, 0.5, 0.7, 1.0);
        assert_eq!((dh, dd), (0.0, 0.0));
    }

    #[test]
    fn original_name_preserved() {
        let term = parse_harmonic("HDSH").unwrap();
        assert_eq!(term.name(), "HDSH");
    }

    #[test]
    fn frequency_3_multiplies_argument() {
        let term = parse_harmonic("HHSH3").unwrap();
        let h = FRAC_PI_2;
        let (dh, _dd) = term.jacobian_equatorial(h, 0.0, 0.0, 1.0);
        assert_eq!(dh, libm::sin(3.0 * h));
    }
}
