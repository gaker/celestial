use celestial_core::Angle;
use celestial_time::JulianDate;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MountType {
    GermanEquatorial,
    ForkEquatorial,
    Altazimuth,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PierSide {
    East,
    West,
    Unknown,
}

impl PierSide {
    pub fn sign(&self) -> f64 {
        match self {
            PierSide::East => 1.0,
            PierSide::West => -1.0,
            PierSide::Unknown => 1.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SiteParams {
    pub latitude: Angle,
    pub longitude: Angle,
    pub temperature: f64,
    pub pressure: f64,
    pub elevation: f64,
    pub humidity: f64,
    pub wavelength: f64,
    pub lapse_rate: f64,
}

#[derive(Debug, Clone)]
pub struct Observation {
    pub catalog_ra: Angle,
    pub catalog_dec: Angle,
    pub observed_ra: Angle,
    pub observed_dec: Angle,
    pub lst: Angle,
    pub commanded_ha: Angle,
    pub actual_ha: Angle,
    pub pier_side: PierSide,
    pub masked: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum IndatOption {
    NoDA,
    AllSky,
    Equinox,
    Equatorial,
    Altaz,
    Gimbal { z: Angle, y: Angle, x: Angle },
    RotatorTelescope,
    RotatorNasmythLeft,
    RotatorNasmythRight,
    RotatorCoudeLeft,
    RotatorCoudeRight,
}

#[derive(Debug, Clone)]
pub struct IndatFile {
    pub site: SiteParams,
    pub options: Vec<IndatOption>,
    pub observations: Vec<Observation>,
    pub mount_type: MountType,
    pub header_lines: Vec<String>,
    pub date: JulianDate,
}

pub fn decode_pier_side(raw_dec_deg: f64) -> (Angle, PierSide) {
    if raw_dec_deg.abs() > 90.0 {
        let sign = raw_dec_deg.signum();
        let dec_sky_deg = sign * (180.0 - raw_dec_deg.abs());
        (Angle::from_degrees(dec_sky_deg), PierSide::West)
    } else {
        (Angle::from_degrees(raw_dec_deg), PierSide::East)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pier_east_positive_dec() {
        let (dec, side) = decode_pier_side(45.0);
        assert_eq!(side, PierSide::East);
        assert_eq!(dec, Angle::from_degrees(45.0));
    }

    #[test]
    fn pier_east_negative_dec() {
        let (dec, side) = decode_pier_side(-30.0);
        assert_eq!(side, PierSide::East);
        assert_eq!(dec, Angle::from_degrees(-30.0));
    }

    #[test]
    fn pier_west_positive_dec() {
        let (dec, side) = decode_pier_side(135.0);
        assert_eq!(side, PierSide::West);
        assert_eq!(dec, Angle::from_degrees(45.0));
    }

    #[test]
    fn pier_west_negative_dec() {
        let (dec, side) = decode_pier_side(-135.0);
        assert_eq!(side, PierSide::West);
        assert_eq!(dec, Angle::from_degrees(-45.0));
    }

    #[test]
    fn pier_east_at_boundary() {
        let (dec, side) = decode_pier_side(90.0);
        assert_eq!(side, PierSide::East);
        assert_eq!(dec, Angle::from_degrees(90.0));
    }

    #[test]
    fn pier_east_zero_dec() {
        let (dec, side) = decode_pier_side(0.0);
        assert_eq!(side, PierSide::East);
        assert_eq!(dec, Angle::from_degrees(0.0));
    }

    #[test]
    fn pier_side_signs() {
        assert_eq!(PierSide::East.sign(), 1.0);
        assert_eq!(PierSide::West.sign(), -1.0);
        assert_eq!(PierSide::Unknown.sign(), 1.0);
    }
}
