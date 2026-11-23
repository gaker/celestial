use celestial_core::angle::Angle;

use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;

// Paper II §2 evaluates celestial axes in degrees. Paper I Table 6 defines
// the other angular units by their size in degrees, so they divide straight
// into degrees rather than rounding once more through radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum AngleUnit {
    Degree,
    Arcmin,
    Arcsec,
    Milliarcsec,
    Radian,
}

impl AngleUnit {
    fn from_cunit(unit: &str) -> Option<Self> {
        match unit {
            "" | "deg" => Some(Self::Degree),
            "arcmin" => Some(Self::Arcmin),
            "arcsec" => Some(Self::Arcsec),
            "mas" => Some(Self::Milliarcsec),
            "rad" => Some(Self::Radian),
            _ => None,
        }
    }

    pub(super) fn to_degrees(self, value: f64) -> f64 {
        match self {
            Self::Degree => value,
            Self::Arcmin => value / 60.0,
            Self::Arcsec => value / 3600.0,
            Self::Milliarcsec => value / 3_600_000.0,
            Self::Radian => Angle::from_radians(value).degrees(),
        }
    }
}

// FITS strings drop trailing blanks; an absent CUNIT is degrees.
pub(super) fn celestial_units(header: &impl KeywordProvider) -> WcsResult<[AngleUnit; 2]> {
    let unit = |key: &str| match header.get_string(key) {
        None => Ok(AngleUnit::Degree),
        Some(unit) => AngleUnit::from_cunit(unit.trim_end()).ok_or_else(|| {
            WcsError::invalid_keyword(key, format!("'{unit}' is not an angular unit"))
        }),
    };
    Ok([unit("CUNIT1")?, unit("CUNIT2")?])
}

pub(super) fn in_degrees(values: [f64; 2], units: [AngleUnit; 2]) -> [f64; 2] {
    [0, 1].map(|i| units[i].to_degrees(values[i]))
}
