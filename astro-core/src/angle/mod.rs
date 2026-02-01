mod core;
mod format;
mod normalize;
mod ops;
mod parse;
#[cfg(feature = "serde")]
mod serde_;
mod validate;

pub use core::Angle;
pub use format::{parse_angle, DmsFmt, HmsFmt, ParsedAngle};
pub use normalize::{clamp_dec, wrap_0_2pi, wrap_pm_pi, NormalizeMode};
pub use parse::{parse_dms, parse_hms, AngleUnits, ParseAngle};
pub use validate::{
    validate_declination, validate_latitude, validate_longitude, validate_right_ascension,
};

pub use core::{arcmin, arcsec, deg, hours, rad};
