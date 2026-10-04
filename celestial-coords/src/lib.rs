pub mod aberration;
pub(crate) mod constants;
pub mod distance;
pub mod eop;
pub mod errors;
pub mod frames;
pub mod lighttime;
pub mod lunar;
pub mod proper_motion;
pub mod solar;
pub mod transforms;

pub use celestial_core::angle::Angle;
pub use distance::Distance;
pub use eop::{EopParameters, EopProvider, EopRecord};
pub use errors::{CoordError, CoordResult};
pub use lighttime::LightTimeCorrection;

pub use frames::{
    CIRSPosition, EclipticCartesian, EclipticPosition, GCRSPosition, GalacticPosition,
    HeliographicCarrington, HeliographicStonyhurst, HourAnglePosition, ICRSPosition, ITRSPosition,
    SelenographicPosition, TIRSPosition, TopocentricPosition,
};

pub use transforms::{CartesianFrame, CoordinateFrame};

pub use celestial_core::{location::Location, matrix::Vector3};
pub use celestial_time::scales::tai::TAI;
pub use celestial_time::scales::tt::TT;
pub use celestial_time::scales::ut1::UT1;
pub use celestial_time::scales::utc::UTC;
pub use celestial_time::{TimeError, TimeResult};
