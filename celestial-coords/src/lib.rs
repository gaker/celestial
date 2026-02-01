pub mod aberration;
pub(crate) mod constants;
pub mod distance;
pub mod eop;
pub mod errors;
pub mod frames;
pub mod lighttime;
pub mod lunar;
pub mod solar;
pub mod transforms;

pub use celestial_core::Angle;
pub use distance::Distance;
pub use eop::{EopManager, EopParameters, EopRecord};
pub use errors::{CoordError, CoordResult};
pub use lighttime::LightTimeCorrection;

pub use frames::{
    CIRSPosition, EclipticCartesian, EclipticPosition, GCRSPosition, GalacticPosition,
    HeliographicCarrington, HeliographicStonyhurst, HourAnglePosition, ICRSPosition, ITRSPosition,
    SelenographicPosition, TIRSPosition, TopocentricPosition,
};

pub use transforms::{CartesianFrame, CoordinateFrame};

pub use celestial_core::{Location, Vector3};
pub use celestial_time::{TimeError, TimeResult, TAI, TT, UT1, UTC};
