pub mod cartesian;

use crate::{frames::ICRSPosition, CoordResult};
use celestial_time::TT;

pub use cartesian::CartesianFrame;

pub trait CoordinateFrame: Sized {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition>;

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self>;
}
