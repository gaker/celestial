pub mod cartesian;

use crate::errors::CoordResult;
use crate::frames::icrs::ICRSPosition;
use celestial_time::scales::tt::TT;

pub trait CoordinateFrame: Sized {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition>;

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self>;
}
