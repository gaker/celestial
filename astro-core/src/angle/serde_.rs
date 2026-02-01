use super::Angle;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

impl Serialize for Angle {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_f64(self.radians())
    }
}

impl<'de> Deserialize<'de> for Angle {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let r = f64::deserialize(d)?;
        Ok(Angle::from_radians(r))
    }
}
