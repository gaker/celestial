use super::normalize::finite_angle;
use super::Angle;
use serde::{de, ser, Deserialize, Deserializer, Serialize, Serializer};

impl Serialize for Angle {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let rad = finite_angle(self.radians(), "serialize Angle").map_err(ser::Error::custom)?;
        s.serialize_f64(rad)
    }
}

impl<'de> Deserialize<'de> for Angle {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let rad = f64::deserialize(d)?;
        let rad = finite_angle(rad, "deserialize Angle").map_err(de::Error::custom)?;
        Ok(Self::from_radians(rad))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::PI;
    use serde::de::value::{Error as ValueError, F64Deserializer};
    use serde::de::IntoDeserializer;

    #[test]
    fn test_round_trips_radians_through_json() {
        let json = serde_json::to_string(&Angle::from_radians(PI)).unwrap();
        assert_eq!(json, "3.141592653589793");
        let angle: Angle = serde_json::from_str(&json).unwrap();
        assert_eq!(angle.radians().to_bits(), PI.to_bits());
    }

    #[test]
    fn test_serialize_rejects_non_finite() {
        for rad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = serde_json::to_string(&Angle::from_radians(rad)).unwrap_err();
            assert_eq!(
                err.to_string(),
                "Math error in serialize Angle (not finite): Angle must be finite"
            );
        }
    }

    #[test]
    fn test_deserialize_rejects_non_finite() {
        for rad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let input: F64Deserializer<ValueError> = rad.into_deserializer();
            assert_eq!(
                Angle::deserialize(input).unwrap_err().to_string(),
                "Math error in deserialize Angle (not finite): Angle must be finite"
            );
        }
    }
}
