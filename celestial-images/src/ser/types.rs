use crate::ser::error::SerError;
use crate::ser::timestamp::SerTimestamp;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum ColorId {
    Mono = 0,
    BayerRggb = 8,
    BayerGrbg = 9,
    BayerGbrg = 10,
    BayerBggr = 11,
    BayerCyym = 16,
    BayerYcmy = 17,
    BayerYmcy = 18,
    BayerMyyc = 19,
    Rgb = 100,
    Bgr = 101,
}

impl ColorId {
    pub(crate) fn from_u32(value: u32) -> Result<Self, SerError> {
        match value {
            0 => Ok(Self::Mono),
            8 => Ok(Self::BayerRggb),
            9 => Ok(Self::BayerGrbg),
            10 => Ok(Self::BayerGbrg),
            11 => Ok(Self::BayerBggr),
            16 => Ok(Self::BayerCyym),
            17 => Ok(Self::BayerYcmy),
            18 => Ok(Self::BayerYmcy),
            19 => Ok(Self::BayerMyyc),
            100 => Ok(Self::Rgb),
            101 => Ok(Self::Bgr),
            _ => Err(SerError::UnsupportedColorFormat(value)),
        }
    }

    pub fn planes(self) -> u32 {
        match self {
            Self::Rgb | Self::Bgr => 3,
            _ => 1,
        }
    }

    pub fn is_bayer(self) -> bool {
        matches!(
            self,
            Self::BayerRggb
                | Self::BayerGrbg
                | Self::BayerGbrg
                | Self::BayerBggr
                | Self::BayerCyym
                | Self::BayerYcmy
                | Self::BayerYmcy
                | Self::BayerMyyc
        )
    }
}

impl fmt::Display for ColorId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Mono => "MONO",
            Self::BayerRggb => "BAYER_RGGB",
            Self::BayerGrbg => "BAYER_GRBG",
            Self::BayerGbrg => "BAYER_GBRG",
            Self::BayerBggr => "BAYER_BGGR",
            Self::BayerCyym => "BAYER_CYYM",
            Self::BayerYcmy => "BAYER_YCMY",
            Self::BayerYmcy => "BAYER_YMCY",
            Self::BayerMyyc => "BAYER_MYYC",
            Self::Rgb => "RGB",
            Self::Bgr => "BGR",
        };
        write!(f, "{}", name)
    }
}

pub type SerFrameId = u32;

#[derive(Debug)]
pub struct SerFrame<'a> {
    pub id: SerFrameId,
    pub data: &'a [u8],
    pub timestamp: Option<SerTimestamp>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn color_id_from_u32_all_variants() {
        assert_eq!(ColorId::from_u32(0).unwrap(), ColorId::Mono);
        assert_eq!(ColorId::from_u32(8).unwrap(), ColorId::BayerRggb);
        assert_eq!(ColorId::from_u32(9).unwrap(), ColorId::BayerGrbg);
        assert_eq!(ColorId::from_u32(10).unwrap(), ColorId::BayerGbrg);
        assert_eq!(ColorId::from_u32(11).unwrap(), ColorId::BayerBggr);
        assert_eq!(ColorId::from_u32(16).unwrap(), ColorId::BayerCyym);
        assert_eq!(ColorId::from_u32(17).unwrap(), ColorId::BayerYcmy);
        assert_eq!(ColorId::from_u32(18).unwrap(), ColorId::BayerYmcy);
        assert_eq!(ColorId::from_u32(19).unwrap(), ColorId::BayerMyyc);
        assert_eq!(ColorId::from_u32(100).unwrap(), ColorId::Rgb);
        assert_eq!(ColorId::from_u32(101).unwrap(), ColorId::Bgr);

        assert!(matches!(
            ColorId::from_u32(999),
            Err(SerError::UnsupportedColorFormat(999))
        ));
    }

    #[test]
    fn color_id_planes() {
        assert_eq!(ColorId::Rgb.planes(), 3);
        assert_eq!(ColorId::Bgr.planes(), 3);
        assert_eq!(ColorId::Mono.planes(), 1);
        assert_eq!(ColorId::BayerRggb.planes(), 1);
        assert_eq!(ColorId::BayerGrbg.planes(), 1);
        assert_eq!(ColorId::BayerGbrg.planes(), 1);
        assert_eq!(ColorId::BayerBggr.planes(), 1);
        assert_eq!(ColorId::BayerCyym.planes(), 1);
        assert_eq!(ColorId::BayerYcmy.planes(), 1);
        assert_eq!(ColorId::BayerYmcy.planes(), 1);
        assert_eq!(ColorId::BayerMyyc.planes(), 1);
    }

    #[test]
    fn color_id_is_bayer() {
        assert!(!ColorId::Mono.is_bayer());
        assert!(!ColorId::Rgb.is_bayer());
        assert!(!ColorId::Bgr.is_bayer());

        assert!(ColorId::BayerRggb.is_bayer());
        assert!(ColorId::BayerGrbg.is_bayer());
        assert!(ColorId::BayerGbrg.is_bayer());
        assert!(ColorId::BayerBggr.is_bayer());
        assert!(ColorId::BayerCyym.is_bayer());
        assert!(ColorId::BayerYcmy.is_bayer());
        assert!(ColorId::BayerYmcy.is_bayer());
        assert!(ColorId::BayerMyyc.is_bayer());
    }

    #[test]
    fn color_id_display_all_variants() {
        assert_eq!(ColorId::Mono.to_string(), "MONO");
        assert_eq!(ColorId::BayerRggb.to_string(), "BAYER_RGGB");
        assert_eq!(ColorId::BayerGrbg.to_string(), "BAYER_GRBG");
        assert_eq!(ColorId::BayerGbrg.to_string(), "BAYER_GBRG");
        assert_eq!(ColorId::BayerBggr.to_string(), "BAYER_BGGR");
        assert_eq!(ColorId::BayerCyym.to_string(), "BAYER_CYYM");
        assert_eq!(ColorId::BayerYcmy.to_string(), "BAYER_YCMY");
        assert_eq!(ColorId::BayerYmcy.to_string(), "BAYER_YMCY");
        assert_eq!(ColorId::BayerMyyc.to_string(), "BAYER_MYYC");
        assert_eq!(ColorId::Rgb.to_string(), "RGB");
        assert_eq!(ColorId::Bgr.to_string(), "BGR");
    }

    #[test]
    fn ser_frame_fields() {
        let data = [1, 2, 3, 4, 5];
        let ts = SerTimestamp::new(1000);
        let frame = SerFrame {
            id: 42,
            data: &data,
            timestamp: Some(ts),
        };
        assert_eq!(frame.id, 42);
        assert_eq!(frame.data, &data);
        assert_eq!(frame.timestamp, Some(ts));
    }

    #[test]
    fn color_id_and_clone() {
        let color = ColorId::BayerRggb;
        let cloned = color;
        assert_eq!(color, cloned);

        let debug_str = format!("{:?}", color);
        assert!(debug_str.contains("BayerRggb"));
    }
}
