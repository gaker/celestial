use crate::core::BitPix;

#[derive(Debug, Clone)]
pub enum PixelData {
    U8(Vec<u8>),
    U16(Vec<u16>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    F32(Vec<f32>),
    F64(Vec<f64>),
}

impl PixelData {
    pub fn len(&self) -> usize {
        match self {
            Self::U8(v) => v.len(),
            Self::U16(v) => v.len(),
            Self::I16(v) => v.len(),
            Self::I32(v) => v.len(),
            Self::F32(v) => v.len(),
            Self::F64(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn bitpix(&self) -> BitPix {
        match self {
            Self::U8(_) => BitPix::U8,
            Self::U16(_) => BitPix::I16,
            Self::I16(_) => BitPix::I16,
            Self::I32(_) => BitPix::I32,
            Self::F32(_) => BitPix::F32,
            Self::F64(_) => BitPix::F64,
        }
    }

    pub fn as_u8(&self) -> Option<&Vec<u8>> {
        match self {
            Self::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u8_mut(&mut self) -> Option<&mut Vec<u8>> {
        match self {
            Self::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u16(&self) -> Option<&Vec<u16>> {
        match self {
            Self::U16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u16_mut(&mut self) -> Option<&mut Vec<u16>> {
        match self {
            Self::U16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i16(&self) -> Option<&Vec<i16>> {
        match self {
            Self::I16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i16_mut(&mut self) -> Option<&mut Vec<i16>> {
        match self {
            Self::I16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i32(&self) -> Option<&Vec<i32>> {
        match self {
            Self::I32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i32_mut(&mut self) -> Option<&mut Vec<i32>> {
        match self {
            Self::I32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f32(&self) -> Option<&Vec<f32>> {
        match self {
            Self::F32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f32_mut(&mut self) -> Option<&mut Vec<f32>> {
        match self {
            Self::F32(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<&Vec<f64>> {
        match self {
            Self::F64(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f64_mut(&mut self) -> Option<&mut Vec<f64>> {
        match self {
            Self::F64(v) => Some(v),
            _ => None,
        }
    }

    pub fn to_f32(&self) -> Vec<f32> {
        match self {
            Self::U8(v) => v.iter().map(|&x| x as f32).collect(),
            Self::U16(v) => v.iter().map(|&x| x as f32).collect(),
            Self::I16(v) => v.iter().map(|&x| x as f32).collect(),
            Self::I32(v) => v.iter().map(|&x| x as f32).collect(),
            Self::F32(v) => v.clone(),
            Self::F64(v) => v.iter().map(|&x| x as f32).collect(),
        }
    }

    pub fn to_f32_normalized(&self) -> Vec<f32> {
        match self {
            Self::U8(v) => v.iter().map(|&x| x as f32 / 255.0).collect(),
            Self::U16(v) => v.iter().map(|&x| x as f32 / 65535.0).collect(),
            Self::I16(v) => v.iter().map(|&x| (x as f32 + 32768.0) / 65535.0).collect(),
            Self::I32(v) => v
                .iter()
                .map(|&x| (x as f64 + 2147483648.0) as f32 / 4294967295.0)
                .collect(),
            Self::F32(v) => v.clone(),
            Self::F64(v) => v.iter().map(|&x| x as f32).collect(),
        }
    }

    pub fn to_f64_normalized(&self) -> Vec<f64> {
        match self {
            Self::U8(v) => v.iter().map(|&x| x as f64 / 255.0).collect(),
            Self::U16(v) => v.iter().map(|&x| x as f64 / 65535.0).collect(),
            Self::I16(v) => v.iter().map(|&x| (x as f64 + 32768.0) / 65535.0).collect(),
            Self::I32(v) => v
                .iter()
                .map(|&x| (x as f64 + 2147483648.0) / 4294967295.0)
                .collect(),
            Self::F32(v) => v.iter().map(|&x| x as f64).collect(),
            Self::F64(v) => v.clone(),
        }
    }

    pub fn convert_to_f32(&mut self) {
        let converted = self.to_f32_normalized();
        *self = Self::F32(converted);
    }
}

impl From<&[u8]> for PixelData {
    fn from(data: &[u8]) -> Self {
        Self::U8(data.to_vec())
    }
}

impl From<&[u16]> for PixelData {
    fn from(data: &[u16]) -> Self {
        Self::U16(data.to_vec())
    }
}

impl From<&[i16]> for PixelData {
    fn from(data: &[i16]) -> Self {
        Self::I16(data.to_vec())
    }
}

impl From<&[i32]> for PixelData {
    fn from(data: &[i32]) -> Self {
        Self::I32(data.to_vec())
    }
}

impl From<&[f32]> for PixelData {
    fn from(data: &[f32]) -> Self {
        Self::F32(data.to_vec())
    }
}

impl From<&[f64]> for PixelData {
    fn from(data: &[f64]) -> Self {
        Self::F64(data.to_vec())
    }
}

impl From<Vec<u8>> for PixelData {
    fn from(data: Vec<u8>) -> Self {
        Self::U8(data)
    }
}

impl From<Vec<u16>> for PixelData {
    fn from(data: Vec<u16>) -> Self {
        Self::U16(data)
    }
}

impl From<Vec<i16>> for PixelData {
    fn from(data: Vec<i16>) -> Self {
        Self::I16(data)
    }
}

impl From<Vec<i32>> for PixelData {
    fn from(data: Vec<i32>) -> Self {
        Self::I32(data)
    }
}

impl From<Vec<f32>> for PixelData {
    fn from(data: Vec<f32>) -> Self {
        Self::F32(data)
    }
}

impl From<Vec<f64>> for PixelData {
    fn from(data: Vec<f64>) -> Self {
        Self::F64(data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn len_and_is_empty_across_variants() {
        assert_eq!(PixelData::U8(vec![]).len(), 0);
        assert!(PixelData::U8(vec![]).is_empty());

        assert_eq!(PixelData::U16(vec![1, 2, 3]).len(), 3);
        assert!(!PixelData::U16(vec![1, 2, 3]).is_empty());

        assert_eq!(PixelData::I16(vec![1; 10]).len(), 10);
        assert_eq!(PixelData::I32(vec![1; 5]).len(), 5);
        assert_eq!(PixelData::F32(vec![1.0; 7]).len(), 7);
        assert_eq!(PixelData::F64(vec![1.0; 2]).len(), 2);
    }

    #[test]
    fn bitpix_maps_variants_correctly() {
        assert_eq!(PixelData::U8(vec![]).bitpix(), BitPix::U8);
        assert_eq!(PixelData::U16(vec![]).bitpix(), BitPix::I16);
        assert_eq!(PixelData::I16(vec![]).bitpix(), BitPix::I16);
        assert_eq!(PixelData::I32(vec![]).bitpix(), BitPix::I32);
        assert_eq!(PixelData::F32(vec![]).bitpix(), BitPix::F32);
        assert_eq!(PixelData::F64(vec![]).bitpix(), BitPix::F64);
    }

    #[test]
    fn as_ref_accessors_return_some_only_for_matching_variant() {
        let pd = PixelData::U8(vec![1, 2, 3]);
        assert_eq!(pd.as_u8(), Some(&vec![1, 2, 3]));
        assert!(pd.as_u16().is_none());
        assert!(pd.as_i16().is_none());
        assert!(pd.as_i32().is_none());
        assert!(pd.as_f32().is_none());
        assert!(pd.as_f64().is_none());

        let pd = PixelData::U16(vec![100, 200]);
        assert_eq!(pd.as_u16(), Some(&vec![100u16, 200]));
        assert!(pd.as_u8().is_none());

        let pd = PixelData::I16(vec![-5, 5]);
        assert_eq!(pd.as_i16(), Some(&vec![-5i16, 5]));

        let pd = PixelData::I32(vec![1000]);
        assert_eq!(pd.as_i32(), Some(&vec![1000i32]));

        let pd = PixelData::F32(vec![1.5]);
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let pd = PixelData::F64(vec![2.5]);
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }

    #[test]
    fn as_mut_accessors_allow_modification() {
        let mut pd = PixelData::U8(vec![1, 2, 3]);
        pd.as_u8_mut().unwrap().push(4);
        assert_eq!(pd.as_u8(), Some(&vec![1u8, 2, 3, 4]));
        assert!(pd.as_u16_mut().is_none());

        let mut pd = PixelData::U16(vec![1]);
        pd.as_u16_mut().unwrap().push(2);
        assert_eq!(pd.as_u16(), Some(&vec![1u16, 2]));
        assert!(pd.as_u8_mut().is_none());

        let mut pd = PixelData::I16(vec![]);
        pd.as_i16_mut().unwrap().push(-1);
        assert_eq!(pd.as_i16(), Some(&vec![-1i16]));

        let mut pd = PixelData::I32(vec![]);
        pd.as_i32_mut().unwrap().push(100);
        assert_eq!(pd.as_i32(), Some(&vec![100i32]));

        let mut pd = PixelData::F32(vec![]);
        pd.as_f32_mut().unwrap().push(1.5);
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let mut pd = PixelData::F64(vec![]);
        pd.as_f64_mut().unwrap().push(2.5);
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }

    #[test]
    fn to_f32_from_every_variant() {
        assert_eq!(PixelData::U8(vec![0, 128, 255]).to_f32(), vec![0.0, 128.0, 255.0]);
        assert_eq!(PixelData::U16(vec![0, 65535]).to_f32(), vec![0.0, 65535.0]);
        assert_eq!(PixelData::I16(vec![-1, 1]).to_f32(), vec![-1.0, 1.0]);
        assert_eq!(PixelData::I32(vec![-42, 42]).to_f32(), vec![-42.0, 42.0]);
        assert_eq!(PixelData::F32(vec![1.5, 2.5]).to_f32(), vec![1.5, 2.5]);
        assert_eq!(PixelData::F64(vec![1.5, 2.5]).to_f32(), vec![1.5f32, 2.5]);
    }

    #[test]
    fn to_f32_normalized_maps_integer_range_to_unit() {
        let normalized = PixelData::U8(vec![0, 255]).to_f32_normalized();
        assert_eq!(normalized, vec![0.0f32, 1.0]);

        let normalized = PixelData::U16(vec![0, 65535]).to_f32_normalized();
        assert_eq!(normalized, vec![0.0f32, 1.0]);

        let normalized = PixelData::I16(vec![-32768, 32767]).to_f32_normalized();
        assert_eq!(normalized[0], 0.0);
        assert!((normalized[1] - 1.0).abs() < 1e-5);

        let normalized = PixelData::I32(vec![i32::MIN, i32::MAX]).to_f32_normalized();
        assert!(normalized[0] >= 0.0 && normalized[0] < 1e-6);
        assert!((normalized[1] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn to_f32_normalized_passes_through_floats() {
        let data = vec![0.0f32, 0.5, 1.0];
        assert_eq!(PixelData::F32(data.clone()).to_f32_normalized(), data);

        let data = vec![0.0f64, 0.5, 1.0];
        let out = PixelData::F64(data).to_f32_normalized();
        assert_eq!(out, vec![0.0f32, 0.5, 1.0]);
    }

    #[test]
    fn to_f64_normalized_behaves_like_f32_variant_but_preserves_precision() {
        assert_eq!(
            PixelData::U8(vec![0, 255]).to_f64_normalized(),
            vec![0.0, 1.0]
        );
        assert_eq!(
            PixelData::U16(vec![0, 65535]).to_f64_normalized(),
            vec![0.0, 1.0]
        );
        let out = PixelData::I16(vec![-32768, 32767]).to_f64_normalized();
        assert_eq!(out[0], 0.0);
        assert!((out[1] - 1.0).abs() < 1e-10);

        let out = PixelData::I32(vec![i32::MIN, i32::MAX]).to_f64_normalized();
        assert!(out[0] >= 0.0 && out[0] < 1e-9);
        assert!((out[1] - 1.0).abs() < 1e-9);

        assert_eq!(
            PixelData::F32(vec![0.5]).to_f64_normalized(),
            vec![0.5f64]
        );
        assert_eq!(
            PixelData::F64(vec![1.5, 2.5]).to_f64_normalized(),
            vec![1.5, 2.5]
        );
    }

    #[test]
    fn convert_to_f32_replaces_variant() {
        let mut pd = PixelData::U8(vec![0, 255]);
        pd.convert_to_f32();
        assert!(matches!(pd, PixelData::F32(_)));
        assert_eq!(pd.as_f32().unwrap(), &vec![0.0f32, 1.0]);

        let mut pd = PixelData::F64(vec![0.25]);
        pd.convert_to_f32();
        assert!(matches!(pd, PixelData::F32(_)));
    }

    #[test]
    fn from_slice_conversions() {
        let pd: PixelData = (&[1u8, 2, 3][..]).into();
        assert_eq!(pd.as_u8(), Some(&vec![1u8, 2, 3]));

        let pd: PixelData = (&[100u16, 200][..]).into();
        assert_eq!(pd.as_u16(), Some(&vec![100u16, 200]));

        let pd: PixelData = (&[-1i16, 1][..]).into();
        assert_eq!(pd.as_i16(), Some(&vec![-1i16, 1]));

        let pd: PixelData = (&[42i32][..]).into();
        assert_eq!(pd.as_i32(), Some(&vec![42i32]));

        let pd: PixelData = (&[1.5f32][..]).into();
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let pd: PixelData = (&[2.5f64][..]).into();
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }

    #[test]
    fn from_vec_conversions() {
        let pd: PixelData = vec![1u8, 2].into();
        assert_eq!(pd.as_u8(), Some(&vec![1u8, 2]));

        let pd: PixelData = vec![100u16].into();
        assert_eq!(pd.as_u16(), Some(&vec![100u16]));

        let pd: PixelData = vec![-1i16].into();
        assert_eq!(pd.as_i16(), Some(&vec![-1i16]));

        let pd: PixelData = vec![42i32].into();
        assert_eq!(pd.as_i32(), Some(&vec![42i32]));

        let pd: PixelData = vec![1.5f32].into();
        assert_eq!(pd.as_f32(), Some(&vec![1.5f32]));

        let pd: PixelData = vec![2.5f64].into();
        assert_eq!(pd.as_f64(), Some(&vec![2.5f64]));
    }
}
