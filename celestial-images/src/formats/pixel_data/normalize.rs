use super::{PixelData, Sample};
use std::ops::Div;

const U32_SPAN: f64 = 4294967295.0;

// The float type a normalized conversion returns.
trait UnitSample: Copy + Div<Output = Self> {
    fn from_f64(value: f64) -> Self;
    fn from_sample<T: Sample>(sample: T) -> Self;
}

impl UnitSample for f32 {
    fn from_f64(value: f64) -> Self {
        value as f32
    }
    fn from_sample<T: Sample>(sample: T) -> Self {
        sample.to_f32()
    }
}

impl UnitSample for f64 {
    fn from_f64(value: f64) -> Self {
        value
    }
    fn from_sample<T: Sample>(sample: T) -> Self {
        sample.to_f64()
    }
}

impl PixelData {
    pub fn to_f32_normalized(&self) -> Vec<f32> {
        self.normalized()
    }

    pub fn to_f64_normalized(&self) -> Vec<f64> {
        self.normalized()
    }

    // Samples scaled to [0, 1]. Integer types up to 32 bits are offset to
    // start at 0 and divided by the width of the type's range; U32 divides
    // in f64 before converting.
    fn normalized<F: UnitSample>(&self) -> Vec<F> {
        match self {
            Self::U8(v) => offset_quotients(v, 0.0, 255.0),
            Self::I8(v) => offset_quotients(v, 128.0, 255.0),
            Self::U16(v) => offset_quotients(v, 0.0, 65535.0),
            Self::I16(v) => offset_quotients(v, 32768.0, 65535.0),
            Self::U32(v) => v
                .iter()
                .map(|&x| F::from_f64(x as f64 / U32_SPAN))
                .collect(),
            Self::I32(v) => offset_quotients(v, 2147483648.0, U32_SPAN),
            Self::U64(v) => data_range_unit(v).into_iter().map(F::from_f64).collect(),
            Self::I64(v) => data_range_unit(v).into_iter().map(F::from_f64).collect(),
            Self::F32(v) => float_unit(v),
            Self::F64(v) => float_unit(v),
        }
    }
}

// (x + offset) / span with the sum taken in f64 and the quotient in `F`.
fn offset_quotients<T: Sample, F: UnitSample>(values: &[T], offset: f64, span: f64) -> Vec<F> {
    let span = F::from_f64(span);
    values
        .iter()
        .map(|&x| F::from_f64(x.to_f64() + offset) / span)
        .collect()
}

// Float data already in [0, 1] converts unchanged; other data is rescaled.
fn float_unit<T: Sample, F: UnitSample>(values: &[T]) -> Vec<F> {
    match unit_rescale(values) {
        None => values.iter().map(|&x| F::from_sample(x)).collect(),
        Some(r) => values
            .iter()
            .map(|&x| F::from_f64(r.apply(x.to_f64())))
            .collect(),
    }
}

// The 64-bit types scale by the data's min and max, because near 2^63 adjacent
// f64 values are 2048 apart and the type range would round ordinary images to
// one value. The offset from the minimum is exact in i128. A constant image
// maps to 0.
fn data_range_unit<T: Copy + Into<i128>>(values: &[T]) -> Vec<f64> {
    let Some((min, max)) = values.iter().map(|&v| v.into()).fold(None, |acc, v| {
        Some(acc.map_or((v, v), |(lo, hi): (i128, i128)| (lo.min(v), hi.max(v))))
    }) else {
        return Vec::new();
    };
    let range = (max - min).max(1) as f64;
    values
        .iter()
        .map(|&v| (v.into() - min) as f64 / range)
        .collect()
}

// The finite min and range of float data, both halved so that max - min stays
// finite for values near f64::MAX.
#[derive(Clone, Copy)]
struct UnitRescale {
    half_min: f64,
    half_range: f64,
}

impl UnitRescale {
    fn apply(self, x: f64) -> f64 {
        (x / 2.0 - self.half_min) / self.half_range
    }
}

// Callers of the normalized conversions expect [0, 1]. Float data with finite
// values outside that range, such as ADU or scaled-integer FITS data, returns
// the finite min and range to rescale by. NaN and infinities do not set them.
// A constant image gets the smallest positive range, so it maps to 0.
fn unit_rescale<T: Sample>(values: &[T]) -> Option<UnitRescale> {
    let (min, max) = values
        .iter()
        .map(|v| v.to_f64())
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (libm::fmin(lo, v), libm::fmax(hi, v))
        });
    if min >= 0.0 && max <= 1.0 {
        return None;
    }
    let half_min = min / 2.0;
    Some(UnitRescale {
        half_min,
        half_range: libm::fmax(max / 2.0 - half_min, f64::MIN_POSITIVE),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_f32_normalized_maps_integer_range_to_unit() {
        let normalized = PixelData::U8(vec![0, 255]).to_f32_normalized();
        assert_eq!(normalized, vec![0.0f32, 1.0]);

        let normalized = PixelData::U16(vec![0, 65535]).to_f32_normalized();
        assert_eq!(normalized, vec![0.0f32, 1.0]);

        let normalized = PixelData::I16(vec![-32768, 32767]).to_f32_normalized();
        assert_eq!(normalized, vec![0.0f32, 1.0]);

        // i32::MAX + 2^31 = 2^32 - 1 rounds to 2^32 in f32, and so does the
        // divisor 4294967295, so the quotient is exactly 1.
        let normalized = PixelData::I32(vec![i32::MIN, i32::MAX]).to_f32_normalized();
        assert_eq!(normalized, vec![0.0f32, 1.0]);
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
        assert_eq!(out, vec![0.0, 1.0]);

        let out = PixelData::I32(vec![i32::MIN, i32::MAX]).to_f64_normalized();
        assert_eq!(out, vec![0.0, 1.0]);

        assert_eq!(PixelData::F32(vec![0.5]).to_f64_normalized(), vec![0.5f64]);
        assert_eq!(
            PixelData::F64(vec![1.5, 2.5]).to_f64_normalized(),
            vec![0.0, 1.0]
        );
    }

    fn bits64(values: &[f64]) -> Vec<u64> {
        values.iter().map(|v| v.to_bits()).collect()
    }

    fn bits32(values: &[f32]) -> Vec<u32> {
        values.iter().map(|v| v.to_bits()).collect()
    }

    #[test]
    fn normalized_conversions_rescale_float_data_outside_unit_range() {
        let expected = [0.0, f64::NAN, 0.5, 1.0];
        let f64s = PixelData::F64(vec![500.0, f64::NAN, 750.0, 1000.0]);
        assert_eq!(bits64(&f64s.to_f64_normalized()), bits64(&expected));
        let f32s = PixelData::F32(vec![-4.0, f32::NAN, -2.0, 0.0]);
        assert_eq!(
            bits32(&f32s.to_f32_normalized()),
            bits32(&[0.0, f32::NAN, 0.5, 1.0])
        );
        assert_eq!(bits64(&f32s.to_f64_normalized()), bits64(&expected));
    }

    #[test]
    fn normalized_conversions_of_extreme_finite_floats_stay_finite() {
        let pd = PixelData::F64(vec![-1e308, 0.0, 1.0, 1e308]);
        assert_eq!(
            bits64(&pd.to_f64_normalized()),
            bits64(&[0.0, 0.5, 0.5, 1.0])
        );
        assert_eq!(
            bits32(&pd.to_f32_normalized()),
            bits32(&[0.0, 0.5, 0.5, 1.0])
        );
    }

    #[test]
    fn normalized_conversions_rescale_64_bit_integers_by_data_range() {
        let base = 1u64 << 63;
        let u64s = PixelData::U64(vec![base + 3, base + 1, base + 5, base + 2]);
        let i64s = PixelData::I64(vec![-2, -4, 0, -3]);
        let expected = [0.5, 0.0, 1.0, 0.25];
        for pd in [u64s, i64s] {
            assert_eq!(bits64(&pd.to_f64_normalized()), bits64(&expected));
            assert_eq!(
                bits32(&pd.to_f32_normalized()),
                bits32(&[0.5, 0.0, 1.0, 0.25])
            );
        }
        let constant = PixelData::I64(vec![i64::MAX; 3]);
        assert_eq!(constant.to_f64_normalized(), vec![0.0; 3]);
        assert_eq!(constant.to_f32_normalized(), vec![0.0f32; 3]);
    }

    #[test]
    fn normalized_conversions_keep_float_data_within_unit_range() {
        let f32s = vec![0.0f32, 0.25, f32::NAN, 1.0, f32::INFINITY];
        let pd = PixelData::F32(f32s.clone());
        assert_eq!(bits32(&pd.to_f32_normalized()), bits32(&f32s));
        let widened: Vec<f64> = f32s.iter().map(|&v| v as f64).collect();
        assert_eq!(bits64(&pd.to_f64_normalized()), bits64(&widened));

        let f64s = vec![0.1, f64::NAN, 0.9, f64::NEG_INFINITY];
        let pd = PixelData::F64(f64s.clone());
        assert_eq!(bits64(&pd.to_f64_normalized()), bits64(&f64s));
    }
}
