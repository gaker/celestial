use std::cmp::Ordering;
use std::ops::{Add, Mul};

// Median of `buf`, which is reordered. `None` when `buf` is empty. Callers pass finite
// values only, for which `partial_cmp` is a total order, so the `Equal` fallback is never
// taken. For an even count the median is the mean of the two middle values.
pub(crate) fn median_in_place<T>(buf: &mut [T]) -> Option<T>
where
    T: Copy + PartialOrd + Add<Output = T> + Mul<Output = T> + From<f32>,
{
    let n = buf.len();
    if n == 0 {
        return None;
    }
    let mid = n / 2;
    let (lower, &mut hi, _) =
        buf.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    if n % 2 == 1 {
        return Some(hi);
    }
    let lo = lower
        .iter()
        .copied()
        .reduce(|a, b| if b > a { b } else { a })?;
    Some((lo + hi) * T::from(0.5))
}

// Median absolute deviation of `data` about `median`. 0 for empty `data`.
pub(super) fn mad(data: &[f64], median: f64) -> f64 {
    let mut deviations: Vec<f64> = data.iter().map(|&v| libm::fabs(v - median)).collect();
    median_in_place(&mut deviations).unwrap_or(0.0)
}

// Sample standard deviation (n - 1 denominator). 0 for fewer than 2 values.
pub(super) fn stddev(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let mean = data.iter().sum::<f64>() / n;
    let var = data.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
    libm::sqrt(var)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_of_odd_count() {
        assert_eq!(median_in_place(&mut [5.0_f64, 1.0, 3.0]), Some(3.0));
    }

    #[test]
    fn median_of_even_count() {
        assert_eq!(median_in_place(&mut [4.0_f64, 1.0, 3.0, 2.0]), Some(2.5));
        assert_eq!(median_in_place(&mut [4.0_f32, 1.0, 3.0, 2.0]), Some(2.5));
    }

    #[test]
    fn mad_and_stddev() {
        assert_eq!(mad(&[1.0, 2.0, 4.0], 2.0), 1.0);
        assert_eq!(mad(&[], 2.0), 0.0);
        assert_eq!(stddev(&[1.0, 3.0]), libm::sqrt(2.0));
        assert_eq!(stddev(&[1.0]), 0.0);
    }

    #[test]
    fn median_of_empty_is_none() {
        assert_eq!(median_in_place::<f64>(&mut []), None);
    }
}
