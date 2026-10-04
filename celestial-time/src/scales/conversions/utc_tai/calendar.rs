use crate::{TimeError, TimeResult};

pub(super) fn check_julian_range(dj: f64) -> TimeResult<()> {
    const DJMIN: f64 = -68569.5;
    const DJMAX: f64 = 1e9;

    if !(DJMIN..=DJMAX).contains(&dj) {
        return Err(TimeError::ConversionError(format!(
            "Julian Date {} out of valid range [{}, {}]",
            dj, DJMIN, DJMAX
        )));
    }
    Ok(())
}

fn nearest_int(a: f64) -> f64 {
    if libm::fabs(a) < 0.5 {
        0.0
    } else if a < 0.0 {
        libm::ceil(a - 0.5)
    } else {
        libm::floor(a + 0.5)
    }
}

// The Julian Day number of the civil day and the time since its midnight.
pub(super) fn day_number_and_fraction(jd1: f64, jd2: f64) -> (i64, f64) {
    let day_int_1 = nearest_int(jd1);
    let day_int_2 = nearest_int(jd2);
    let (carry, kahan) = sum_fractions([jd1 - day_int_1, jd2 - day_int_2]);
    let (adjust, fraction) = normalize_fraction(kahan);
    (
        day_int_1 as i64 + day_int_2 as i64 + carry + adjust,
        fraction,
    )
}

struct Kahan {
    sum: f64,
    correction: f64,
}

// Starting the sum at 0.5 moves the day boundary from noon to midnight.
fn sum_fractions(fractions: [f64; 2]) -> (i64, Kahan) {
    let mut carry = 0;
    let mut sum = 0.5;
    let mut correction = 0.0;
    for frac in fractions {
        let temp = sum + frac;
        correction += if libm::fabs(sum) >= libm::fabs(frac) {
            (sum - temp) + frac
        } else {
            (frac - temp) + sum
        };
        sum = temp;
        if sum >= 1.0 {
            carry += 1;
            sum -= 1.0;
        }
    }
    (carry, Kahan { sum, correction })
}

// Brings the compensated fraction into [0, 1), returning the days moved.
fn normalize_fraction(kahan: Kahan) -> (i64, f64) {
    let fraction = kahan.sum + kahan.correction;
    let kahan = Kahan {
        sum: kahan.sum,
        correction: fraction - kahan.sum,
    };
    let (borrow, kahan, fraction) = if fraction < 0.0 {
        borrow_day(kahan)
    } else {
        (0, kahan, fraction)
    };
    let (carry, fraction) = carry_day(kahan, fraction);
    (carry - borrow, fraction)
}

fn borrow_day(kahan: Kahan) -> (i64, Kahan, f64) {
    let sum = kahan.sum + 1.0;
    let correction = kahan.correction + ((1.0 - sum) + kahan.sum);
    let fraction = sum + correction;
    let kahan = Kahan {
        sum,
        correction: fraction - sum,
    };
    (1, kahan, fraction)
}

// A fraction within a few ulp below 1.0 rounds up into the next day.
fn carry_day(kahan: Kahan, fraction: f64) -> (i64, f64) {
    const DBL_EPSILON: f64 = 2.220446049250313e-16;
    if (fraction - 1.0) < -DBL_EPSILON / 4.0 {
        return (0, fraction);
    }
    let sum = kahan.sum - 1.0;
    let correction = kahan.correction + ((kahan.sum - sum) - 1.0);
    let fraction = sum + correction;
    if (-DBL_EPSILON / 2.0) < fraction {
        (1, if fraction > 0.0 { fraction } else { 0.0 })
    } else {
        (0, fraction)
    }
}

pub(super) fn gregorian_date(jd: i64) -> (i32, i32, i32) {
    let mut l = jd + 68569_i64;
    let n = (4_i64 * l) / 146097_i64;
    l -= (146097_i64 * n + 3_i64) / 4_i64;
    let i = (4000_i64 * (l + 1_i64)) / 1461001_i64;
    l -= (1461_i64 * i) / 4_i64 - 31_i64;
    let k = (80_i64 * l) / 2447_i64;
    let day = (l - (2447_i64 * k) / 80_i64) as i32;
    let l_final = k / 11_i64;
    let month = (k + 2_i64 - 12_i64 * l_final) as i32;
    let year = (100_i64 * (n - 49_i64) + i + l_final) as i32;
    (year, month, day)
}
