use celestial_core::constants::HALF_PI;

use super::*;

fn check(x: [f64; LANES]) {
    let (sin, cos) = sincos(x);
    for i in 0..LANES {
        let (s, c) = libm::sincos(x[i]);
        assert_eq!(
            (sin[i].to_bits(), cos[i].to_bits()),
            (s.to_bits(), c.to_bits()),
            "sincos({:e}), bits {:016x}",
            x[i],
            x[i].to_bits()
        );
    }
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn signed(&mut self, x: f64) -> f64 {
        if self.next() & 1 == 0 {
            x
        } else {
            -x
        }
    }
}

#[test]
fn magnitudes_from_tiny_to_past_the_medium_range() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..100_000 {
        check([0; LANES].map(|_| {
            let x = libm::exp2(rng.unit() * 90.0 - 66.0);
            rng.signed(x)
        }));
    }
}

#[test]
fn arguments_up_to_the_moon_series_reach() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for _ in 0..100_000 {
        check([0; LANES].map(|_| {
            let x = rng.unit() * 4.0e6;
            rng.signed(x)
        }));
    }
}

// Arguments a few ulps from a multiple of π/2 lose the most bits when it is
// subtracted, which is what the second and third rounds are for.
#[test]
fn arguments_near_multiples_of_half_pi() {
    let mut rng = Rng(0x1234_5678_9abc_def1);
    for _ in 0..100_000 {
        check([0; LANES].map(|_| {
            let k = (rng.next() % 2_000_000) as f64;
            let ulps = (rng.next() % 4096) as i64 - 2048;
            let x = f64::from_bits(((k * HALF_PI).to_bits() as i64 + ulps) as u64);
            rng.signed(x)
        }));
    }
}

#[test]
fn either_side_of_every_branch() {
    let words = [
        0x3e46_a09e,
        0x3fe9_21fb,
        0x3ff9_21fb,
        0x4002_d97c,
        0x4009_21fb,
        0x400f_6a7a,
        0x4012_d97c,
        0x4015_fdbc,
        0x4019_21fb,
        0x401c_463b,
        0x4139_21fb,
    ];
    for word in words {
        for start in [word << 32, (word + 1) << 32] {
            for offset in (0..4096).step_by(LANES) {
                let bits = |i: usize| start - 2048 + offset + i as u64;
                check(std::array::from_fn(|i| f64::from_bits(bits(i))));
                check(std::array::from_fn(|i| -f64::from_bits(bits(i))));
            }
        }
    }
}

#[test]
fn special_values() {
    check([
        0.0,
        -0.0,
        5e-324,
        -1e-310,
        f64::MIN_POSITIVE,
        1e-9,
        0.5,
        -1.0,
    ]);
    check([
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        1e300,
        -3.7e6,
        1.7e6,
        1.6e6,
        2.0,
    ]);
}

// The doubles nearest multiples of π/2 lose the most bits to cancellation, and
// some of them need libm's third round.
#[test]
fn multiples_of_half_pi_through_the_medium_range() {
    for k in (1..1 << 20).step_by(LANES / 2) {
        check(std::array::from_fn(|i| {
            let x = (k + i / 2) as f64 * HALF_PI;
            if i % 2 == 0 {
                x
            } else {
                -x
            }
        }));
    }
}

// Arguments that lose exactly 50 bits of x's exponent in the second round, the
// fewest that make libm run a third.
#[test]
fn arguments_at_the_edge_of_the_third_round() {
    let x = [
        0x4113_75a8_67ad_0070,
        0x411a_817f_1bf0_40a8,
        0x411e_65c6_2ac8_9f31,
        0x4122_5a0b_0ea2_7597,
        0x4132_6c54_07bc_843d,
        0x4138_7ccd_7df6_daff,
    ]
    .map(f64::from_bits);
    for sign in [1.0, -1.0] {
        check(std::array::from_fn(|i| sign * x[i % x.len()]));
    }
}
