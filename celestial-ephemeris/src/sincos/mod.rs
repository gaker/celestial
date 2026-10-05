// The sines and cosines of eight arguments at once. This is libm's sincos,
// rem_pio2, k_sin and k_cos worked lane by lane: each lane runs every branch
// and keeps the result of the one libm takes for its argument, so the bits
// match libm::sincos. Arguments past 2²⁰·π/2, infinities and NaNs still go to
// libm.
//
// Those libm files are ports of FreeBSD's msun, which comes from Sun's
// fdlibm, and carry these notices:
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
//
// Optimized by Bruce D. Evans.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
//
// libm is Copyright (c) 2018 Jorge Aparicio, and the musl libc it is ported
// from is Copyright © 2005-2020 Rich Felker, et al., both under this license:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

use std::array::from_fn;

use wide::{f64x8, u64x8};

#[cfg(test)]
mod tests;

const LANES: usize = 8;

// Every constant is a whole vector: a scalar operand would be splatted at run
// time, through memory.
const fn splat(bits: u64) -> f64x8 {
    f64x8::splat(f64::from_bits(bits))
}

const S1: f64x8 = splat(0xbfc5_5555_5555_5549);
const S2: f64x8 = splat(0x3f81_1111_1110_f8a6);
const S3: f64x8 = splat(0xbf2a_01a0_19c1_61d5);
const S4: f64x8 = splat(0x3ec7_1de3_57b1_fe7d);
const S5: f64x8 = splat(0xbe5a_e5e6_8a2b_9ceb);
const S6: f64x8 = splat(0x3de5_d93a_5acf_d57c);
const C1: f64x8 = splat(0x3fa5_5555_5555_554c);
const C2: f64x8 = splat(0xbf56_c16c_16c1_5177);
const C3: f64x8 = splat(0x3efa_01a0_19cb_1590);
const C4: f64x8 = splat(0xbe92_7e4f_809c_52ad);
const C5: f64x8 = splat(0x3e21_ee9e_bdb4_b1c4);
const C6: f64x8 = splat(0xbda8_fae9_be88_38d4);
const INV_PIO2: f64x8 = splat(0x3fe4_5f30_6dc9_c883);
const PIO2_1: f64x8 = splat(0x3ff9_21fb_5440_0000);
const PIO2_1T: f64x8 = splat(0x3dd0_b461_1a62_6331);
const PIO2_2: f64x8 = splat(0x3dd0_b461_1a60_0000);
const PIO2_2T: f64x8 = splat(0x3ba3_198a_2e03_7073);
const PIO2_3: f64x8 = splat(0x3ba3_198a_2e00_0000);
const PIO2_3T: f64x8 = splat(0x397b_839a_2520_49c1);
const TO_INT: f64x8 = splat(0x4338_0000_0000_0000);
const HALF: f64x8 = f64x8::splat(0.5);
const EXPONENT: f64x8 = f64x8::splat(f64::INFINITY);
const SECOND_ROUND: f64x8 = splat(0x3ef0_0000_0000_0000);
const THIRD_ROUND: f64x8 = splat(0x3ce0_0000_0000_0000);

// libm branches on the high word of |x|. Ordering non-negative doubles by
// their bits orders them by value, so "high word < w" is "|x| < starting(w)",
// and each branch here is a compare on |x|.
const fn starting(word: u64) -> f64x8 {
    splat(word << 32)
}

const TINY: f64x8 = starting(0x3e46_a09e);
const UNREDUCED_END: f64x8 = starting(0x3fe9_21fc);
const STEPS: [f64x8; 4] = [
    starting(0x3fe9_21fc),
    starting(0x4002_d97d),
    starting(0x400f_6a7b),
    starting(0x4015_fdbd),
];
const STEPPED_END: f64x8 = starting(0x401c_463c);
const MEDIUM_END: f64x8 = starting(0x4139_21fb);

// High words close to multiples of π/2, which libm moves from the steps of
// π/2 to rounding x·2/π.
const NEAR: [(f64x8, f64x8); 4] = [
    word(0x3ff9_21fb),
    word(0x4009_21fb),
    word(0x4012_d97c),
    word(0x4019_21fb),
];

const fn word(word: u64) -> (f64x8, f64x8) {
    (starting(word), starting(word + 1))
}

// Calls `add` with each term, its argument's rate, and the sine and cosine of
// its argument, in order. The sines and cosines are worked out LANES terms at
// a time; the padding lanes past the last term are never read.
pub(crate) fn for_each<T>(
    terms: &[T],
    argument: impl Fn(&T) -> (f64, f64),
    mut add: impl FnMut(&T, f64, f64, f64),
) {
    for block in terms.chunks(LANES) {
        let args: [(f64, f64); LANES] = from_fn(|i| block.get(i).map_or((0.0, 0.0), &argument));
        let (sin, cos) = sincos(args.map(|(arg, _)| arg));
        for (i, term) in block.iter().enumerate() {
            add(term, args[i].1, sin[i], cos[i]);
        }
    }
}

fn sincos(x: [f64; LANES]) -> ([f64; LANES], [f64; LANES]) {
    let (sin, cos, unhandled) = lanes(f64x8::new(x));
    let (mut sin, mut cos) = (sin.to_array(), cos.to_array());
    if unhandled.any() {
        fall_back(&x, unhandled, &mut sin, &mut cos);
    }
    (sin, cos)
}

#[cold]
fn fall_back(x: &[f64; LANES], unhandled: f64x8, sin: &mut [f64; LANES], cos: &mut [f64; LANES]) {
    for (i, bits) in unhandled.to_bits().to_array().into_iter().enumerate() {
        if bits != 0 {
            (sin[i], cos[i]) = libm::sincos(x[i]);
        }
    }
}

#[inline(always)]
fn lanes(x: f64x8) -> (f64x8, f64x8, f64x8) {
    let ax = x.abs();
    let (y0, y1, n) = reduce(x, ax);
    let s = k_sin(y0, y1, ax.simd_lt(UNREDUCED_END));
    let c = k_cos(y0, y1);
    let tiny = ax.simd_lt(TINY);
    let (sin, cos) = place(n, tiny.select(x, s), tiny.select(f64x8::ONE, c));
    let unhandled = ax.simd_ge(MEDIUM_END) | x.is_nan();
    (sin, cos, unhandled)
}

// x − n·π/2 as y0 + y1, and the bits of n's low end.
#[inline(always)]
fn reduce(x: f64x8, ax: f64x8) -> (f64x8, f64x8, u64x8) {
    let (f_n, stepped) = quotient(x, ax);
    let r1 = x - f_n * PIO2_1;
    let w1 = f_n * PIO2_1T;
    let y01 = r1 - w1;
    // libm runs another round when the last lost more than 16, then 49, bits
    // of x's exponent to cancellation.
    let unit = ax & EXPONENT;
    let second = not(stepped) & y01.abs().simd_lt(unit * SECOND_ROUND);
    let (r2, w2, y02) = refine(r1, f_n, PIO2_2, PIO2_2T);
    let third = second & y02.abs().simd_lt(unit * THIRD_ROUND);
    let (r3, w3, y03) = refine(r2, f_n, PIO2_3, PIO2_3T);
    let r = third.select(r3, second.select(r2, r1));
    let w = third.select(w3, second.select(w2, w1));
    let y0 = third.select(y03, second.select(y02, y01));
    (y0, (r - y0) - w, (f_n + TO_INT).to_bits())
}

// n: whole steps of π/2 up to 9π/4, except near multiples of π/2, and x·2/π
// rounded to the nearest integer past that.
#[inline(always)]
fn quotient(x: f64x8, ax: f64x8) -> (f64x8, f64x8) {
    let mut k = f64x8::ZERO;
    for step in STEPS {
        k += ax.simd_ge(step) & f64x8::ONE;
    }
    let k = x.simd_lt(f64x8::ZERO).select(f64x8::ZERO - k, k);
    let near = NEAR.iter().fold(f64x8::ZERO, |m, &(lo, hi)| {
        m | (ax.simd_ge(lo) & ax.simd_lt(hi))
    });
    let stepped = ax.simd_lt(STEPPED_END) & not(near);
    let rounded = (x * INV_PIO2 + TO_INT) - TO_INT;
    (stepped.select(k, rounded), stepped)
}

#[inline(always)]
fn refine(t: f64x8, f_n: f64x8, hi: f64x8, lo: f64x8) -> (f64x8, f64x8, f64x8) {
    let w = f_n * hi;
    let r = t - w;
    let w = f_n * lo - ((t - r) - w);
    (r, w, r - w)
}

#[inline(always)]
fn not(mask: f64x8) -> f64x8 {
    f64x8::from_bits(!mask.to_bits())
}

#[inline(always)]
fn k_sin(x: f64x8, y: f64x8, unreduced: f64x8) -> f64x8 {
    let z = x * x;
    let w = z * z;
    let r = S2 + z * (S3 + z * S4) + z * w * (S5 + z * S6);
    let v = z * x;
    let reduced = x - ((z * (HALF * y - v * r) - y) - v * S1);
    unreduced.select(x + v * (S1 + z * r), reduced)
}

#[inline(always)]
fn k_cos(x: f64x8, y: f64x8) -> f64x8 {
    let z = x * x;
    let w = z * z;
    let r = z * (C1 + z * (C2 + z * C3)) + w * w * (C4 + z * (C5 + z * C6));
    let hz = HALF * z;
    let w = f64x8::ONE - hz;
    w + (((f64x8::ONE - w) - hz) + (z * r - x * y))
}

const ONE: u64x8 = u64x8::splat(1);
const TWO: u64x8 = u64x8::splat(2);

// The sine and cosine of y + n·π/2 from those of y.
#[inline(always)]
fn place(n: u64x8, s: f64x8, c: f64x8) -> (f64x8, f64x8) {
    let odd = f64x8::from_bits(u64x8::ZERO - (n & ONE));
    let (sin, cos) = (odd.select(c, s), odd.select(s, c));
    let flip_sin = (n & TWO) << 62;
    let flip_cos = ((n + ONE) & TWO) << 62;
    (
        f64x8::from_bits(sin.to_bits() ^ flip_sin),
        f64x8::from_bits(cos.to_bits() ^ flip_cos),
    )
}
