use celestial_core::constants::HALF_PI;

use super::fit::{Amplitudes, Arguments, Poly5};
use crate::lunar_coefficients::moon::{
    MainTerm, PertBlock, PertTerm, MAIN_DISTANCE, MAIN_LATITUDE, MAIN_LONGITUDE, PERT_DISTANCE,
    PERT_LATITUDE, PERT_LONGITUDE,
};
use crate::sincos;

// A term under one fit: its amplitude, and its argument as a polynomial in t.
// Neither depends on the epoch, so each fit works them out once.
struct Term {
    amplitude: f64,
    arg: Poly5,
}

struct Coordinate {
    main: Vec<Term>,
    pert: Vec<(usize, Vec<Term>)>,
}

// The longitude, latitude and distance series under one fit.
pub(super) struct Series([Coordinate; 3]);

impl Series {
    pub(super) fn new(args: &Arguments, amplitudes: &Amplitudes) -> Self {
        let coordinate = |main: &[MainTerm], pert: &[PertBlock], is_distance: bool| Coordinate {
            main: main
                .iter()
                .map(|term| main_term(term, args, amplitudes, is_distance))
                .collect(),
            pert: pert
                .iter()
                .map(|block| (usize::from(block.power), pert_terms(block.terms, args)))
                .collect(),
        };
        Self([
            coordinate(MAIN_LONGITUDE, PERT_LONGITUDE, false),
            coordinate(MAIN_LATITUDE, PERT_LATITUDE, false),
            coordinate(MAIN_DISTANCE, PERT_DISTANCE, true),
        ])
    }

    // The longitude, latitude and distance series, then their rates per
    // century. Each sum runs in the order of the authors' code so the results
    // match it.
    pub(super) fn sums<const RATES: bool>(&self, t: &Poly5) -> [f64; 6] {
        let mut v = [0.0; 6];
        for (iv, coordinate) in self.0.iter().enumerate() {
            let mut sum = (0.0, 0.0);
            main_sum::<RATES>(&coordinate.main, t, &mut sum);
            for (power, terms) in &coordinate.pert {
                pert_sum::<RATES>(terms, *power, t, &mut sum);
            }
            (v[iv], v[iv + 3]) = sum;
        }
        v
    }
}

// The authors' code starts the argument at the phase and adds each power's
// coefficient times tᵏ. With t⁰ = 1 the phase can join the constant
// coefficient here and give the same sums.
fn main_term(term: &MainTerm, args: &Arguments, amplitudes: &Amplitudes, distance: bool) -> Term {
    let phase = if distance { HALF_PI } else { 0.0 };
    let mut arg = [0.0; 5];
    for (k, a) in arg.iter_mut().enumerate() {
        *a = (0..4)
            .map(|i| term.delaunay[i] as f64 * args.del[i][k])
            .sum();
    }
    arg[0] += phase;
    Term {
        amplitude: amplitudes.corrected(&term.coeffs, distance),
        arg,
    }
}

fn pert_terms(terms: &[PertTerm], args: &Arguments) -> Vec<Term> {
    let term = |term: &PertTerm| Term {
        amplitude: term.amplitude,
        arg: [0, 1, 2, 3, 4].map(|k| pert_coefficient(term, args, k)),
    };
    terms.iter().map(term).collect()
}

// The phase joins the constant part of the argument before the multiples of
// the arguments are added, as in the authors' code.
fn pert_coefficient(term: &PertTerm, args: &Arguments, k: usize) -> f64 {
    let mut arg = if k == 0 { term.phase } else { 0.0 };
    for (&m, d) in term.multipliers[..4].iter().zip(&args.del) {
        arg += m as f64 * d[k];
    }
    for (&m, p) in term.multipliers[4..12].iter().zip(&args.p) {
        arg += m as f64 * p[k];
    }
    arg + term.multipliers[12] as f64 * args.zeta[k]
}

fn main_sum<const RATES: bool>(terms: &[Term], t: &Poly5, (val, deriv): &mut (f64, f64)) {
    sincos::for_each(
        terms,
        |term| argument::<RATES>(&term.arg, t),
        |term, yp, sin_y, cos_y| {
            *val += term.amplitude * sin_y;
            if RATES {
                *deriv += term.amplitude * yp * cos_y;
            }
        },
    );
}

fn pert_sum<const RATES: bool>(
    terms: &[Term],
    power: usize,
    t: &Poly5,
    (val, deriv): &mut (f64, f64),
) {
    sincos::for_each(
        terms,
        |term| argument::<RATES>(&term.arg, t),
        |term, yp, sin_y, cos_y| {
            let x = term.amplitude;
            *val += x * t[power] * sin_y;
            if RATES {
                let xp = amplitude_rate(x, power, t);
                *deriv = *deriv + xp * sin_y + x * t[power] * yp * cos_y;
            }
        },
    );
}

fn argument<const RATES: bool>(arg: &Poly5, t: &Poly5) -> (f64, f64) {
    let mut y = arg[0];
    let mut yp = 0.0;
    for k in 1..5 {
        y += arg[k] * t[k];
        if RATES {
            yp += (k as f64) * arg[k] * t[k - 1];
        }
    }
    (y, yp)
}

// The rate of the amplitude x tⁿ.
fn amplitude_rate(x: f64, n: usize, t: &Poly5) -> f64 {
    if n > 0 {
        (n as f64) * x * t[n - 1]
    } else {
        0.0
    }
}
