use crate::planetary_coefficients::{Term, TimeBlock};
use crate::sincos;

// The series of the elements a, λ, k, h, q and p.
pub(super) type Series = [&'static [TimeBlock]; 6];

// The arguments at J2000 (rad) and their rates (rad per millennium).
const LAMBDA0: [f64; 17] = [
    4.402608631669, // Mercury
    3.176134461576, // Venus
    1.753470369433, // Earth-Moon Barycenter
    6.203500014141, // Mars
    4.09136000305,  // Vesta
    1.713740719173, // Iris
    5.598641292287, // Bamberga
    2.805136360408, // Ceres
    2.32698973462,  // Pallas
    0.599546107035, // Jupiter
    0.874018510107, // Saturn
    5.481225395663, // Uranus
    5.311897933164, // Neptune
    0.0,            // mu = (n5 - n6) t / 880 has no constant term
    5.19846640063,  // Moon D
    1.62790513602,  // Moon F
    2.35555563875,  // Moon l
];

const LAMBDA_DOT: [f64; 17] = [
    26087.90314068555,  // Mercury
    10213.28554743445,  // Venus
    6283.075850353215,  // Earth-Moon Barycenter
    3340.612434145457,  // Mars
    1731.170452721855,  // Vesta
    1704.450855027201,  // Iris
    1428.948917844273,  // Bamberga
    1364.75651362999,   // Ceres
    1361.923207632842,  // Pallas
    529.690961562325,   // Jupiter
    213.299086108488,   // Saturn
    74.781659030778,    // Uranus
    38.132972226125,    // Neptune
    0.3595362285049309, // mu
    77713.7714481804,   // Moon D
    84334.6615717837,   // Moon F
    83286.9142477147,   // Moon l
];

// The elements at t millennia from J2000 and their rates per millennium. The
// rates are computed even when only the values are wanted: next to the sines
// they cost about 3%, which doesn't pay for a second copy of the loop.
pub(super) fn elements(series: &Series, t: f64) -> ([f64; 6], [f64; 6]) {
    let mut lambdas = [0.0; 17];
    for ((l, l0), n) in lambdas.iter_mut().zip(LAMBDA0).zip(LAMBDA_DOT) {
        *l = l0 + n * t;
    }
    let mut values = [0.0; 6];
    let mut rates = [0.0; 6];
    for ((value, rate), blocks) in values.iter_mut().zip(&mut rates).zip(series) {
        (*value, *rate) = variable(blocks, t, &lambdas);
    }
    (values, rates)
}

fn variable(blocks: &[TimeBlock], t: f64, lambdas: &[f64; 17]) -> (f64, f64) {
    let mut value = 0.0;
    let mut rate = 0.0;
    for block in blocks {
        let (sum, sum_rate) = periodic(block.terms, lambdas);
        let (tn, tn_rate) = power(t, block.power);
        value += sum * tn;
        rate += sum_rate * tn + sum * tn_rate;
    }
    (value, rate)
}

// t^n and its derivative n·t^(n−1).
fn power(t: f64, n: u8) -> (f64, f64) {
    let mut previous = 0.0;
    let mut tn = 1.0;
    for _ in 0..n {
        previous = tn;
        tn *= t;
    }
    (tn, f64::from(n) * previous)
}

fn periodic(terms: &[Term], lambdas: &[f64; 17]) -> (f64, f64) {
    let mut sum = 0.0;
    let mut rate = 0.0;
    sincos::for_each(
        terms,
        |term| argument(term, lambdas),
        |term, arg_rate, sin_arg, cos_arg| {
            sum += term.s * sin_arg + term.c * cos_arg;
            rate += arg_rate * (term.s * cos_arg - term.c * sin_arg);
        },
    );
    (sum, rate)
}

// The padding slots add 0 · λ = ±0, which leaves the sum unchanged: it starts
// at +0 and never becomes −0, so this matches skipping them.
fn argument(term: &Term, lambdas: &[f64; 17]) -> (f64, f64) {
    let mut arg = 0.0;
    let mut rate = 0.0;
    for (&m, &i) in term.mult.iter().zip(&term.index) {
        let m = f64::from(m);
        arg += m * lambdas[usize::from(i)];
        rate += m * LAMBDA_DOT[usize::from(i)];
    }
    (arg, rate)
}
