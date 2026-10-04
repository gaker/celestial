// Terms of the s + XY/2 series, in arcseconds. SP holds the polynomial coefficients and
// S0..S4 the periodic terms that multiply t^0..t^4. Multipliers are for the arguments
// l, l', F, D, Ω, L_Ve, L_E, p_A.

#[derive(Clone, Copy)]
pub(super) struct SeriesTerm {
    pub(super) coeffs: [i8; 8],
    pub(super) sine: f64,
    pub(super) cosine: f64,
}

pub(super) const SP: [f64; 6] = [
    94.00e-6,
    3808.65e-6,
    -122.68e-6,
    -72574.11e-6,
    27.98e-6,
    15.62e-6,
];

pub(super) const S0: [SeriesTerm; 33] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: -2640.73e-6,
        cosine: 0.39e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: -63.53e-6,
        cosine: 0.02e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 3, 0, 0, 0],
        sine: -11.75e-6,
        cosine: -0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 1, 0, 0, 0],
        sine: -11.21e-6,
        cosine: -0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 2, 0, 0, 0],
        sine: 4.57e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 3, 0, 0, 0],
        sine: -2.02e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 1, 0, 0, 0],
        sine: -1.98e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 3, 0, 0, 0],
        sine: 1.72e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 0, 0, 1, 0, 0, 0],
        sine: 1.41e-6,
        cosine: 0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 0, 0, -1, 0, 0, 0],
        sine: 1.26e-6,
        cosine: 0.01e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, -1, 0, 0, 0],
        sine: 0.63e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, 1, 0, 0, 0],
        sine: 0.63e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 3, 0, 0, 0],
        sine: -0.46e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 1, 0, 0, 0],
        sine: -0.45e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 4, -4, 4, 0, 0, 0],
        sine: -0.36e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 1, -1, 1, -8, 12, 0],
        sine: 0.24e-6,
        cosine: 0.12e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 0, 0, 0, 0],
        sine: -0.32e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 2, 0, 0, 0],
        sine: -0.28e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 3, 0, 0, 0],
        sine: -0.27e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 1, 0, 0, 0],
        sine: -0.26e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 0, 0, 0, 0],
        sine: 0.21e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, -2, 2, -3, 0, 0, 0],
        sine: -0.19e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, -2, 2, -1, 0, 0, 0],
        sine: -0.18e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 0, 8, -13, -1],
        sine: 0.10e-6,
        cosine: -0.05e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 2, 0, 0, 0, 0],
        sine: -0.15e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, -2, 0, -1, 0, 0, 0],
        sine: 0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 2, 0, 0, 0],
        sine: 0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, -2, 1, 0, 0, 0],
        sine: -0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, -2, -1, 0, 0, 0],
        sine: -0.14e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 4, -2, 4, 0, 0, 0],
        sine: -0.13e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 4, 0, 0, 0],
        sine: 0.11e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, 0, -3, 0, 0, 0],
        sine: -0.11e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, 0, -1, 0, 0, 0],
        sine: -0.11e-6,
        cosine: 0.00e-6,
    },
];

pub(super) const S1: [SeriesTerm; 3] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: -0.07e-6,
        cosine: 3.57e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: 1.73e-6,
        cosine: -0.03e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 3, 0, 0, 0],
        sine: 0.00e-6,
        cosine: 0.48e-6,
    },
];

pub(super) const S2: [SeriesTerm; 25] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: 743.52e-6,
        cosine: -0.17e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 2, 0, 0, 0],
        sine: 56.91e-6,
        cosine: 0.06e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 2, 0, 0, 0],
        sine: 9.84e-6,
        cosine: -0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: -8.85e-6,
        cosine: 0.01e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 0, 0, 0, 0, 0, 0],
        sine: -6.38e-6,
        cosine: -0.05e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, 0, 0, 0, 0],
        sine: -3.07e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, 2, -2, 2, 0, 0, 0],
        sine: 2.23e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 1, 0, 0, 0],
        sine: 1.67e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 2, 0, 0, 0],
        sine: 1.30e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 1, -2, 2, -2, 0, 0, 0],
        sine: 0.93e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, -2, 0, 0, 0, 0],
        sine: 0.68e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 1, 0, 0, 0],
        sine: -0.55e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, 0, -2, 0, 0, 0],
        sine: 0.53e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 2, 0, 0, 0, 0],
        sine: -0.27e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, 1, 0, 0, 0],
        sine: -0.27e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, -2, -2, -2, 0, 0, 0],
        sine: -0.26e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 0, 0, -1, 0, 0, 0],
        sine: -0.25e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, 0, 1, 0, 0, 0],
        sine: 0.22e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, 0, -2, 0, 0, 0, 0],
        sine: -0.21e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, -2, 0, -1, 0, 0, 0],
        sine: 0.20e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 2, 2, 0, 0, 0],
        sine: 0.17e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, 2, 0, 2, 0, 0, 0],
        sine: 0.13e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [2, 0, 0, 0, 0, 0, 0, 0],
        sine: -0.13e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [1, 0, 2, -2, 2, 0, 0, 0],
        sine: -0.12e-6,
        cosine: 0.00e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 0, 0, 0, 0],
        sine: -0.11e-6,
        cosine: 0.00e-6,
    },
];

pub(super) const S3: [SeriesTerm; 4] = [
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
        sine: 0.30e-6,
        cosine: -23.42e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, -2, 2, 0, 0, 0],
        sine: -0.03e-6,
        cosine: -1.46e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 2, 0, 2, 0, 0, 0],
        sine: -0.01e-6,
        cosine: -0.25e-6,
    },
    SeriesTerm {
        coeffs: [0, 0, 0, 0, 2, 0, 0, 0],
        sine: 0.00e-6,
        cosine: 0.23e-6,
    },
];

pub(super) const S4: [SeriesTerm; 1] = [SeriesTerm {
    coeffs: [0, 0, 0, 0, 1, 0, 0, 0],
    sine: -0.26e-6,
    cosine: -0.01e-6,
}];
