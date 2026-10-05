use std::sync::OnceLock;

use celestial_core::constants::{ARCSEC_PER_RAD, DEG_TO_RAD, PI};

use super::frame::{self, Frame};
use super::series::Series;

pub(super) type Poly5 = [f64; 5];

const AM: f64 = 0.074801329;
const ALPHA: f64 = 0.002571881;
const DTASM: f64 = (2.0 * ALPHA) / (3.0 * AM);
const XA: f64 = (2.0 * ALPHA) / 3.0;
const DPREC: f64 = -0.29965;

// Partial derivatives of the node (W2) and perigee (W3) rates with respect to
// the fitted constants.
const BP: [Poly5; 2] = [
    [
        0.311079095,
        -0.4482398e-2,
        -0.1102485e-2,
        0.1056062e-2,
        0.50928e-4,
    ],
    [
        -0.103837907,
        0.668287e-3,
        -0.1298072e-2,
        -0.178028e-3,
        -0.37342e-4,
    ],
];

// One fit's corrections to the constants of the solution, in arcseconds and
// arcseconds per century.
struct Adjustments {
    dw1: [f64; 3],
    dw2: [f64; 2],
    dw3: [f64; 2],
    deart: [f64; 2],
    dperi: f64,
    dgam: f64,
    de: f64,
    dep: f64,
    secular: [Poly5; 3],
    frame: (f64, f64),
}

const LLR: Adjustments = Adjustments {
    dw1: [-0.10525, -0.32311, -0.03794],
    dw2: [0.16826, 0.08017],
    dw3: [-0.10760, -0.04317],
    deart: [-0.04012, 0.01442],
    dperi: -0.04854,
    dgam: 0.00069,
    de: 0.00005,
    dep: 0.00226,
    secular: [[0.0; 5]; 3],
    frame: frame::LLR,
};

const DE405: Adjustments = Adjustments {
    dw1: [-0.07008, -0.35106, -0.03743],
    dw2: [0.20794, 0.08017],
    dw3: [-0.07215, -0.04317],
    deart: [-0.00033, 0.00732],
    dperi: -0.00749,
    dgam: 0.00085,
    de: -0.00006,
    dep: 0.00224,
    secular: [
        [0.0, 0.0, 0.0, -0.00018865, -0.00001024],
        [0.0, 0.0, 0.00470602, -0.00025213, 0.0],
        [0.0, 0.0, -0.00261070, -0.00010712, 0.0],
    ],
    frame: frame::DE405,
};

// The fitted constants that the B coefficients of the main problem are
// derivatives with respect to: Δν/ν, Δe, Δγ, Δn′/ν and Δe′.
pub(super) struct Amplitudes {
    delnu: f64,
    dele: f64,
    delg: f64,
    delnp: f64,
    delep: f64,
}

impl Amplitudes {
    fn new(adj: &Adjustments, w11: f64) -> Self {
        Self {
            delnu: (0.55604 + adj.dw1[1]) / ARCSEC_PER_RAD / w11,
            dele: (0.01789 + adj.de) / ARCSEC_PER_RAD,
            delg: (-0.08066 + adj.dgam) / ARCSEC_PER_RAD,
            delnp: (-0.06424 + adj.deart[1]) / ARCSEC_PER_RAD / w11,
            delep: (-0.12879 + adj.dep) / ARCSEC_PER_RAD,
        }
    }

    pub(super) fn corrected(&self, coeffs: &[f64; 7], is_distance: bool) -> f64 {
        let [a, b1, b2, b3, b4, b5, _] = *coeffs;
        let tgv = b1 + DTASM * b5;
        let a = if is_distance {
            a - 2.0 * a * self.delnu / 3.0
        } else {
            a
        };
        a + tgv * (self.delnp - AM * self.delnu) + b2 * self.delg + b3 * self.dele + b4 * self.delep
    }
}

// The polynomials in t of the arguments the series are built from: the
// Delaunay arguments, the planets' mean longitudes and ζ.
pub(super) struct Arguments {
    pub(super) del: [Poly5; 4],
    pub(super) p: [Poly5; 8],
    pub(super) zeta: Poly5,
}

// Everything in the solution that depends on the fit but not on the epoch.
// Each is built on first use and kept for the life of the program.
pub(super) struct Fit {
    pub(super) w1: Poly5,
    pub(super) series: Series,
    pub(super) frame: Frame,
}

impl Fit {
    pub(super) fn llr() -> &'static Self {
        static FIT: OnceLock<Fit> = OnceLock::new();
        FIT.get_or_init(|| Self::new(&LLR))
    }

    pub(super) fn de405() -> &'static Self {
        static FIT: OnceLock<Fit> = OnceLock::new();
        FIT.get_or_init(|| Self::new(&DE405))
    }

    fn new(adj: &Adjustments) -> Self {
        let w = moon_longitudes(adj);
        let mut zeta = w[0];
        zeta[1] += (5029.0966 + DPREC) / ARCSEC_PER_RAD;
        let args = Arguments {
            del: delaunay(&w, &earth_moon_barycenter(adj), &perihelion(adj)),
            p: planets(),
            zeta,
        };
        Self {
            w1: w[0],
            series: Series::new(&args, &Amplitudes::new(adj, w[0][1])),
            frame: Frame::new(adj.frame),
        }
    }
}

fn dms(deg: i32, min: i32, sec: f64) -> f64 {
    (deg as f64 + min as f64 / 60.0 + sec / 3600.0) * DEG_TO_RAD
}

// Mean longitudes of the Moon (W1), its node (W2) and its perigee (W3).
fn moon_longitudes(adj: &Adjustments) -> [Poly5; 3] {
    let mut w = [
        [
            dms(218, 18, 59.95571 + adj.dw1[0]),
            (1732559343.73604 + adj.dw1[1]) / ARCSEC_PER_RAD,
            (-6.8084 + adj.dw1[2]) / ARCSEC_PER_RAD,
            0.66040e-2 / ARCSEC_PER_RAD,
            -0.31690e-4 / ARCSEC_PER_RAD,
        ],
        [
            dms(83, 21, 11.67475 + adj.dw2[0]),
            (14643420.3171 + adj.dw2[1]) / ARCSEC_PER_RAD,
            -38.2631 / ARCSEC_PER_RAD,
            -0.45047e-1 / ARCSEC_PER_RAD,
            0.21301e-3 / ARCSEC_PER_RAD,
        ],
        [
            dms(125, 2, 40.39816 + adj.dw3[0]),
            (-6967919.5383 + adj.dw3[1]) / ARCSEC_PER_RAD,
            6.3590 / ARCSEC_PER_RAD,
            0.76250e-2 / ARCSEC_PER_RAD,
            -0.35860e-4 / ARCSEC_PER_RAD,
        ],
    ];
    add_secular_terms(&mut w, &adj.secular);
    correct_node_and_perigee_rates(&mut w, adj);
    w
}

fn add_secular_terms(w: &mut [Poly5; 3], secular: &[Poly5; 3]) {
    for (wi, si) in w.iter_mut().zip(secular) {
        for (x, s) in wi.iter_mut().zip(si) {
            *x += s / ARCSEC_PER_RAD;
        }
    }
}

fn correct_node_and_perigee_rates(w: &mut [Poly5; 3], adj: &Adjustments) {
    for (k, bp) in BP.iter().enumerate() {
        let x = w[k + 1][1] / w[0][1];
        let y = AM * bp[0] + XA * bp[4];
        let d1 = x - y;
        let d2 = w[0][1] * bp[1];
        let d3 = w[0][1] * bp[2];
        let d4 = w[0][1] * bp[3];
        let d5 = y / AM;
        let cw = d1 * adj.dw1[1] + d5 * adj.deart[1] + d2 * adj.dgam + d3 * adj.de + d4 * adj.dep;
        w[k + 1][1] += cw / ARCSEC_PER_RAD;
    }
}

fn earth_moon_barycenter(adj: &Adjustments) -> Poly5 {
    [
        dms(100, 27, 59.13885 + adj.deart[0]),
        (129597742.29300 + adj.deart[1]) / ARCSEC_PER_RAD,
        -0.020200 / ARCSEC_PER_RAD,
        0.90000e-5 / ARCSEC_PER_RAD,
        0.15000e-6 / ARCSEC_PER_RAD,
    ]
}

fn perihelion(adj: &Adjustments) -> Poly5 {
    [
        dms(102, 56, 14.45766 + adj.dperi),
        1161.24342 / ARCSEC_PER_RAD,
        0.529265 / ARCSEC_PER_RAD,
        -0.11814e-3 / ARCSEC_PER_RAD,
        0.11379e-4 / ARCSEC_PER_RAD,
    ]
}

// The Delaunay arguments D, F, l and l′.
fn delaunay(w: &[Poly5; 3], eart: &Poly5, peri: &Poly5) -> [Poly5; 4] {
    let mut del = [[0.0; 5]; 4];
    for i in 0..5 {
        del[0][i] = w[0][i] - eart[i];
        del[1][i] = w[0][i] - w[2][i];
        del[2][i] = w[0][i] - w[1][i];
        del[3][i] = eart[i] - peri[i];
    }
    del[0][0] += PI;
    del
}

// Mean longitudes of Mercury through Neptune, with the Earth-Moon barycenter
// in third place.
fn planets() -> [Poly5; 8] {
    let longitude =
        |deg, min, sec, rate: f64| [dms(deg, min, sec), rate / ARCSEC_PER_RAD, 0.0, 0.0, 0.0];
    [
        longitude(252, 15, 3.216919, 538101628.66888),
        longitude(181, 58, 44.758419, 210664136.45777),
        longitude(100, 27, 59.138850, 129597742.29300),
        longitude(355, 26, 3.642778, 68905077.65936),
        longitude(34, 21, 5.379392, 10925660.57335),
        longitude(50, 4, 38.902495, 4399609.33632),
        longitude(314, 3, 4.354234, 1542482.57845),
        longitude(304, 20, 56.808371, 786547.89700),
    ]
}
