use celestial_core::constants::TWOPI;
use celestial_core::errors::{AstroError, AstroResult, MathErrorKind};

const KEPLER_MAX_ITERATIONS: usize = 32;

// An ellipse given by the elements a, λ, k, h, q and p, solved for the
// eccentric longitude F.
pub(super) struct Orbit {
    elements: [f64; 6],
    xfi: f64,
    xki: f64,
    u: f64,
    sin_f: f64,
    cos_f: f64,
    // k cos F + h sin F and k sin F − h cos F.
    z3: (f64, f64),
}

impl Orbit {
    pub(super) fn new(elements: [f64; 6]) -> AstroResult<Self> {
        let [_, lambda, k, h, q, p] = elements;
        let xfi = libm::sqrt(1.0 - k * k - h * h);
        let xki = libm::sqrt(1.0 - q * q - p * p);
        let ex = libm::sqrt(k * k + h * h);
        check_elliptic(ex, xki)?;
        let (sin_f, cos_f) = libm::sincos(solve_kepler(mean_longitude(lambda), k, h, ex)?);
        Ok(Self {
            elements,
            xfi,
            xki,
            u: 1.0 / (1.0 + xfi),
            sin_f,
            cos_f,
            z3: (k * cos_f + h * sin_f, k * sin_f - h * cos_f),
        })
    }

    pub(super) fn position(&self) -> [f64; 3] {
        let [a, _, _, _, q, p] = self.elements;
        let (bx, by) = self.in_plane();
        let rsa = 1.0 - self.z3.0;
        let (zto_real, zto_imag) = (bx / rsa, by / rsa);
        let xm = p * zto_real - q * zto_imag;
        let xr = a * rsa;
        [
            xr * (zto_real - 2.0 * p * xm),
            xr * (zto_imag + 2.0 * q * xm),
            -2.0 * xr * self.xki * xm,
        ]
    }

    // The derivative of the position, given the elements' rates, in the rates'
    // unit of time.
    pub(super) fn velocity(&self, rates: &[f64; 6]) -> [f64; 3] {
        let [a, _, _, _, q, p] = self.elements;
        let [da, _, _, _, dq, dp] = *rates;
        let (bx, by) = self.in_plane();
        let (dbx, dby) = self.in_plane_rates(rates);
        let (x1, y1) = (a * bx, a * by);
        let (dx1, dy1) = (da * bx + a * dbx, da * by + a * dby);
        let m = p * x1 - q * y1;
        let dm = dp * x1 + p * dx1 - dq * y1 - q * dy1;
        let dxki = -(q * dq + p * dp) / self.xki;
        [
            dx1 - 2.0 * (dp * m + p * dm),
            dy1 + 2.0 * (dq * m + q * dm),
            -2.0 * (dxki * m + self.xki * dm),
        ]
    }

    // The coordinates in the orbital plane, over a.
    fn in_plane(&self) -> (f64, f64) {
        let [_, _, k, h, ..] = self.elements;
        let w = self.z3.1;
        (
            -k + self.cos_f + self.u * h * w,
            -h + self.sin_f - self.u * k * w,
        )
    }

    fn in_plane_rates(&self, rates: &[f64; 6]) -> (f64, f64) {
        let [_, _, k, h, ..] = self.elements;
        let [_, dlambda, dk, dh, ..] = *rates;
        let (z3r, w) = self.z3;
        let df = (dlambda + dk * self.sin_f - dh * self.cos_f) / (1.0 - z3r);
        let dw = dk * self.sin_f - dh * self.cos_f + df * z3r;
        let du = self.u * self.u * (k * dk + h * dh) / self.xfi;
        (
            -dk - self.sin_f * df + du * h * w + self.u * (dh * w + h * dw),
            -dh + self.cos_f * df - du * k * w - self.u * (dk * w + k * dw),
        )
    }
}

fn check_elliptic(ex: f64, xki: f64) -> AstroResult<()> {
    if ex < 1.0 && xki > 0.0 {
        return Ok(());
    }
    Err(AstroError::math_error(
        "VSOP2013 elements",
        MathErrorKind::OutOfRange,
        "eccentricity or inclination is outside the elliptic domain",
    ))
}

fn mean_longitude(lambda: f64) -> f64 {
    let l = lambda % TWOPI;
    if l < 0.0 {
        l + TWOPI
    } else {
        l
    }
}

// Solves λ = F − k sin F + h cos F for F.
fn solve_kepler(gl: f64, k: f64, h: f64, ex: f64) -> AstroResult<f64> {
    let gm = gl - libm::atan2(h, k);
    let ex3 = ex * ex * ex;
    let mut f = gl
        + (ex - 0.125 * ex3) * libm::sin(gm)
        + 0.5 * (ex * ex) * libm::sin(2.0 * gm)
        + 0.375 * ex3 * libm::sin(3.0 * gm);
    for _ in 0..KEPLER_MAX_ITERATIONS {
        let (sin_f, cos_f) = libm::sincos(f);
        let dl = gl - f + (k * sin_f - h * cos_f);
        f += dl / (1.0 - (k * cos_f + h * sin_f));
        if libm::fabs(dl) < 1e-15 {
            return Ok(f);
        }
    }
    Err(AstroError::calculation_error(
        "VSOP2013 Kepler equation",
        &format!("no convergence after {} iterations", KEPLER_MAX_ITERATIONS),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn orbit(lambda: f64, k: f64) -> AstroResult<Orbit> {
        Orbit::new([1.0, lambda, k, 0.0, 0.0, 0.0])
    }

    #[test]
    fn position_of_ten_digit_elements() {
        // Pluto's elements at J2000 from the full series, rounded to ten
        // digits. The unrounded elements give -9.8753625435, -27.9588613710
        // and 5.8504463318 AU, and the rounding moves each by up to 2e-9 AU.
        let orbit = Orbit::new([
            39.2648542648,
            4.1726045776,
            -0.1758641167,
            -0.1701234143,
            -0.0517015914,
            0.1398654514,
        ])
        .unwrap();
        assert_eq!(
            orbit.position(),
            [-9.875362545449178, -27.95886136998148, 5.850446333112353]
        );
    }

    #[test]
    fn circular_orbit_velocity_is_the_mean_motion_along_the_track() {
        let rates = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0];
        for lambda in [0.3, 2.0, 4.5] {
            let (sin_l, cos_l) = libm::sincos(lambda);
            let velocity = orbit(lambda, 0.0).unwrap().velocity(&rates);
            assert_eq!(velocity, [-2.0 * sin_l, 2.0 * cos_l, 0.0]);
        }
    }

    #[test]
    fn non_converging_kepler_equation_is_an_error() {
        let result = orbit(f64::NAN, 0.1);
        assert!(
            matches!(result, Err(AstroError::CalculationError { .. })),
            "{:?}",
            result.map(|o| o.position())
        );
    }

    #[test]
    fn hyperbolic_elements_are_an_error() {
        let result = orbit(1.0, 1.5);
        assert!(
            matches!(
                result,
                Err(AstroError::MathError {
                    kind: MathErrorKind::OutOfRange,
                    ..
                })
            ),
            "{:?}",
            result.map(|o| o.position())
        );
    }
}
