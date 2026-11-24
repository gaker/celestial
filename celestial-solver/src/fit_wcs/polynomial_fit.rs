//! Weighted least-squares polynomials from pixel offsets to tangent-plane positions,
//! with the order chosen by leave-one-out cross-validation.
//!
//! The order, from 1 (TAN only) to the highest allowed, is the one with the least
//! prediction error sum of squares, PRESS (Allen 1974, Technometrics 16, 125): the sum
//! over the pairs of the weighted squared residual of each pair under the fit to the
//! other pairs. For weighted least squares that residual is e_i / (1 - h_i), where e_i
//! is the pair's residual under the fit to all pairs and h_i = w_i x_iᵀ (Xᵀ W X)⁻¹ x_i
//! is its leverage, so no refit is needed. An order whose fit interpolates a pair
//! (h_i = 1) predicts nothing there and is not chosen. An order with at least as many
//! terms per axis as pairs interpolates every pair, and rounding can leave each h_i a
//! little below 1, so it is not scored. Of two orders with equal PRESS the lower is
//! kept.

use super::cd_matrix::check_cd_determinant;
use super::normal_equations::{monomials, terms_up_to, NormalEquations};
use super::solution::InverseCd;
use super::tangent_projection::{ProjectedPair, TangentProjection};
use crate::error::SolveError;

// The polynomial of one order for each tangent-plane axis, in degrees, with
// coefficients of the monomials of (u / scale.0, v / scale.1).
pub(super) struct Polynomial {
    pub(super) order: u32,
    terms: Vec<(u32, u32)>,
    scale: (f64, f64),
    xi: Vec<f64>,
    eta: Vec<f64>,
    inverse_cd: InverseCd,
}

impl Polynomial {
    // The (xi, eta) coefficients of the monomial u^p v^q, in degrees, 0 when the
    // polynomial has no such term.
    pub(super) fn coefficients(&self, term: (u32, u32)) -> (f64, f64) {
        (
            self.coefficient(&self.xi, term),
            self.coefficient(&self.eta, term),
        )
    }

    // The coefficients of u^p v^q mapped through the inverse CD matrix, in pixels.
    pub(super) fn pixel_coefficients(&self, term: (u32, u32)) -> (f64, f64) {
        let (xi, eta) = self.coefficients(term);
        self.inverse_cd.apply(xi, eta)
    }

    fn coefficient(&self, axis: &[f64], term: (u32, u32)) -> f64 {
        let at = self.terms.iter().position(|&t| t == term);
        at.map_or(0.0, |i| axis[i] / scale_of(term, self.scale))
    }

    pub(super) fn cd(&self) -> [f64; 4] {
        let c = |axis: &[f64], term| self.coefficient(axis, term);
        [
            c(&self.xi, (1, 0)),
            c(&self.xi, (0, 1)),
            c(&self.eta, (1, 0)),
            c(&self.eta, (0, 1)),
        ]
    }

    fn row(&self, u: f64, v: f64) -> Vec<f64> {
        monomials(&self.terms, u / self.scale.0, v / self.scale.1, self.order)
    }

    // The residual in pixels of a tangent-plane row: the tangent-plane residual mapped
    // through the inverse CD matrix.
    pub(super) fn pixel_error(&self, row: &(f64, f64, f64, f64)) -> f64 {
        let x = self.row(row.0, row.1);
        let dot = |c: &[f64]| c.iter().zip(&x).map(|(a, b)| a * b).sum::<f64>();
        let (du, dv) = self
            .inverse_cd
            .apply(row.2 - dot(&self.xi), row.3 - dot(&self.eta));
        libm::hypot(du, dv)
    }
}

// The fit to the pairs `mask` keeps, about the tangent point of `tp`, with the order up
// to `max_order` that has the least PRESS. An error when the TAN fit itself fails.
pub(super) fn fit_best_order(
    tp: &TangentProjection<'_>,
    mask: &[bool],
    max_order: u32,
    scale: (f64, f64),
) -> Result<Polynomial, SolveError> {
    let (mut best, normal) = fit_order(tp, mask, 1, scale)?;
    if max_order < 2 {
        return Ok(best);
    }
    let mut best_press = press(&best, &normal, tp, mask).unwrap_or(f64::INFINITY);
    let pairs = mask.iter().filter(|&&on| on).count();
    let scored = (2..=max_order).take_while(|&order| terms_up_to(order).len() < pairs);
    for order in scored {
        let Ok((fit, normal)) = fit_order(tp, mask, order, scale) else {
            continue;
        };
        let Some(score) = press(&fit, &normal, tp, mask) else {
            continue;
        };
        log::debug!("SIP order {order}: PRESS {score:.6e} (best {best_press:.6e})");
        if score < best_press {
            (best, best_press) = (fit, score);
        }
    }
    Ok(best)
}

// The factor by which the coefficient of u^p v^q fitted to (u / scale.0, v / scale.1)
// exceeds the coefficient of u^p v^q: scale.0^p scale.1^q.
pub(super) fn scale_of((p, q): (u32, u32), scale: (f64, f64)) -> f64 {
    libm::pow(scale.0, f64::from(p)) * libm::pow(scale.1, f64::from(q))
}

// The weighted least-squares polynomial of `order`, and its normal equations.
pub(super) fn fit_order(
    tp: &TangentProjection<'_>,
    mask: &[bool],
    order: u32,
    scale: (f64, f64),
) -> Result<(Polynomial, NormalEquations), SolveError> {
    let terms = terms_up_to(order);
    let mut normal = NormalEquations::new(terms.len());
    for (row, w) in active_rows(tp, mask) {
        let x = monomials(&terms, row.0 / scale.0, row.1 / scale.1, order);
        normal.add(&x, w, row.2, row.3);
    }
    let (xi, eta) = normal.solve()?;
    let mut fit = Polynomial {
        order,
        terms,
        scale,
        xi,
        eta,
        inverse_cd: InverseCd::new(0.0, 0.0, 0.0, 0.0),
    };
    let [a, b, c, d] = fit.cd();
    check_cd_determinant(a * d - b * c)?;
    fit.inverse_cd = InverseCd::new(a, b, c, d);
    Ok((fit, normal))
}

// PRESS of `fit` over the pairs `mask` keeps, in square pixels. `None` when the normal
// matrix has no inverse or a pair has leverage 1 or more.
fn press(
    fit: &Polynomial,
    normal: &NormalEquations,
    tp: &TangentProjection<'_>,
    mask: &[bool],
) -> Option<f64> {
    let inverse = normal.inverse().ok()?;
    let mut sum = 0.0;
    for (row, w) in active_rows(tp, mask) {
        let x = fit.row(row.0, row.1);
        let h = w * quadratic_form(&inverse, &x);
        if h.is_nan() || h >= 1.0 {
            return None;
        }
        let e = fit.pixel_error(row) / (1.0 - h);
        sum += w * e * e;
    }
    Some(sum)
}

// The rows of `tp` that `mask` keeps, with their weights.
fn active_rows<'a>(
    tp: &'a TangentProjection<'_>,
    mask: &'a [bool],
) -> impl Iterator<Item = (&'a ProjectedPair, f64)> {
    tp.rows
        .iter()
        .zip(&tp.weights)
        .zip(mask)
        .filter(|(_, &on)| on)
        .map(|((row, &w), _)| (row, w))
}

// x^T M x for the symmetric matrix `m`, given as its rows.
fn quadratic_form(m: &[Vec<f64>], x: &[f64]) -> f64 {
    m.iter()
        .zip(x)
        .map(|(r, xi)| xi * r.iter().zip(x).map(|(a, b)| a * b).sum::<f64>())
        .sum()
}

pub(super) fn pixel_errors(
    fit: &Polynomial,
    tp: &TangentProjection<'_>,
    mask: &[bool],
) -> Vec<f64> {
    let on = tp.rows.iter().zip(mask);
    on.map(|(row, &on)| if on { fit.pixel_error(row) } else { 0.0 })
        .collect()
}
