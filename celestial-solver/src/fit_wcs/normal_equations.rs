//! Least-squares fits of two quantities as polynomials in two variables, through the
//! normal equations they share.

use crate::error::SolveError;

// Every (p, q) with p + q at most `order`, ordered by p and then by q.
pub(super) fn terms_up_to(order: u32) -> Vec<(u32, u32)> {
    (0..=order)
        .flat_map(|p| (0..=order - p).map(move |q| (p, q)))
        .collect()
}

// u^p v^q for each (p, q) of `terms`, whose p and q are at most `order`.
pub(super) fn monomials(terms: &[(u32, u32)], u: f64, v: f64, order: u32) -> Vec<f64> {
    let powers = |x: f64| {
        let mut p = vec![1.0_f64; order as usize + 1];
        for k in 1..p.len() {
            p[k] = p[k - 1] * x;
        }
        p
    };
    let (pu, pv) = (powers(u), powers(v));
    terms
        .iter()
        .map(|&(p, q)| pu[p as usize] * pv[q as usize])
        .collect()
}

// Normal equations shared by two right-hand sides.
pub(super) struct NormalEquations {
    n: usize,
    ata: Vec<f64>,
    atb: [Vec<f64>; 2],
}

impl NormalEquations {
    pub(super) fn new(n: usize) -> Self {
        Self {
            n,
            ata: vec![0.0; n * n],
            atb: [vec![0.0; n], vec![0.0; n]],
        }
    }

    // Adds the design row `row` with weight `w` and the targets `a` and `b`.
    pub(super) fn add(&mut self, row: &[f64], w: f64, a: f64, b: f64) {
        for (i, &ri) in row.iter().enumerate() {
            let wi = w * ri;
            for (j, &rj) in row.iter().enumerate() {
                self.ata[i * self.n + j] += wi * rj;
            }
            self.atb[0][i] += wi * a;
            self.atb[1][i] += wi * b;
        }
    }

    // The coefficients for the two targets. An error when the system is singular.
    pub(super) fn solve(&self) -> Result<(Vec<f64>, Vec<f64>), SolveError> {
        let [mut a, mut b] = self.atb.clone();
        solve_normal_equations(&mut self.ata.clone(), &mut a, self.n)?;
        solve_normal_equations(&mut self.ata.clone(), &mut b, self.n)?;
        Ok((a, b))
    }

    // The inverse of the normal matrix as a list of its columns, which equal its rows
    // since the matrix is symmetric. An error when it is singular.
    pub(super) fn inverse(&self) -> Result<Vec<Vec<f64>>, SolveError> {
        (0..self.n)
            .map(|j| {
                let mut unit = vec![0.0; self.n];
                unit[j] = 1.0;
                solve_normal_equations(&mut self.ata.clone(), &mut unit, self.n)?;
                Ok(unit)
            })
            .collect()
    }
}

// Gaussian elimination with partial pivoting. Rounding in the elimination leaves
// errors of about n times the machine epsilon times the largest entry (Golub and Van
// Loan, Matrix Computations, sec. 3.3), so a pivot no larger than that is
// indistinguishable from 0 and the matrix is singular to working precision.
// The solution replaces `b`.
pub(super) fn solve_normal_equations(
    a: &mut [f64],
    b: &mut [f64],
    n: usize,
) -> Result<(), SolveError> {
    let largest = a.iter().fold(0.0_f64, |m, &x| m.max(libm::fabs(x)));
    let singular_below = n as f64 * f64::EPSILON * largest;
    for col in 0..n {
        let (pivot_row, pivot_val) = largest_in_column(a, n, col);
        if pivot_val.is_nan() || pivot_val <= singular_below {
            return Err(SolveError::fit_failed(format!(
                "singular normal equations at col {col}"
            )));
        }
        if pivot_row != col {
            for k in 0..n {
                a.swap(col * n + k, pivot_row * n + k);
            }
            b.swap(col, pivot_row);
        }
        eliminate_below(a, b, n, col);
    }
    back_substitute(a, b, n);
    Ok(())
}

// The row at or below `col` with the largest entry in size in column `col`, the first
// of them on a tie, and that size.
fn largest_in_column(a: &[f64], n: usize, col: usize) -> (usize, f64) {
    let mut pivot_row = col;
    let mut pivot_val = libm::fabs(a[col * n + col]);
    for row in (col + 1)..n {
        let v = libm::fabs(a[row * n + col]);
        if v > pivot_val {
            pivot_val = v;
            pivot_row = row;
        }
    }
    (pivot_row, pivot_val)
}

// Subtracts multiples of row `col` from the rows below it to zero column `col` there.
fn eliminate_below(a: &mut [f64], b: &mut [f64], n: usize, col: usize) {
    let diag = a[col * n + col];
    for row in (col + 1)..n {
        let factor = a[row * n + col] / diag;
        for k in col..n {
            let above = a[col * n + k];
            a[row * n + k] -= factor * above;
        }
        b[row] -= factor * b[col];
    }
}

// Solves the upper-triangular system left by the elimination, in place in `b`.
fn back_substitute(a: &[f64], b: &mut [f64], n: usize) {
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= a[i * n + j] * b[j];
        }
        b[i] = sum / a[i * n + i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_normal_equations_identity_returns_b() {
        // 3x3 identity matrix, b = [1, 2, 3]
        let mut a = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let mut b = vec![1.0, 2.0, 3.0];
        solve_normal_equations(&mut a, &mut b, 3).unwrap();
        assert_eq!(b, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn solve_normal_equations_diagonal_returns_componentwise_divide() {
        // diag(2, 4, 8), b = [2, 8, 24] -> x = [1, 2, 3]
        let mut a = vec![2.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 8.0];
        let mut b = vec![2.0, 8.0, 24.0];
        solve_normal_equations(&mut a, &mut b, 3).unwrap();
        assert_eq!(b, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn solve_normal_equations_needs_pivoting_still_works() {
        // first-column pivot is zero, so row swap is required
        let mut a = vec![0.0, 1.0, 2.0, 0.0];
        let mut b = vec![3.0, 4.0];
        // x_0 = 4/2 = 2, x_1 = 3/1 = 3
        solve_normal_equations(&mut a, &mut b, 2).unwrap();
        assert_eq!(b, vec![2.0, 3.0]);
    }

    #[test]
    fn solve_normal_equations_singular_errors() {
        // rank-1 matrix: row 1 = 2x row 0
        let mut a = vec![1.0, 2.0, 2.0, 4.0];
        let mut b = vec![1.0, 2.0];
        assert!(solve_normal_equations(&mut a, &mut b, 2).is_err());
    }
}
