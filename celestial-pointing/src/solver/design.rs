use crate::observation::Observation;
use crate::terms::Term;
use nalgebra::{DMatrix, DVector};

pub(crate) fn build_residuals(observations: &[&Observation]) -> DVector<f64> {
    let n = observations.len();
    let mut b = DVector::zeros(2 * n);
    for (i, obs) in observations.iter().enumerate() {
        b[2 * i] = (obs.actual_ha - obs.commanded_ha).arcseconds();
        b[2 * i + 1] = (obs.observed_dec - obs.catalog_dec).arcseconds();
    }
    b
}

pub(super) fn build_design_matrix(
    observations: &[&Observation],
    terms: &[Box<dyn Term>],
    lat: f64,
) -> DMatrix<f64> {
    let n = observations.len();
    let m = terms.len();
    let mut a = DMatrix::zeros(2 * n, m);
    for (i, obs) in observations.iter().enumerate() {
        let h = obs.commanded_ha.radians();
        let dec = obs.catalog_dec.radians();
        let pier = obs.pier_side.sign();
        for (j, term) in terms.iter().enumerate() {
            let (jh, jd) = term.jacobian_equatorial(h, dec, lat, pier);
            a[(2 * i, j)] = jh;
            a[(2 * i + 1, j)] = jd;
        }
    }
    a
}

pub(super) fn cos_dec_per_obs(observations: &[&Observation]) -> Vec<f64> {
    observations
        .iter()
        .map(|o| libm::cos(o.catalog_dec.radians()))
        .collect()
}

pub(super) fn scale_ha_rows_by_cos_dec(b: &mut DVector<f64>, cos_dec: &[f64]) {
    for (i, &c) in cos_dec.iter().enumerate() {
        b[2 * i] *= c;
    }
}

pub(super) fn scale_a_ha_rows_by_cos_dec(a: &DMatrix<f64>, cos_dec: &[f64]) -> DMatrix<f64> {
    let mut out = a.clone();
    let cols = out.ncols();
    for (i, &c) in cos_dec.iter().enumerate() {
        for j in 0..cols {
            out[(2 * i, j)] *= c;
        }
    }
    out
}

pub(super) fn extract_columns(a: &DMatrix<f64>, cols: &[usize]) -> DMatrix<f64> {
    let rows = a.nrows();
    let m = cols.len();
    let mut out = DMatrix::zeros(rows, m);
    for (j, &col) in cols.iter().enumerate() {
        for i in 0..rows {
            out[(i, j)] = a[(i, col)];
        }
    }
    out
}

pub(super) fn subtract_all_contributions(
    b: &mut DVector<f64>,
    a: &DMatrix<f64>,
    coefficients: &[f64],
) {
    for (idx, &coeff) in coefficients.iter().enumerate() {
        if coeff == 0.0 {
            continue;
        }
        for row in 0..a.nrows() {
            b[row] -= a[(row, idx)] * coeff;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::terms::create_term;
    use crate::test_support::ObsBuilder;

    fn obs(commanded_ha_arcsec: f64, actual_ha_arcsec: f64, dec_deg: f64, obs_dec_deg: f64) -> Observation {
        ObsBuilder::new()
            .commanded_ha_arcsec(commanded_ha_arcsec)
            .actual_ha_arcsec(actual_ha_arcsec)
            .catalog_dec_deg(dec_deg)
            .observed_dec_deg(obs_dec_deg)
            .build()
    }

    // --- build_residuals ---------------------------------------------

    #[test]
    fn build_residuals_pairs_ha_and_dec_per_observation() {
        let o1 = obs(0.0, 100.0, 30.0, 30.5);
        let o2 = obs(0.0, -50.0, 45.0, 45.0);
        let observations: Vec<&Observation> = vec![&o1, &o2];
        let b = build_residuals(&observations);
        assert_eq!(b.len(), 4);
        assert!((b[0] - 100.0).abs() < 1e-9);
        assert!((b[1] - (0.5 * 3600.0)).abs() < 1e-6);
        assert!((b[2] - (-50.0)).abs() < 1e-9);
        assert!((b[3] - 0.0).abs() < 1e-9);
    }

    #[test]
    fn build_residuals_empty_input_returns_empty_vector() {
        let observations: Vec<&Observation> = Vec::new();
        let b = build_residuals(&observations);
        assert_eq!(b.len(), 0);
    }

    // --- build_design_matrix -----------------------------------------

    #[test]
    fn build_design_matrix_has_correct_shape() {
        let o = obs(0.0, 0.0, 30.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![
            create_term("IH").unwrap(),
            create_term("ID").unwrap(),
        ];
        let a = build_design_matrix(&observations, &terms, 0.0);
        assert_eq!(a.nrows(), 6); // 2 * n_obs
        assert_eq!(a.ncols(), 2); // n_terms
    }

    // IH affects HA only — its dec-row (every 2*i+1) column should be zero.
    // ID affects Dec only — its HA-row (every 2*i) column should be zero.
    #[test]
    fn build_design_matrix_ih_only_affects_ha_rows() {
        let o = obs(0.0, 0.0, 30.0, 30.0);
        let observations: Vec<&Observation> = vec![&o];
        let terms = vec![create_term("IH").unwrap(), create_term("ID").unwrap()];
        let a = build_design_matrix(&observations, &terms, 0.0);
        // IH column = 0: HA row nonzero, Dec row zero
        assert!(a[(0, 0)].abs() > 0.0);
        assert_eq!(a[(1, 0)], 0.0);
        // ID column = 1: HA row zero, Dec row nonzero
        assert_eq!(a[(0, 1)], 0.0);
        assert!(a[(1, 1)].abs() > 0.0);
    }

    #[test]
    fn build_design_matrix_zero_obs_zero_terms_yields_empty() {
        let observations: Vec<&Observation> = Vec::new();
        let terms: Vec<Box<dyn crate::terms::Term>> = Vec::new();
        let a = build_design_matrix(&observations, &terms, 0.0);
        assert_eq!(a.nrows(), 0);
        assert_eq!(a.ncols(), 0);
    }

    // --- cos_dec_per_obs ----------------------------------------------

    #[test]
    fn cos_dec_per_obs_returns_cosine_for_each_observation() {
        let o0 = obs(0.0, 0.0, 0.0, 0.0);
        let o60 = obs(0.0, 0.0, 60.0, 60.0);
        let o90 = obs(0.0, 0.0, 90.0, 90.0);
        let observations: Vec<&Observation> = vec![&o0, &o60, &o90];
        let cos = cos_dec_per_obs(&observations);
        assert_eq!(cos.len(), 3);
        assert!((cos[0] - 1.0).abs() < 1e-12);
        assert!((cos[1] - 0.5).abs() < 1e-12);
        assert!(cos[2].abs() < 1e-12); // cos(90°) ≈ 0
    }

    #[test]
    fn cos_dec_per_obs_empty_returns_empty() {
        let observations: Vec<&Observation> = Vec::new();
        assert!(cos_dec_per_obs(&observations).is_empty());
    }

    // --- scale_ha_rows_by_cos_dec ------------------------------------

    #[test]
    fn scale_ha_rows_by_cos_dec_only_touches_ha_rows() {
        let mut b = DVector::from_vec(vec![10.0, 20.0, 30.0, 40.0]);
        let cos_dec = vec![0.5, 2.0];
        scale_ha_rows_by_cos_dec(&mut b, &cos_dec);
        assert_eq!(b[0], 5.0); // HA row 0 scaled by 0.5
        assert_eq!(b[1], 20.0); // Dec row 0 unchanged
        assert_eq!(b[2], 60.0); // HA row 1 scaled by 2.0
        assert_eq!(b[3], 40.0); // Dec row 1 unchanged
    }

    // --- scale_a_ha_rows_by_cos_dec ----------------------------------

    #[test]
    fn scale_a_ha_rows_by_cos_dec_clones_and_scales_ha_rows() {
        let mut a = DMatrix::zeros(4, 2);
        for i in 0..4 {
            for j in 0..2 {
                a[(i, j)] = (i * 10 + j) as f64;
            }
        }
        let cos_dec = vec![0.5, 2.0];
        let scaled = scale_a_ha_rows_by_cos_dec(&a, &cos_dec);
        // Original unchanged (clone semantics)
        assert_eq!(a[(0, 0)], 0.0);
        // HA rows scaled
        assert_eq!(scaled[(0, 0)], 0.0); // 0 * 0.5
        assert_eq!(scaled[(0, 1)], 0.5); // 1 * 0.5
        assert_eq!(scaled[(2, 0)], 40.0); // 20 * 2.0
        // Dec rows untouched
        assert_eq!(scaled[(1, 0)], 10.0);
        assert_eq!(scaled[(3, 1)], 31.0);
    }

    // --- extract_columns ----------------------------------------------

    #[test]
    fn extract_columns_picks_named_columns_in_order() {
        let mut a = DMatrix::zeros(2, 4);
        for i in 0..2 {
            for j in 0..4 {
                a[(i, j)] = (i * 10 + j) as f64;
            }
        }
        let extracted = extract_columns(&a, &[2, 0]);
        assert_eq!(extracted.nrows(), 2);
        assert_eq!(extracted.ncols(), 2);
        assert_eq!(extracted[(0, 0)], 2.0); // col 2 first
        assert_eq!(extracted[(0, 1)], 0.0); // col 0 second
        assert_eq!(extracted[(1, 0)], 12.0);
        assert_eq!(extracted[(1, 1)], 10.0);
    }

    #[test]
    fn extract_columns_with_empty_cols_yields_zero_width() {
        let a = DMatrix::zeros(3, 2);
        let extracted = extract_columns(&a, &[]);
        assert_eq!(extracted.nrows(), 3);
        assert_eq!(extracted.ncols(), 0);
    }

    // --- subtract_all_contributions ----------------------------------

    #[test]
    fn subtract_all_contributions_subtracts_column_coeff_products() {
        let mut b = DVector::from_vec(vec![100.0, 200.0]);
        let mut a = DMatrix::zeros(2, 2);
        a[(0, 0)] = 1.0;
        a[(1, 0)] = 2.0;
        a[(0, 1)] = 3.0;
        a[(1, 1)] = 4.0;
        let coeffs = vec![10.0, 5.0];
        subtract_all_contributions(&mut b, &a, &coeffs);
        // b[0] = 100 - 1*10 - 3*5 = 75
        // b[1] = 200 - 2*10 - 4*5 = 160
        assert_eq!(b[0], 75.0);
        assert_eq!(b[1], 160.0);
    }

    // Zero coefficients are skipped — important optimization, document it.
    #[test]
    fn subtract_all_contributions_skips_zero_coefficients() {
        let mut b = DVector::from_vec(vec![100.0]);
        let mut a = DMatrix::zeros(1, 2);
        a[(0, 0)] = 999.0; // would wreak havoc if not skipped
        a[(0, 1)] = 1.0;
        let coeffs = vec![0.0, 5.0];
        subtract_all_contributions(&mut b, &a, &coeffs);
        assert_eq!(b[0], 95.0);
    }

    #[test]
    fn subtract_all_contributions_empty_coeffs_leaves_b_unchanged() {
        let mut b = DVector::from_vec(vec![1.0, 2.0]);
        let a = DMatrix::zeros(2, 0);
        subtract_all_contributions(&mut b, &a, &[]);
        assert_eq!(b[0], 1.0);
        assert_eq!(b[1], 2.0);
    }
}
