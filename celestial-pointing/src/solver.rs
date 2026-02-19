use crate::error::{Error, Result};
use crate::observation::Observation;
use crate::terms::Term;
use nalgebra::{DMatrix, DVector};

#[derive(Clone)]
pub struct FitResult {
    pub coefficients: Vec<f64>,
    pub sigma: Vec<f64>,
    pub sky_rms: f64,
    pub term_names: Vec<String>,
}

pub fn fit_model(
    observations: &[&Observation],
    terms: &[Box<dyn Term>],
    fixed: &[bool],
    coefficients: &[f64],
    latitude: f64,
) -> Result<FitResult> {
    let free_count = fixed.iter().filter(|&&f| !f).count();
    if free_count == 0 && !terms.is_empty() {
        return Err(Error::Fit("all terms are fixed".into()));
    }
    if terms.is_empty() {
        return Err(Error::Fit("no terms to fit".into()));
    }
    if observations.len() < free_count {
        return Err(Error::Fit("insufficient observations".into()));
    }

    let free_indices: Vec<usize> = fixed
        .iter()
        .enumerate()
        .filter(|(_, &f)| !f)
        .map(|(i, _)| i)
        .collect();
    let fixed_indices: Vec<usize> = fixed
        .iter()
        .enumerate()
        .filter(|(_, &f)| f)
        .map(|(i, _)| i)
        .collect();

    let mut b = build_residuals(observations);
    let w = build_weights(observations);
    let a_full = build_design_matrix(observations, terms, latitude);

    subtract_fixed_contributions(&mut b, &a_full, coefficients, &fixed_indices);

    let a_free = extract_columns(&a_full, &free_indices);
    let free_coeffs = solve_weighted(&a_free, &b, &w)?;
    let free_residuals = &b - &a_free * &free_coeffs;

    let mut all_coeffs = vec![0.0; terms.len()];
    for (fi, &idx) in free_indices.iter().enumerate() {
        all_coeffs[idx] = free_coeffs[fi];
    }

    let full_residuals = build_residuals(observations);
    let full_a = &a_full;
    let all_coeffs_dv = DVector::from_vec(all_coeffs.clone());
    let actual_residuals = &full_residuals - full_a * &all_coeffs_dv;

    let sigma = compute_sigma_free(&a_free, &free_residuals, &w, &free_indices, terms.len());
    let sky_rms = compute_sky_rms(&actual_residuals, observations);
    let term_names = terms.iter().map(|t| t.name().to_string()).collect();

    Ok(FitResult {
        coefficients: all_coeffs,
        sigma,
        sky_rms,
        term_names,
    })
}

fn subtract_fixed_contributions(
    b: &mut DVector<f64>,
    a: &DMatrix<f64>,
    coefficients: &[f64],
    fixed_indices: &[usize],
) {
    for &idx in fixed_indices {
        let coeff = coefficients[idx];
        for row in 0..a.nrows() {
            b[row] -= a[(row, idx)] * coeff;
        }
    }
}

fn extract_columns(a: &DMatrix<f64>, cols: &[usize]) -> DMatrix<f64> {
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

pub fn build_residuals(observations: &[&Observation]) -> DVector<f64> {
    let n = observations.len();
    let mut b = DVector::zeros(2 * n);
    for (i, obs) in observations.iter().enumerate() {
        b[2 * i] = (obs.actual_ha - obs.commanded_ha).arcseconds();
        b[2 * i + 1] = (obs.observed_dec - obs.catalog_dec).arcseconds();
    }
    b
}

fn build_weights(observations: &[&Observation]) -> DVector<f64> {
    let n = observations.len();
    let mut w = DVector::zeros(2 * n);
    for (i, obs) in observations.iter().enumerate() {
        let cos_dec = libm::cos(obs.catalog_dec.radians());
        w[2 * i] = cos_dec * cos_dec;
        w[2 * i + 1] = 1.0;
    }
    w
}

fn build_design_matrix(
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

fn solve_weighted(a: &DMatrix<f64>, b: &DVector<f64>, w: &DVector<f64>) -> Result<DVector<f64>> {
    let sqrt_w = w.map(libm::sqrt);
    let rows = a.nrows();
    let cols = a.ncols();
    let a_w = DMatrix::from_fn(rows, cols, |i, j| a[(i, j)] * sqrt_w[i]);
    let b_w = DVector::from_fn(rows, |i, _| b[i] * sqrt_w[i]);
    let svd = a_w.svd(true, true);
    svd.solve(&b_w, 1e-10)
        .map_err(|e| Error::Fit(format!("SVD solve failed: {}", e)))
}

fn compute_sigma_free(
    a_free: &DMatrix<f64>,
    residuals: &DVector<f64>,
    w: &DVector<f64>,
    free_indices: &[usize],
    total_terms: usize,
) -> Vec<f64> {
    let n = a_free.nrows();
    let m = a_free.ncols();
    let dof = n.saturating_sub(m).max(1);
    let sqrt_w = w.map(libm::sqrt);
    let a_w = DMatrix::from_fn(n, m, |i, j| a_free[(i, j)] * sqrt_w[i]);
    let r_w = DVector::from_fn(n, |i, _| residuals[i] * sqrt_w[i]);
    let s2 = r_w.dot(&r_w) / dof as f64;
    let ata = a_w.transpose() * &a_w;
    let free_sigma = match ata.try_inverse() {
        Some(inv) => (0..m)
            .map(|j| libm::sqrt((s2 * inv[(j, j)]).abs()))
            .collect::<Vec<_>>(),
        None => vec![f64::NAN; m],
    };
    let mut sigma = vec![0.0; total_terms];
    for (fi, &idx) in free_indices.iter().enumerate() {
        sigma[idx] = free_sigma[fi];
    }
    sigma
}

pub fn compute_sky_rms(residuals: &DVector<f64>, observations: &[&Observation]) -> f64 {
    let n = observations.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum_sq = 0.0;
    for i in 0..n {
        let dh = residuals[2 * i];
        let dd = residuals[2 * i + 1];
        let cos_dec = libm::cos(observations[i].catalog_dec.radians());
        let dx = dh * cos_dec;
        sum_sq += dx * dx + dd * dd;
    }
    libm::sqrt(sum_sq / n as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::PierSide;
    use crate::terms::create_term;
    use celestial_core::Angle;

    fn make_obs(cmd_ha_arcsec: f64, act_ha_arcsec: f64, dec_deg: f64) -> Observation {
        Observation {
            catalog_ra: Angle::from_hours(0.0),
            catalog_dec: Angle::from_degrees(dec_deg),
            observed_ra: Angle::from_hours(0.0),
            observed_dec: Angle::from_degrees(dec_deg),
            lst: Angle::from_hours(0.0),
            commanded_ha: Angle::from_arcseconds(cmd_ha_arcsec),
            actual_ha: Angle::from_arcseconds(act_ha_arcsec),
            pier_side: PierSide::East,
            masked: false,
        }
    }

    #[test]
    fn fit_ih_recovers_known_coefficient() {
        let obs = [
            make_obs(0.0, 100.0, 30.0),
            make_obs(0.0, 100.0, 45.0),
            make_obs(0.0, 100.0, 60.0),
        ];
        let refs: Vec<&Observation> = obs.iter().collect();
        let terms: Vec<Box<dyn Term>> = vec![create_term("IH").unwrap()];
        let fixed = [false];
        let coeffs = [0.0];
        let result = fit_model(&refs, &terms, &fixed, &coeffs, 0.7).unwrap();
        assert_eq!(result.term_names, vec!["IH"]);
        assert!((result.coefficients[0] - (-100.0)).abs() < 1e-6);
    }

    #[test]
    fn fit_insufficient_observations() {
        let obs = [make_obs(0.0, 100.0, 30.0)];
        let refs: Vec<&Observation> = obs.iter().collect();
        let terms: Vec<Box<dyn Term>> =
            vec![create_term("IH").unwrap(), create_term("ID").unwrap()];
        let fixed = [false, false];
        let coeffs = [0.0, 0.0];
        let result = fit_model(&refs, &terms, &fixed, &coeffs, 0.7);
        assert!(result.is_err());
    }

    #[test]
    fn fit_no_terms() {
        let obs = [make_obs(0.0, 100.0, 30.0)];
        let refs: Vec<&Observation> = obs.iter().collect();
        let terms: Vec<Box<dyn Term>> = vec![];
        let fixed: [bool; 0] = [];
        let coeffs: [f64; 0] = [];
        let result = fit_model(&refs, &terms, &fixed, &coeffs, 0.7);
        assert!(result.is_err());
    }

    #[test]
    fn sky_rms_known_residuals() {
        let obs = [make_obs(0.0, 0.0, 0.0), make_obs(0.0, 0.0, 0.0)];
        let refs: Vec<&Observation> = obs.iter().collect();
        let n = obs.len();
        let mut residuals = DVector::zeros(2 * n);
        residuals[0] = 3.0;
        residuals[1] = 4.0;
        residuals[2] = 3.0;
        residuals[3] = 4.0;
        let rms = compute_sky_rms(&residuals, &refs);
        assert_eq!(rms, 5.0);
    }
}
