use crate::error::{Error, Result};
use crate::observation::PierSide;
use crate::terms::{create_term, Term};
use celestial_core::Angle;

#[derive(Default)]
pub struct PointingModel {
    terms: Vec<Box<dyn Term>>,
    coefficients: Vec<f64>,
    fixed: Vec<bool>,
    parallel: Vec<bool>,
}

impl PointingModel {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(&mut self, name: &str) -> Result<()> {
        let term = create_term(name)?;
        if self.terms.iter().any(|t| t.name() == term.name()) {
            return Ok(());
        }
        self.terms.push(term);
        self.coefficients.push(0.0);
        self.fixed.push(false);
        self.parallel.push(false);
        Ok(())
    }

    pub fn remove_term(&mut self, name: &str) {
        if let Some(idx) = self.terms.iter().position(|t| t.name() == name) {
            self.terms.remove(idx);
            self.coefficients.remove(idx);
            self.fixed.remove(idx);
            self.parallel.remove(idx);
        }
    }

    pub fn remove_all(&mut self) {
        self.terms.clear();
        self.coefficients.clear();
        self.fixed.clear();
        self.parallel.clear();
    }

    pub fn fix_term(&mut self, name: &str) -> bool {
        if let Some(idx) = self.terms.iter().position(|t| t.name() == name) {
            self.fixed[idx] = true;
            return true;
        }
        false
    }

    pub fn fix_all(&mut self) {
        self.fixed.iter_mut().for_each(|f| *f = true);
    }

    pub fn unfix_term(&mut self, name: &str) -> bool {
        if let Some(idx) = self.terms.iter().position(|t| t.name() == name) {
            self.fixed[idx] = false;
            return true;
        }
        false
    }

    pub fn unfix_all(&mut self) {
        self.fixed.iter_mut().for_each(|f| *f = false);
    }

    pub fn is_fixed(&self, idx: usize) -> bool {
        self.fixed.get(idx).copied().unwrap_or(false)
    }

    pub fn fixed_flags(&self) -> &[bool] {
        &self.fixed
    }

    pub fn set_parallel(&mut self, name: &str) -> bool {
        if let Some(idx) = self.terms.iter().position(|t| t.name() == name) {
            self.parallel[idx] = true;
            return true;
        }
        false
    }

    pub fn set_chained(&mut self, name: &str) -> bool {
        if let Some(idx) = self.terms.iter().position(|t| t.name() == name) {
            self.parallel[idx] = false;
            return true;
        }
        false
    }

    pub fn set_all_parallel(&mut self) {
        self.parallel.iter_mut().for_each(|p| *p = true);
    }

    pub fn set_all_chained(&mut self) {
        self.parallel.iter_mut().for_each(|p| *p = false);
    }

    pub fn is_parallel(&self, idx: usize) -> bool {
        self.parallel.get(idx).copied().unwrap_or(true)
    }

    pub fn parallel_flags(&self) -> &[bool] {
        &self.parallel
    }

    pub fn zero_coefficients(&mut self) {
        self.coefficients.iter_mut().for_each(|c| *c = 0.0);
    }

    pub fn term_count(&self) -> usize {
        self.terms.len()
    }

    pub fn term_names(&self) -> Vec<&str> {
        self.terms.iter().map(|t| t.name()).collect()
    }

    pub fn terms(&self) -> &[Box<dyn Term>] {
        &self.terms
    }

    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    pub fn set_coefficients(&mut self, coeffs: &[f64]) -> Result<()> {
        if coeffs.len() != self.terms.len() {
            return Err(Error::Fit(format!(
                "coefficient count {} does not match term count {}",
                coeffs.len(),
                self.terms.len()
            )));
        }
        self.coefficients.copy_from_slice(coeffs);
        Ok(())
    }

    pub fn apply_equatorial(&self, h: f64, dec: f64, lat: f64, pier: f64) -> (f64, f64) {
        let mut dh = 0.0;
        let mut ddec = 0.0;
        for (term, &coeff) in self.terms.iter().zip(self.coefficients.iter()) {
            let (jh, jd) = term.jacobian_equatorial(h, dec, lat, pier);
            dh += coeff * jh;
            ddec += coeff * jd;
        }
        (dh, ddec)
    }

    pub fn apply_altaz(&self, az: f64, el: f64, lat: f64) -> (f64, f64) {
        let mut daz = 0.0;
        let mut del = 0.0;
        for (term, &coeff) in self.terms.iter().zip(self.coefficients.iter()) {
            let (ja, je) = term.jacobian_altaz(az, el, lat);
            daz += coeff * ja;
            del += coeff * je;
        }
        (daz, del)
    }

    pub fn apply_equatorial_chained(&self, h: f64, dec: f64, lat: f64, pier: f64) -> (f64, f64) {
        let mut h_corr = h;
        let mut dec_corr = dec;

        for (i, term) in self.terms.iter().enumerate() {
            if !self.parallel[i] {
                let (jh, jd) = term.jacobian_equatorial(h_corr, dec_corr, lat, pier);
                h_corr += self.coefficients[i] * jh;
                dec_corr += self.coefficients[i] * jd;
            }
        }

        let mut dh = 0.0;
        let mut ddec = 0.0;
        for (i, term) in self.terms.iter().enumerate() {
            if self.parallel[i] {
                let (jh, jd) = term.jacobian_equatorial(h_corr, dec_corr, lat, pier);
                dh += self.coefficients[i] * jh;
                ddec += self.coefficients[i] * jd;
            }
        }

        (h_corr + dh - h, dec_corr + ddec - dec)
    }

    pub fn target_to_command(
        &self,
        ra: Angle,
        dec: Angle,
        lst: Angle,
        lat: Angle,
        pier: PierSide,
    ) -> (Angle, Angle) {
        let ha = lst - ra;
        let (dh, dd) =
            self.apply_equatorial(ha.radians(), dec.radians(), lat.radians(), pier.sign());
        let cmd_ha = ha + Angle::from_arcseconds(dh);
        let cmd_dec = dec + Angle::from_arcseconds(dd);
        let cmd_ra = lst - cmd_ha;
        (cmd_ra, cmd_dec)
    }

    pub fn command_to_target(
        &self,
        ra_encoder: Angle,
        dec_encoder: Angle,
        lst: Angle,
        lat: Angle,
        pier: PierSide,
    ) -> (Angle, Angle) {
        let ha = lst - ra_encoder;
        let (dh, dd) = self.apply_equatorial(
            ha.radians(),
            dec_encoder.radians(),
            lat.radians(),
            pier.sign(),
        );
        let true_ha = ha - Angle::from_arcseconds(dh);
        let true_dec = dec_encoder - Angle::from_arcseconds(dd);
        let true_ra = lst - true_ha;
        (true_ra, true_dec)
    }

    pub fn predict_breakdown(
        &self,
        h: f64,
        dec: f64,
        lat: f64,
        pier: f64,
    ) -> Vec<(String, f64, f64)> {
        self.terms
            .iter()
            .zip(self.coefficients.iter())
            .map(|(term, &coeff)| {
                let (jh, jd) = term.jacobian_equatorial(h, dec, lat, pier);
                (term.name().to_string(), coeff * jh, coeff * jd)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_4;

    #[test]
    fn apply_equatorial_ih_id() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.add_term("ID").unwrap();
        model.set_coefficients(&[10.0, 20.0]).unwrap();

        let (dh, ddec) = model.apply_equatorial(FRAC_PI_4, 0.5, 0.7, 1.0);
        assert_eq!(dh, -10.0);
        assert_eq!(ddec, -20.0);
    }

    #[test]
    fn apply_equatorial_id_west_pier() {
        let mut model = PointingModel::new();
        model.add_term("ID").unwrap();
        model.set_coefficients(&[20.0]).unwrap();

        let (_, ddec) = model.apply_equatorial(0.0, 0.0, 0.0, -1.0);
        assert_eq!(ddec, 20.0 * 1.0);
    }

    #[test]
    fn add_remove_terms() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.add_term("ID").unwrap();
        model.add_term("CH").unwrap();
        assert_eq!(model.term_count(), 3);
        assert_eq!(model.term_names(), vec!["IH", "ID", "CH"]);

        model.remove_term("ID");
        assert_eq!(model.term_count(), 2);
        assert_eq!(model.term_names(), vec!["IH", "CH"]);
        assert_eq!(model.coefficients().len(), 2);
    }

    #[test]
    fn remove_all_clears_model() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.add_term("ID").unwrap();
        model.remove_all();
        assert_eq!(model.term_count(), 0);
        assert_eq!(model.coefficients().len(), 0);
    }

    #[test]
    fn remove_nonexistent_is_noop() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.remove_term("ZZZZ");
        assert_eq!(model.term_count(), 1);
    }

    #[test]
    fn set_coefficients_wrong_length() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        let result = model.set_coefficients(&[1.0, 2.0]);
        assert!(result.is_err());
    }

    #[test]
    fn empty_model_returns_zero_correction() {
        let model = PointingModel::new();
        let (dh, ddec) = model.apply_equatorial(1.0, 0.5, 0.7, 1.0);
        assert_eq!(dh, 0.0);
        assert_eq!(ddec, 0.0);
    }

    #[test]
    fn add_unknown_term_returns_error() {
        let mut model = PointingModel::new();
        let result = model.add_term("ZZZZ");
        assert!(result.is_err());
    }

    // Sign conventions anchored to the fit residual definition:
    //   apply_equatorial(true_pos) = encoder − true_pos       (the model predicts the residual)
    // Therefore:
    //   target_to_command:  encoder_cmd = target + apply(target)
    //   command_to_target:  true_pos    = encoder − apply(encoder)
    // TPOINT manual section 7.2 equivalents:
    //   "ADD corrections to encoder to get true sky"     ⇒ correction_manual = −apply
    //   "SUBTRACT corrections from target to command"    ⇒ command = target + apply
    //
    // The tests below are the *only* place this convention is locked down. If you
    // change the residual sign in `build_residuals` or flip any per-term Jacobian
    // sign, expect at least one of these to fail.

    #[test]
    fn target_to_command_matches_residual_convention() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.set_coefficients(&[-100.0]).unwrap();
        let lst = Angle::from_hours(6.0);
        let lat = Angle::from_degrees(39.0);
        let target_ra = Angle::from_hours(6.0);
        let target_dec = Angle::from_degrees(30.0);
        let ha = lst - target_ra;
        let (dh, _) = model.apply_equatorial(
            ha.radians(),
            target_dec.radians(),
            lat.radians(),
            PierSide::East.sign(),
        );
        assert!((dh - 100.0).abs() < 1e-9, "apply.dh = {} for IH=-100", dh);
        let (cmd_ra, _) = model.target_to_command(target_ra, target_dec, lst, lat, PierSide::East);
        let cmd_ha = lst - cmd_ra;
        let actual_offset = (cmd_ha - ha).arcseconds();
        assert!(
            (actual_offset - dh).abs() < 1e-6,
            "target_to_command should produce cmd_ha = ha + apply.dh; offset={}, expected {}",
            actual_offset,
            dh,
        );
    }

    #[test]
    fn command_to_target_matches_residual_convention() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.set_coefficients(&[-100.0]).unwrap();
        let lst = Angle::from_hours(6.0);
        let lat = Angle::from_degrees(39.0);
        let target_ra = Angle::from_hours(6.0);
        let target_dec = Angle::from_degrees(30.0);
        let ha = lst - target_ra;
        let (dh, _) = model.apply_equatorial(
            ha.radians(),
            target_dec.radians(),
            lat.radians(),
            PierSide::East.sign(),
        );
        assert!((dh - 100.0).abs() < 1e-9, "apply.dh = {} for IH=-100", dh);
        let encoder_ra = lst - (ha + Angle::from_arcseconds(dh));
        let encoder_dec = target_dec;
        let (true_ra, true_dec) =
            model.command_to_target(encoder_ra, encoder_dec, lst, lat, PierSide::East);
        let dra = (true_ra - target_ra).arcseconds();
        let ddec = (true_dec - target_dec).arcseconds();
        assert!(
            dra.abs() < 1e-6 && ddec.abs() < 1e-6,
            "command_to_target failed: dra={}, ddec={}",
            dra,
            ddec,
        );
    }

    #[test]
    fn target_command_target_roundtrip() {
        let mut model = PointingModel::new();
        model.add_term("IH").unwrap();
        model.add_term("ID").unwrap();
        model.set_coefficients(&[-100.0, 50.0]).unwrap();
        let lst = Angle::from_hours(6.0);
        let lat = Angle::from_degrees(39.0);
        let target_ra = Angle::from_hours(6.0);
        let target_dec = Angle::from_degrees(30.0);
        let (cmd_ra, cmd_dec) =
            model.target_to_command(target_ra, target_dec, lst, lat, PierSide::East);
        let (back_ra, back_dec) =
            model.command_to_target(cmd_ra, cmd_dec, lst, lat, PierSide::East);
        let dra = (back_ra - target_ra).arcseconds();
        let ddec = (back_dec - target_dec).arcseconds();
        assert!(
            dra.abs() < 1e-6 && ddec.abs() < 1e-6,
            "round-trip failed: dra={}, ddec={}",
            dra,
            ddec,
        );
    }
}
