use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Mvet;

impl Command for Mvet {
    fn name(&self) -> &str {
        "MVET"
    }
    
    fn description(&self) -> &str {
        "Find and optionally remove weak terms"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse(
                "MVET requires a significance threshold".into(),
            ));
        }
        let threshold: f64 = args[0]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid threshold: {}", e)))?;
        let remove = args.get(1).is_some_and(|a| a.eq_ignore_ascii_case("R"));

        let fit = session
            .last_fit
            .as_ref()
            .ok_or_else(|| Error::Fit("no fit results available (run FIT first)".into()))?;

        let mut weak: Vec<(String, f64, f64, f64)> = Vec::new();
        for (i, name) in fit.term_names.iter().enumerate() {
            let coeff = fit.coefficients[i];
            let sigma = fit.sigma[i];
            if sigma > 0.0 {
                let significance = (coeff / sigma).abs();
                if significance < threshold {
                    weak.push((name.clone(), coeff, sigma, significance));
                }
            }
        }

        if weak.is_empty() {
            return Ok(CommandOutput::Text(format!(
                "No weak terms (all significance >= {:.1})",
                threshold
            )));
        }

        let mut output = format!("Weak terms (significance < {:.1}):\n", threshold);
        for (name, coeff, sigma, sig) in &weak {
            output += &format!(
                "  {}:  coeff={:.1}  sigma={:.1}  sig={:.2}\n",
                name, coeff, sigma, sig
            );
        }

        if remove {
            for (name, _, _, _) in &weak {
                session.model.remove_term(name);
            }
            session.last_fit = None;
            output += &format!("\nRemoved {} terms", weak.len());
        } else {
            output += &format!("\nUse MVET {:.1} R to remove", threshold);
        }

        Ok(CommandOutput::Text(output))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::FitResult;
    use crate::test_support::FitResultBuilder;

    fn fit(term_names: Vec<&str>, coeffs: Vec<f64>, sigma: Vec<f64>) -> FitResult {
        FitResultBuilder::new()
            .coefficients(coeffs)
            .sigma(sigma)
            .sky_rms(1.0)
            .popn_sd(1.0)
            .term_names(term_names)
            .build()
    }

    fn session_with_fit(terms: &[&str], coeffs: Vec<f64>, sigma: Vec<f64>) -> Session {
        let mut s = Session::new();
        for t in terms {
            s.model.add_term(t).unwrap();
        }
        s.last_fit = Some(fit(terms.to_vec(), coeffs, sigma));
        s
    }

    fn text(out: CommandOutput) -> String {
        match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    fn parse_msg(err: &Error) -> &str {
        match err {
            Error::Parse(m) => m,
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Mvet.name(), "MVET");
        assert_eq!(Mvet.description(), "Find and optionally remove weak terms");
    }

    #[test]
    fn no_args_errors() {
        let mut s = Session::new();
        let err = Mvet.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("MVET requires"));
    }

    #[test]
    fn unparseable_threshold_errors() {
        let mut s = Session::new();
        let err = Mvet.execute(&mut s, &["banana"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid threshold"));
    }

    #[test]
    fn missing_fit_errors() {
        let mut s = Session::new();
        let err = Mvet.execute(&mut s, &["1.0"]).unwrap_err();
        match err {
            Error::Fit(m) => assert!(m.contains("no fit results")),
            other => panic!("expected Fit, got {:?}", other),
        }
    }

    #[test]
    fn no_weak_terms_reports_clean() {
        // 100/1 = 100, well above threshold 2.0
        let mut s = session_with_fit(&["IH"], vec![100.0], vec![1.0]);
        let body = text(Mvet.execute(&mut s, &["2.0"]).unwrap());
        assert!(body.contains("No weak terms (all significance >= 2.0)"));
        assert_eq!(s.model.term_count(), 1, "no removal without R flag");
        assert!(s.last_fit.is_some(), "last_fit preserved when no removal");
    }

    #[test]
    fn weak_terms_listed_without_r_flag_does_not_remove() {
        // IH: 1/1 = 1.0 (weak vs threshold 2.0). ID: 100/1 = 100 (strong).
        let mut s = session_with_fit(&["IH", "ID"], vec![1.0, 100.0], vec![1.0, 1.0]);
        let body = text(Mvet.execute(&mut s, &["2.0"]).unwrap());
        assert!(body.contains("Weak terms (significance < 2.0)"));
        assert!(body.contains("IH:"));
        assert!(!body.contains("ID:"));
        assert!(body.contains("Use MVET 2.0 R to remove"));
        assert_eq!(s.model.term_names(), vec!["IH", "ID"]);
        assert!(s.last_fit.is_some());
    }

    #[test]
    fn r_flag_removes_weak_terms_and_clears_last_fit() {
        let mut s = session_with_fit(
            &["IH", "ID", "CH"],
            vec![1.0, 100.0, 0.5],
            vec![1.0, 1.0, 1.0],
        );
        let body = text(Mvet.execute(&mut s, &["2.0", "R"]).unwrap());
        assert!(body.contains("Removed 2 terms"));
        assert_eq!(s.model.term_names(), vec!["ID"]);
        assert!(s.last_fit.is_none());
    }

    #[test]
    fn r_flag_is_case_insensitive() {
        for token in ["R", "r"] {
            let mut s = session_with_fit(&["IH"], vec![0.5], vec![1.0]);
            text(Mvet.execute(&mut s, &["2.0", token]).unwrap());
            assert_eq!(s.model.term_count(), 0, "token {:?}", token);
        }
    }

    // Non-R second arg is ignored — the boolean is just `eq_ignore_ascii_case("R")`.
    #[test]
    fn non_r_second_arg_does_not_trigger_removal() {
        let mut s = session_with_fit(&["IH"], vec![0.5], vec![1.0]);
        let body = text(Mvet.execute(&mut s, &["2.0", "remove"]).unwrap());
        assert!(body.contains("Use MVET"));
        assert_eq!(s.model.term_count(), 1);
    }

    // Terms with sigma <= 0 are skipped entirely (fixed/degenerate parameters).
    #[test]
    fn zero_sigma_terms_are_skipped() {
        let mut s = session_with_fit(&["IH", "ID"], vec![0.5, 0.5], vec![0.0, 1.0]);
        let body = text(Mvet.execute(&mut s, &["2.0"]).unwrap());
        assert!(!body.contains("IH:"));
        assert!(body.contains("ID:"));
    }

    #[test]
    fn negative_coefficient_uses_absolute_significance() {
        // -100/1 = -100; |sig| = 100. Not weak vs threshold 2.0.
        let mut s = session_with_fit(&["IH"], vec![-100.0], vec![1.0]);
        let body = text(Mvet.execute(&mut s, &["2.0"]).unwrap());
        assert!(body.contains("No weak terms"));
    }

    #[test]
    fn empty_model_with_empty_fit_reports_no_weak_terms() {
        let mut s = Session::new();
        s.last_fit = Some(fit(vec![], vec![], vec![]));
        let body = text(Mvet.execute(&mut s, &["2.0"]).unwrap());
        assert!(body.contains("No weak terms"));
    }
}
