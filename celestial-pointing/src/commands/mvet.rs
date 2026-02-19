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
