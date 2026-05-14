use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Mask;
pub struct Unmask;

impl Command for Mask {
    fn name(&self) -> &str {
        "MASK"
    }
    
    fn description(&self) -> &str {
        "Mask observations (exclude from fit)"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("MASK requires observation numbers".into()));
        }
        let indices = parse_obs_indices(args, session.observations.len())?;
        let mut count = 0;
        for idx in &indices {
            if !session.observations[*idx].masked {
                session.observations[*idx].masked = true;
                count += 1;
            }
        }
        Ok(CommandOutput::Text(format!(
            "Masked {} observations",
            count
        )))
    }
}

impl Command for Unmask {
    fn name(&self) -> &str {
        "UNMASK"
    }
    fn description(&self) -> &str {
        "Unmask observations (include in fit)"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse(
                "UNMASK requires observation numbers or ALL".into(),
            ));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            let count = session.observations.iter().filter(|o| o.masked).count();
            for obs in &mut session.observations {
                obs.masked = false;
            }
            return Ok(CommandOutput::Text(format!(
                "Unmasked {} observations",
                count
            )));
        }
        let indices = parse_obs_indices(args, session.observations.len())?;
        let mut count = 0;
        for idx in &indices {
            if session.observations[*idx].masked {
                session.observations[*idx].masked = false;
                count += 1;
            }
        }
        Ok(CommandOutput::Text(format!(
            "Unmasked {} observations",
            count
        )))
    }
}

fn parse_obs_indices(args: &[&str], total: usize) -> Result<Vec<usize>> {
    let mut indices = Vec::new();
    for arg in args {
        if arg.contains('-') {
            let parts: Vec<&str> = arg.splitn(2, '-').collect();
            let start: usize = parts[0]
                .parse()
                .map_err(|e| Error::Parse(format!("invalid range start: {}", e)))?;
            let end: usize = parts[1]
                .parse()
                .map_err(|e| Error::Parse(format!("invalid range end: {}", e)))?;
            if start < 1 || end < 1 || start > total || end > total {
                return Err(Error::Parse(format!(
                    "range {}-{} out of bounds (1-{})",
                    start, end, total
                )));
            }
            for i in start..=end {
                indices.push(i - 1);
            }
        } else {
            let num: usize = arg
                .parse()
                .map_err(|e| Error::Parse(format!("invalid observation number: {}", e)))?;
            if num < 1 || num > total {
                return Err(Error::Parse(format!(
                    "observation {} out of bounds (1-{})",
                    num, total
                )));
            }
            indices.push(num - 1);
        }
    }
    Ok(indices)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::obs;

    fn session_with(n: usize) -> Session {
        let mut s = Session::new();
        s.observations = (0..n).map(|_| obs()).collect();
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

    fn masked_flags(s: &Session) -> Vec<bool> {
        s.observations.iter().map(|o| o.masked).collect()
    }

    #[test]
    fn metadata() {
        assert_eq!(Mask.name(), "MASK");
        assert_eq!(Mask.description(), "Mask observations (exclude from fit)");
        assert_eq!(Unmask.name(), "UNMASK");
        assert_eq!(
            Unmask.description(),
            "Unmask observations (include in fit)"
        );
    }

    #[test]
    fn mask_no_args_errors() {
        let mut s = session_with(3);
        let err = Mask.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("MASK requires"));
        assert_eq!(masked_flags(&s), vec![false; 3]);
    }

    #[test]
    fn unmask_no_args_errors() {
        let mut s = session_with(3);
        let err = Unmask.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("UNMASK requires"));
    }

    #[test]
    fn mask_single_observation() {
        let mut s = session_with(3);
        let body = text(Mask.execute(&mut s, &["2"]).unwrap());
        assert_eq!(body, "Masked 1 observations");
        assert_eq!(masked_flags(&s), vec![false, true, false]);
    }

    #[test]
    fn mask_range_inclusive() {
        let mut s = session_with(5);
        let body = text(Mask.execute(&mut s, &["2-4"]).unwrap());
        assert_eq!(body, "Masked 3 observations");
        assert_eq!(masked_flags(&s), vec![false, true, true, true, false]);
    }

    #[test]
    fn mask_mixed_single_and_range() {
        let mut s = session_with(5);
        text(Mask.execute(&mut s, &["1", "3-4"]).unwrap());
        assert_eq!(masked_flags(&s), vec![true, false, true, true, false]);
    }

    #[test]
    fn mask_singleton_range() {
        let mut s = session_with(3);
        let body = text(Mask.execute(&mut s, &["2-2"]).unwrap());
        assert_eq!(body, "Masked 1 observations");
        assert_eq!(masked_flags(&s), vec![false, true, false]);
    }

    // Count reflects "newly masked", not "indices touched". Re-masking an
    // already-masked observation does not increment the counter.
    #[test]
    fn mask_already_masked_does_not_increment_count() {
        let mut s = session_with(3);
        s.observations[1].masked = true;
        let body = text(Mask.execute(&mut s, &["1", "2", "3"]).unwrap());
        assert_eq!(body, "Masked 2 observations");
        assert_eq!(masked_flags(&s), vec![true, true, true]);
    }

    #[test]
    fn unmask_single_observation() {
        let mut s = session_with(3);
        for o in &mut s.observations {
            o.masked = true;
        }
        let body = text(Unmask.execute(&mut s, &["2"]).unwrap());
        assert_eq!(body, "Unmasked 1 observations");
        assert_eq!(masked_flags(&s), vec![true, false, true]);
    }

    #[test]
    fn unmask_already_unmasked_does_not_increment_count() {
        let mut s = session_with(3);
        s.observations[0].masked = true;
        let body = text(Unmask.execute(&mut s, &["1", "2", "3"]).unwrap());
        assert_eq!(body, "Unmasked 1 observations");
        assert_eq!(masked_flags(&s), vec![false; 3]);
    }

    #[test]
    fn unmask_all_clears_every_mask_case_insensitive() {
        for token in ["ALL", "all", "All"] {
            let mut s = session_with(4);
            s.observations[0].masked = true;
            s.observations[2].masked = true;
            let body = text(Unmask.execute(&mut s, &[token]).unwrap());
            assert_eq!(body, "Unmasked 2 observations", "token {:?}", token);
            assert_eq!(masked_flags(&s), vec![false; 4]);
        }
    }

    #[test]
    fn unmask_all_with_no_masked_observations_reports_zero() {
        let mut s = session_with(3);
        let body = text(Unmask.execute(&mut s, &["ALL"]).unwrap());
        assert_eq!(body, "Unmasked 0 observations");
    }

    #[test]
    fn parse_obs_indices_zero_is_out_of_bounds() {
        let mut s = session_with(3);
        let err = Mask.execute(&mut s, &["0"]).unwrap_err();
        assert!(parse_msg(&err).contains("out of bounds"));
    }

    #[test]
    fn parse_obs_indices_above_total_is_out_of_bounds() {
        let mut s = session_with(3);
        let err = Mask.execute(&mut s, &["4"]).unwrap_err();
        assert!(parse_msg(&err).contains("out of bounds"));
    }

    #[test]
    fn parse_obs_indices_unparseable_single_errors() {
        let mut s = session_with(3);
        let err = Mask.execute(&mut s, &["abc"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid observation number"));
    }

    #[test]
    fn parse_obs_indices_unparseable_range_endpoints_error() {
        let mut s = session_with(3);
        let err = Mask.execute(&mut s, &["abc-2"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid range start"));
        let err = Mask.execute(&mut s, &["1-xyz"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid range end"));
    }

    #[test]
    fn parse_obs_indices_range_out_of_bounds() {
        let mut s = session_with(3);
        let err = Mask.execute(&mut s, &["0-2"]).unwrap_err();
        assert!(parse_msg(&err).contains("out of bounds"));
        let err = Mask.execute(&mut s, &["1-5"]).unwrap_err();
        assert!(parse_msg(&err).contains("out of bounds"));
    }

    // A reversed range like "4-2" silently iterates over an empty range,
    // producing zero index pushes — neither an error nor any state change.
    // Pinning this so anyone tightening the bounds check sees it.
    #[test]
    fn reversed_range_is_silent_noop() {
        let mut s = session_with(5);
        let body = text(Mask.execute(&mut s, &["4-2"]).unwrap());
        assert_eq!(body, "Masked 0 observations");
        assert_eq!(masked_flags(&s), vec![false; 5]);
    }

    // Parse errors happen before any mutation. Earlier valid args don't get
    // applied if a later one fails.
    #[test]
    fn mask_parse_failure_leaves_observations_untouched() {
        let mut s = session_with(3);
        let _ = Mask.execute(&mut s, &["1", "999"]).unwrap_err();
        assert_eq!(masked_flags(&s), vec![false; 3]);
    }

    #[test]
    fn unmask_named_indices_path_runs_for_non_all_first_arg() {
        let mut s = session_with(3);
        for o in &mut s.observations {
            o.masked = true;
        }
        // First arg is not "ALL", so falls through to the index path.
        text(Unmask.execute(&mut s, &["1", "3"]).unwrap());
        assert_eq!(masked_flags(&s), vec![false, true, false]);
    }
}
