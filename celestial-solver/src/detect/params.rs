use super::DetectionParams;
use crate::error::SolveError;

// Above this the structure filter is 2^10 + 1 = 1025 px wide, with a Gaussian kernel of
// 1027 taps per pass: far wider than any star and tens of seconds per frame. At 32 the
// filter-size shift overflows.
const MAX_STRUCTURE_LAYERS: u32 = 10;

pub(super) fn check_params(p: &DetectionParams) -> Result<(), SolveError> {
    let non_negative = |v: f64| v.is_finite() && v >= 0.0;
    let checks = [
        (
            "structure_layers",
            p.structure_layers <= MAX_STRUCTURE_LAYERS,
            "at most 10",
        ),
        (
            "min_structure_size",
            p.min_structure_size <= p.max_structure_size,
            "at most max_structure_size",
        ),
        (
            "min_coverage",
            (0.0..=1.0).contains(&p.min_coverage),
            "from 0 to 1",
        ),
        ("min_snr", non_negative(p.min_snr), "finite and at least 0"),
        (
            "min_spacing",
            non_negative(p.min_spacing),
            "finite and at least 0",
        ),
        // +inf is allowed and means that no pixel counts as saturated.
        (
            "saturation_adu",
            p.saturation_adu.is_none_or(|v| v >= 0.0),
            "at least 0, or +inf for no limit",
        ),
    ];
    match checks.iter().find(|(_, ok, _)| !ok) {
        Some(&(name, _, requirement)) => Err(SolveError::InvalidParameter {
            name,
            reason: format!("must be {requirement}"),
        }),
        None => Ok(()),
    }
}
