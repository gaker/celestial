// The TAN parameters and residuals of the model the report prints.

use super::super::solution::{StarResidual, WcsSolution};
use crate::solve::SolveResult;

// The TAN parameters of the projection the report prints, with the residuals of the
// fitted pairs under the whole model: those of the SIP model, which its A and B
// polynomials need, or the TAN WCS when there is no SIP.
pub(super) fn projection_model(result: &SolveResult) -> WcsSolution {
    let Some(sip) = &result.sip else {
        return result.wcs.clone();
    };
    WcsSolution {
        crval1: sip.crval1,
        crval2: sip.crval2,
        cd1_1: sip.cd1_1,
        cd1_2: sip.cd1_2,
        cd2_1: sip.cd2_1,
        cd2_2: sip.cd2_2,
        residuals: sip_residuals(result),
        ..result.wcs.clone()
    }
}

// The residuals of the pairs `result.wcs` was fitted on, with each pixel predicted by
// the FITS WCS that includes the SIP distortion. A pair with no predicted pixel is left
// out.
fn sip_residuals(result: &SolveResult) -> Vec<StarResidual> {
    let residual = |r: &StarResidual| {
        let (x, y) = result.sky_to_pixel(r.ra_deg, r.dec_deg).ok()?;
        let (err_x, err_y) = (r.px_x - x, r.px_y - y);
        Some(StarResidual {
            err_x,
            err_y,
            err_px: libm::hypot(err_x, err_y),
            ..r.clone()
        })
    };
    result.wcs.residuals.iter().filter_map(residual).collect()
}
