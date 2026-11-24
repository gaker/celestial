//! The FITS WCS of a solution: TAN, or TAN-SIP with the inverse polynomials.

use celestial_wcs::builder::{Wcs, WcsBuilder};
use celestial_wcs::distortion::sip::SipDistortion;
use celestial_wcs::distortion::DistortionModel;
use celestial_wcs::spherical::Projection;

use crate::error::SolveError;
use crate::pixel_frame::detection_to_fits;

use super::sip::SipSolution;
use super::sip_inverse::{fit_inverse, inverse_order};
use super::solution::WcsSolution;

// `wcs` in FITS pixel coordinates, or with `sip` the TAN parameters of `sip` with its
// distortion and AP and BP fitted to it. An error when the parameters do not build a
// WCS or the inverse fit fails.
pub(crate) fn build_fits_wcs(
    wcs: &WcsSolution,
    sip: Option<&SipSolution>,
) -> Result<Wcs, SolveError> {
    let crpix = detection_to_fits((wcs.crpix1, wcs.crpix2));
    let (crval, cd) = match sip {
        Some(s) => (
            (s.crval1, s.crval2),
            [[s.cd1_1, s.cd1_2], [s.cd2_1, s.cd2_2]],
        ),
        None => (
            (wcs.crval1, wcs.crval2),
            [[wcs.cd1_1, wcs.cd1_2], [wcs.cd2_1, wcs.cd2_2]],
        ),
    };
    let mut builder = WcsBuilder::new()
        .crpix(crpix.0, crpix.1)
        .crval(crval.0, crval.1)
        .cd_matrix(cd)
        .projection(Projection::tan());
    if let Some(sip) = sip {
        builder = builder.distortion(DistortionModel::Sip(sip_distortion(sip, crpix)?));
    }
    Ok(builder.build()?)
}

fn sip_distortion(sip: &SipSolution, crpix: (f64, f64)) -> Result<SipDistortion, SolveError> {
    let mut distortion = SipDistortion::new([crpix.0, crpix.1], sip.a_order, sip.b_order);
    for (&(p, q), &c) in &sip.a_coeffs {
        distortion.set_a(p, q, c)?;
    }
    for (&(p, q), &c) in &sip.b_coeffs {
        distortion.set_b(p, q, c)?;
    }
    let order = inverse_order(sip);
    let (ap, bp) = fit_inverse(sip, crpix, order).map_err(|e| {
        SolveError::fit_failed(format!("inverse SIP polynomials of order {order}: {e}"))
    })?;
    distortion.set_inverse_order(order, order)?;
    for ((p, q), c) in ap {
        distortion.set_ap(p, q, c)?;
    }
    for ((p, q), c) in bp {
        distortion.set_bp(p, q, c)?;
    }
    Ok(distortion)
}
