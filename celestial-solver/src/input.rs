//! Checks of the inputs of the public entry points. A value that is not finite or is
//! out of its range is `SolveError::InvalidParameter` before any work, so it cannot
//! reach the catalog search or a fit as NaN.

use celestial_time::julian::JulianDate;

use crate::detect::DetectedStar;
use crate::error::SolveError;
use crate::fit_wcs::solution::WcsSolution;
use crate::match_field::{MatchParams, StarPair};

fn invalid(name: &'static str, reason: String) -> SolveError {
    SolveError::InvalidParameter { name, reason }
}

pub(crate) fn finite(name: &'static str, value: f64) -> Result<(), SolveError> {
    if value.is_finite() {
        return Ok(());
    }
    Err(invalid(name, format!("must be finite, got {value}")))
}

pub(crate) fn positive(name: &'static str, value: f64) -> Result<(), SolveError> {
    if value.is_finite() && value > 0.0 {
        return Ok(());
    }
    Err(invalid(
        name,
        format!("must be finite and positive, got {value}"),
    ))
}

pub(crate) fn positive_if_set(name: &'static str, value: Option<f64>) -> Result<(), SolveError> {
    value.map_or(Ok(()), |v| positive(name, v))
}

// A cone search radius in degrees: positive, and at most 180 since a cone cannot be
// larger than the sphere.
pub(crate) fn search_radius(name: &'static str, value: Option<f64>) -> Result<(), SolveError> {
    let Some(radius) = value else {
        return Ok(());
    };
    positive(name, radius)?;
    if radius > 180.0 {
        return Err(invalid(name, format!("must be at most 180, got {radius}")));
    }
    Ok(())
}

pub(crate) fn declination(name: &'static str, value: f64) -> Result<(), SolveError> {
    if (-90.0..=90.0).contains(&value) {
        return Ok(());
    }
    Err(invalid(
        name,
        format!("must be from -90 to 90, got {value}"),
    ))
}

pub(crate) fn epoch(epoch: JulianDate) -> Result<(), SolveError> {
    let (jd1, jd2) = (epoch.jd1(), epoch.jd2());
    if jd1.is_finite() && jd2.is_finite() {
        return Ok(());
    }
    Err(invalid(
        "epoch",
        format!("must be finite, got ({jd1}, {jd2})"),
    ))
}

pub(crate) fn stars(stars: &[DetectedStar]) -> Result<(), SolveError> {
    let finite = |s: &DetectedStar| [s.x, s.y, s.snr].iter().all(|v| v.is_finite());
    match stars.iter().position(|s| !finite(s)) {
        None => Ok(()),
        Some(i) => {
            let s = &stars[i];
            let reason = format!(
                "star {i} is not finite: pixel ({}, {}), snr {}",
                s.x, s.y, s.snr
            );
            Err(invalid("stars", reason))
        }
    }
}

pub(crate) fn pairs(pairs: &[StarPair]) -> Result<(), SolveError> {
    let finite = |p: &StarPair| {
        let values = [p.px_x, p.px_y, p.ra_deg, p.dec_deg, p.snr];
        values.iter().all(|v| v.is_finite())
    };
    match pairs.iter().position(|p| !finite(p)) {
        None => Ok(()),
        Some(i) => {
            let p = &pairs[i];
            let reason = format!(
                "pair {i} is not finite: pixel ({}, {}), sky ({}, {}), snr {}",
                p.px_x, p.px_y, p.ra_deg, p.dec_deg, p.snr
            );
            Err(invalid("pairs", reason))
        }
    }
}

// A quad is a star and three of its nearest neighbours, so fewer than 4 stars or 3
// neighbours form none, and every match would fail.
pub(crate) fn match_params(params: &MatchParams) -> Result<(), SolveError> {
    if params.max_stars < 4 {
        let reason = format!(
            "must be at least 4, the stars of one quad, got {}",
            params.max_stars
        );
        return Err(invalid("max_stars", reason));
    }
    if params.k_neighbors < 3 {
        let reason = format!(
            "must be at least 3, the neighbours a quad adds to its first star, got {}",
            params.k_neighbors
        );
        return Err(invalid("k_neighbors", reason));
    }
    search_radius("search_radius_deg", params.search_radius_deg)
}

// The CRPIX, CRVAL and CD matrix of `wcs` are finite, CRVAL2 is a declination, and the CD
// matrix has an inverse.
pub(crate) fn wcs(name: &'static str, wcs: &WcsSolution) -> Result<(), SolveError> {
    let fields = [
        ("crpix1", wcs.crpix1),
        ("crpix2", wcs.crpix2),
        ("crval1", wcs.crval1),
        ("crval2", wcs.crval2),
        ("cd1_1", wcs.cd1_1),
        ("cd1_2", wcs.cd1_2),
        ("cd2_1", wcs.cd2_1),
        ("cd2_2", wcs.cd2_2),
    ];
    if let Some((field, value)) = fields.iter().find(|(_, v)| !v.is_finite()) {
        return Err(invalid(
            name,
            format!("{field} must be finite, got {value}"),
        ));
    }
    if !(-90.0..=90.0).contains(&wcs.crval2) {
        let reason = format!("crval2 must be from -90 to 90, got {}", wcs.crval2);
        return Err(invalid(name, reason));
    }
    let det = wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1;
    if det == 0.0 {
        let reason = format!("the CD matrix must be invertible, got determinant {det}");
        return Err(invalid(name, reason));
    }
    Ok(())
}

// A SIP fit needs order 2 or above, and the FITS reader of celestial_wcs reads A_ORDER
// and B_ORDER only up to 9, so a header with a higher order could not be read back.
pub(crate) fn sip_order(name: &'static str, sip_order: Option<u32>) -> Result<(), SolveError> {
    match sip_order {
        Some(order) if !(MIN_SIP_ORDER..=MAX_SIP_ORDER).contains(&order) => {
            let reason = format!("must be from {MIN_SIP_ORDER} to {MAX_SIP_ORDER}, got {order}");
            Err(invalid(name, reason))
        }
        _ => Ok(()),
    }
}

const MIN_SIP_ORDER: u32 = 2;
const MAX_SIP_ORDER: u32 = 9;
