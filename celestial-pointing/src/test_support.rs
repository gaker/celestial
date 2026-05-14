//! Shared fixtures for unit tests across the crate.
//!
//! Provides builders for `Observation` and `FitResult` so individual test
//! modules don't each copy-paste the boilerplate. Compiled only under `cfg(test)`.

use crate::observation::{Observation, PierSide};
use crate::solver::{FitResult, IterReport, ObsDiagnostic, RankInfo};
use celestial_core::Angle;

/// Build an `Observation` with overridable fields.
///
/// Defaults: every angular value is zero, pier east, unmasked. Tests override
/// only the fields they care about.
#[derive(Clone)]
pub(crate) struct ObsBuilder {
    catalog_ra: Angle,
    catalog_dec: Angle,
    observed_ra: Angle,
    observed_dec: Angle,
    lst: Angle,
    commanded_ha: Angle,
    actual_ha: Angle,
    pier_side: PierSide,
    masked: bool,
}

impl ObsBuilder {
    pub(crate) fn new() -> Self {
        Self {
            catalog_ra: Angle::from_hours(0.0),
            catalog_dec: Angle::from_degrees(0.0),
            observed_ra: Angle::from_hours(0.0),
            observed_dec: Angle::from_degrees(0.0),
            lst: Angle::from_hours(0.0),
            commanded_ha: Angle::from_arcseconds(0.0),
            actual_ha: Angle::from_arcseconds(0.0),
            pier_side: PierSide::East,
            masked: false,
        }
    }

    pub(crate) fn catalog_ra_hours(mut self, h: f64) -> Self {
        self.catalog_ra = Angle::from_hours(h);
        self
    }

    pub(crate) fn catalog_dec_deg(mut self, d: f64) -> Self {
        self.catalog_dec = Angle::from_degrees(d);
        self
    }

    pub(crate) fn observed_dec_deg(mut self, d: f64) -> Self {
        self.observed_dec = Angle::from_degrees(d);
        self
    }

    pub(crate) fn lst_hours(mut self, h: f64) -> Self {
        self.lst = Angle::from_hours(h);
        self
    }

    pub(crate) fn commanded_ha_arcsec(mut self, s: f64) -> Self {
        self.commanded_ha = Angle::from_arcseconds(s);
        self
    }

    pub(crate) fn actual_ha_arcsec(mut self, s: f64) -> Self {
        self.actual_ha = Angle::from_arcseconds(s);
        self
    }

    pub(crate) fn pier(mut self, side: PierSide) -> Self {
        self.pier_side = side;
        self
    }

    pub(crate) fn masked(mut self, m: bool) -> Self {
        self.masked = m;
        self
    }

    pub(crate) fn build(self) -> Observation {
        Observation {
            catalog_ra: self.catalog_ra,
            catalog_dec: self.catalog_dec,
            observed_ra: self.observed_ra,
            observed_dec: self.observed_dec,
            lst: self.lst,
            commanded_ha: self.commanded_ha,
            actual_ha: self.actual_ha,
            pier_side: self.pier_side,
            masked: self.masked,
        }
    }
}

/// Shorthand: a fully default observation.
pub(crate) fn obs() -> Observation {
    ObsBuilder::new().build()
}

/// Shorthand: a single-coefficient observation with `actual_ha = +offset"` at the given dec.
///
/// The canonical "IH should fit to `-offset`" smoke fixture.
pub(crate) fn obs_with_ha_offset(actual_ha_arcsec: f64, dec_deg: f64) -> Observation {
    ObsBuilder::new()
        .catalog_dec_deg(dec_deg)
        .observed_dec_deg(dec_deg)
        .actual_ha_arcsec(actual_ha_arcsec)
        .build()
}

/// Build a `FitResult` with overridable fields. Sensible defaults so tests
/// only specify what they care about.
#[derive(Clone)]
pub(crate) struct FitResultBuilder {
    coefficients: Vec<f64>,
    change: Option<Vec<f64>>,
    sigma: Option<Vec<f64>>,
    sky_rms: f64,
    popn_sd: f64,
    term_names: Vec<String>,
    rank_info: RankInfo,
    diagnostics: Vec<ObsDiagnostic>,
    iter_report: IterReport,
}

impl FitResultBuilder {
    pub(crate) fn new() -> Self {
        Self {
            coefficients: Vec::new(),
            change: None,
            sigma: None,
            sky_rms: 0.0,
            popn_sd: 0.0,
            term_names: Vec::new(),
            rank_info: RankInfo::default(),
            diagnostics: Vec::new(),
            iter_report: IterReport::default(),
        }
    }

    pub(crate) fn coefficients(mut self, c: Vec<f64>) -> Self {
        self.coefficients = c;
        self
    }

    pub(crate) fn sigma(mut self, s: Vec<f64>) -> Self {
        self.sigma = Some(s);
        self
    }

    pub(crate) fn sky_rms(mut self, r: f64) -> Self {
        self.sky_rms = r;
        self
    }

    pub(crate) fn popn_sd(mut self, p: f64) -> Self {
        self.popn_sd = p;
        self
    }

    pub(crate) fn term_names<I, S>(mut self, names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.term_names = names.into_iter().map(Into::into).collect();
        self
    }

    pub(crate) fn rank_info(mut self, info: RankInfo) -> Self {
        self.rank_info = info;
        self
    }

    pub(crate) fn diagnostics(mut self, d: Vec<ObsDiagnostic>) -> Self {
        self.diagnostics = d;
        self
    }

    pub(crate) fn build(self) -> FitResult {
        let n = self.coefficients.len();
        FitResult {
            change: self.change.unwrap_or_else(|| vec![0.0; n]),
            sigma: self.sigma.unwrap_or_else(|| vec![0.1; n]),
            coefficients: self.coefficients,
            sky_rms: self.sky_rms,
            popn_sd: self.popn_sd,
            term_names: self.term_names,
            rank_info: self.rank_info,
            diagnostics: self.diagnostics,
            iter_report: self.iter_report,
        }
    }
}

/// One-liner `ObsDiagnostic` with the given sky residual.
pub(crate) fn diag(residual_sky: f64) -> ObsDiagnostic {
    ObsDiagnostic {
        residual_sky,
        leverage: 0.0,
        studentized: 0.0,
    }
}
