use std::collections::HashMap;

use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::spherical::Projection;

use super::keyword::WcsKeyword;

// Paper II §2.5: the native point (φ0, θ0) whose celestial coordinates CRVAL
// gives, and, when asked for, its projection, which is taken off (x, y) so
// that the point falls at CRPIX.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Fiducial {
    phi_0: f64,
    theta_0: f64,
    offset: Option<(f64, f64)>,
}

impl Fiducial {
    pub(super) fn of(projection: &Projection) -> Self {
        let (phi_0, theta_0) = projection.native_reference();
        Self {
            phi_0,
            theta_0,
            offset: None,
        }
    }

    // PV1_1 and PV1_2 on the longitude axis replace the projection's (φ0, θ0),
    // and a PV1_0 other than zero asks for the offset.
    pub(super) fn from_pv(projection: &Projection, pv: &HashMap<(u8, u8), f64>) -> WcsResult<Self> {
        let default = Self::of(projection);
        let phi_0 = pv.get(&(1, 1)).copied().unwrap_or(default.phi_0);
        let theta_0 = pv.get(&(1, 2)).copied().unwrap_or(default.theta_0);
        let offset = match pv.get(&(1, 0)).copied() {
            None | Some(0.0) => None,
            Some(flag) if flag.is_nan() => {
                return Err(WcsError::invalid_parameter("PV1_0 is NaN"));
            }
            Some(_) => {
                let native = NativeCoord::new(phi_0, theta_0);
                let at = projection.project(native)?;
                Some((at.x_deg(), at.y_deg()))
            }
        };
        Ok(Self {
            phi_0,
            theta_0,
            offset,
        })
    }

    pub(super) fn phi_0(&self) -> f64 {
        self.phi_0
    }

    pub(super) fn theta_0(&self) -> f64 {
        self.theta_0
    }

    pub(super) fn add_offset(&self, inter: IntermediateCoord) -> IntermediateCoord {
        match self.offset {
            Some((x0, y0)) => IntermediateCoord::new(inter.x_deg() + x0, inter.y_deg() + y0),
            None => inter,
        }
    }

    pub(super) fn remove_offset(&self, plane: IntermediateCoord) -> IntermediateCoord {
        match self.offset {
            Some((x0, y0)) => IntermediateCoord::new(plane.x_deg() - x0, plane.y_deg() - y0),
            None => plane,
        }
    }

    pub(super) fn keywords(&self, projection: &Projection) -> Vec<WcsKeyword> {
        let mut keywords = Vec::new();
        if self.offset.is_some() {
            keywords.push(WcsKeyword::real("PV1_0", 1.0));
        }
        if (self.phi_0, self.theta_0) != projection.native_reference() {
            keywords.push(WcsKeyword::real("PV1_1", self.phi_0));
            keywords.push(WcsKeyword::real("PV1_2", self.theta_0));
        }
        keywords
    }
}
