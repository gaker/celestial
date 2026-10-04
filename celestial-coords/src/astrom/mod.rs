pub mod earth;

use crate::aberration::{
    aberrate, apply_light_deflection, compute_earth_state, observer_motion,
    remove_light_deflection, unaberrate,
};
use crate::distance::Distance;
use crate::eop::record::EopParameters;
use crate::errors::{CoordError, CoordResult};
use crate::frames::cirs::CIRSPosition;
use crate::frames::gcrs::GCRSPosition;
use crate::frames::icrs::ICRSPosition;
use celestial_core::cio::{gcrs_to_cirs_matrix, CioSolution};
use celestial_core::matrix::{RotationMatrix3, Vector3};
use celestial_core::precession::PrecessionIAU2006;
use celestial_core::utils::jd_to_centuries;
use celestial_time::scales::tt::TT;
use celestial_time::transforms::nutation::NutationCalculator;

// Everything about the Earth's position, motion and orientation that the
// ICRS, GCRS and CIRS conversions need at one epoch, computed once.
#[derive(Debug, Clone, PartialEq)]
pub struct Astrom {
    epoch: TT,
    sun_to_earth: Vector3,
    sun_distance_au: f64,
    velocity_c: Vector3,
    bm1: f64,
    gcrs_to_cirs: RotationMatrix3,
}

impl Astrom {
    pub fn new(epoch: &TT) -> CoordResult<Self> {
        Self::with_matrix(epoch, celestial_to_intermediate(epoch, None)?)
    }

    pub fn with_eop(epoch: &TT, eop: &EopParameters) -> CoordResult<Self> {
        Self::with_matrix(epoch, celestial_to_intermediate(epoch, Some(eop))?)
    }

    fn with_matrix(epoch: &TT, gcrs_to_cirs: RotationMatrix3) -> CoordResult<Self> {
        let earth = compute_earth_state(epoch)?;
        let (velocity_c, bm1) = observer_motion(earth.barycentric_velocity);
        Ok(Self {
            epoch: *epoch,
            sun_to_earth: earth.heliocentric_position.normalize()?,
            sun_distance_au: earth.heliocentric_position.magnitude(),
            velocity_c,
            bm1,
            gcrs_to_cirs,
        })
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn icrs_to_gcrs(&self, icrs: &ICRSPosition) -> CoordResult<GCRSPosition> {
        let mut gcrs = GCRSPosition::from_unit_vector(self.apparent(icrs)?, self.epoch)?;
        if let Some(distance) = icrs.distance() {
            gcrs.set_distance(distance);
        }
        Ok(gcrs)
    }

    pub fn gcrs_to_icrs(&self, gcrs: &GCRSPosition) -> CoordResult<ICRSPosition> {
        check_epoch(self.epoch, gcrs.epoch())?;
        let mut icrs = ICRSPosition::from_unit_vector(self.natural(gcrs.unit_vector())?)?;
        if let Some(distance) = gcrs.distance() {
            icrs.set_distance(distance);
        }
        Ok(icrs)
    }

    pub fn gcrs_to_cirs(&self, gcrs: &GCRSPosition) -> CoordResult<CIRSPosition> {
        check_epoch(self.epoch, gcrs.epoch())?;
        self.cirs_from_gcrs_vector(gcrs.unit_vector(), gcrs.distance())
    }

    pub fn cirs_to_gcrs(&self, cirs: &CIRSPosition) -> CoordResult<GCRSPosition> {
        check_epoch(self.epoch, cirs.epoch())?;
        let p = self.gcrs_to_cirs.transpose() * cirs.unit_vector();
        let mut gcrs = GCRSPosition::from_unit_vector(p, self.epoch)?;
        if let Some(distance) = cirs.distance() {
            gcrs.set_distance(distance);
        }
        Ok(gcrs)
    }

    pub fn icrs_to_cirs(&self, icrs: &ICRSPosition) -> CoordResult<CIRSPosition> {
        self.cirs_from_gcrs_vector(self.apparent(icrs)?, icrs.distance())
    }

    pub fn cirs_to_icrs(&self, cirs: &CIRSPosition) -> CoordResult<ICRSPosition> {
        check_epoch(self.epoch, cirs.epoch())?;
        let p = self.gcrs_to_cirs.transpose() * cirs.unit_vector();
        let mut icrs = ICRSPosition::from_unit_vector(self.natural(p)?)?;
        if let Some(distance) = cirs.distance() {
            icrs.set_distance(distance);
        }
        Ok(icrs)
    }

    // ICRS catalogue direction to the GCRS direction seen from the geocentre:
    // light deflection by the Sun, then annual aberration.
    fn apparent(&self, icrs: &ICRSPosition) -> CoordResult<Vector3> {
        let p = icrs.unit_vector().normalize()?;
        let deflected = apply_light_deflection(p, self.sun_to_earth, self.sun_distance_au);
        Ok(aberrate(
            deflected,
            self.velocity_c,
            self.sun_distance_au,
            self.bm1,
        ))
    }

    fn natural(&self, apparent: Vector3) -> CoordResult<Vector3> {
        let deflected = unaberrate(apparent, self.velocity_c, self.sun_distance_au, self.bm1)?;
        remove_light_deflection(deflected, self.sun_to_earth, self.sun_distance_au)
    }

    fn cirs_from_gcrs_vector(
        &self,
        p: Vector3,
        distance: Option<Distance>,
    ) -> CoordResult<CIRSPosition> {
        let mut cirs = CIRSPosition::from_unit_vector(self.gcrs_to_cirs * p, self.epoch)?;
        if let Some(distance) = distance {
            cirs.set_distance(distance);
        }
        Ok(cirs)
    }
}

fn check_epoch(context: TT, epoch: TT) -> CoordResult<()> {
    if epoch == context {
        return Ok(());
    }
    Err(CoordError::invalid_coordinate(format!(
        "position epoch {epoch:?} does not match the context epoch {context:?}"
    )))
}

// IAU 2006/2000A GCRS-to-CIRS matrix. With EOP, the observed celestial pole
// offsets dX and dY are added to the model CIP; s stays the model value.
fn celestial_to_intermediate(
    epoch: &TT,
    eop: Option<&EopParameters>,
) -> CoordResult<RotationMatrix3> {
    let jd = epoch.to_julian_date();
    let t = jd_to_centuries(jd.jd1(), jd.jd2());
    let nutation = epoch.nutation_iau2006a()?;
    let npb = PrecessionIAU2006::new().npb_matrix_iau2006a(
        t,
        nutation.nutation_longitude(),
        nutation.nutation_obliquity(),
    );
    let cio = CioSolution::calculate(&npb, t)?;
    let (x, y) = match eop {
        Some(eop) => (
            eop.corrected_cip_x(cio.cip.x),
            eop.corrected_cip_y(cio.cip.y),
        ),
        None => (cio.cip.x, cio.cip.y),
    };
    Ok(gcrs_to_cirs_matrix(x, y, cio.s)?)
}
