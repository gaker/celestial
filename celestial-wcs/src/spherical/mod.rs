use crate::common::{require_finite, require_on_sphere};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::WcsResult;
use crate::trig::sincosd;

mod code;
mod conic;
mod cylindrical;
mod pole;
mod polyconic;
mod pseudocylindrical;
mod quadcube;
mod rotation;
mod zenithal;

pub use rotation::SphericalRotation;

use conic::{cod_cone, coe_cone, coo_cone, cop_cone, CodCone, CoeCone, CooCone, CopCone};
use conic::{deproject_cod, deproject_coe, deproject_coo, deproject_cop};
use conic::{project_cod, project_coe, project_coo, project_cop};
use cylindrical::{deproject_car, deproject_cea, deproject_cyp, deproject_mer};
use cylindrical::{project_car, project_cea, project_cyp, project_mer};
use polyconic::{bonne_y0, deproject_bon, deproject_pco, project_bon, project_pco};
use pseudocylindrical::{deproject_ait, deproject_mol, deproject_par, deproject_sfl};
use pseudocylindrical::{project_ait, project_mol, project_par, project_sfl};
use quadcube::{deproject_csc, deproject_qsc, deproject_tsc};
use quadcube::{project_csc, project_qsc, project_tsc};
use zenithal::{air_k, projection_point};
use zenithal::{deproject_air, deproject_arc, deproject_azp, deproject_sin, deproject_stg};
use zenithal::{deproject_szp, deproject_tan, deproject_zea, deproject_zpn};
use zenithal::{project_air, project_arc, project_azp, project_sin, project_stg};
use zenithal::{project_szp, project_tan, project_zea, project_zpn};

#[derive(Debug, Clone)]
pub struct Projection(Kind);

// The parameters as given, and what each projection derives from them once.
// The constructors cannot fail, so a parameter with no valid derivation keeps
// its error for each call that needs it.
#[derive(Debug, Clone)]
enum Kind {
    Tan,
    Sin {
        xi: f64,
        eta: f64,
    },
    Arc,
    Stg,
    Zea,
    Azp {
        mu: f64,
        gamma: f64,
        sincos_gamma: (f64, f64),
    },
    Szp {
        mu: f64,
        phi_c: f64,
        theta_c: f64,
        point: WcsResult<(f64, f64, f64)>,
    },
    Zpn {
        coeffs: Vec<f64>,
    },
    Air {
        theta_b: f64,
        k: WcsResult<f64>,
    },
    Car,
    Mer,
    Cea {
        lambda: f64,
    },
    Cyp {
        mu: f64,
        lambda: f64,
    },
    Sfl,
    Par,
    Mol,
    Ait,
    Cop {
        theta_a: f64,
        eta: f64,
        cone: WcsResult<CopCone>,
    },
    Coe {
        theta_a: f64,
        eta: f64,
        cone: WcsResult<CoeCone>,
    },
    Cod {
        theta_a: f64,
        eta: f64,
        cone: WcsResult<CodCone>,
    },
    Coo {
        theta_a: f64,
        eta: f64,
        cone: WcsResult<CooCone>,
    },
    Bon {
        theta_1: f64,
        y0: WcsResult<f64>,
    },
    Pco,
    Tsc,
    Csc,
    Qsc,
}

impl Projection {
    pub fn tan() -> Self {
        Self(Kind::Tan)
    }

    pub fn sin() -> Self {
        Self(Kind::Sin { xi: 0.0, eta: 0.0 })
    }

    pub fn sin_with_params(xi: f64, eta: f64) -> Self {
        Self(Kind::Sin { xi, eta })
    }

    pub fn arc() -> Self {
        Self(Kind::Arc)
    }

    pub fn stg() -> Self {
        Self(Kind::Stg)
    }

    pub fn zea() -> Self {
        Self(Kind::Zea)
    }

    pub fn azp(mu: f64, gamma: f64) -> Self {
        Self(Kind::Azp {
            mu,
            gamma,
            sincos_gamma: sincosd(gamma),
        })
    }

    pub fn szp(mu: f64, phi_c: f64, theta_c: f64) -> Self {
        Self(Kind::Szp {
            mu,
            phi_c,
            theta_c,
            point: projection_point(mu, phi_c, theta_c),
        })
    }

    // A zero at the top of eq. 68's sum adds nothing, so trimming the trailing
    // zeros makes equal polynomials one projection and keeps a header's 30
    // defaulted terms out of every evaluation.
    pub fn zpn(mut coeffs: Vec<f64>) -> Self {
        while coeffs.last() == Some(&0.0) {
            coeffs.pop();
        }
        Self(Kind::Zpn { coeffs })
    }

    pub fn air(theta_b: f64) -> Self {
        Self(Kind::Air {
            theta_b,
            k: air_k(theta_b),
        })
    }

    pub fn car() -> Self {
        Self(Kind::Car)
    }

    pub fn mer() -> Self {
        Self(Kind::Mer)
    }

    pub fn cea() -> Self {
        Self(Kind::Cea { lambda: 1.0 })
    }

    pub fn cea_with_lambda(lambda: f64) -> Self {
        Self(Kind::Cea { lambda })
    }

    pub fn cyp(mu: f64, lambda: f64) -> Self {
        Self(Kind::Cyp { mu, lambda })
    }

    pub fn sfl() -> Self {
        Self(Kind::Sfl)
    }

    pub fn par() -> Self {
        Self(Kind::Par)
    }

    pub fn mol() -> Self {
        Self(Kind::Mol)
    }

    pub fn ait() -> Self {
        Self(Kind::Ait)
    }

    pub fn cop(theta_a: f64, eta: f64) -> Self {
        Self(Kind::Cop {
            theta_a,
            eta,
            cone: cop_cone(theta_a, eta),
        })
    }

    pub fn coe(theta_a: f64, eta: f64) -> Self {
        Self(Kind::Coe {
            theta_a,
            eta,
            cone: coe_cone(theta_a, eta),
        })
    }

    pub fn cod(theta_a: f64, eta: f64) -> Self {
        Self(Kind::Cod {
            theta_a,
            eta,
            cone: cod_cone(theta_a, eta),
        })
    }

    pub fn coo(theta_a: f64, eta: f64) -> Self {
        Self(Kind::Coo {
            theta_a,
            eta,
            cone: coo_cone(theta_a, eta),
        })
    }

    pub fn bon(theta_1: f64) -> Self {
        Self(Kind::Bon {
            theta_1,
            y0: bonne_y0(theta_1),
        })
    }

    pub fn pco() -> Self {
        Self(Kind::Pco)
    }

    pub fn tsc() -> Self {
        Self(Kind::Tsc)
    }

    pub fn csc() -> Self {
        Self(Kind::Csc)
    }

    pub fn qsc() -> Self {
        Self(Kind::Qsc)
    }

    pub fn native_reference(&self) -> (f64, f64) {
        match &self.0 {
            Kind::Tan
            | Kind::Sin { .. }
            | Kind::Arc
            | Kind::Stg
            | Kind::Zea
            | Kind::Azp { .. }
            | Kind::Szp { .. }
            | Kind::Zpn { .. }
            | Kind::Air { .. } => (0.0, 90.0),
            Kind::Car | Kind::Mer | Kind::Cea { .. } | Kind::Cyp { .. } => (0.0, 0.0),
            Kind::Sfl | Kind::Par | Kind::Mol | Kind::Ait => (0.0, 0.0),
            Kind::Cop { theta_a, .. }
            | Kind::Coe { theta_a, .. }
            | Kind::Cod { theta_a, .. }
            | Kind::Coo { theta_a, .. } => (0.0, *theta_a),
            Kind::Bon { .. } => (0.0, 0.0),
            Kind::Pco => (0.0, 0.0),
            Kind::Tsc | Kind::Csc | Kind::Qsc => (0.0, 0.0),
        }
    }

    pub fn project(&self, native: NativeCoord) -> WcsResult<IntermediateCoord> {
        require_on_sphere(native.phi_deg(), native.theta_deg())?;
        match &self.0 {
            Kind::Tan => project_tan(native),
            Kind::Sin { xi, eta } => project_sin(native, *xi, *eta),
            Kind::Arc => project_arc(native),
            Kind::Stg => project_stg(native),
            Kind::Zea => project_zea(native),
            Kind::Azp {
                mu, sincos_gamma, ..
            } => project_azp(native, *mu, *sincos_gamma),
            Kind::Szp { point, .. } => project_szp(native, point),
            Kind::Zpn { coeffs } => project_zpn(native, coeffs),
            Kind::Air { k, .. } => project_air(native, k),
            Kind::Car => project_car(native),
            Kind::Mer => project_mer(native),
            Kind::Cea { lambda } => project_cea(native, *lambda),
            Kind::Cyp { mu, lambda } => project_cyp(native, *mu, *lambda),
            Kind::Sfl => project_sfl(native),
            Kind::Par => project_par(native),
            Kind::Mol => project_mol(native),
            Kind::Ait => project_ait(native),
            Kind::Cop { theta_a, cone, .. } => project_cop(native, *theta_a, cone),
            Kind::Coe { cone, .. } => project_coe(native, cone),
            Kind::Cod { theta_a, cone, .. } => project_cod(native, *theta_a, cone),
            Kind::Coo { theta_a, cone, .. } => project_coo(native, *theta_a, cone),
            Kind::Bon { theta_1, y0 } => project_bon(native, *theta_1, y0),
            Kind::Pco => project_pco(native),
            Kind::Tsc => project_tsc(native),
            Kind::Csc => project_csc(native),
            Kind::Qsc => project_qsc(native),
        }
    }

    pub fn deproject(&self, inter: IntermediateCoord) -> WcsResult<NativeCoord> {
        require_finite(&[inter.x_deg(), inter.y_deg()], "intermediate coordinates")?;
        match &self.0 {
            Kind::Tan => deproject_tan(inter),
            Kind::Sin { xi, eta } => deproject_sin(inter, *xi, *eta),
            Kind::Arc => deproject_arc(inter),
            Kind::Stg => deproject_stg(inter),
            Kind::Zea => deproject_zea(inter),
            Kind::Azp {
                mu, sincos_gamma, ..
            } => deproject_azp(inter, *mu, *sincos_gamma),
            Kind::Szp { point, .. } => deproject_szp(inter, point),
            Kind::Zpn { coeffs } => deproject_zpn(inter, coeffs),
            Kind::Air { k, .. } => deproject_air(inter, k),
            Kind::Car => deproject_car(inter),
            Kind::Mer => deproject_mer(inter),
            Kind::Cea { lambda } => deproject_cea(inter, *lambda),
            Kind::Cyp { mu, lambda } => deproject_cyp(inter, *mu, *lambda),
            Kind::Sfl => deproject_sfl(inter),
            Kind::Par => deproject_par(inter),
            Kind::Mol => deproject_mol(inter),
            Kind::Ait => deproject_ait(inter),
            Kind::Cop { theta_a, cone, .. } => deproject_cop(inter, *theta_a, cone),
            Kind::Coe { theta_a, cone, .. } => deproject_coe(inter, *theta_a, cone),
            Kind::Cod { theta_a, cone, .. } => deproject_cod(inter, *theta_a, cone),
            Kind::Coo { theta_a, cone, .. } => deproject_coo(inter, *theta_a, cone),
            Kind::Bon { theta_1, y0 } => deproject_bon(inter, *theta_1, y0),
            Kind::Pco => deproject_pco(inter),
            Kind::Tsc => deproject_tsc(inter),
            Kind::Csc => deproject_csc(inter),
            Kind::Qsc => deproject_qsc(inter),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinate::CelestialCoord;
    use crate::error::WcsError;
    use celestial_core::constants::SQRT2;

    fn every_projection() -> [Projection; 26] {
        [
            Projection::tan(),
            Projection::sin(),
            Projection::arc(),
            Projection::stg(),
            Projection::zea(),
            Projection::azp(2.0, 30.0),
            Projection::szp(2.0, 45.0, 60.0),
            Projection::zpn(vec![0.0, 1.0]),
            Projection::air(45.0),
            Projection::car(),
            Projection::mer(),
            Projection::cea_with_lambda(0.5),
            Projection::cyp(1.0, 2.0),
            Projection::sfl(),
            Projection::par(),
            Projection::mol(),
            Projection::ait(),
            Projection::cop(45.0, 0.0),
            Projection::coe(45.0, 0.0),
            Projection::cod(45.0, 0.0),
            Projection::coo(45.0, 0.0),
            Projection::bon(45.0),
            Projection::pco(),
            Projection::tsc(),
            Projection::csc(),
            Projection::qsc(),
        ]
    }

    // COE and COO bring the origin back an ulp below θa = 45°, a fifth and a
    // third of the half-ulp rounding of their inverse steps.
    #[test]
    fn test_all_projections_map_reference_to_origin_and_back() {
        for proj in &every_projection() {
            let (phi0, theta0) = proj.native_reference();
            let inter = proj.project(NativeCoord::new(phi0, theta0)).unwrap();
            assert_eq!((inter.x_deg(), inter.y_deg()), (0.0, 0.0), "{proj:?}");
            let native = proj.deproject(IntermediateCoord::new(0.0, 0.0)).unwrap();
            let back = match proj.code() {
                "COE" | "COO" => (0.0, 44.99999999999999),
                _ => (phi0, theta0),
            };
            assert_eq!((native.phi_deg(), native.theta_deg()), back, "{proj:?}");
        }
    }

    // Paper II §2.1 puts latitude on [−90°, 90°], and no equation has an
    // image for a NaN or infinite coordinate.
    const OFF_THE_SPHERE: [(f64, f64); 6] = [
        (f64::NAN, 37.5),
        (37.5, f64::NAN),
        (f64::INFINITY, 37.5),
        (37.5, f64::NEG_INFINITY),
        (37.5, 100.0),
        (37.5, -100.0),
    ];

    #[test]
    fn test_projections_reject_points_off_the_sphere() {
        let off_the_plane = [
            (f64::NAN, 10.0),
            (10.0, f64::NAN),
            (f64::INFINITY, 10.0),
            (10.0, f64::NEG_INFINITY),
        ];
        for proj in every_projection() {
            for (phi, theta) in OFF_THE_SPHERE {
                let result = proj.project(NativeCoord::new(phi, theta));
                assert!(result.is_err(), "{proj:?} ({phi}, {theta}): {result:?}");
            }
            for (x, y) in off_the_plane {
                let result = proj.deproject(IntermediateCoord::new(x, y));
                assert!(result.is_err(), "{proj:?} ({x}, {y}): {result:?}");
            }
        }
    }

    #[test]
    fn test_rotation_rejects_points_off_the_sphere() {
        let rotation = SphericalRotation::new(123.456, 52.5, 180.0);
        for (lon, lat) in OFF_THE_SPHERE {
            let celestial = rotation.native_to_celestial(NativeCoord::new(lon, lat));
            assert!(celestial.is_err(), "native ({lon}, {lat}): {celestial:?}");
            let native = rotation.celestial_to_native(CelestialCoord::new(lon, lat));
            assert!(native.is_err(), "celestial ({lon}, {lat}): {native:?}");
        }
    }

    #[test]
    fn test_deprojection_rejects_points_beyond_the_boundary() {
        // Off the end of a cylinder, past the edge of a bounded map, behind a
        // cone's apex: each comes back with φ or θ off the native sphere.
        let cases = [
            (Projection::car(), (200.0, 10.0)),
            (Projection::car(), (10.0, 200.0)),
            (Projection::mer(), (200.0, 10.0)),
            (Projection::cea(), (200.0, 10.0)),
            (Projection::cyp(1.0, 1.0), (10.0, 200.0)),
            (Projection::sfl(), (180.0, 10.0)),
            (Projection::sfl(), (10.0, 200.0)),
            (Projection::par(), (150.0, 60.0)),
            (Projection::mol(), (180.0, 10.0)),
            (Projection::ait(), (150.0, 60.0)),
            (Projection::cop(45.0, 25.0), (0.0, 57.29577951308232)),
            (Projection::cod(45.0, 25.0), (0.0, 57.29577951308232)),
            (Projection::cod(45.0, 25.0), (500.0, -500.0)),
            (Projection::coo(45.0, 25.0), (0.0, 57.29577951308232)),
            (Projection::bon(45.0), (10.0, 200.0)),
            (Projection::pco(), (180.0, 10.0)),
            (Projection::arc(), (0.0, -200.0)),
        ];
        for (proj, (x, y)) in cases {
            let err = proj.deproject(IntermediateCoord::new(x, y)).unwrap_err();
            let out = matches!(err, WcsError::OutOfBounds { .. });
            assert!(out, "{proj:?} ({x}, {y}): {err:?}");
        }
    }

    #[test]
    fn test_deprojection_keeps_points_on_the_boundary() {
        // Rounding can put the forward of the seam or a pole just past the
        // edge; it is still on the boundary and deprojects. The poles of SFL,
        // PAR, MOL and BON do not deproject at all yet, so they take the seam only.
        let with_poles = [
            Projection::car(),
            Projection::cea_with_lambda(0.75),
            Projection::cyp(1.0, SQRT2 / 2.0),
            Projection::ait(),
            Projection::cop(45.0, 25.0),
            Projection::coe(-45.0, 10.0),
            Projection::cod(45.0, 25.0),
            Projection::coo(45.0, 25.0),
            Projection::pco(),
            Projection::arc(),
        ];
        let seam_only = [
            Projection::sfl(),
            Projection::par(),
            Projection::mol(),
            Projection::bon(45.0),
        ];
        let keeps = |proj: &Projection, phi: f64, theta: f64| {
            let native = NativeCoord::new(phi, theta);
            if let Ok(inter) = proj.project(native) {
                let back = proj.deproject(inter);
                assert!(back.is_ok(), "{proj:?} ({phi}, {theta}): {back:?}");
            }
        };
        let halves = || (-179..=179).map(|t| t as f64 / 2.0);
        for proj in with_poles.iter().chain(&seam_only) {
            for theta in halves() {
                keeps(proj, -180.0, theta);
                keeps(proj, 180.0, theta);
            }
        }
        for proj in &with_poles {
            for phi in halves() {
                keeps(proj, phi, 90.0);
                keeps(proj, phi, -90.0);
            }
        }
    }
}
