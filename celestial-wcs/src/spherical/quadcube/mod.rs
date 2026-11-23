use crate::common::native_coord;
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::WcsResult;
use crate::trig::{atan2d, sincosd};

mod csc;
mod qsc;
mod tsc;

pub(super) use csc::{deproject_csc, project_csc};
pub(super) use qsc::{deproject_qsc, project_qsc};
pub(super) use tsc::{deproject_tsc, project_tsc};

// Paper II §5.6: face 0 is centred on θ = 90°, faces 1 to 4 on the equator
// at φ = 0°, 90°, 180° and 270°, and face 5 on θ = −90°.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Face {
    Top,
    Equator0,
    Equator90,
    Equator180,
    Equator270,
    Bottom,
}

impl Face {
    const ALL: [Self; 6] = [
        Self::Top,
        Self::Equator0,
        Self::Equator90,
        Self::Equator180,
        Self::Equator270,
        Self::Bottom,
    ];

    // Face 4 is centred at φ = −90°, left of face 1 in the forward layout.
    fn centre(self) -> (f64, f64) {
        match self {
            Self::Top => (0.0, 90.0),
            Self::Equator0 => (0.0, 0.0),
            Self::Equator90 => (90.0, 0.0),
            Self::Equator180 => (180.0, 0.0),
            Self::Equator270 => (-90.0, 0.0),
            Self::Bottom => (0.0, -90.0),
        }
    }

    // Table 3: (ξ, η, ζ) on this face, the inverse of direction_cosines.
    fn plane_cosines(self, (l, m, n): (f64, f64, f64)) -> (f64, f64, f64) {
        match self {
            Self::Top => (m, -l, n),
            Self::Equator0 => (m, n, l),
            Self::Equator90 => (-l, n, m),
            Self::Equator180 => (-m, n, -l),
            Self::Equator270 => (l, n, -m),
            Self::Bottom => (m, l, -n),
        }
    }

    fn direction_cosines(self, xi: f64, eta: f64, zeta: f64) -> (f64, f64, f64) {
        match self {
            Self::Top => (-eta, xi, zeta),
            Self::Equator0 => (zeta, xi, eta),
            Self::Equator90 => (-xi, zeta, eta),
            Self::Equator180 => (-zeta, -xi, eta),
            Self::Equator270 => (xi, -zeta, eta),
            Self::Bottom => (eta, xi, -zeta),
        }
    }

    // Paper II leaves (φ, θ) to be read from (l, m, n). Near a polar face
    // centre n rounds to 1 and sin⁻¹ n loses θ, so θ is arg(√(l² + m²), n).
    // At a pole φ has no value and 0 stands for it; l can be −0 there, which
    // would give 180°.
    fn native(self, xi: f64, eta: f64, zeta: f64) -> WcsResult<NativeCoord> {
        let (l, m, n) = self.direction_cosines(xi, eta, zeta);
        let phi = if l == 0.0 && m == 0.0 {
            0.0
        } else {
            atan2d(m, l)
        };
        native_coord(phi, atan2d(n, libm::hypot(l, m)))
    }

    fn native_from_plane(self, chi: f64, psi: f64) -> WcsResult<NativeCoord> {
        let zeta = 1.0 / libm::sqrt(1.0 + chi * chi + psi * psi);
        self.native(chi * zeta, psi * zeta, zeta)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct QuadcubeFace {
    face: Face,
    xi: f64,
    eta: f64,
    zeta: f64,
}

impl QuadcubeFace {
    fn new(face: Face, (xi, eta, zeta): (f64, f64, f64)) -> Self {
        Self {
            face,
            xi,
            eta,
            zeta,
        }
    }

    fn plane(&self) -> (f64, f64) {
        (self.xi / self.zeta, self.eta / self.zeta)
    }

    fn intermediate(&self, x: f64, y: f64) -> IntermediateCoord {
        let (phi_c, theta_c) = self.face.centre();
        IntermediateCoord::new(phi_c + x, theta_c + y)
    }
}

// The face with the largest direction cosine ζ, which is at least 1/√3.
fn front_face(native: NativeCoord) -> QuadcubeFace {
    select_quadcube_face(native.phi_deg(), native.theta_deg())
}

fn select_quadcube_face(phi: f64, theta: f64) -> QuadcubeFace {
    let (sin_theta, cos_theta) = sincosd(theta);
    let (sin_phi, cos_phi) = sincosd(phi);
    let lmn = (cos_theta * cos_phi, cos_theta * sin_phi, sin_theta);
    let mut best = QuadcubeFace::new(Face::Top, (0.0, 0.0, f64::NEG_INFINITY));
    for face in Face::ALL {
        let candidate = QuadcubeFace::new(face, face.plane_cosines(lmn));
        if candidate.zeta > best.zeta {
            best = candidate;
        }
    }
    best
}

fn quadcube_face_from_xy(x: f64, y: f64) -> (Face, f64, f64) {
    let x = forward_layout_x(x);
    let face = if y > 45.0 {
        Face::Top
    } else if y < -45.0 {
        Face::Bottom
    } else if (-45.0..45.0).contains(&x) {
        Face::Equator0
    } else if (45.0..135.0).contains(&x) {
        Face::Equator90
    } else if !(-135.0..135.0).contains(&x) {
        Face::Equator180
    } else {
        Face::Equator270
    };
    let (phi_c, theta_c) = face.centre();
    (face, x - phi_c, y - theta_c)
}

// Faces 2, 3 and 4 may sit on either side of face 1, or be split between the
// two, so the equatorial band spans |x| ≤ 315° (Paper II §5.6). A 360° shift
// takes any such layout onto the forward projection's, x in [−135°, 225°).
// Wider x lies on no face and is left for the face bounds check to reject.
fn forward_layout_x(x: f64) -> f64 {
    if (225.0..=315.0).contains(&x) {
        x - 360.0
    } else if (-315.0..-135.0).contains(&x) {
        x + 360.0
    } else {
        x
    }
}

#[cfg(test)]
mod tests {
    use crate::common::assert_pinned;
    use crate::coordinate::IntermediateCoord;
    use crate::spherical::Projection;

    // Per-projection native_reference checks, both ways, are covered by
    // spherical::tests::test_all_projections_map_reference_to_origin_and_back.

    #[test]
    fn test_quadcube_deproject_accepts_either_layout() {
        // Faces 2, 3 and 4 may sit on either side of face 1 (Paper II §5.6), so
        // a point and its copy 360° along the equatorial band are the same.
        let shifted = [(280.0, -80.0), (-200.0, 160.0), (-260.0, 100.0)];
        for proj in [Projection::tsc(), Projection::csc(), Projection::qsc()] {
            for (x, x_forward) in shifted {
                let got = proj.deproject(IntermediateCoord::new(x, 10.0)).unwrap();
                let expected = proj
                    .deproject(IntermediateCoord::new(x_forward, 10.0))
                    .unwrap();
                assert_eq!(got, expected, "{proj:?} x = {x}");
            }
            let beyond = proj.deproject(IntermediateCoord::new(400.0, 10.0));
            assert!(beyond.is_err(), "{proj:?} x = 400");
        }
    }

    // The polar faces put each pole at (0°, ±90°) for any φ. Coming back, φ
    // has no value there and 0 stands for it.
    #[test]
    fn test_quadcube_poles_are_exact() {
        let rows = [
            ((37.5, 90.0), (0.0, 90.0), (0.0, 90.0)),
            ((37.5, -90.0), (0.0, -90.0), (0.0, -90.0)),
        ];
        for proj in [Projection::tsc(), Projection::csc(), Projection::qsc()] {
            assert_pinned(&proj, &rows);
        }
    }
}
