use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

use super::{front_face, quadcube_face_from_xy};

pub(crate) fn project_csc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let face = front_face(native);
    let (chi, psi) = face.plane();
    let x_face = csc_forward_poly(chi, psi);
    let y_face = csc_forward_poly(psi, chi);
    Ok(face.intermediate(45.0 * x_face, 45.0 * y_face))
}

fn csc_forward_poly(chi: f64, psi: f64) -> f64 {
    const GAMMA_STAR: f64 = 1.37484847732;
    const M: f64 = 0.004869491981;
    const GAMMA: f64 = -0.13161671474;
    const OMEGA1: f64 = -0.159596235474;
    const C00: f64 = 0.141189631152;
    const C10: f64 = 0.0809701286525;
    const C01: f64 = -0.281528535557;
    const C20: f64 = -0.178251207466;
    const C11: f64 = 0.15384112876;
    const C02: f64 = 0.106959469314;
    const D0: f64 = 0.0759196200467;
    const D1: f64 = -0.0217762490699;

    let chi2 = chi * chi;
    let psi2 = psi * psi;

    let c_poly =
        C00 + C10 * chi2 + C01 * psi2 + C20 * chi2 * chi2 + C11 * chi2 * psi2 + C02 * psi2 * psi2;

    let d_poly = D0 + D1 * chi2;

    let term1 = chi * GAMMA_STAR + chi * chi2 * (1.0 - GAMMA_STAR);
    let term2 = chi * psi2 * (1.0 - chi2) * (GAMMA + (M - GAMMA) * chi2 + (1.0 - psi2) * c_poly);
    let term3 = chi * chi2 * (1.0 - chi2) * (OMEGA1 - (1.0 - chi2) * d_poly);

    term1 + term2 + term3
}

pub(crate) fn deproject_csc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (face, x_rel, y_rel) = quadcube_face_from_xy(inter.x_deg(), inter.y_deg());

    let x_norm = x_rel / 45.0;
    let y_norm = y_rel / 45.0;

    if x_norm.abs() > 1.0 || y_norm.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "CSC deprojection: point outside cube face",
        ));
    }

    let chi = csc_inverse_poly(x_norm, y_norm);
    let psi = csc_inverse_poly(y_norm, x_norm);
    face.native_from_plane(chi, psi)
}

fn csc_inverse_poly(x: f64, y: f64) -> f64 {
    // CSC inverse polynomial coefficients from Paper II (Calabretta & Greisen 2002)
    // Section 5.6.2, page 24. P[i][j] corresponds to P_ij where the polynomial is:
    // f(X,Y) = X + X(1 - X²) Σᵢⱼ PᵢⱼX²ⁱY²ʲ
    // The matrix is organized as P[i][j] where i is X² power and j is Y² power.
    const P: [[f64; 7]; 7] = [
        // i=0: P_00, P_01, P_02, P_03, P_04, P_05, P_06
        [
            -0.27292696,
            -0.02819452,
            0.27058160,
            -0.60441560,
            0.93412077,
            -0.63915306,
            0.14381585,
        ],
        // i=1: P_10, P_11, P_12, P_13, P_14, P_15
        [
            -0.07629969,
            -0.01471565,
            -0.56800938,
            1.50880086,
            -1.41601920,
            0.52032238,
            0.0,
        ],
        // i=2: P_20, P_21, P_22, P_23, P_24
        [
            -0.22797056,
            0.48051509,
            0.30803317,
            -0.93678576,
            0.33887446,
            0.0,
            0.0,
        ],
        // i=3: P_30, P_31, P_32, P_33
        [
            0.54852384,
            -1.74114454,
            0.98938102,
            0.08693841,
            0.0,
            0.0,
            0.0,
        ],
        // i=4: P_40, P_41, P_42
        [-0.62930065, 1.71547508, -0.83180469, 0.0, 0.0, 0.0, 0.0],
        // i=5: P_50, P_51
        [0.25795794, -0.53022337, 0.0, 0.0, 0.0, 0.0, 0.0],
        // i=6: P_60
        [0.02584375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    let x2 = x * x;
    let y2 = y * y;

    let mut sum = 0.0;
    let mut x_pow = 1.0;
    for (i, p_row) in P.iter().enumerate() {
        let mut y_pow = 1.0;
        for p_val in p_row.iter().take(7 - i) {
            sum += p_val * x_pow * y_pow;
            y_pow *= y2;
        }
        x_pow *= x2;
    }

    x + x * (1.0 - x2) * sum
}

#[cfg(test)]
mod tests {
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Eqs. 172 and 175 are separate fits (Paper II §5.6.2), neither the
    // inverse of the other, so each direction is pinned against its own
    // polynomial evaluated to 60 digits with the printed coefficients. A round
    // trip misses by 38″ on the sky at (20°, −75°). One row on each face, then
    // the seam and the edge between faces 3 and 4, which goes to face 3. The
    // projected values are within 2.77 ulp and the deprojected ones within
    // 1.68 ulp, each within the summed half-ulp rounding of the steps as
    // evaluated, the coefficients' rounding to binary among them.
    #[test]
    fn test_csc_pinned_both_ways() {
        let rows = [
            (
                (-30.0, 15.0),
                (-31.33571373148643, 18.068687049929647),
                (-30.0011301109495, 14.998471910659005),
            ),
            (
                (120.0, 15.0),
                (121.33571373148644, 18.068687049929647),
                (120.0011301109495, 14.998471910659005),
            ),
            (
                (180.0, -30.0),
                (180.0, -31.287623083766523),
                (180.0, -29.999896037503916),
            ),
            (
                (-90.0, -40.0),
                (-90.0, -40.55123726629974),
                (-90.0, -39.99806563906088),
            ),
            (
                (-135.0, 80.0),
                (-7.661615045826882, 97.66161504582689),
                (-135.00000000000006, 80.0033276256917),
            ),
            (
                (20.0, -75.0),
                (5.645215910238105, -74.84412250466025),
                (19.964995487867185, -74.99433052454935),
            ),
            (
                (-135.0, -30.0),
                (225.0, -39.36550667963484),
                (-135.0, -29.99714296090563),
            ),
        ];
        assert_pinned(&Projection::csc(), &rows);
    }
}
