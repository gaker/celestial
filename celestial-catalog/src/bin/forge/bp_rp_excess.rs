// Corrected BP and RP flux excess factor C* from Riello et al. 2021, A&A 649, A3,
// "Gaia Early Data Release 3: Photometric content and validation". Equation and table numbers
// are those of arXiv:2012.01916v1.
//
// Eq. 6 (Sec. 6): C* = C - f(G_BP - G_RP), where C is phot_bp_rp_excess_factor and f is the
// expected excess at that colour for sources with good photometry.

// Table 2: coefficients a_0 to a_3 of f(x) = sum a_i x^i, x = G_BP - G_RP, for x < 0.5,
// 0.5 <= x < 4.0 and x >= 4.0. The blue range is quadratic and the red range is linear.
const BLUE: [f64; 4] = [1.154360, 0.033772, 0.032277, 0.0];
const MIDDLE: [f64; 4] = [1.162004, 0.011464, 0.049255, -0.005879];
const RED: [f64; 4] = [1.057572, 0.140537, 0.0, 0.0];

// Sec. 6: the fit is valid for -1.0 <= G_BP - G_RP <= 7.0.
const MIN_COLOUR: f64 = -1.0;
const MAX_COLOUR: f64 = 7.0;

// Eq. 18 (Sec. 9.4): sigma_C*(G) = c0 + c1 G^m, the 1-sigma scatter of C* for well-behaved
// isolated stellar sources.
const SIGMA_C0: f64 = 0.0059898;
const SIGMA_C1: f64 = 8.817481e-12;
const SIGMA_M: f64 = 7.618399;

// Sec. 9.4 defines the filter |C*| < N sigma_C* and says to apply it only for G > 4, because
// saturation still dominates at brighter magnitudes. The paper shows 3-sigma and 5-sigma cuts.
// The flag marks a source without removing it, so it uses the tighter 3-sigma cut.
const MIN_G: f64 = 4.0;
const N_SIGMA: f64 = 3.0;

// Returns f(bp_rp), or None outside the colour range of the fit.
fn expected_excess(bp_rp: f64) -> Option<f64> {
    if !(MIN_COLOUR..=MAX_COLOUR).contains(&bp_rp) {
        return None;
    }
    let a = if bp_rp < 0.5 {
        BLUE
    } else if bp_rp < 4.0 {
        MIDDLE
    } else {
        RED
    };
    let x = bp_rp;
    Some(a[0] + a[1] * x + a[2] * x * x + a[3] * x * x * x)
}

fn sigma_corrected_excess(g: f64) -> f64 {
    SIGMA_C0 + SIGMA_C1 * libm::pow(g, SIGMA_M)
}

// True when |C*| >= 3 sigma_C*(G). A source at G <= 4 or outside the colour range of the fit
// has no defined test and returns false.
pub(crate) fn is_excess_suspect(excess: f64, bp_rp: f64, g: f64) -> bool {
    let Some(expected) = expected_excess(bp_rp) else {
        return false;
    };
    let corrected = excess - expected;
    g > MIN_G && libm::fabs(corrected) >= N_SIGMA * sigma_corrected_excess(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Each value is a_0 + a_1 x + a_2 x^2 + a_3 x^3 with the Table 2 coefficients, worked by
    // hand: at x = 1, 1.162004 + 0.011464 + 0.049255 - 0.005879 = 1.216844.
    #[test]
    fn test_expected_excess_at_hand_worked_colours() {
        assert_eq!(expected_excess(0.0), Some(1.154360));
        assert_eq!(expected_excess(1.0), Some(1.216844));
        assert_eq!(expected_excess(5.0), Some(1.760257));
        // 1.154360 - 0.033772 + 0.032277 = 1.152865, one ulp above in binary.
        assert_eq!(expected_excess(-1.0), Some(1.1528650000000003));
    }

    #[test]
    fn test_expected_excess_is_undefined_outside_the_fit() {
        assert_eq!(expected_excess(-1.01), None);
        assert_eq!(expected_excess(7.01), None);
    }

    #[test]
    fn test_sigma_at_g_1_is_c0_plus_c1() {
        assert_eq!(sigma_corrected_excess(1.0), 0.0059898 + 8.817481e-12);
    }

    // At G = 10, sigma = 0.0059898 + 8.817481e-12 * 10^7.618399 = 0.0063560 to 7 digits, so
    // the 3-sigma limit is 0.019068. At x = 1, f = 1.216844.
    #[test]
    fn test_flag_at_three_sigma() {
        assert_eq!(sigma_corrected_excess(10.0), 0.006356021243144275);
        assert!(!is_excess_suspect(1.216844 + 0.0190, 1.0, 10.0));
        assert!(is_excess_suspect(1.216844 + 0.0191, 1.0, 10.0));
        assert!(is_excess_suspect(1.216844 - 0.0191, 1.0, 10.0));
    }

    #[test]
    fn test_no_flag_at_g_4_or_brighter() {
        assert!(!is_excess_suspect(2.0, 1.0, 4.0));
        assert!(is_excess_suspect(2.0, 1.0, 4.001));
    }
}
