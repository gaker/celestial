use crate::distortion::{DistortionModel, SipDistortion, TpvDistortion};
use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;
use crate::spherical::Projection;

use super::fiducial::Fiducial;
use super::keyword::WcsKeyword;

// The distortion a header's projection code and CTYPE suffix name. TPV is a
// projection code of its own, TAN with PVi_m as polynomial terms; SIP is the
// "-SIP" suffix on any code.
pub(super) fn header_distortion(
    header: &impl KeywordProvider,
    code: &str,
    suffix: Option<&str>,
    crpix: [f64; 2],
) -> WcsResult<Option<DistortionModel>> {
    match (code, suffix) {
        ("TPV", None) => Ok(Some(DistortionModel::Tpv(Box::new(read_tpv(header)?)))),
        (_, None) => Ok(None),
        (_, Some("SIP")) => Ok(Some(DistortionModel::Sip(read_sip(header, crpix)?))),
        (_, Some(other)) => Err(WcsError::unsupported_projection(format!("{code}-{other}"))),
    }
}

// TPV convention: absent terms are 0, except PV1_1 and PV2_1, which are 1.
fn read_tpv(header: &impl KeywordProvider) -> WcsResult<TpvDistortion> {
    let mut tpv = TpvDistortion::new();
    for m in 0..40 {
        if let Some(value) = header.get_float(&format!("PV1_{m}")) {
            tpv.set_pv(1, m, value)?;
        }
        if let Some(value) = header.get_float(&format!("PV2_{m}")) {
            tpv.set_pv(2, m, value)?;
        }
    }
    Ok(tpv)
}

// SIP convention §2: A_ORDER and B_ORDER are 2 to 9, absent terms are 0, and
// the reverse AP/BP polynomials are optional.
fn read_sip(header: &impl KeywordProvider, crpix: [f64; 2]) -> WcsResult<SipDistortion> {
    let a_order = require_order(header, "A_ORDER")?;
    let b_order = require_order(header, "B_ORDER")?;
    let mut sip = SipDistortion::new(crpix, a_order, b_order);
    read_terms(header, "A", a_order, |p, q, v| sip.set_a(p, q, v))?;
    read_terms(header, "B", b_order, |p, q, v| sip.set_b(p, q, v))?;
    match (order(header, "AP_ORDER", 0)?, order(header, "BP_ORDER", 0)?) {
        (Some(ap_order), Some(bp_order)) => {
            sip.set_inverse_order(ap_order, bp_order)?;
            read_terms(header, "AP", ap_order, |p, q, v| sip.set_ap(p, q, v))?;
            read_terms(header, "BP", bp_order, |p, q, v| sip.set_bp(p, q, v))?;
        }
        (None, None) => {}
        _ => return Err(WcsError::missing_keyword("AP_ORDER and BP_ORDER together")),
    }
    Ok(sip)
}

fn require_order(header: &impl KeywordProvider, key: &str) -> WcsResult<u32> {
    order(header, key, 2)?.ok_or_else(|| WcsError::missing_keyword(key))
}

// The convention bounds A_ORDER and B_ORDER to 2-9 and says nothing of the
// reverse orders; the upper bound applies to them too, which also caps the
// number of keywords looked up.
fn order(header: &impl KeywordProvider, key: &str, lowest: u32) -> WcsResult<Option<u32>> {
    let Some(value) = header.get_int(key) else {
        return Ok(None);
    };
    match u32::try_from(value) {
        Ok(order) if (lowest..=9).contains(&order) => Ok(Some(order)),
        _ => Err(WcsError::invalid_keyword(
            key,
            format!("{value} is outside {lowest} to 9"),
        )),
    }
}

fn read_terms(
    header: &impl KeywordProvider,
    name: &str,
    order: u32,
    mut set: impl FnMut(u32, u32, f64) -> WcsResult<()>,
) -> WcsResult<()> {
    for p in 0..=order {
        for q in 0..=order - p {
            if let Some(value) = header.get_float(&format!("{name}_{p}_{q}")) {
                set(p, q, value)?;
            }
        }
    }
    Ok(())
}

// Only what a header can carry is built: SIP measures (u, v) from CRPIX (SIP
// convention §2), and TPV is TAN with every PV1_m taken as a polynomial term,
// which leaves no PV1_0..2 for a moved fiducial point.
pub(super) fn check_distortion(
    distortion: Option<&DistortionModel>,
    crpix: [f64; 2],
    projection: &Projection,
    fiducial: &Fiducial,
) -> WcsResult<()> {
    let problem = match distortion {
        Some(DistortionModel::Sip(sip)) if sip.crpix() != crpix => {
            format!("SIP is centred on {:?}, not CRPIX {crpix:?}", sip.crpix())
        }
        Some(DistortionModel::Tpv(_)) if *projection != Projection::tan() => {
            "TPV distortion needs the TAN projection".to_string()
        }
        Some(DistortionModel::Tpv(_)) if *fiducial != Fiducial::of(projection) => {
            "TPV distortion cannot move the fiducial point".to_string()
        }
        _ => return Ok(()),
    };
    Err(WcsError::invalid_parameter(problem))
}

pub(super) fn ctype_code(proj_code: &str, distortion: Option<&DistortionModel>) -> String {
    match distortion {
        Some(DistortionModel::Sip(_)) => format!("{proj_code}-SIP"),
        Some(DistortionModel::Tpv(_)) => "TPV".to_string(),
        _ => proj_code.to_string(),
    }
}

pub(super) fn distortion_keywords(distortion: Option<&DistortionModel>) -> Vec<WcsKeyword> {
    match distortion {
        Some(DistortionModel::Sip(sip)) => sip_keywords(sip),
        Some(DistortionModel::Tpv(tpv)) => tpv_keywords(tpv),
        None => Vec::new(),
    }
}

fn sip_keywords(sip: &SipDistortion) -> Vec<WcsKeyword> {
    let mut keywords = Vec::new();
    for (name, order, terms) in sip.polynomials() {
        let Some(order) = order else {
            continue;
        };
        keywords.push(WcsKeyword::integer(
            format!("{name}_ORDER"),
            i64::from(order),
        ));
        keywords.extend(
            terms
                .iter()
                .map(|&((p, q), value)| WcsKeyword::real(format!("{name}_{p}_{q}"), value)),
        );
    }
    keywords
}

// Terms left at the defaults read_tpv starts from are left out.
fn tpv_keywords(tpv: &TpvDistortion) -> Vec<WcsKeyword> {
    let identity = TpvDistortion::new();
    let mut keywords = Vec::new();
    for ((axis, terms), (_, defaults)) in tpv.axes().into_iter().zip(identity.axes()) {
        for (m, (&value, &default)) in terms.iter().zip(defaults).enumerate() {
            if value != default {
                keywords.push(WcsKeyword::real(format!("{axis}_{m}"), value));
            }
        }
    }
    keywords
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::keyword::keyword_map;
    use crate::builder::{Wcs, WcsBuilder};
    use crate::coordinate::PixelCoord;
    use crate::distortion::{DistortionModel, SipDistortion, TpvDistortion};
    use crate::header::KeywordMap;

    fn tan_header(code: &str, floats: &[(&str, f64)], ints: &[(&str, i64)]) -> KeywordMap {
        let mut h = KeywordMap::new();
        h.set_string("CTYPE1", format!("RA---{code}"))
            .set_string("CTYPE2", format!("DEC--{code}"))
            .set_float("CRPIX1", 512.5)
            .set_float("CRPIX2", 512.5)
            .set_float("CRVAL1", 123.456)
            .set_float("CRVAL2", 37.5)
            .set_float("CD1_1", -1e-3)
            .set_float("CD2_2", 1e-3);
        for &(key, value) in floats {
            h.set_float(key, value);
        }
        for &(key, value) in ints {
            h.set_int(key, value);
        }
        h
    }

    fn header_wcs(header: &KeywordMap) -> WcsResult<Wcs> {
        WcsBuilder::from_header(header)?.build()
    }

    fn tan_with(distortion: DistortionModel) -> Wcs {
        WcsBuilder::new()
            .crpix(512.5, 512.5)
            .crval(123.456, 37.5)
            .cd_matrix([[-1e-3, 0.0], [0.0, 1e-3]])
            .proj_code("TAN")
            .distortion(distortion)
            .build()
            .unwrap()
    }

    #[track_caller]
    fn assert_same_mapping(got: &Wcs, expected: &Wcs) {
        for (x, y) in [(612.5, 412.5), (12.0, 1000.0), (900.25, 30.75)] {
            let pixel = PixelCoord::new(x, y);
            let world = expected.pixel_to_celestial(pixel).unwrap();
            assert_eq!(got.pixel_to_celestial(pixel).unwrap(), world, "({x}, {y})");
            let back = expected.celestial_to_pixel(world).unwrap();
            assert_eq!(got.celestial_to_pixel(world).unwrap(), back, "{world:?}");
        }
    }

    #[test]
    fn test_tan_sip_header_reads_the_polynomials() {
        let forward = [("A_2_0", 1e-6), ("A_1_1", -2e-7), ("B_0_2", -5e-7)];
        let reverse = [("AP_2_0", -1e-6), ("BP_0_2", 5e-7)];
        let orders = [("A_ORDER", 2), ("B_ORDER", 2)];
        let reverse_orders = [("AP_ORDER", 2), ("BP_ORDER", 2)];
        let mut sip = SipDistortion::new([512.5, 512.5], 2, 2);
        sip.set_a(2, 0, 1e-6).unwrap();
        sip.set_a(1, 1, -2e-7).unwrap();
        sip.set_b(0, 2, -5e-7).unwrap();
        let header = tan_header("TAN-SIP", &forward, &orders);
        assert_same_mapping(
            &header_wcs(&header).unwrap(),
            &tan_with(DistortionModel::Sip(sip.clone())),
        );

        sip.set_inverse_order(2, 2).unwrap();
        sip.set_ap(2, 0, -1e-6).unwrap();
        sip.set_bp(0, 2, 5e-7).unwrap();
        let floats = [forward.as_slice(), &reverse].concat();
        let header = tan_header("TAN-SIP", &floats, &[orders, reverse_orders].concat());
        assert_same_mapping(
            &header_wcs(&header).unwrap(),
            &tan_with(DistortionModel::Sip(sip)),
        );
    }

    #[test]
    fn test_tpv_header_reads_pv_as_distortion_terms() {
        // TPV convention: absent PV1_1 and PV2_1 are 1. PV1_0 and PV1_3 are
        // terms here, not the fiducial offset flag and LONPOLE of Paper II.
        let terms = [
            ("PV1_0", 1e-4),
            ("PV1_3", 2e-3),
            ("PV1_4", 1e-3),
            ("PV2_2", 1e-3),
            ("PV2_39", 1e-9),
        ];
        let mut tpv = TpvDistortion::new();
        tpv.set_pv(1, 0, 1e-4).unwrap();
        tpv.set_pv(1, 3, 2e-3).unwrap();
        tpv.set_pv(1, 4, 1e-3).unwrap();
        tpv.set_pv(2, 2, 1e-3).unwrap();
        tpv.set_pv(2, 39, 1e-9).unwrap();
        let got = header_wcs(&tan_header("TPV", &terms, &[])).unwrap();
        assert_same_mapping(&got, &tan_with(DistortionModel::Tpv(Box::new(tpv))));
    }

    #[test]
    fn test_malformed_distortion_headers_are_rejected() {
        let terms = [("A_2_0", 1e-6)];
        let mut mixed = tan_header("TAN-SIP", &terms, &[("A_ORDER", 2), ("B_ORDER", 2)]);
        mixed.set_string("CTYPE2", "DEC--TAN");
        let cases = [
            (
                "no A_ORDER",
                tan_header("TAN-SIP", &terms, &[("B_ORDER", 2)]),
            ),
            (
                "A_ORDER 10",
                tan_header("TAN-SIP", &terms, &[("A_ORDER", 10), ("B_ORDER", 2)]),
            ),
            (
                "AP_ORDER alone",
                tan_header(
                    "TAN-SIP",
                    &terms,
                    &[("A_ORDER", 2), ("B_ORDER", 2), ("AP_ORDER", 2)],
                ),
            ),
            ("unknown suffix", tan_header("TAN-XYZ", &[], &[])),
            ("suffix on one axis", mixed),
        ];
        for (label, header) in cases {
            assert!(header_wcs(&header).is_err(), "{label}");
        }
    }

    #[test]
    fn test_to_keywords_writes_the_distortion() {
        let mut sip = SipDistortion::new([512.5, 512.5], 3, 2);
        sip.set_a(3, 0, 1e-9).unwrap();
        sip.set_a(1, 1, -2e-7).unwrap();
        sip.set_b(0, 2, -5e-7).unwrap();
        sip.set_inverse_order(2, 2).unwrap();
        sip.set_ap(2, 0, -1e-6).unwrap();
        sip.set_bp(0, 2, 5e-7).unwrap();
        // PV1_1 away from its default of 1 must be written to survive.
        let mut tpv = TpvDistortion::new();
        tpv.set_pv(1, 1, 0.999).unwrap();
        tpv.set_pv(1, 4, 1e-3).unwrap();
        tpv.set_pv(2, 0, 1e-4).unwrap();
        tpv.set_pv(2, 39, 1e-9).unwrap();
        let cases = [
            ("RA---TAN-SIP", DistortionModel::Sip(sip)),
            ("RA---TPV", DistortionModel::Tpv(Box::new(tpv))),
        ];
        for (ctype, distortion) in cases {
            let wcs = tan_with(distortion);
            let header = keyword_map(wcs.to_keywords());
            assert_eq!(header.get_string("CTYPE1").as_deref(), Some(ctype));
            assert_eq!(header.get_string("CUNIT1").as_deref(), Some("deg"));
            assert_eq!(header.get_string("CUNIT2").as_deref(), Some("deg"));
            assert_same_mapping(&header_wcs(&header).unwrap(), &wcs);
        }
    }

    // AP_ORDER and BP_ORDER are optional (SIP §2); without the reverse
    // polynomials none of their keywords is written.
    #[test]
    fn test_sip_without_reverse_polynomials_writes_a_and_b_only() {
        let mut sip = SipDistortion::new([512.5, 512.5], 2, 2);
        sip.set_a(2, 0, 1e-6).unwrap();
        sip.set_b(0, 2, -1e-6).unwrap();
        let wcs = tan_with(DistortionModel::Sip(sip));
        let names: Vec<String> = wcs.to_keywords().into_iter().map(|k| k.name).collect();
        assert!(names.iter().any(|n| n == "A_ORDER"), "{names:?}");
        assert!(
            !names
                .iter()
                .any(|n| n.starts_with("AP_") || n.starts_with("BP_")),
            "{names:?}"
        );
        let header = keyword_map(wcs.to_keywords());
        assert_same_mapping(&header_wcs(&header).unwrap(), &wcs);
    }

    #[test]
    fn test_build_rejects_distortion_a_header_cannot_carry() {
        let base = || {
            WcsBuilder::new()
                .crpix(512.5, 512.5)
                .crval(123.456, 37.5)
                .cd_matrix([[-1e-3, 0.0], [0.0, 1e-3]])
        };
        let tpv = || DistortionModel::Tpv(Box::default());
        let off_centre = DistortionModel::Sip(SipDistortion::new([500.0, 500.0], 2, 2));
        let cases = [
            (
                "SIP away from CRPIX",
                base().proj_code("TAN").distortion(off_centre),
            ),
            ("TPV on SIN", base().proj_code("SIN").distortion(tpv())),
            (
                "TPV with a moved fiducial point",
                base().proj_code("TAN").pv(1, 2, 60.0).distortion(tpv()),
            ),
        ];
        for (label, builder) in cases {
            assert!(builder.build().is_err(), "{label}");
        }
    }
}
