// A WCS saved to XISF and read back from the solution properties alone must
// map pixels to the sky exactly as the FITS header it came from.
use super::Image;
use crate::core::errors::ImageError;
use crate::fits::header::keywords::Keyword;
use crate::formats::pixel_data::PixelData;
use crate::xisf::reader::XisfFile;
use crate::xisf::solution::{celestial_wcs, solution_keywords};
use crate::xisf::writer::XisfWriter;
use celestial_wcs::builder::Wcs;
use celestial_wcs::coordinate::PixelCoord;

// Corners, edges and interior points of a 1024 x 768 image, in FITS pixels.
const PIXELS: [(f64, f64); 8] = [
    (0.5, 0.5),
    (1.0, 1.0),
    (1024.0, 1.0),
    (1.0, 768.0),
    (1024.0, 768.0),
    (1024.5, 768.5),
    (512.5, 384.5),
    (300.25, 700.75),
];

fn header(code: &str, crval: (f64, f64), extra: Vec<Keyword>) -> Vec<Keyword> {
    let mut keywords = vec![
        Keyword::string("CTYPE1", format!("RA---{code}").as_str()),
        Keyword::string("CTYPE2", format!("DEC--{code}").as_str()),
        Keyword::real("CRPIX1", 512.3),
        Keyword::real("CRPIX2", 380.7),
        Keyword::real("CRVAL1", crval.0),
        Keyword::real("CRVAL2", crval.1),
        Keyword::real("CD1_1", -1.1e-3),
        Keyword::real("CD1_2", 2.3e-5),
        Keyword::real("CD2_1", -1.7e-5),
        Keyword::real("CD2_2", 1.05e-3),
    ];
    for keyword in extra {
        keywords.retain(|k| k.name != keyword.name);
        keywords.push(keyword);
    }
    keywords
}

// Saves an image with the keywords, then builds a WCS from the solution
// properties read back through the public reader.
fn saved_and_read_back(keywords: &[Keyword]) -> Result<Wcs, ImageError> {
    let tmp = tempfile::Builder::new().suffix(".xisf").tempfile().unwrap();
    let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
    img.keywords = keywords.to_vec();
    img.save(tmp.path())?;
    let properties = XisfFile::open(tmp.path()).unwrap().properties(0).unwrap();
    let decoded = solution_keywords(&properties).expect("the saved solution decodes");
    Ok(celestial_wcs(&decoded).unwrap().unwrap())
}

fn assert_same_mapping(original: &Wcs, decoded: &Wcs, case: &str) {
    for (x, y) in PIXELS {
        let pixel = PixelCoord::new(x, y);
        let want = original.pixel_to_celestial(pixel).unwrap();
        assert_eq!(
            decoded.pixel_to_celestial(pixel).unwrap(),
            want,
            "{case} ({x}, {y})"
        );
        assert_eq!(
            decoded.celestial_to_pixel(want).unwrap(),
            original.celestial_to_pixel(want).unwrap(),
            "{case} ({x}, {y}) inverse"
        );
    }
}

fn assert_round_trip(case: &str, keywords: Vec<Keyword>) {
    let original = celestial_wcs(&keywords).unwrap().unwrap();
    let decoded = saved_and_read_back(&keywords).unwrap();
    assert_same_mapping(&original, &decoded, case);
}

#[test]
fn table_16_projections_round_trip_exactly() {
    for code in ["TAN", "STG", "SIN", "ZEA", "CAR", "MER", "AIT"] {
        assert_round_trip(code, header(code, (37.5, 30.0), vec![]));
    }
}

#[test]
fn car_with_a_moved_fiducial_point_round_trips_exactly() {
    let extra = vec![
        Keyword::real("PV1_1", 5.0),
        Keyword::real("PV1_2", 20.0),
        Keyword::real("LONPOLE", 10.0),
    ];
    assert_round_trip("CAR PV1_1 PV1_2", header("CAR", (37.5, 30.0), extra));
}

#[test]
fn southern_fields_round_trip_exactly() {
    for code in ["TAN", "CAR", "AIT"] {
        assert_round_trip(code, header(code, (250.0, -65.25), vec![]));
    }
}

#[test]
fn integer_valued_keywords_round_trip_exactly() {
    let keywords = vec![
        Keyword::string("CTYPE1", "RA---TAN"),
        Keyword::string("CTYPE2", "DEC--TAN"),
        Keyword::integer("CRPIX1", 512),
        Keyword::integer("CRPIX2", 384),
        Keyword::integer("CRVAL1", 150),
        Keyword::integer("CRVAL2", -30),
        Keyword::integer("CD1_1", -1),
        Keyword::integer("CD1_2", 0),
        Keyword::integer("CD2_1", 0),
        Keyword::integer("CD2_2", 1),
    ];
    assert_round_trip("integer TAN", keywords);
}

// The properties hold the linear part of a TAN-SIP solution.
#[test]
fn tan_sip_round_trips_its_linear_part_exactly() {
    let sip = vec![
        Keyword::string("CTYPE1", "RA---TAN-SIP"),
        Keyword::string("CTYPE2", "DEC--TAN-SIP"),
        Keyword::integer("A_ORDER", 2),
        Keyword::integer("B_ORDER", 2),
        Keyword::real("A_2_0", 2.5e-6),
        Keyword::real("A_1_1", -1.0e-6),
        Keyword::real("B_0_2", 3.0e-6),
    ];
    let linear = celestial_wcs(&header("TAN", (37.5, 30.0), vec![]))
        .unwrap()
        .unwrap();
    let decoded = saved_and_read_back(&header("TAN", (37.5, 30.0), sip)).unwrap();
    assert_same_mapping(&linear, &decoded, "TAN-SIP linear part");
}

// These projections have no Table 16 name, so the file holds no solution
// properties, and the FITSKeyword elements carry the WCS back.
#[test]
fn projections_outside_table_16_round_trip_through_the_keywords() {
    let quad_cube = header("TSC", (37.5, 30.0), vec![]);
    let conic = |code| header(code, (37.5, 30.0), vec![Keyword::real("PV2_1", 45.0)]);
    for (case, keywords) in [
        ("ARC", header("ARC", (37.5, 30.0), vec![])),
        ("CEA", header("CEA", (37.5, 30.0), vec![])),
        ("MOL", header("MOL", (37.5, 30.0), vec![])),
        ("COP", conic("COP")),
        ("COE", conic("COE")),
        ("TSC", quad_cube),
    ] {
        let tmp = tempfile::Builder::new().suffix(".xisf").tempfile().unwrap();
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.keywords = keywords.clone();
        img.save(tmp.path()).unwrap();
        let properties = XisfFile::open(tmp.path()).unwrap().properties(0).unwrap();
        assert_eq!(properties, vec![], "{case}");

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.keywords, keywords, "{case}");
        let original = celestial_wcs(&keywords).unwrap().unwrap();
        let rebuilt = celestial_wcs(&restored.keywords).unwrap().unwrap();
        assert_same_mapping(&original, &rebuilt, case);
    }
}

// XISF stores CRPIX - 0.5 and the reader adds 0.5 back. FITS to XISF holding
// only the properties, then back to FITS, must give the same CRPIX.
#[test]
fn crpix_round_trips_through_the_reference_image_coordinates() {
    let dir = tempfile::tempdir().unwrap();
    for crpix in [0.25, 0.5, 1.0, -3.7, 1024.5, 2048.25] {
        let extra = vec![
            Keyword::real("CRPIX1", crpix),
            Keyword::real("CRPIX2", crpix),
        ];
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.keywords = header("TAN", (37.5, 30.0), extra);
        img.save(dir.path().join("in.fits")).unwrap();
        let from_fits = Image::open(dir.path().join("in.fits")).unwrap();
        from_fits.save(dir.path().join("full.xisf")).unwrap();

        let mut xisf = XisfFile::open(dir.path().join("full.xisf")).unwrap();
        let mut writer = XisfWriter::create(dir.path().join("props.xisf")).unwrap();
        writer.add_image(&[0u8; 4], 2, 2, 1).unwrap();
        writer.add_properties(xisf.properties(0).unwrap());
        writer.write().unwrap();
        Image::open(dir.path().join("props.xisf"))
            .unwrap()
            .save(dir.path().join("out.fits"))
            .unwrap();

        let out = Image::open(dir.path().join("out.fits")).unwrap();
        assert_eq!(out.get("CRPIX1").as_f64(), Some(crpix), "{crpix}");
        assert_eq!(out.get("CRPIX2").as_f64(), Some(crpix), "{crpix}");
    }
}
