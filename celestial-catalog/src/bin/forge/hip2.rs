// Records of hip2.dat, the Hipparcos new reduction (van Leeuwen 2007, CDS catalogue I/311),
// and their propagation from the Hipparcos epoch J1991.25 to the Gaia DR3 epoch.

use crate::gaia::GAIA_DR3_EPOCH;
use crate::intermediate;
use celestial_catalog::query::catalog::{
    StarRecord, FLAG_HAS_PARALLAX, FLAG_HAS_PROPER_MOTION, FLAG_SOURCE_HIPPARCOS,
};
use celestial_catalog::query::healpix::reduce_ra_deg;
use celestial_coords::errors::CoordResult;
use celestial_coords::proper_motion::CatalogStar;
use celestial_core::angle::Angle;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tt::TT;
use std::str::FromStr;

const HIPPARCOS_EPOCH: f64 = 1991.25;

// Byte ranges (start inclusive, end exclusive) of the hip2.dat fields read here.
const HIP: (usize, usize) = (0, 6);
const RA_RAD: (usize, usize) = (15, 28);
const DEC_RAD: (usize, usize) = (29, 42);
const PARALLAX_MAS: (usize, usize) = (43, 50);
const PMRA_MAS_YR: (usize, usize) = (51, 59);
const PMDEC_MAS_YR: (usize, usize) = (60, 68);
const HPMAG: (usize, usize) = (129, 136);
const B_V: (usize, usize) = (152, 158);

// A hip2.dat record: its J1991.25 astrometry and its place and motion at J2016.0.
pub(crate) struct Hip2Star {
    pub(crate) hip: u32,
    pub(crate) hpmag: f32,
    pub(crate) b_v: Option<f64>,
    pub(crate) at_hipparcos_epoch: CatalogStar,
    pub(crate) at_gaia_epoch: Place,
}

// RA and Dec in degrees, RA in [0, 360); proper motion in mas/yr, pmra including cos(Dec).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Place {
    pub(crate) ra: f64,
    pub(crate) dec: f64,
    pub(crate) pmra: f64,
    pub(crate) pmdec: f64,
}

impl Hip2Star {
    // Returns None when a field other than B-V is missing or does not parse, or the record
    // cannot be propagated, for example a Dec beyond the pole. Every hip2.dat record has a
    // parallax and a proper motion, so a field that fails to parse marks a corrupt line.
    pub(crate) fn parse(line: &str) -> Option<Self> {
        let bytes = line.as_bytes();
        // hip2.dat has no radial velocities, so the propagation uses zero.
        let star = CatalogStar {
            ra: Angle::from_radians(finite(bytes, RA_RAD)?),
            dec: Angle::from_radians(finite(bytes, DEC_RAD)?),
            pm_ra_cos_dec_mas_per_year: finite(bytes, PMRA_MAS_YR)?,
            pm_dec_mas_per_year: finite(bytes, PMDEC_MAS_YR)?,
            parallax_mas: finite(bytes, PARALLAX_MAS)?,
            radial_velocity_km_per_s: 0.0,
        };
        Some(Self {
            hip: field(bytes, HIP)?,
            hpmag: field::<f32>(bytes, HPMAG).filter(|m| m.is_finite())?,
            b_v: finite(bytes, B_V),
            at_hipparcos_epoch: star,
            at_gaia_epoch: propagate(&star).ok()?,
        })
    }

    // The record keeps the measured parallax. CatalogStar::propagate raises a small or
    // negative parallax to a floor, and that value is not a measurement.
    pub(crate) fn record(&self, source_id: i64) -> StarRecord {
        let (measured, place) = (&self.at_hipparcos_epoch, self.at_gaia_epoch);
        let mut record = intermediate::zeroed_record();
        record.source_id = source_id;
        (record.ra, record.dec) = (place.ra, place.dec);
        (record.pmra, record.pmdec) = (place.pmra, place.pmdec);
        record.parallax = measured.parallax_mas;
        record.mag = self.hpmag;
        record.flags = FLAG_SOURCE_HIPPARCOS;
        if measured.pm_ra_cos_dec_mas_per_year != 0.0 || measured.pm_dec_mas_per_year != 0.0 {
            record.flags |= FLAG_HAS_PROPER_MOTION;
        }
        if measured.parallax_mas != 0.0 {
            record.flags |= FLAG_HAS_PARALLAX;
        }
        record
    }
}

fn propagate(star: &CatalogStar) -> CoordResult<Place> {
    let moved = star.propagate(&epoch(HIPPARCOS_EPOCH), &epoch(GAIA_DR3_EPOCH))?;
    Ok(Place {
        ra: reduce_ra_deg(moved.ra.degrees()),
        dec: moved.dec.degrees(),
        pmra: moved.pm_ra_cos_dec_mas_per_year,
        pmdec: moved.pm_dec_mas_per_year,
    })
}

// A Julian epoch, read as TT.
fn epoch(julian_year: f64) -> TT {
    TT::from_julian_date(JulianDate::from_julian_year(julian_year))
}

fn finite(bytes: &[u8], range: (usize, usize)) -> Option<f64> {
    field::<f64>(bytes, range).filter(|v| v.is_finite())
}

fn field<T: FromStr>(bytes: &[u8], (start, end): (usize, usize)) -> Option<T> {
    let text = std::str::from_utf8(bytes.get(start..end)?).ok()?;
    text.trim().parse().ok()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    // A hip2.dat line with the fields this parser reads, in their byte columns, and B-V blank.
    pub(crate) fn hip2_line(hip: u32, astrometry: [&str; 5], hpmag: &str) -> String {
        let mut line = vec![b' '; 171];
        let mut put = |(start, end): (usize, usize), text: &str| {
            assert!(
                text.len() <= end - start,
                "{text:?} is wider than its field"
            );
            line[end - text.len()..end].copy_from_slice(text.as_bytes());
        };
        put(HIP, &hip.to_string());
        let [ra, dec, plx, pmra, pmdec] = astrometry;
        put(RA_RAD, ra);
        put(DEC_RAD, dec);
        put(PARALLAX_MAS, plx);
        put(PMRA_MAS_YR, pmra);
        put(PMDEC_MAS_YR, pmdec);
        put(HPMAG, hpmag);
        String::from_utf8(line).unwrap()
    }

    fn place(astrometry: [&str; 5]) -> Place {
        Hip2Star::parse(&hip2_line(1, astrometry, "5.0000"))
            .unwrap()
            .at_gaia_epoch
    }

    // hip2.dat fields of real stars, recovered exactly from the J1991.25 places behind the
    // Hipparcos records of the shipped catalog. The expected J2016.0 places are the output of
    // CatalogStar::propagate; erfa.pmsafe agrees to within 4 ulp in radians.
    #[test]
    fn test_real_stars_propagate_as_catalog_star() {
        let cases = [
            // HIP 32349, Sirius
            (
                [
                    "1.7678185359",
                    "-0.2916993748",
                    "379.21",
                    "-546.01",
                    "-1223.07",
                ],
                (101.28462149217447, -16.721551631891643),
                (-546.0340532982693, -1223.0592238146976),
            ),
            // HIP 91262, Vega
            (
                ["4.8735545728", "0.6768909262", "130.23", "200.94", "286.23"],
                (279.23588048446805, 38.784961079827326),
                (200.94554530484754, 286.2261062366564),
            ),
            // HIP 71683, alpha Centauri A
            (
                [
                    "3.8383352142",
                    "-1.0617735850",
                    "754.81",
                    "-3679.25",
                    "473.67",
                ],
                (219.86850784424806, -60.83187873405317),
                (-3678.8734575708913, 476.57996333912325),
            ),
            // HIP 5336, mu Cassiopeiae
            (
                [
                    "0.2976456257",
                    "0.9586086146",
                    "132.38",
                    "3422.23",
                    "-1598.93",
                ],
                (17.09476918322385, 54.91322830106515),
                (3421.293969616139, -1600.9300512312352),
            ),
            // HIP 11767, Polaris
            (
                ["0.6622851337", "1.5579531082", "7.54", "44.48", "-11.85"],
                (37.9699512436537, 89.26405624686883),
                (44.47507250741573, -11.868480329463171),
            ),
        ];
        for (astrometry, (ra, dec), (pmra, pmdec)) in cases {
            let expected = Place {
                ra,
                dec,
                pmra,
                pmdec,
            };
            assert_eq!(place(astrometry), expected, "{astrometry:?}");
        }
    }

    // 10"/yr toward the pole from Dec 89.99 moves 0.069 degrees in 24.75 years, so the star
    // crosses the pole and comes out on the far side at RA 1 rad + 180 degrees.
    #[test]
    fn test_star_near_the_pole_crosses_it() {
        let moved = place(["1.0000000000", "1.5706217940", "0.00", "0.00", "10000.00"]);
        let expected = Place {
            ra: 237.29577951308232,
            dec: 89.94125002552929,
            pmra: 7.846930551687332e-13,
            pmdec: -9999.98560217442,
        };
        assert_eq!(moved, expected);
    }

    #[test]
    fn test_ra_past_360_reduces_into_range() {
        let moved = place(["6.2831850000", "0.0000000000", "5.00", "1000.00", "0.00"]);
        assert_eq!((moved.ra, moved.dec), (0.006857399873122728, 0.0));
    }

    #[test]
    fn test_unparseable_parallax_or_motion_rejects_the_line() {
        let good = ["1.0000000000", "0.5000000000", "10.00", "10.00", "-10.00"];
        assert!(Hip2Star::parse(&hip2_line(1, good, "5.0000")).is_some());
        for field in 2..5 {
            let mut bad = good;
            bad[field] = "1O.00";
            let line = hip2_line(1, bad, "5.0000");
            assert!(Hip2Star::parse(&line).is_none(), "{bad:?}");
        }
    }

    #[test]
    fn test_record_keeps_measured_parallax() {
        let line = hip2_line(
            7,
            ["1.0000000000", "0.5000000000", "0.00", "3.00", "0.00"],
            "5.0000",
        );
        let record = Hip2Star::parse(&line).unwrap().record(-7);
        assert_eq!(
            (record.source_id, record.parallax, record.mag),
            (-7, 0.0, 5.0)
        );
        assert_eq!(record.flags, FLAG_SOURCE_HIPPARCOS | FLAG_HAS_PROPER_MOTION);
    }
}
