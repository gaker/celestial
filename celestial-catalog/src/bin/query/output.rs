//! Table, JSON and CSV output of cone search results.

use celestial_catalog::query::catalog::FLAG_SOURCE_HIPPARCOS;
use celestial_catalog::query::cone::ConeSearchResult;
use celestial_core::angle::{Angle, DmsFmt, HmsFmt};

#[derive(serde::Serialize)]
struct JsonStar {
    source_id: i64,
    ra_deg: f64,
    dec_deg: f64,
    mag: f32,
    distance_deg: f64,
    source: &'static str,
}

pub(crate) fn print_table(results: &[ConeSearchResult], raw: bool) {
    for (i, result) in results.iter().enumerate() {
        println!("{}", table_row(i + 1, result, raw));
    }
    if results.is_empty() {
        println!("No stars found matching the search criteria.");
    } else {
        println!("\nTotal results: {}", results.len());
    }
}

fn table_row(number: usize, result: &ConeSearchResult, raw: bool) -> String {
    let (ra, dec) = if raw {
        (
            format!("{:.6}°", result.ra_deg),
            format!("{:+.6}°", result.dec_deg),
        )
    } else {
        (
            HmsFmt { frac_digits: 4 }.fmt(Angle::from_degrees(result.ra_deg)),
            DmsFmt { frac_digits: 4 }.fmt(Angle::from_degrees(result.dec_deg)),
        )
    };
    format!(
        "{:4}: {:>20} RA={} Dec={} Mag={:5.2} Dist={:.4}° Source={}",
        number,
        result.star.source_id,
        ra,
        dec,
        result.star.mag,
        result.distance_deg,
        source_name(result.star.flags)
    )
}

pub(crate) fn print_json(results: &[ConeSearchResult]) -> anyhow::Result<()> {
    let stars: Vec<JsonStar> = results
        .iter()
        .map(|r| JsonStar {
            source_id: r.star.source_id,
            ra_deg: r.ra_deg,
            dec_deg: r.dec_deg,
            mag: r.star.mag,
            distance_deg: r.distance_deg,
            source: source_name(r.star.flags),
        })
        .collect();
    println!("{}", serde_json::to_string_pretty(&stars)?);
    Ok(())
}

pub(crate) fn print_csv(results: &[ConeSearchResult]) {
    println!("source_id,ra_deg,dec_deg,mag,distance_deg,source");
    for r in results {
        println!(
            "{},{},{},{},{},{}",
            r.star.source_id,
            r.ra_deg,
            r.dec_deg,
            r.star.mag,
            r.distance_deg,
            source_name(r.star.flags)
        );
    }
}

fn source_name(flags: u16) -> &'static str {
    if (flags & FLAG_SOURCE_HIPPARCOS) != 0 {
        "HIPPARCOS"
    } else {
        "GAIA"
    }
}
