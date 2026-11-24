use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::cone::{cone_search, ConeSearchParams, ConeSearchResult};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let Some(path) = std::env::args().nth(1) else {
        return Err("Usage: cone_search <catalog.bin>".into());
    };

    let catalog = Catalog::open(&path)?;
    println!("{}", catalog.header());

    let params = ConeSearchParams {
        ra_deg: 83.633,
        dec_deg: -5.375,
        radius_deg: 0.5,
        max_mag: Some(10.0),
        max_results: Some(20),
        epoch: None,
    };

    let results = cone_search(&catalog, &params)?;
    print_results(&params, &results);
    Ok(())
}

fn print_results(params: &ConeSearchParams, results: &[ConeSearchResult]) {
    println!(
        "\n{} stars within {:.1}° of ({:.3}, {:.3}):\n",
        results.len(),
        params.radius_deg,
        params.ra_deg,
        params.dec_deg,
    );

    for r in results {
        println!(
            "  {:>20}  RA {:.6}°  Dec {:+.6}°  mag {:.2}  dist {:.4}°",
            r.star.source_id, r.ra_deg, r.dec_deg, r.star.mag, r.distance_deg,
        );
    }
}
