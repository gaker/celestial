use anyhow::Result;

use celestial_images::formats::image::Image;
use celestial_solver::fit_wcs::display::SolveDisplay;
use celestial_solver::metadata::time::jd_from_image;
use celestial_solver::observation::geodetic::geodetic_from_image;
use celestial_solver::solve::SolveResult;
use celestial_time::scales::utc::UTC;

pub(crate) fn report(img: &Image, result: &SolveResult) -> Result<String> {
    let software = format!("Celestial Image Solver v{}", env!("CARGO_PKG_VERSION"));
    let mut display = SolveDisplay::new(result)
        .with_creation_time(UTC::now()?)
        .with_software(&software);
    if let Some(start) = observation_time(img) {
        display = display.with_observation_time(start);
    }
    if let Some(geo) = geodetic_from_image(img) {
        display = display.with_geodetic(geo.lon_deg, geo.lat_deg, geo.alt_m);
    }
    Ok(display.to_string())
}

// The epoch the solve used, which `jd_from_image` reads from MJD-OBS or DATE-OBS.
fn observation_time(img: &Image) -> Option<UTC> {
    match jd_from_image(img) {
        Ok(jd) => Some(UTC::from_julian_date(jd)),
        Err(e) => {
            log::warn!("the report omits the observation time: {e}");
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::keywords::Keyword;
    use celestial_images::formats::pixel_data::PixelData;
    use celestial_time::julian::JulianDate;

    // MJD-OBS sets the epoch of the solve, so the report prints that instant, including
    // when DATE-OBS is a date alone or names another time.
    #[test]
    fn observation_time_is_the_epoch_of_the_solve() {
        for date_obs in ["1998-04-20", "1998-04-20T18:00:00"] {
            let mut img = Image::new(PixelData::F32(vec![0.0; 4]), [2, 2]);
            img.set_keyword(Keyword::string("DATE-OBS", date_obs));
            img.set_keyword(Keyword::real("MJD-OBS", 50923.7765625));
            let epoch = UTC::from_julian_date(JulianDate::new(2_400_000.5, 50923.7765625));
            assert_eq!(observation_time(&img), Some(epoch), "{date_obs}");
        }
    }

    #[test]
    fn observation_time_without_mjd_obs_is_date_obs_and_time_obs() {
        let mut img = Image::new(PixelData::F32(vec![0.0; 4]), [2, 2]);
        img.set_keyword(Keyword::string("DATE-OBS", "1998-04-20"));
        img.set_keyword(Keyword::string("TIME-OBS", "18:38:15"));
        let expected = "1998-04-20T18:38:15".parse::<UTC>().unwrap();
        assert_eq!(observation_time(&img), Some(expected));
    }
}
