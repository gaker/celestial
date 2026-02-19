use crate::error::{Error, Result};
use crate::observation::{
    decode_pier_side, IndatFile, IndatOption, MountType, Observation, PierSide, SiteParams,
};
use celestial_core::Angle;
use celestial_time::JulianDate;

pub fn parse_indat(content: &str) -> Result<IndatFile> {
    let mut header_lines = Vec::new();
    let mut options = Vec::new();
    let mut mount_type = MountType::GermanEquatorial;
    let mut site_and_date = None;
    let mut observations = Vec::new();

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if site_and_date.is_some() {
            observations.push(parse_observation_line(trimmed)?);
        } else {
            classify_line(
                trimmed,
                &mut header_lines,
                &mut options,
                &mut mount_type,
                &mut site_and_date,
            )?;
        }
    }

    let (site, date) = site_and_date.ok_or_else(|| Error::Parse("no site line found".into()))?;
    Ok(IndatFile {
        site,
        options,
        observations,
        mount_type,
        header_lines,
        date,
    })
}

fn classify_line(
    line: &str,
    headers: &mut Vec<String>,
    options: &mut Vec<IndatOption>,
    mount: &mut MountType,
    site: &mut Option<(SiteParams, JulianDate)>,
) -> Result<()> {
    if line.starts_with('!') {
        headers.push(line.to_string());
    } else if line.starts_with(':') {
        let (opt, maybe_mount) = parse_option(line)?;
        options.push(opt);
        if let Some(m) = maybe_mount {
            *mount = m;
        }
    } else if is_site_line(line) {
        *site = Some(parse_site_line(line)?);
    } else {
        headers.push(line.to_string());
    }
    Ok(())
}

fn is_site_line(line: &str) -> bool {
    let first = line.trim().as_bytes().first().copied().unwrap_or(0);
    first == b'+' || first == b'-' || first.is_ascii_digit()
}

fn parse_option(line: &str) -> Result<(IndatOption, Option<MountType>)> {
    let keyword = line.trim_start_matches(':').trim().to_uppercase();
    match keyword.as_str() {
        "NODA" => Ok((IndatOption::NoDA, None)),
        "ALLSKY" => Ok((IndatOption::AllSky, None)),
        "EQUINOX" => Ok((IndatOption::Equinox, None)),
        "EQUAT" => Ok((IndatOption::Equatorial, Some(MountType::GermanEquatorial))),
        "ALTAZ" => Ok((IndatOption::Altaz, Some(MountType::Altazimuth))),
        "ROTTEL" => Ok((IndatOption::RotatorTelescope, None)),
        "ROTNL" => Ok((IndatOption::RotatorNasmythLeft, None)),
        "ROTNR" => Ok((IndatOption::RotatorNasmythRight, None)),
        "ROTCL" => Ok((IndatOption::RotatorCoudeLeft, None)),
        "ROTCR" => Ok((IndatOption::RotatorCoudeRight, None)),
        s if s.starts_with("GIMBAL") => parse_gimbal(s),
        _ => Err(Error::Parse(format!("unknown option: {}", keyword))),
    }
}

fn parse_gimbal(s: &str) -> Result<(IndatOption, Option<MountType>)> {
    let parts: Vec<&str> = s.split_whitespace().collect();
    if parts.len() != 4 {
        return Err(Error::Parse(format!(
            "GIMBAL requires 3 angles, got: {}",
            s
        )));
    }
    let z = parse_f64(parts[1], "gimbal z")?;
    let y = parse_f64(parts[2], "gimbal y")?;
    let x = parse_f64(parts[3], "gimbal x")?;
    Ok((
        IndatOption::Gimbal {
            z: Angle::from_degrees(z),
            y: Angle::from_degrees(y),
            x: Angle::from_degrees(x),
        },
        None,
    ))
}

fn parse_site_line(line: &str) -> Result<(SiteParams, JulianDate)> {
    let p: Vec<&str> = line.split_whitespace().collect();
    if p.len() < 12 {
        return Err(Error::Parse(format!(
            "site line needs 12 fields, got {}",
            p.len()
        )));
    }
    let lat = parse_dms_latitude(
        parse_f64(p[0], "lat_d")?,
        parse_f64(p[1], "lat_m")?,
        parse_f64(p[2], "lat_s")?,
    );
    let date = build_julian_date(&p[3..6])?;
    let site = SiteParams {
        latitude: lat,
        longitude: Angle::from_degrees(0.0),
        temperature: parse_f64(p[6], "temp")?,
        pressure: parse_f64(p[7], "pressure")?,
        elevation: parse_f64(p[8], "elevation")?,
        humidity: parse_f64(p[9], "humidity")?,
        wavelength: parse_f64(p[10], "wavelength")?,
        lapse_rate: parse_f64(p[11], "lapse_rate")?,
    };
    Ok((site, date))
}

fn build_julian_date(parts: &[&str]) -> Result<JulianDate> {
    let year: i32 = parts[0]
        .parse()
        .map_err(|e| Error::Parse(format!("year: {}", e)))?;
    let month: u8 = parts[1]
        .parse()
        .map_err(|e| Error::Parse(format!("month: {}", e)))?;
    let day: u8 = parts[2]
        .parse()
        .map_err(|e| Error::Parse(format!("day: {}", e)))?;
    Ok(JulianDate::from_calendar(year, month, day, 0, 0, 0.0))
}

fn parse_dms_latitude(d: f64, m: f64, s: f64) -> Angle {
    let sign = if d < 0.0 || (d == 0.0 && d.is_sign_negative()) {
        -1.0
    } else {
        1.0
    };
    let deg = d.abs() + m / 60.0 + s / 3600.0;
    Angle::from_degrees(sign * deg)
}

fn parse_observation_line(line: &str) -> Result<Observation> {
    let p: Vec<&str> = line.split_whitespace().collect();
    if p.len() < 14 {
        return Err(Error::Parse(format!(
            "obs line needs 14 fields, got {}",
            p.len()
        )));
    }
    let catalog_ra = parse_ra(&p[0..3])?;
    let catalog_dec = parse_dec_as_angle(&p[3..6])?;
    let tel_ra = parse_ra(&p[6..9])?;
    let raw_tel_dec_deg = parse_dec_raw(&p[9..12])?;
    let lst = parse_lst(&p[12..14])?;

    let (observed_dec, pier_side) = decode_pier_side(raw_tel_dec_deg);
    let observed_ra = compute_observed_ra(tel_ra, &pier_side);
    let commanded_ha = (lst - catalog_ra).wrapped();
    let actual_ha = (lst - observed_ra).wrapped();

    Ok(Observation {
        catalog_ra,
        catalog_dec,
        observed_ra,
        observed_dec,
        lst,
        commanded_ha,
        actual_ha,
        pier_side,
        masked: false,
    })
}

fn compute_observed_ra(tel_ra: Angle, pier_side: &PierSide) -> Angle {
    match pier_side {
        PierSide::West => (tel_ra + Angle::from_hours(12.0)).normalized(),
        _ => tel_ra,
    }
}

fn parse_ra(parts: &[&str]) -> Result<Angle> {
    let h = parse_f64(parts[0], "ra_h")?;
    let m = parse_f64(parts[1], "ra_m")?;
    let s = parse_f64(parts[2], "ra_s")?;
    Ok(Angle::from_hours(h + m / 60.0 + s / 3600.0))
}

fn parse_dec_raw(parts: &[&str]) -> Result<f64> {
    let d = parse_f64(parts[0], "dec_d")?;
    let m = parse_f64(parts[1], "dec_m")?;
    let s = parse_f64(parts[2], "dec_s")?;
    let sign = if d < 0.0 || (d == 0.0 && parts[0].starts_with('-')) {
        -1.0
    } else {
        1.0
    };
    Ok(sign * (d.abs() + m / 60.0 + s / 3600.0))
}

fn parse_dec_as_angle(parts: &[&str]) -> Result<Angle> {
    Ok(Angle::from_degrees(parse_dec_raw(parts)?))
}

fn parse_lst(parts: &[&str]) -> Result<Angle> {
    let h = parse_f64(parts[0], "lst_h")?;
    let m = parse_f64(parts[1], "lst_m")?;
    Ok(Angle::from_hours(h + m / 60.0))
}

fn parse_lst_hms(parts: &[&str]) -> Result<Angle> {
    let h = parse_f64(parts[0], "lst_h")?;
    let m = parse_f64(parts[1], "lst_m")?;
    let s = parse_f64(parts[2], "lst_s")?;
    Ok(Angle::from_hours(h + m / 60.0 + s / 3600.0))
}

pub fn parse_coordinates(args: &[&str]) -> Result<(Angle, Angle)> {
    match args.len() {
        6 => {
            let ra = parse_ra(&args[0..3])?;
            let dec = parse_dec_as_angle(&args[3..6])?;
            Ok((ra, dec))
        }
        2 => {
            let ra_hours = parse_f64(args[0], "ra_hours")?;
            let dec_deg = parse_f64(args[1], "dec_degrees")?;
            Ok((Angle::from_hours(ra_hours), Angle::from_degrees(dec_deg)))
        }
        _ => Err(Error::Parse(
            "expected 6 args (h m s d m s) or 2 args (decimal_hours decimal_degrees)".into(),
        )),
    }
}

pub fn parse_lst_args(args: &[&str]) -> Result<Angle> {
    match args.len() {
        3 => parse_lst_hms(args),
        1 => {
            let hours = parse_f64(args[0], "lst_hours")?;
            Ok(Angle::from_hours(hours))
        }
        _ => Err(Error::Parse(
            "expected 3 args (h m s) or 1 arg (decimal_hours)".into(),
        )),
    }
}

fn parse_f64(s: &str, field: &str) -> Result<f64> {
    s.parse::<f64>()
        .map_err(|e| Error::Parse(format!("{}: {}", field, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_DAT: &str = "\
ASCOM Mount
:NODA
:EQUAT
+39 00 26 2024 7 14 29.20 987.00 231.65  0.94 0.5500 0.0065
21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130
23 46 02.2988 +77 38 38.725 11 26 17.6308 +104 03 28.734  16 24.711";

    #[test]
    fn parse_header_lines() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        assert_eq!(indat.header_lines.len(), 1);
        assert_eq!(indat.header_lines[0], "ASCOM Mount");
    }

    #[test]
    fn parse_options_count_and_values() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        assert_eq!(indat.options.len(), 2);
        assert_eq!(indat.options[0], IndatOption::NoDA);
        assert_eq!(indat.options[1], IndatOption::Equatorial);
    }

    #[test]
    fn parse_mount_type() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        assert_eq!(indat.mount_type, MountType::GermanEquatorial);
    }

    #[test]
    fn parse_observation_count() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        assert_eq!(indat.observations.len(), 2);
    }

    #[test]
    fn parse_site_latitude() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let expected = Angle::from_degrees(39.0 + 0.0 / 60.0 + 26.0 / 3600.0);
        assert_eq!(indat.site.latitude, expected);
    }

    #[test]
    fn parse_site_conditions() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        assert_eq!(indat.site.temperature, 29.20);
        assert_eq!(indat.site.pressure, 987.00);
        assert_eq!(indat.site.elevation, 231.65);
        assert_eq!(indat.site.humidity, 0.94);
        assert_eq!(indat.site.wavelength, 0.5500);
        assert_eq!(indat.site.lapse_rate, 0.0065);
        assert_eq!(indat.site.longitude, Angle::from_degrees(0.0));
    }

    #[test]
    fn parse_site_date() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let expected = JulianDate::from_calendar(2024, 7, 14, 0, 0, 0.0);
        assert_eq!(indat.date, expected);
    }

    #[test]
    fn first_obs_catalog_ra() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let expected = Angle::from_hours(21.0 + 43.0 / 60.0 + 18.4460 / 3600.0);
        assert_eq!(obs.catalog_ra, expected);
    }

    #[test]
    fn first_obs_catalog_dec() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let expected = Angle::from_degrees(72.0 + 29.0 / 60.0 + 8.368 / 3600.0);
        assert_eq!(obs.catalog_dec, expected);
    }

    #[test]
    fn first_obs_pier_side_west() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        assert_eq!(obs.pier_side, PierSide::West);
    }

    #[test]
    fn second_obs_pier_side_west() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[1];
        assert_eq!(obs.pier_side, PierSide::West);
    }

    #[test]
    fn first_obs_lst() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let expected = Angle::from_hours(16.0 + 23.130 / 60.0);
        assert_eq!(obs.lst, expected);
    }

    #[test]
    fn first_obs_observed_ra_includes_12h_flip() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let tel_ra = Angle::from_hours(9.0 + 28.0 / 60.0 + 59.9527 / 3600.0);
        let expected = (tel_ra + Angle::from_hours(12.0)).normalized();
        assert_eq!(obs.observed_ra, expected);
    }

    #[test]
    fn first_obs_observed_dec_decoded() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let raw_deg = 109.0 + 20.0 / 60.0 + 6.469 / 3600.0;
        let (expected_dec, _) = decode_pier_side(raw_deg);
        assert_eq!(obs.observed_dec, expected_dec);
    }

    #[test]
    fn first_obs_commanded_ha() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let expected = obs.lst - obs.catalog_ra;
        assert_eq!(obs.commanded_ha, expected);
    }

    #[test]
    fn first_obs_actual_ha() {
        let indat = parse_indat(SIMPLE_DAT).unwrap();
        let obs = &indat.observations[0];
        let expected = obs.lst - obs.observed_ra;
        assert_eq!(obs.actual_ha, expected);
    }

    #[test]
    fn parse_dms_latitude_positive() {
        let lat = parse_dms_latitude(39.0, 0.0, 26.0);
        assert_eq!(lat, Angle::from_degrees(39.0 + 26.0 / 3600.0));
    }

    #[test]
    fn parse_dms_latitude_negative() {
        let lat = parse_dms_latitude(-33.0, 15.0, 30.0);
        let expected = Angle::from_degrees(-(33.0 + 15.0 / 60.0 + 30.0 / 3600.0));
        assert_eq!(lat, expected);
    }

    #[test]
    fn option_case_insensitive() {
        let (opt, _) = parse_option(":noda").unwrap();
        assert_eq!(opt, IndatOption::NoDA);
    }

    #[test]
    fn option_altaz_sets_mount() {
        let (opt, mount) = parse_option(":ALTAZ").unwrap();
        assert_eq!(opt, IndatOption::Altaz);
        assert_eq!(mount, Some(MountType::Altazimuth));
    }

    #[test]
    fn empty_content_errors() {
        let result = parse_indat("");
        assert!(result.is_err());
    }

    #[test]
    fn no_site_line_errors() {
        let result = parse_indat("!comment\n:NODA\n");
        assert!(result.is_err());
    }

    #[test]
    fn parse_coordinates_sexagesimal() {
        let args = vec!["12", "30", "00", "+45", "00", "00"];
        let (ra, dec) = parse_coordinates(&args).unwrap();
        assert_eq!(ra, Angle::from_hours(12.0 + 30.0 / 60.0));
        assert_eq!(dec, Angle::from_degrees(45.0));
    }

    #[test]
    fn parse_coordinates_sexagesimal_negative_dec() {
        let args = vec!["6", "0", "0", "-30", "15", "30"];
        let (ra, dec) = parse_coordinates(&args).unwrap();
        assert_eq!(ra, Angle::from_hours(6.0));
        let expected_dec = Angle::from_degrees(-(30.0 + 15.0 / 60.0 + 30.0 / 3600.0));
        assert_eq!(dec, expected_dec);
    }

    #[test]
    fn parse_coordinates_decimal() {
        let args = vec!["12.5", "45.0"];
        let (ra, dec) = parse_coordinates(&args).unwrap();
        assert_eq!(ra, Angle::from_hours(12.5));
        assert_eq!(dec, Angle::from_degrees(45.0));
    }

    #[test]
    fn parse_coordinates_wrong_arg_count() {
        let args = vec!["12", "30", "00", "+45"];
        assert!(parse_coordinates(&args).is_err());
    }

    #[test]
    fn parse_coordinates_zero_args() {
        let args: Vec<&str> = vec![];
        assert!(parse_coordinates(&args).is_err());
    }

    #[test]
    fn parse_lst_args_hms() {
        let args = vec!["14", "30", "0"];
        let lst = parse_lst_args(&args).unwrap();
        assert_eq!(lst, Angle::from_hours(14.5));
    }

    #[test]
    fn parse_lst_args_decimal() {
        let args = vec!["14.5"];
        let lst = parse_lst_args(&args).unwrap();
        assert_eq!(lst, Angle::from_hours(14.5));
    }

    #[test]
    fn parse_lst_args_wrong_count() {
        let args = vec!["14", "30"];
        assert!(parse_lst_args(&args).is_err());
    }
}
