use celestial_coords::eop::record::EopRecord;
use celestial_coords::frames::{CIRSPosition, ICRSPosition};
use celestial_coords::transforms::CoordinateFrame;
use celestial_coords::{Distance, EopProvider, Location};
use celestial_time::{tt_from_calendar, ToTAI, ToUTC};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- Setup: observer, time, EOP ---

    // McDonald Observatory, Texas
    let observer = Location::from_degrees(30.6714, -104.0225, 2070.0)?;

    // 2023-06-15 03:00:00 TT (nighttime in Texas)
    let tt = tt_from_calendar(2023, 6, 15, 3, 0, 0.0);
    let utc = tt.to_tai()?.to_utc()?;
    println!("Epoch: TT JD = {:.6}", tt.to_julian_date().to_f64());
    println!("       UTC JD = {:.6}", utc.to_julian_date().to_f64());
    println!(
        "Observer: McDonald Observatory ({:.4}°N, {:.4}°W, {:.0}m)\n",
        observer.latitude_angle().degrees(),
        -observer.longitude_angle().degrees(),
        2070.0
    );

    // Realistic EOP data around this date
    let eop_records = vec![
        EopRecord::new(60108.0, 0.183, 0.343, -0.0298, 0.00071)?.with_cip_offsets(0.198, -0.102)?,
        EopRecord::new(60109.0, 0.184, 0.341, -0.0305, 0.00073)?.with_cip_offsets(0.201, -0.105)?,
        EopRecord::new(60110.0, 0.185, 0.339, -0.0312, 0.00075)?.with_cip_offsets(0.204, -0.108)?,
        EopRecord::new(60111.0, 0.186, 0.337, -0.0319, 0.00077)?.with_cip_offsets(0.207, -0.111)?,
        EopRecord::new(60112.0, 0.187, 0.335, -0.0326, 0.00079)?.with_cip_offsets(0.210, -0.114)?,
    ];

    let provider = EopProvider::from_records(eop_records)?;

    // --- Vega: ICRS catalog position ---

    let vega = ICRSPosition::from_hours_degrees(18.61564, 38.78369)?;
    println!("=== Vega ===");
    println!(
        "ICRS:  RA = {:.5}h  Dec = {:+.5}°",
        vega.ra().hours(),
        vega.dec().degrees()
    );

    // ICRS → CIRS (applies frame bias, precession, nutation, aberration, light deflection)
    let cirs = CIRSPosition::from_icrs(&vega, &tt)?;
    println!(
        "CIRS:  RA = {:.5}h  Dec = {:+.5}°",
        cirs.ra().hours(),
        cirs.dec().degrees()
    );

    // CIRS → Hour Angle / Declination
    let mjd = tt.to_julian_date().to_f64() - 2400000.5;
    let mut eop = provider.get(mjd)?;
    eop.compute_s_prime();
    let delta_t = eop.ut1_utc; // UT1-UTC in seconds

    let ha_pos = cirs.to_hour_angle(&observer, -delta_t)?;
    println!(
        "HA/Dec: HA = {:.5}h  Dec = {:+.5}°",
        ha_pos.hour_angle().hours(),
        ha_pos.declination().degrees()
    );

    // Hour Angle → Topocentric (Az/El)
    let topo = ha_pos.to_topocentric()?;
    println!(
        "Topo:  Az = {:.4}°  El = {:.4}°",
        topo.azimuth().degrees(),
        topo.elevation().degrees()
    );
    println!("       Air mass = {:.3}", topo.air_mass());

    // Atmospheric refraction at standard conditions
    let refraction = topo.atmospheric_refraction(1013.25, 15.0, 0.3, 0.55);
    println!("       Refraction = {:.2}\"", refraction.degrees() * 3600.0);
    println!();

    // --- CIRS → TIRS → ITRS (full EOP chain) ---

    println!("=== Full IAU 2000/2006 chain ===");
    let tirs = cirs.to_tirs(&eop)?;
    println!("TIRS:  ({:.9}, {:.9}, {:.9})", tirs.x(), tirs.y(), tirs.z());

    let itrs = tirs.to_itrs(&tt, &eop)?;
    println!("ITRS:  ({:.9}, {:.9}, {:.9})", itrs.x(), itrs.y(), itrs.z());
    println!();

    // --- Roundtrip precision test: ICRS → CIRS → ICRS ---

    let roundtrip = cirs.to_icrs(&tt)?;
    let ra_diff_mas = (roundtrip.ra().degrees() - vega.ra().degrees()).abs() * 3600.0 * 1000.0;
    let dec_diff_mas = (roundtrip.dec().degrees() - vega.dec().degrees()).abs() * 3600.0 * 1000.0;
    println!("=== Roundtrip precision (ICRS → CIRS → ICRS) ===");
    println!("  ΔRA  = {:.6} mas", ra_diff_mas);
    println!("  ΔDec = {:.6} mas", dec_diff_mas);
    println!();

    // --- Sirius: with distance ---

    let sirius = ICRSPosition::from_hours_degrees(6.75248, -16.71612)?;
    let sirius = sirius.with_distance_value(Distance::from_parsecs(2.637)?);
    println!(
        "=== Sirius (d = {:.3} pc = {:.2} ly) ===",
        sirius.distance().unwrap().parsecs(),
        sirius.distance().unwrap().light_years()
    );

    let cirs = CIRSPosition::from_icrs(&sirius, &tt)?;
    let ha_pos = cirs.to_hour_angle(&observer, -delta_t)?;
    let topo = ha_pos.to_topocentric()?;
    println!(
        "ICRS:  RA = {:.5}h  Dec = {:+.5}°",
        sirius.ra().hours(),
        sirius.dec().degrees()
    );
    println!(
        "Topo:  Az = {:.4}°  El = {:.4}°",
        topo.azimuth().degrees(),
        topo.elevation().degrees()
    );

    if topo.elevation().degrees() < 0.0 {
        println!("       (below horizon)");
    } else {
        println!("       Air mass = {:.3}", topo.air_mass());
    }
    println!();

    // --- Angular separation ---

    let betelgeuse = ICRSPosition::from_hours_degrees(5.91953, 7.40706)?;
    let rigel = ICRSPosition::from_hours_degrees(5.24230, -8.20164)?;
    let sep = betelgeuse.angular_separation(&rigel);
    println!("=== Angular separation ===");
    println!("Betelgeuse ↔ Rigel = {:.4}°", sep.degrees());

    Ok(())
}

trait WithDistanceValue {
    fn with_distance_value(self, distance: Distance) -> Self;
}

impl WithDistanceValue for ICRSPosition {
    fn with_distance_value(mut self, distance: Distance) -> Self {
        self.set_distance(distance);
        self
    }
}
