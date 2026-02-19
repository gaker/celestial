use celestial_coords::frames::{EclipticPosition, GalacticPosition, ICRSPosition};
use celestial_coords::transforms::CoordinateFrame;
use celestial_coords::Distance;
use celestial_time::tt_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tt = tt_from_calendar(2024, 6, 21, 12, 0, 0.0); // Summer solstice 2024

    // --- Galactic coordinates ---
    // The galactic frame is fixed (IAU 1958 definition, refined by Hipparcos).
    // No epoch dependence — the rotation matrix from ICRS to galactic is constant.

    println!("=== Galactic Coordinate System ===\n");

    // Sagittarius A* (galactic center)
    let sgr_a = ICRSPosition::from_hours_degrees(17.76112, -29.00781)?;
    let gal = sgr_a.to_galactic(&tt)?;
    println!("Sgr A* (Galactic Center):");
    println!(
        "  ICRS:     RA = {:.5}h  Dec = {:+.5}°",
        sgr_a.ra().hours(),
        sgr_a.dec().degrees()
    );
    println!(
        "  Galactic: l = {:.4}°  b = {:+.4}°",
        gal.longitude().degrees(),
        gal.latitude().degrees()
    );
    println!(
        "  Near plane? {}  In bulge? {}\n",
        gal.is_near_galactic_plane(),
        gal.is_in_galactic_bulge()
    );

    // Polaris — far from the galactic plane
    let polaris = ICRSPosition::from_hours_degrees(2.53030, 89.26411)?;
    let gal = polaris.to_galactic(&tt)?;
    println!("Polaris:");
    println!(
        "  ICRS:     RA = {:.5}h  Dec = {:+.5}°",
        polaris.ra().hours(),
        polaris.dec().degrees()
    );
    println!(
        "  Galactic: l = {:.4}°  b = {:+.4}°",
        gal.longitude().degrees(),
        gal.latitude().degrees()
    );
    println!("  Near pole? {}\n", gal.is_near_galactic_pole());

    // Roundtrip: galactic → ICRS
    let gc = GalacticPosition::galactic_center();
    let icrs = gc.to_icrs(&tt)?;
    println!("Galactic center reference point:");
    println!(
        "  l = {:.1}°, b = {:.1}° → RA = {:.5}h, Dec = {:+.5}°\n",
        gc.longitude().degrees(),
        gc.latitude().degrees(),
        icrs.ra().hours(),
        icrs.dec().degrees()
    );

    // Notable galactic reference points
    let ngp = GalacticPosition::north_galactic_pole();
    let ngp_icrs = ngp.to_icrs(&tt)?;
    println!("North Galactic Pole:");
    println!(
        "  l = {:.1}°, b = {:+.1}° → RA = {:.5}h, Dec = {:+.5}°",
        ngp.longitude().degrees(),
        ngp.latitude().degrees(),
        ngp_icrs.ra().hours(),
        ngp_icrs.dec().degrees()
    );

    let anti = GalacticPosition::galactic_anticenter();
    let anti_icrs = anti.to_icrs(&tt)?;
    println!("Galactic Anticenter:");
    println!(
        "  l = {:.1}°, b = {:+.1}° → RA = {:.5}h, Dec = {:+.5}°\n",
        anti.longitude().degrees(),
        anti.latitude().degrees(),
        anti_icrs.ra().hours(),
        anti_icrs.dec().degrees()
    );

    // --- Ecliptic coordinates ---
    // Ecliptic longitude/latitude relative to the ecliptic plane.
    // Epoch-dependent: the ecliptic precesses with time.

    println!("=== Ecliptic Coordinate System ===\n");

    // Objects near the ecliptic tend to be solar system bodies.
    // Stars far from the ecliptic have large |β|.

    let vega = ICRSPosition::from_hours_degrees(18.61564, 38.78369)?;
    let ecl = vega.to_ecliptic(&tt)?;
    println!("Vega:");
    println!(
        "  ICRS:     RA = {:.5}h  Dec = {:+.5}°",
        vega.ra().hours(),
        vega.dec().degrees()
    );
    println!(
        "  Ecliptic: λ = {:.4}°  β = {:+.4}°",
        ecl.lambda().degrees(),
        ecl.beta().degrees()
    );
    println!("  Near ecliptic? {}", ecl.is_near_ecliptic_plane());
    println!(
        "  Mean obliquity = {:.6}°\n",
        ecl.mean_obliquity().degrees()
    );

    // Aldebaran — near the ecliptic (zodiac star)
    let aldebaran = ICRSPosition::from_hours_degrees(4.59868, 16.50930)?;
    let ecl = aldebaran.to_ecliptic(&tt)?;
    println!("Aldebaran (in Taurus, near ecliptic):");
    println!(
        "  ICRS:     RA = {:.5}h  Dec = {:+.5}°",
        aldebaran.ra().hours(),
        aldebaran.dec().degrees()
    );
    println!(
        "  Ecliptic: λ = {:.4}°  β = {:+.4}°",
        ecl.lambda().degrees(),
        ecl.beta().degrees()
    );
    println!("  Near ecliptic? {}\n", ecl.is_near_ecliptic_plane());

    // Ecliptic reference points
    let ve = EclipticPosition::vernal_equinox(tt);
    let ve_icrs = ve.to_icrs(&tt)?;
    println!("Vernal equinox (λ=0°, β=0°):");
    println!(
        "  → RA = {:.5}h, Dec = {:+.5}°",
        ve_icrs.ra().hours(),
        ve_icrs.dec().degrees()
    );

    let ss = EclipticPosition::summer_solstice(tt);
    let ss_icrs = ss.to_icrs(&tt)?;
    println!("Summer solstice (λ=90°, β=0°):");
    println!(
        "  → RA = {:.5}h, Dec = {:+.5}°\n",
        ss_icrs.ra().hours(),
        ss_icrs.dec().degrees()
    );

    // --- Angular separation across frames ---

    let m31 = ICRSPosition::from_hours_degrees(0.71222, 41.26917)?;
    let m31 = with_distance(m31, Distance::from_parsecs(778_000.0)?);
    let m31_gal = m31.to_galactic(&tt)?;

    println!("=== Distances ===\n");
    println!("M31 (Andromeda):");
    println!(
        "  ICRS:     RA = {:.5}h  Dec = {:+.5}°",
        m31.ra().hours(),
        m31.dec().degrees()
    );
    println!(
        "  Galactic: l = {:.4}°  b = {:+.4}°",
        m31_gal.longitude().degrees(),
        m31_gal.latitude().degrees()
    );
    println!(
        "  Distance: {:.0} kpc = {:.2} Mly = {:.0} AU",
        m31.distance().unwrap().parsecs() / 1000.0,
        m31.distance().unwrap().light_years() / 1e6,
        m31.distance().unwrap().au()
    );
    println!(
        "  Distance modulus: {:.2} mag",
        m31.distance().unwrap().distance_modulus()
    );
    println!(
        "  Local Group member? {}",
        m31.distance().unwrap().is_local_group()
    );

    let sep = sgr_a.angular_separation(&m31);
    println!("  Sgr A* ↔ M31 = {:.2}°", sep.degrees());

    Ok(())
}

fn with_distance(mut pos: ICRSPosition, d: Distance) -> ICRSPosition {
    pos.set_distance(d);
    pos
}
