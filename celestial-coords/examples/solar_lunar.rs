use celestial_coords::frames::{
    HeliographicCarrington, HeliographicStonyhurst, SelenographicPosition,
};
use celestial_coords::lunar::compute_lunar_orientation;
use celestial_coords::solar::{
    carrington_rotation_number, compute_solar_orientation, sun_earth_distance,
};
use celestial_time::tt_from_calendar;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- Solar orientation ---
    // B0: heliographic latitude of the sub-solar point (disk center tilt)
    // L0: Carrington longitude of central meridian
    // P:  position angle of the solar rotation axis

    let solstice = tt_from_calendar(2024, 6, 21, 12, 0, 0.0);
    let equinox = tt_from_calendar(2024, 3, 20, 3, 6, 0.0);

    println!("=== Solar Orientation ===\n");

    for (label, epoch) in [
        ("Summer solstice 2024", &solstice),
        ("Vernal equinox 2024", &equinox),
    ] {
        let orient = compute_solar_orientation(epoch);
        let dist = sun_earth_distance(epoch);
        let cr = carrington_rotation_number(epoch);

        println!("{label}:");
        println!(
            "  B0 = {:+.4}°  (disk center heliographic latitude)",
            orient.b0.degrees()
        );
        println!(
            "  L0 = {:.4}°   (central meridian Carrington longitude)",
            orient.l0.degrees()
        );
        println!(
            "  P  = {:+.4}°  (solar north pole position angle)",
            orient.p.degrees()
        );
        println!("  Sun-Earth = {:.6} AU", dist);
        println!("  Carrington rotation #{cr}\n");
    }

    // --- Heliographic coordinates ---
    // Two systems: Stonyhurst (Earth-fixed meridian) and Carrington (rotating with Sun).
    // Active regions are commonly reported in both.

    println!("=== Heliographic Coordinates ===\n");

    // A sunspot near disk center
    let spot_stonyhurst = HeliographicStonyhurst::from_degrees(15.0, -10.0)?;
    let spot_carrington = spot_stonyhurst.to_carrington(&solstice)?;

    println!("Sunspot at Stonyhurst (15°N, 10°W):");
    println!(
        "  Stonyhurst: lat = {:.1}°, lon = {:.1}°",
        spot_stonyhurst.latitude().degrees(),
        spot_stonyhurst.longitude().degrees()
    );
    println!(
        "  Carrington: lat = {:.1}°, lon = {:.2}°",
        spot_carrington.latitude().degrees(),
        spot_carrington.longitude().degrees()
    );

    // Roundtrip
    let back = spot_carrington.to_stonyhurst(&solstice)?;
    println!(
        "  Roundtrip:  lat = {:.1}°, lon = {:.1}°\n",
        back.latitude().degrees(),
        back.longitude().degrees()
    );

    // Disk center in Stonyhurst is always (B0, 0) by definition
    let center = HeliographicStonyhurst::disk_center(&solstice);
    println!("Disk center at solstice:");
    println!(
        "  lat = {:+.4}°, lon = {:.4}°\n",
        center.latitude().degrees(),
        center.longitude().degrees()
    );

    // Carrington rotation number (fractional)
    let cr = HeliographicCarrington::carrington_rotation_number(&solstice);
    println!("Carrington rotation number: {:.3}", cr);
    println!("  (integer part = rotation count, fraction = phase within rotation)\n");

    // --- Lunar orientation and libration ---
    // The Moon's apparent tilt changes due to optical libration.
    // Longitude libration: ±7.9° (due to eccentric orbit)
    // Latitude libration: ±6.7° (due to orbital inclination)

    println!("=== Lunar Orientation ===\n");

    let dates = [
        (
            "New Moon ~2024-01-11",
            tt_from_calendar(2024, 1, 11, 12, 0, 0.0),
        ),
        (
            "Full Moon ~2024-01-25",
            tt_from_calendar(2024, 1, 25, 17, 0, 0.0),
        ),
        (
            "New Moon ~2024-02-09",
            tt_from_calendar(2024, 2, 9, 23, 0, 0.0),
        ),
        (
            "Full Moon ~2024-06-22",
            tt_from_calendar(2024, 6, 22, 1, 0, 0.0),
        ),
    ];

    for (label, epoch) in &dates {
        let orient = compute_lunar_orientation(epoch);

        println!("{label}:");
        println!(
            "  Optical libration: lon = {:+.3}°, lat = {:+.3}°",
            orient.optical_libration.longitude.degrees(),
            orient.optical_libration.latitude.degrees()
        );
        println!("  Position angle: {:+.3}°", orient.position_angle.degrees());
        println!();
    }

    // --- Selenographic coordinates ---
    // Coordinates on the lunar surface. Used for crater positions, landing sites, etc.

    println!("=== Selenographic Coordinates ===\n");

    // Apollo 11 landing site: Sea of Tranquility
    let apollo11 = SelenographicPosition::from_degrees(0.6744, 23.4730)?;
    println!("Apollo 11 landing site (Tranquility Base):");
    println!(
        "  lat = {:.4}°, lon = {:.4}°",
        apollo11.latitude().degrees(),
        apollo11.longitude().degrees()
    );
    println!(
        "  Visible from Earth? {}\n",
        apollo11.is_visible_from_earth(&solstice)
    );

    // Tycho crater
    let tycho = SelenographicPosition::from_degrees(-43.31, -11.36)?;
    println!("Tycho crater:");
    println!(
        "  lat = {:.2}°, lon = {:.2}°",
        tycho.latitude().degrees(),
        tycho.longitude().degrees()
    );
    println!(
        "  Visible from Earth? {}\n",
        tycho.is_visible_from_earth(&solstice)
    );

    // South Pole-Aitken Basin (far side)
    let spa = SelenographicPosition::from_degrees(-53.0, 169.0)?;
    println!("South Pole-Aitken Basin (far side):");
    println!(
        "  lat = {:.1}°, lon = {:.1}°",
        spa.latitude().degrees(),
        spa.longitude().degrees()
    );
    println!(
        "  Visible from Earth? {}\n",
        spa.is_visible_from_earth(&solstice)
    );

    // Sub-Earth point
    let sub_earth = SelenographicPosition::sub_earth_point(&solstice)?;
    println!("Sub-Earth point at solstice:");
    println!(
        "  lat = {:+.3}°, lon = {:+.3}°",
        sub_earth.latitude().degrees(),
        sub_earth.longitude().degrees()
    );

    // Angular separation between Apollo 11 and Tycho
    let sep = apollo11.angular_separation(&tycho);
    println!(
        "\nApollo 11 ↔ Tycho = {:.2}° on lunar surface",
        sep.degrees()
    );

    // Reference points
    let nearside = SelenographicPosition::nearside_center();
    let farside = SelenographicPosition::farside_center();
    println!(
        "\nNearside center: ({:.0}°, {:.0}°)",
        nearside.latitude().degrees(),
        nearside.longitude().degrees()
    );
    println!(
        "Farside center:  ({:.0}°, {:.0}°)",
        farside.latitude().degrees(),
        farside.longitude().degrees()
    );

    Ok(())
}
