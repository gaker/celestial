use celestial_core::constants::DAYS_PER_JULIAN_CENTURY;

// Days from J2000.0. libm's fmod loops once per bit of x/y, and the fundamental
// arguments grow with |t|, so the cost can drift with epoch; five centuries either
// side shows how much.
pub const EPOCHS: [(&str, f64); 4] = [
    ("j2000", 0.0),
    ("2026", 0.2675 * DAYS_PER_JULIAN_CENTURY),
    ("plus_5c", 5.0 * DAYS_PER_JULIAN_CENTURY),
    ("minus_5c", -5.0 * DAYS_PER_JULIAN_CENTURY),
];
