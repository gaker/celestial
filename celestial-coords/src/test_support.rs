// Meeus prints his worked examples to a fixed number of decimals; compare at that precision.
pub(crate) fn rounded(x: f64, places: i32) -> f64 {
    let scale = libm::pow(10.0, places as f64);
    libm::round(x * scale) / scale
}
