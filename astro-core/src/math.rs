#[inline]
pub fn fmod(x: f64, y: f64) -> f64 {
    libm::fmod(x, y)
}

#[inline]
pub fn vincenty_angular_separation(
    sin_lat1: f64,
    cos_lat1: f64,
    sin_lat2: f64,
    cos_lat2: f64,
    delta_lon: f64,
) -> f64 {
    let (sin_delta_lon, cos_delta_lon) = delta_lon.sin_cos();

    let num = ((cos_lat2 * sin_delta_lon).powi(2)
        + (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_delta_lon).powi(2))
    .sqrt();
    let den = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lon;

    num.atan2(den)
}
