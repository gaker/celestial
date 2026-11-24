// Synthetic star images shared by the detection and solve tests.

// A circular Gaussian of height 1 and width `sigma` centred on `(cx, cy)`, evaluated at
// each pixel of a `w` x `h` image in row-major order.
pub(crate) fn gaussian_profile(
    w: usize,
    h: usize,
    (cx, cy): (f64, f64),
    sigma: f64,
) -> impl Iterator<Item = f64> {
    (0..h).flat_map(move |y| {
        (0..w).map(move |x| {
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            libm::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma))
        })
    })
}

// Adds a circular Gaussian of height `peak` and width `sigma` centred on `centre`.
pub(crate) fn add_gaussian(
    image: &mut [f32],
    w: usize,
    h: usize,
    centre: (f64, f64),
    peak: f32,
    sigma: f64,
) {
    for (v, g) in image.iter_mut().zip(gaussian_profile(w, h, centre, sigma)) {
        *v += peak * g as f32;
    }
}

// A `w` x `h` image of constant `bg` with one Gaussian star added.
pub(crate) fn gaussian_image(
    w: usize,
    h: usize,
    centre: (f64, f64),
    peak: f32,
    sigma: f64,
    bg: f32,
) -> Vec<f32> {
    let mut image = vec![bg; w * h];
    add_gaussian(&mut image, w, h, centre, peak, sigma);
    image
}

// A 128 x 128 u8 frame at 10 ADU with about 2% of pixels at 11, the quantised sky of a
// low-gain camera. Each star is a 3 x 3 pyramid: `peak` at the centre, halfway between
// 10 and `peak` on the ring.
pub(crate) fn quantised_frame(stars: &[(usize, usize, u8)]) -> Vec<u8> {
    let mut state = 3_u64;
    let mut image: Vec<u8> = (0..128 * 128)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            if (state >> 33).is_multiple_of(50) {
                11
            } else {
                10
            }
        })
        .collect();
    for &(cx, cy, peak) in stars {
        for y in cy - 1..=cy + 1 {
            for x in cx - 1..=cx + 1 {
                image[y * 128 + x] = 10 + (peak - 10) / 2;
            }
        }
        image[cy * 128 + cx] = peak;
    }
    image
}

// Sample standard deviation of `image - 10` from integer sums. Every intermediate the
// float code forms is exact for these frames, so the result matches it bit for bit.
pub(crate) fn sample_sigma_about_10(image: &[u8]) -> f64 {
    let n = image.len() as i64;
    let (s, q) = image.iter().fold((0_i64, 0_i64), |(s, q), &v| {
        let r = v as i64 - 10;
        (s + r, q + r * r)
    });
    let var_sum = (n * q - s * s) as f64 / n as f64;
    libm::sqrt(var_sum / (n - 1) as f64)
}
