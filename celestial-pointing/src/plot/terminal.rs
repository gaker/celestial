use textplots::{Chart, Plot, Shape};

pub fn scatter_terminal(
    points: &[(f64, f64)],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> String {
    if points.is_empty() {
        return format!("{title}\n  (no data)\n");
    }
    let f32_pts = to_f32_points(points);
    let (xmin, xmax) = f32_extent(f32_pts.iter().map(|p| p.0));
    let chart_body = render_chart(&f32_pts, xmin, xmax);
    format!("{title}\n  {y_label} vs {x_label}\n{chart_body}")
}

pub fn histogram_terminal(values: &[f64], title: &str, label: &str) -> String {
    if values.is_empty() {
        return format!("{title}\n  (no data)\n");
    }
    let (bins, bin_width, min_val) = bin_values(values, 20);
    let max_count = bins.iter().copied().max().unwrap_or(1);
    let (mean, sigma) = mean_sigma(values);
    let mut out = format!("{title}\n\n");
    for (i, &count) in bins.iter().enumerate() {
        let edge = min_val + i as f64 * bin_width;
        let bar_len = if max_count > 0 {
            (count * 40) / max_count
        } else {
            0
        };
        let bar: String = "\u{2588}".repeat(bar_len);
        out.push_str(&format!("  {edge:>8.1} \u{2502}{bar}\n"));
    }
    out.push_str(&format!(
        "\n  {label}  mean: {mean:.1}\"  sigma: {sigma:.1}\"\n"
    ));
    out
}

pub fn xy_plot_terminal(
    points: &[(f64, f64)],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> String {
    if points.is_empty() {
        return format!("{title}\n  (no data)\n");
    }
    let f32_pts = to_f32_points(points);
    let (xmin, xmax) = f32_extent(f32_pts.iter().map(|p| p.0));
    let chart_body = render_chart(&f32_pts, xmin, xmax);
    format!("{title}\n  {y_label} vs {x_label}\n{chart_body}")
}

fn render_chart(pts: &[(f32, f32)], xmin: f32, xmax: f32) -> String {
    let shape = Shape::Points(pts);
    let mut chart = Chart::new(80, 24, xmin, xmax);
    let rendered = chart.lineplot(&shape);
    rendered.axis();
    rendered.figures();
    format!("{rendered}")
}

fn to_f32_points(points: &[(f64, f64)]) -> Vec<(f32, f32)> {
    points.iter().map(|&(x, y)| (x as f32, y as f32)).collect()
}

fn f32_extent(iter: impl Iterator<Item = f32>) -> (f32, f32) {
    let (lo, hi) = iter.fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), v| {
        (lo.min(v), hi.max(v))
    });
    if (hi - lo).abs() < 1e-6 {
        (lo - 1.0, hi + 1.0)
    } else {
        (lo, hi)
    }
}

fn bin_values(values: &[f64], n_bins: usize) -> (Vec<usize>, f64, f64) {
    let (min_val, max_val) = f64_extent(values.iter().copied());
    let range = (max_val - min_val).max(1e-10);
    let bin_width = range / n_bins as f64;
    let mut bins = vec![0usize; n_bins];
    for &v in values {
        let idx = libm::floor((v - min_val) / bin_width) as usize;
        bins[idx.min(n_bins - 1)] += 1;
    }
    (bins, bin_width, min_val)
}

fn f64_extent(iter: impl Iterator<Item = f64>) -> (f64, f64) {
    iter.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
        (lo.min(v), hi.max(v))
    })
}

fn mean_sigma(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n;
    (mean, libm::sqrt(variance))
}
