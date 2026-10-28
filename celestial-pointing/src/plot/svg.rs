use rizzma::artist::Rgba;
use rizzma::core::AxisBelow;
use rizzma::figure::{Bins, HistOptions, QuiverAngles, QuiverOptions, QuiverUnits};
use rizzma::{Axes, Figure};
use std::path::Path;

type PlotResult = std::result::Result<(), Box<dyn std::error::Error>>;

const TITLE_SIZE_PX: f64 = 19.0;
const MARKER_DIAMETER_PT: f64 = 6.0 * 72.0 / 100.0;

pub fn scatter_svg(
    points: &[(f64, f64)],
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
) -> PlotResult {
    let (x_range, y_range) = padded_ranges(points);
    let mut fig = new_figure();
    let ax = fig.add_subplot(1, 1, 1);
    setup_axes(ax, title, x_label, y_label, &x_range, &y_range);
    draw_crosshairs(ax);
    draw_points(ax, points);
    fig.save_svg(path)?;
    Ok(())
}

pub fn histogram_svg(values: &[f64], path: &Path, title: &str, x_label: &str) -> PlotResult {
    if values.is_empty() {
        return Ok(());
    }
    let mut fig = new_figure();
    let ax = fig.add_subplot(1, 1, 1);
    let hist = ax.hist_with(
        values,
        &HistOptions::new().bins(Bins::Count(20)).color(Rgba::BLUE),
    );
    let (Some(&x_min), Some(&x_max)) = (hist.edges.first(), hist.edges.last()) else {
        return Ok(());
    };
    let max_count = hist.counts.iter().copied().fold(1.0, f64::max);
    setup_axes(
        ax,
        title,
        x_label,
        "Count",
        &(x_min, x_max),
        &(0.0, max_count + 1.0),
    );
    // Keep the grid behind the bars.
    ax.set_axisbelow(AxisBelow::Below);
    fig.save_svg(path)?;
    Ok(())
}

pub fn vector_map_svg(
    positions: &[(f64, f64)],
    vectors: &[(f64, f64)],
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    scale: f64,
) -> PlotResult {
    let (x_range, y_range) = padded_ranges(positions);
    let mut fig = new_figure();
    let ax = fig.add_subplot(1, 1, 1);
    setup_axes(ax, title, x_label, y_label, &x_range, &y_range);
    draw_points(ax, positions);
    draw_vectors(ax, positions, vectors, scale);
    fig.save_svg(path)?;
    Ok(())
}

fn new_figure() -> Figure {
    Figure::new(8.0, 6.0).with_dpi(100.0)
}

fn padded_ranges(points: &[(f64, f64)]) -> ((f64, f64), (f64, f64)) {
    if points.is_empty() {
        return ((-1.0, 1.0), (-1.0, 1.0));
    }
    let (mut x_min, mut x_max) = extent(points.iter().map(|p| p.0));
    let (mut y_min, mut y_max) = extent(points.iter().map(|p| p.1));
    let x_pad = (x_max - x_min).abs() * 0.1 + 1e-6;
    let y_pad = (y_max - y_min).abs() * 0.1 + 1e-6;
    x_min -= x_pad;
    x_max += x_pad;
    y_min -= y_pad;
    y_max += y_pad;
    ((x_min, x_max), (y_min, y_max))
}

fn extent(iter: impl Iterator<Item = f64>) -> (f64, f64) {
    iter.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
        (lo.min(v), hi.max(v))
    })
}

fn setup_axes(
    ax: &mut Axes,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_range: &(f64, f64),
    y_range: &(f64, f64),
) {
    ax.set_title(title)
        .set_title_size(TITLE_SIZE_PX)
        .set_xlabel(x_label)
        .set_ylabel(y_label)
        .set_xlim(x_range.0, x_range.1)
        .set_ylim(y_range.0, y_range.1)
        .grid(true);
}

fn draw_crosshairs(ax: &mut Axes) {
    let color = Rgba::BLACK.with_alpha(0.3);
    ax.axhline_with(0.0, color, 1.0);
    ax.axvline_with(0.0, color, 1.0);
}

fn draw_points(ax: &mut Axes, points: &[(f64, f64)]) {
    let xs: Vec<f64> = points.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = points.iter().map(|p| p.1).collect();
    ax.scatter(&xs, &ys)
        .set_facecolors(vec![Rgba::BLUE])
        .set_sizes(vec![MARKER_DIAMETER_PT]);
}

/// Each residual becomes an arrow from its position to `position + scale * vector`, in data units.
fn draw_vectors(ax: &mut Axes, positions: &[(f64, f64)], vectors: &[(f64, f64)], scale: f64) {
    if scale <= 0.0 || !scale.is_finite() {
        return;
    }
    let (xs, ys): (Vec<f64>, Vec<f64>) = positions.iter().copied().unzip();
    let (us, vs): (Vec<f64>, Vec<f64>) = vectors.iter().copied().unzip();
    ax.quiver_with(
        &xs,
        &ys,
        &us,
        &vs,
        &QuiverOptions::new()
            .scale_units(QuiverUnits::Xy)
            .angles(QuiverAngles::Xy)
            .scale(1.0 / scale)
            .color(Rgba::RED),
    );
}
