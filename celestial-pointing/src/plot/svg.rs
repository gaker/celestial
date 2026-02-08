use plotters::prelude::*;
use std::path::Path;

type PlotResult = std::result::Result<(), Box<dyn std::error::Error>>;

pub fn scatter_svg(
    points: &[(f64, f64)],
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
) -> PlotResult {
    let (x_range, y_range) = padded_ranges(points);
    let root = SVGBackend::new(path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = build_chart(&root, title, x_label, y_label, &x_range, &y_range)?;
    draw_crosshairs(&mut chart, &x_range, &y_range)?;
    draw_points(&mut chart, points)?;
    root.present()?;
    Ok(())
}

pub fn histogram_svg(
    values: &[f64],
    path: &Path,
    title: &str,
    x_label: &str,
) -> PlotResult {
    if values.is_empty() {
        return Ok(());
    }
    let (bins, bin_width, min_val) = bin_values(values, 20);
    let max_count = bins.iter().copied().max().unwrap_or(1);
    let root = SVGBackend::new(path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let x_max = min_val + 20.0 * bin_width;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_val..x_max, 0u32..(max_count + 1))?;
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc("Count")
        .draw()?;
    draw_histogram_bars(&mut chart, &bins, min_val, bin_width)?;
    root.present()?;
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
    let root = SVGBackend::new(path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = build_chart(&root, title, x_label, y_label, &x_range, &y_range)?;
    draw_vectors(&mut chart, positions, vectors, scale)?;
    root.present()?;
    Ok(())
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

fn build_chart<'a, DB: DrawingBackend + 'a>(
    area: &'a DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_range: &(f64, f64),
    y_range: &(f64, f64),
) -> std::result::Result<
    ChartContext<'a, DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>,
    Box<dyn std::error::Error>,
>
where
    DB::ErrorType: 'static,
{
    let mut chart = ChartBuilder::on(area)
        .caption(title, ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_range.0..x_range.1, y_range.0..y_range.1)?;
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()?;
    Ok(chart)
}

fn draw_crosshairs<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>,
    x_range: &(f64, f64),
    y_range: &(f64, f64),
) -> PlotResult
where
    DB::ErrorType: 'static,
{
    let style = BLACK.mix(0.3).stroke_width(1);
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(x_range.0, 0.0), (x_range.1, 0.0)],
        style,
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, y_range.0), (0.0, y_range.1)],
        style,
    )))?;
    Ok(())
}

fn draw_points<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>,
    points: &[(f64, f64)],
) -> PlotResult
where
    DB::ErrorType: 'static,
{
    chart.draw_series(
        points
            .iter()
            .map(|&(x, y)| Circle::new((x, y), 3, BLUE.filled())),
    )?;
    Ok(())
}

fn draw_histogram_bars<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordu32>>,
    bins: &[u32],
    min_val: f64,
    bin_width: f64,
) -> PlotResult
where
    DB::ErrorType: 'static,
{
    chart.draw_series(bins.iter().enumerate().map(|(i, &count)| {
        let x0 = min_val + i as f64 * bin_width;
        let x1 = x0 + bin_width;
        Rectangle::new([(x0, 0), (x1, count)], BLUE.filled())
    }))?;
    Ok(())
}

fn draw_vectors<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>,
    positions: &[(f64, f64)],
    vectors: &[(f64, f64)],
    scale: f64,
) -> PlotResult
where
    DB::ErrorType: 'static,
{
    for (&(px, py), &(vx, vy)) in positions.iter().zip(vectors.iter()) {
        chart.draw_series(std::iter::once(Circle::new(
            (px, py),
            3,
            BLUE.filled(),
        )))?;
        let ex = px + vx * scale;
        let ey = py + vy * scale;
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(px, py), (ex, ey)],
            RED.stroke_width(1),
        )))?;
    }
    Ok(())
}

fn bin_values(values: &[f64], n_bins: usize) -> (Vec<u32>, f64, f64) {
    let (min_val, max_val) = extent(values.iter().copied());
    let range = (max_val - min_val).max(1e-10);
    let bin_width = range / n_bins as f64;
    let mut bins = vec![0u32; n_bins];
    for &v in values {
        let idx = ((v - min_val) / bin_width).floor() as usize;
        bins[idx.min(n_bins - 1)] += 1;
    }
    (bins, bin_width, min_val)
}
