mod common;

use std::hint::black_box;

use celestial_core::constants::J2000_JD;
use celestial_core::matrix::RotationMatrix3;
use celestial_core::nutation::NutationIAU2006A;
use celestial_core::precession::PrecessionIAU2006;
use celestial_core::utils::jd_to_centuries;
use celestial_core::{
    cio::{gcrs_to_cirs_matrix, CioLocator, CioSolution},
    errors::AstroResult,
};
use common::EPOCHS;
use criterion::{criterion_group, criterion_main, Criterion};

fn npb_at(jd2: f64) -> AstroResult<(f64, RotationMatrix3)> {
    let t = jd_to_centuries(J2000_JD, jd2);
    let nutation = NutationIAU2006A::new().compute(J2000_JD, jd2)?;
    let npb =
        PrecessionIAU2006::new().npb_matrix_iau2006a(t, nutation.delta_psi, nutation.delta_eps);
    Ok((t, npb))
}

// The same chain celestial-coords runs for GCRS to CIRS, so the CIO step can be
// read as a share of the whole.
fn gcrs_to_cirs(jd2: f64) -> AstroResult<RotationMatrix3> {
    let (t, npb) = npb_at(jd2)?;
    let cio = CioSolution::calculate(&npb, t)?;
    gcrs_to_cirs_matrix(cio.cip.x, cio.cip.y, cio.s)
}

fn bench_solution(c: &mut Criterion) {
    let mut group = c.benchmark_group("cio_solution_calculate");
    for (name, jd2) in EPOCHS {
        let (t, npb) = npb_at(jd2).expect("epoch inside the model range");
        CioSolution::calculate(&npb, t).expect("CIP inside the model range");
        group.bench_function(name, |b| {
            b.iter(|| CioSolution::calculate(black_box(&npb), black_box(t)))
        });
    }
    group.finish();
}

fn bench_locator(c: &mut Criterion) {
    let mut group = c.benchmark_group("cio_locator_iau2006a");
    for (name, jd2) in EPOCHS {
        let (t, npb) = npb_at(jd2).expect("epoch inside the model range");
        let (x, y) = (npb.elements()[2][0], npb.elements()[2][1]);
        group.bench_function(name, |b| {
            b.iter(|| CioLocator::iau2006a(black_box(t)).calculate(black_box(x), black_box(y)))
        });
    }
    group.finish();
}

fn bench_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("cio_gcrs_to_cirs_pipeline");
    for (name, jd2) in EPOCHS {
        let _ = gcrs_to_cirs(jd2).expect("epoch inside the model range");
        group.bench_function(name, |b| b.iter(|| gcrs_to_cirs(black_box(jd2))));
    }
    group.finish();
}

criterion_group!(benches, bench_solution, bench_locator, bench_pipeline);
criterion_main!(benches);
