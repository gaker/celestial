mod common;

use std::hint::black_box;

use celestial_core::constants::J2000_JD;
use celestial_core::nutation::NutationIAU2006A;
use celestial_core::precession::PrecessionIAU2006;
use celestial_core::utils::jd_to_centuries;
use common::EPOCHS;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_compute(c: &mut Criterion) {
    let model = PrecessionIAU2006::new();
    let mut group = c.benchmark_group("precession_iau2006_compute");
    for (name, jd2) in EPOCHS {
        model
            .compute(J2000_JD, jd2)
            .expect("epoch inside the model range");
        group.bench_function(name, |b| {
            b.iter(|| model.compute(black_box(J2000_JD), black_box(jd2)))
        });
    }
    group.finish();
}

fn bench_npb_matrix(c: &mut Criterion) {
    let model = PrecessionIAU2006::new();
    let mut group = c.benchmark_group("precession_npb_matrix_iau2006a");
    for (name, jd2) in EPOCHS {
        let t = jd_to_centuries(J2000_JD, jd2);
        let nutation = NutationIAU2006A::new()
            .compute(J2000_JD, jd2)
            .expect("epoch inside the model range");
        let (dpsi, deps) = (nutation.delta_psi, nutation.delta_eps);
        group.bench_function(name, |b| {
            b.iter(|| model.npb_matrix_iau2006a(black_box(t), black_box(dpsi), black_box(deps)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_compute, bench_npb_matrix);
criterion_main!(benches);
