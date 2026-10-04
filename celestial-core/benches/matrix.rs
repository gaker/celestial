use std::hint::black_box;

use celestial_core::constants::{PI, TWOPI};
use celestial_core::errors::AstroResult;
use celestial_core::matrix::{RotationMatrix3, Vector3};
use criterion::{criterion_group, criterion_main, Criterion};

const BATCH: usize = 1024;

fn rotation(z1: f64, x: f64, z2: f64) -> RotationMatrix3 {
    let mut m = RotationMatrix3::identity();
    m.rotate_z(z1);
    m.rotate_x(x);
    m.rotate_z(z2);
    m
}

// Golden-angle spiral: evenly spread unit vectors without a random-number crate.
fn unit_vectors(n: usize) -> Vec<Vector3> {
    let golden = PI * (3.0 - libm::sqrt(5.0));
    (0..n)
        .map(|i| {
            let ra = libm::fmod(i as f64 * golden, TWOPI);
            let dec = libm::asin(2.0 * (i as f64 + 0.5) / n as f64 - 1.0);
            Vector3::from_spherical(ra, dec)
        })
        .collect()
}

fn bench_matrix_ops(c: &mut Criterion) {
    let a = rotation(0.3, -0.4, 1.1);
    let b = rotation(-2.0, 0.7, 0.05);
    let v = [0.6, -0.48, 0.64];
    let mut group = c.benchmark_group("matrix");
    group.bench_function("multiply", |bn| {
        bn.iter(|| black_box(&a).multiply(black_box(&b)))
    });
    group.bench_function("apply_to_vector", |bn| {
        bn.iter(|| black_box(&a).apply_to_vector(black_box(v)))
    });
    group.bench_function("transpose", |bn| bn.iter(|| black_box(&a).transpose()));
    group.finish();
}

fn bench_matrix_rotations(c: &mut Criterion) {
    let a = rotation(0.3, -0.4, 1.1);
    let mut group = c.benchmark_group("matrix");
    group.bench_function("rotate_x", |bn| {
        bn.iter(|| {
            let mut m = black_box(a);
            m.rotate_x(black_box(0.25));
            m
        })
    });
    group.bench_function("rotate_z", |bn| {
        bn.iter(|| {
            let mut m = black_box(a);
            m.rotate_z(black_box(0.25));
            m
        })
    });
    group.finish();
}

fn bench_vector_ops(c: &mut Criterion) {
    let u = Vector3::new(0.6, -0.48, 0.64);
    let w = Vector3::new(-0.3, 0.9, 0.1);
    let mut group = c.benchmark_group("vector");
    group.bench_function("add", |bn| bn.iter(|| black_box(u) + black_box(w)));
    group.bench_function("dot", |bn| bn.iter(|| black_box(u).dot(black_box(&w))));
    group.bench_function("cross", |bn| bn.iter(|| black_box(u).cross(black_box(&w))));
    group.bench_function("normalize", |bn| bn.iter(|| black_box(w).normalize()));
    group.finish();
}

// Per-call timings of nanosecond methods are mostly call overhead; a loop over many
// vectors shows whether inlining lets the compiler hoist and vectorise.
fn bench_matrix_batch(c: &mut Criterion) {
    let m = rotation(0.3, -0.4, 1.1);
    let vectors = unit_vectors(BATCH);
    let mut out = vec![Vector3::zeros(); BATCH];
    c.bench_function("batch_1024/matrix_times_vector", |bn| {
        bn.iter(|| {
            let m = black_box(&m);
            for (o, v) in out.iter_mut().zip(black_box(&vectors)) {
                *o = m * *v;
            }
        })
    });
}

fn bench_vector_batch(c: &mut Criterion) {
    let vectors = unit_vectors(BATCH);
    let mut out: Vec<AstroResult<Vector3>> = (1..BATCH).map(|_| Ok(Vector3::zeros())).collect();
    c.bench_function("batch_1024/cross_then_normalize", |bn| {
        bn.iter(|| {
            for (o, pair) in out.iter_mut().zip(black_box(&vectors).windows(2)) {
                *o = pair[0].cross(&pair[1]).normalize();
            }
        })
    });
    assert!(
        out.iter().all(|r| r.is_ok()),
        "spiral neighbours are never parallel"
    );
}

criterion_group!(
    benches,
    bench_matrix_ops,
    bench_matrix_rotations,
    bench_vector_ops,
    bench_matrix_batch,
    bench_vector_batch
);
criterion_main!(benches);
