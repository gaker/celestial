mod common;

use std::hint::black_box;

use celestial_core::constants::J2000_JD;
use celestial_core::errors::AstroResult;
use celestial_core::nutation::{
    NutationIAU2000A, NutationIAU2000B, NutationIAU2006A, NutationResult,
};
use common::EPOCHS;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_model<F>(c: &mut Criterion, group_name: &str, compute: F)
where
    F: Fn(f64, f64) -> AstroResult<NutationResult>,
{
    let mut group = c.benchmark_group(group_name);
    for (name, jd2) in EPOCHS {
        compute(J2000_JD, jd2).expect("epoch inside the model range");
        group.bench_function(name, |b| {
            b.iter(|| compute(black_box(J2000_JD), black_box(jd2)))
        });
    }
    group.finish();
}

fn nutation(c: &mut Criterion) {
    let iau2000a = NutationIAU2000A::new();
    bench_model(c, "nutation_iau2000a", |jd1, jd2| {
        iau2000a.compute(jd1, jd2)
    });
    let iau2000b = NutationIAU2000B::new();
    bench_model(c, "nutation_iau2000b", |jd1, jd2| {
        iau2000b.compute(jd1, jd2)
    });
    let iau2006a = NutationIAU2006A::new();
    bench_model(c, "nutation_iau2006a", |jd1, jd2| {
        iau2006a.compute(jd1, jd2)
    });
}

criterion_group!(benches, nutation);
criterion_main!(benches);
