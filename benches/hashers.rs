use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fasthash::CityHasher;
use fnv::FnvHashMap;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;

fn criterion_benchmark(c: &mut Criterion) {
    let input_motif = black_box("ABCDEF");

    let mut fxMap = FxHashMap::default();
    let mut fnvMap = FnvHashMap::default();
    let mut cityMap = HashMap::<&str, i32, BuildHasherDefault<CityHasher>>::default();

    c.bench_function("inserting into fxmap", |b| {
        b.iter(|| fxMap.insert(input_motif, 0))
    });

    c.bench_function("inserting into fnvmap", |b| {
        b.iter(|| fnvMap.insert(input_motif, 0))
    });

    c.bench_function("inserting into citymap", |b| {
        b.iter(|| cityMap.insert(input_motif, 0))
    });

    c.bench_function("getting out of fxmap", |b| {
        b.iter(|| fxMap.get(input_motif))
    });

    c.bench_function("getting out of fnvmap", |b| {
        b.iter(|| fnvMap.get(input_motif))
    });

    c.bench_function("getting out of citymap", |b| {
        b.iter(|| cityMap.get(input_motif))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
