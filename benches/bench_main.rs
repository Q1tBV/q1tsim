use criterion::criterion_main;

mod benchmarks;

criterion_main!(
    benchmarks::randomwalk::benches,
    benchmarks::manybits::benches
);
