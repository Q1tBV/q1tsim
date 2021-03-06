use criterion::{criterion_group, Criterion};
use q1tsim::circuit::{QuStateRepr, Circuit};

use rand_core::SeedableRng;

fn build_manybits_circuit(nr_bits: usize, nr_gates: usize) -> q1tsim::error::Result<Circuit>
{
    let mut circuit = Circuit::new(nr_bits, 2);

    for i in 0..nr_gates
    {
        circuit.h(i % nr_bits)?;
        circuit.cx(i % nr_bits, (i+1) % nr_bits)?;
    }

    circuit.measure(0, 0)?;
    circuit.measure(nr_bits-1, 1)?;

    Ok(circuit)
}

fn manybits(nr_bits: usize, nr_gates: usize, nr_shots: usize, stabilizer: bool)
{
    let mut rng = rand_hc::Hc128Rng::seed_from_u64(0x1f67a51423cd2615);

    let mut circuit = build_manybits_circuit(nr_bits, nr_gates).expect("Failed to build circuit");
    if stabilizer
    {
        circuit.execute_with_rng(nr_shots, &mut rng).expect("Failed to execute circuit");
    }
    else
    {
        let q_state = QuStateRepr::vector(nr_bits, nr_shots);
        circuit.execute_with(nr_shots, &mut rng, q_state).expect("Failed to execute circuit");
    }
}

fn bench_manybits(c: &mut Criterion)
{
    c.bench_function("mb vec 10", |b| b.iter(|| manybits(10, 100, 1_000, false)));
    c.bench_function("mb vec 11", |b| b.iter(|| manybits(11, 100, 1_000, false)));
    c.bench_function("mb vec 12", |b| b.iter(|| manybits(12, 100, 1_000, false)));

    c.bench_function("mb stab 10", |b| b.iter(|| manybits(10, 100, 1_000, true)));
    c.bench_function("mb stab 11", |b| b.iter(|| manybits(11, 100, 1_000, true)));
    c.bench_function("mb stab 12", |b| b.iter(|| manybits(12, 100, 1_000, true)));
    c.bench_function("mb stab 20", |b| b.iter(|| manybits(20, 100, 1_000, true)));
    c.bench_function("mb stab 35", |b| b.iter(|| manybits(35, 100, 1_000, true)));
    c.bench_function("mb stab 50", |b| b.iter(|| manybits(50, 100, 1_000, true)));
    c.bench_function("mb stab 100", |b| b.iter(|| manybits(100, 100, 1_000, true)));
}

criterion_group!(benches, bench_manybits);
