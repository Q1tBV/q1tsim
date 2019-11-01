use criterion::{criterion_group, Criterion};
use q1tsim::{declare_controlled, declare_controlled_cost, declare_controlled_type,
    declare_controlled_impl, declare_controlled_impl_gate, declare_controlled_latex,
    declare_controlled_qasm, declare_controlled_square};
use q1tsim::circuit::Circuit;
use q1tsim::gates::{CCX};

use rand_core::SeedableRng;

declare_controlled!(C3X, CCX);
declare_controlled!(C4X, C3X);
declare_controlled!(C5X, C4X);
declare_controlled!(C6X, C5X);
declare_controlled!(C7X, C6X);
declare_controlled!(C8X, C7X);
declare_controlled!(C9X, C8X);

fn add_cnx(circuit: &mut Circuit, i: usize, nr_pos_bits: usize)
    -> q1tsim::error::Result<()>
{
    let n = nr_pos_bits+1-i;
    let bits: Vec<usize> = (0..=n).collect();
    match n
    {
        9 => circuit.add_gate(C9X::new(), &bits),
        8 => circuit.add_gate(C8X::new(), &bits),
        7 => circuit.add_gate(C7X::new(), &bits),
        6 => circuit.add_gate(C6X::new(), &bits),
        5 => circuit.add_gate(C5X::new(), &bits),
        4 => circuit.add_gate(C4X::new(), &bits),
        3 => circuit.add_gate(C3X::new(), &bits),
        2 => circuit.add_gate(CCX::new(), &bits),
        _ => panic!("Unable to add C{}X gate", n)
    }
}

fn build_randomwalk_circuit(nr_pos_bits: usize, measure: bool) -> q1tsim::error::Result<Circuit>
{
    let mut circuit = Circuit::new(nr_pos_bits+1, nr_pos_bits);

    for _ in 0..((1 << nr_pos_bits) - 1)
    {
        circuit.h(0)?;

        for i in 1..nr_pos_bits
        {
            add_cnx(&mut circuit, i, nr_pos_bits)?;
        }
        circuit.cx(0, 1)?;

        if measure
        {
            circuit.measure(0, 0)?;
        }
    }

    for i in 0..nr_pos_bits
    {
        circuit.measure(i+1, i)?;
    }

    Ok(circuit)
}

fn quantum_random_walk(nr_pos_bits: usize, nr_shots: usize, measure: bool)
{
    let mut rng = rand_hc::Hc128Rng::seed_from_u64(0x1f67a51423cd2615);

    let mut circuit = build_randomwalk_circuit(nr_pos_bits, measure).expect("Failed to build circuit");
    circuit.execute_with_rng(nr_shots, &mut rng).expect("Failed to execute circuit");
}

fn bench_randomwalk(c: &mut Criterion)
{
    c.bench_function("qrw 6 1,000", |b| b.iter(|| quantum_random_walk(6, 1_000, false)));
    c.bench_function("qrw 6 1,000,000", |b| b.iter(|| quantum_random_walk(6, 1_000_000, false)));
    c.bench_function("crw 4 10", |b| b.iter(|| quantum_random_walk(4, 10, true)));
    c.bench_function("crw 4 100", |b| b.iter(|| quantum_random_walk(4, 100, true)));
}

criterion_group!(benches, bench_randomwalk);
