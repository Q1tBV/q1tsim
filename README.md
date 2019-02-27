# [q1tsim](https://github.com/Q1tBV/q1tsim/)

[![Build Status](https://travis-ci.org/Q1tBV/q1tsim.svg?branch=master)](https://travis-ci.org/Q1tBV/q1tsim)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![crates.io](http://meritbadge.herokuapp.com/q1tsim)](https://crates.io/crates/q1tsim)
[![Released API docs](https://docs.rs/q1tsim/badge.svg)](https://docs.rs/q1tsim)

A simple, efficient, quantum computer simulator.

Overview
========

q1tsim is a simulator library for a quantum computer, written in Rust. Its goal
is to be an easy to use, efficient simulator for the development and testing of
quantum algorithms.

Features
========
* Easy implementation and simulation of quantum circuits
* Supports the creation of arbitrary quantum gates
* Most common quantum gates already included
* Measurement in `X`, `Y`, or `Z` basis.
* Creation of histograms of measurement results over multiple runs
* Operations conditional on classical values
* Export of circuits to Open QASM and c-QASM

Usage
=====
To use q1tsim in your Rust application, add the following to your `Cargo.toml` file:

```toml
[dependencies]
q1tsim = "0.2"
```

As an example, here is a 3-qubit quantum Fourier transform of the |000〉quantum
state:
```
extern crate q1tsim;

use q1tsim::{circuit, gates};

fn main()
{
    // The number of times this circuit os evaulated
    let nr_runs = 8192;

    // Create a quantum circuit with 3 quantum bits and 3 classical (measurement)
    // bits, that is evaluated `nr_runs` times. The circuit starts by default
    // with all quantum bits in the |0〉state, so in this case |000〉.
    let mut circuit = circuit::Circuit::new(3, 3, nr_runs);

    // Set up a 3-qubit quantum Fourier transform
    // There is no predefined method on Circuit that implements a controlled
    // `S` or `T` gate, so we use the `add_gate()` method for those.
    circuit.h(2);
    circuit.add_gate(gates::CS::new(), &[1, 2]);
    circuit.add_gate(gates::CT::new(), &[0, 2]);
    circuit.h(1);
    circuit.add_gate(gates::CS::new(), &[0, 1]);
    circuit.h(0);
    circuit.add_gate(gates::Swap::new(), &[0, 2]);

    // Measure all quantum bits in the Pauli `Z` basis
    circuit.measure_all(&[0, 1, 2]);

    // Actually calculate the resulting quantum state and perform the measurements
    circuit.execute();

    // And print the results.
    let hist = circuit.histogram_string();
    for (bits, count) in hist
    {
        println!("{}: {}", bits, count);
    }
}
```
The result should be a more or less equal distribution over the eight possible
states (000, 001, ..., 111).

Read the complete source code documentation on [docs.rs](https://docs.rs/q1tsim).
