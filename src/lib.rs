// Copyright 2019 Q1t BV
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! A simple, efficient, quantum computer simulator.
//!
//! Overview
//! ========
//!
//! q1tsim is a simulator library for a quantum computer, written in Rust. Its goal
//! is to be an easy to use, efficient simulator for the development and testing of
//! quantum algorithms.
//!
//! Features
//! ========
//! * Easy implementation and simulation of quantum circuits
//! * Supports the creation of arbitrary quantum gates
//! * Most common quantum gates already included
//! * Measurement in `X`, `Y`, or `Z` basis
//! * Possibility of measurement without affecting the quantum state
//! * Creation of histograms of measurement results over multiple runs
//! * Operations conditional on classical values
//! * Export of circuits to Open QASM and c-QASM for running your programs on other computers or simulators
//! * Export of circuits to LaTeX, for drawing pictures of your circuit
//! * Efficient simulation of stabilizer circuits
//!
//! Usage
//! =====
//! To use q1tsim in your Rust application, add the following to your `Cargo.toml` file:
//!
//! ```toml
//! [dependencies]
//! q1tsim = "0.3"
//! ```
//!
//! As an example, here is a 3-qubit quantum Fourier transform of the |000⟩ quantum
//! state:
//! ```
//! use q1tsim::circuit::Circuit;
//! use q1tsim::gates::{CS, CT, Swap};
//!
//! fn main()
//! {
//!     // The number of times this circuit is evaluated
//!     let nr_runs = 8192;
//!
//!     // Create a quantum circuit with 3 quantum bits and 3 classical (measurement)
//!     // bits. The circuit starts by default with all quantum bits in the |0⟩
//!     // state, so in this case |000⟩.
//!     let mut circuit = Circuit::new(3, 3);
//!
//!     // Set up a 3-qubit quantum Fourier transform
//!     // There is no predefined method on Circuit that implements a controlled
//!     // `S` or `T` gate, so we use the `add_gate()` method for those.
//!     circuit.h(2);
//!     circuit.add_gate(CS::new(), &[1, 2]);
//!     circuit.add_gate(CT::new(), &[0, 2]);
//!     circuit.h(1);
//!     circuit.add_gate(CS::new(), &[0, 1]);
//!     circuit.h(0);
//!     circuit.add_gate(Swap::new(), &[0, 2]);
//!
//!     // Measure all quantum bits in the Pauli `Z` basis
//!     circuit.measure_all(&[0, 1, 2]);
//!
//!     // Actually calculate the resulting quantum state, and perform the
//!     // measurements, averaging over `nr_runs` runs.
//!     circuit.execute(nr_runs);
//!
//!     // And print the results.
//!     let hist = circuit.histogram_string().unwrap();
//!     for (bits, count) in hist.iter()
//!     {
//!         println!("{}: {}", bits, count);
//!     }
//! }
//! ```
//! The result should be a more or less equal distribution over the eight possible
//! states (000, 001, ..., 111).
//!
//! Creating a circuit
//! ==================
//! Struct [Circuit](circuit/struct.Circuit.html) is the main structure used
//! in creating a quantum program. The basic layout of a program to create
//! a quantum circuit, execute it, and collect the results, is as follows:
//! ```
//! use q1tsim::circuit::Circuit;
//!
//! // Create a new circuit with `nr_qbits` quantum bits and `nr_cbits`
//! // classical bits
//! let nr_qbits = 2;
//! let nr_cbits = 2;
//! let mut circuit = Circuit::new(nr_qbits, nr_cbits);
//!
//! // Add operations on the circuit. In this case, a Hadamard transform on the
//! // first bit, followwed by a CNOT gate with the first bit as control and
//! // the second bit as target.
//! circuit.h(0);
//! circuit.cx(0, 1);
//!
//! // Add a measurement of the resulting quantum state. This measures the first
//! // qbit into classical bit 0, and the second qbit into classical bit 1.
//! circuit.measure_all(&[0, 1]);
//!
//! // Now execute the circuit, averaging measurements over `nr_runs` runs
//! // of the circuit.
//! let nr_runs = 1024;
//! circuit.execute(nr_runs);
//!
//! // And finally collect the results. The `histogram_vec()` method returns a
//! // vector with at each index `i` the number if times the measurement returned
//! // `i` in the classical register.
//! let hist = circuit.histogram_vec();
//! ```
//! Since version 0.3, many of the methods on `Circuit` will return a `Result`,
//! possibly containing an error code (e.g. if invalid bit numbers are used).
//! Checking the result of each modification of the circuit quickly becomes
//! tedious, so the [circuit](macro.circuit.html) macro was added that can make
//! multiple method calls and immediately returns on the first error encountered
//! (or returns `Ok(())` if all calls were successful). With this, the previous
//! program can be written as
//! ```rust
//! # #[macro_use] extern crate q1tsim; fn main() {
//! let nr_qbits = 2;
//! let nr_cbits = 2;
//! let mut circuit = circuit!(nr_qbits, nr_cbits, {
//!     h(0);
//!     cx(0, 1);
//!     measure_all(&[0, 1]);
//! }).expect("Failed to build circuit");
//!
//! let nr_runs = 1024;
//! circuit.execute(nr_runs);
//! let hist = circuit.histogram_vec();
//! # }
//!```
//!
//! Custom gates
//! ============
//! Using the [Circuit::add_gate()](circuit/struct.Circuit.html#method.add_gate)
//! method, arbitrary gates can be added to a circuit. You can define your own
//! custom gates by implementing the [Gate](gates/trait.Gate.html) trait. To
//! implement this trait, the type should implement at least the
//! [description()](gates/trait.Gate.html#tymethod.description),
//! [nr_affected_bits()](gates/trait.Gate.html#tymethod.nr_affected_bits),
//! and [matrix()](gates/trait.Gate.html#tymethod.matrix) methods. The
//! `description()` method should return a short textual identifier or label
//! for the gate, while `nr_affected_bits()` returns the number of qubits on
//! which the gate operates. The `matrix()` method should return a matrix of size
//! `2`<sup>`n`</sup>`×2`<sup>`n`</sup>, where `n` is the number of affected bits,
//! that describes the unitary transformation that the gate implements. An example
//! of a simple custom gate that rotates the `|01⟩` and `|10⟩` components of
//! a pair of qubits, is given below:
//! ```
//! use ndarray::array;
//! use q1tsim::ExportGate;
//!
//! #[derive(ExportGate)]
//! struct Mix
//! {
//!    alpha: f64
//! }
//!
//! impl q1tsim::gates::Gate for Mix
//! {
//!     fn description(&self) -> &str { "M" }
//!     fn nr_affected_bits(&self) -> usize { 2 }
//!     fn matrix(&self) -> q1tsim::cmatrix::CMatrix
//!     {
//!         let o = q1tsim::cmatrix::COMPLEX_ONE;
//!         let z = q1tsim::cmatrix::COMPLEX_ZERO;
//!         let c = self.alpha.cos() * o;
//!         let s = self.alpha.sin() * o;
//!         array![
//!             [o, z,  z, z],
//!             [z, c, -s, z],
//!             [z, s,  c, z],
//!             [z, z,  z, o]
//!         ]
//!     }
//! }
//! ```
//! Types implementing the `Gate` trait may optionally also implement the
//! [apply()](gates/trait.Gate.html#method.apply) family of methods if a more
//! optimal implementation than simply multiplying by its associated matrix can
//! be found.
//!
//! Exporting gates and circuits
//! ============================
//! The discerning reader may have notices the `#[derive(ExportGate)]` statement
//! on the custom gate in the listing above. This makes the type use the default
//! implementations of the export functions for a gate. Currently, there are
//! three traits for exporting a gate:
//! - [OpenQasm](export/trait.OpenQasm.html) for exporting a gate to OpenQasm code.
//! - [CQasm](export/trait.CQasm.html) for exporting a gate to c-Qasm code.
//! - [Latex](export/trait.Latex.html) for exporting a gate to LaTeX.
//!
//! You can use the default implementation for each of these traits by deriving
//! them, e.g.
//! ```
//! use q1tsim::OpenQasm;
//! use q1tsim::gates::Gate;
//!
//! #[derive(OpenQasm)]
//! struct MySpecialGate {}
//!
//! impl Gate for MySpecialGate {
//!     /* ... */
//!     # fn description(&self) -> &str { "" }
//!     # fn nr_affected_bits(&self) -> usize { 0 }
//!     # fn matrix(&self) -> q1tsim::cmatrix::CMatrix { q1tsim::cmatrix::CMatrix::zeros((0,0)) }
//! }
//! ```
//! The default implementations for OpenQasm and CQasm simply return an error,
//! since there is no way [^no_qasm] to know how to encode a custom gate
//! in these formats. The default implementation for the LaTeX export simply draws
//! a rectangular box with the gate description inside. As seen before, if you
//! want to use default definitions for all export traits, derive from `ExportGate`.
//!
//! Note that to use a gate type in a circuit, it must be exportable, so an
//! implementation for the export traits must be defined for your custom type,
//! either through deriving or by providing your own implementation.
//!
//! [^no_qasm]: No reasonable way at least. Technically, we could take the matrix
//! for the gate, decompose it into primitive gates, and export the corresponding
//! code.
//!
//! Stabilizer circuits
//! ===================
//! Stabilizer circuits are circuits that can be expressed entirely in terms of
//! the Clifford gates `H`, `S`, and `CX`, and qubit measurements. Since version
//! 0.4.0, `q1tsim` can simulate these circuits much more efficiently than
//! general circuits. If you have a custom gate type that can be represented
//! in terms of Clifford gates, and wish to use it with the stabilizer backend,
//! you should override the default implementations of the
//! [is_stabilizer()](gates/trait.Gate.html#method.is_stabilizer) and
//! [conjugate()](gates/trait.Gate.html#method.conjugate) methods. As an example,
//! the implementation for a hypothetical `HX` gate that first performs a Hadamard
//! transform, followed by an `X` gate, could look like
//! ```
//! use q1tsim::stabilizer::PauliOp;
//! use q1tsim::gates::Gate;
//!
//! struct HX {}
//!
//! impl Gate for HX {
//!     # fn description(&self) -> &str { "" }
//!     # fn nr_affected_bits(&self) -> usize { 0 }
//!     # fn matrix(&self) -> q1tsim::cmatrix::CMatrix { q1tsim::cmatrix::CMatrix::zeros((0,0)) }
//!     fn is_stabilizer(&self) -> bool
//!     {
//!         true
//!     }
//!
//!     fn conjugate(&self, ops: &mut [PauliOp]) -> q1tsim::error::Result<bool>
//!     {
//!         let (op, sign) = match ops[0]
//!             {
//!                 PauliOp::I => (PauliOp::I, false),
//!                 PauliOp::Z => (PauliOp::X, true),
//!                 PauliOp::X => (PauliOp::Z, false),
//!                 PauliOp::Y => (PauliOp::Y, false)
//!             };
//!         ops[0] = op;
//!         Ok(sign)
//!     }
//!
//!     /* ... */
//! }
//! ```

#[macro_use] extern crate ndarray;
#[cfg(test)] #[macro_use] extern crate matches;

#[macro_use] pub mod cmatrix;
#[macro_use] pub mod gates;
pub mod circuit;
pub mod error;
pub mod export;
pub mod permutation;
pub mod qustate;
pub mod vectorstate;
pub mod stabilizer;

mod idhash;
mod support;
#[cfg(test)] mod stats;

pub use q1tsim_derive::*;
