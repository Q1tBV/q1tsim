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

use crate::export::{CircuitGate, CQasm, OpenQasm};

/// Basis in which to perform measurements
#[derive(Clone, Copy)]
pub enum Basis
{
    /// Pauli `X` basis
    X,
    /// Pauli `Y` basis
    Y,
    /// Pauli `Z` basis
    Z
}

/// A single operation in a circuit
enum CircuitOp
{
    /// Apply a gate to the state
    Gate(Box<CircuitGate>, Vec<usize>),
    /// Conditionally apply a gate, depending on classical bits
    ConditionalGate(Vec<usize>, u64, Box<CircuitGate>, Vec<usize>),
    /// Reset a qubit to |0〉
    Reset(usize),
    /// Reset the quantum state to |00...0〉
    ResetAll,
    /// Measure a qubit in a certain basis
    Measure(usize, usize, Basis),
    /// Measure all qubits
    MeasureAll(Vec<usize>, Basis),
    /// Measure a single qubit in a certain basis without affecting state
    Peek(usize, usize, Basis),
    /// Measure all qubits in a certain basis without affecting state
    PeekAll(Vec<usize>, Basis),
    /// Prevent gate reordering on the associated bits across the barrier
    Barrier(Vec<usize>)
}

/// A quantum circuit
///
/// Struct Circuit represents a quantum circuit, holding a quantum state and the
/// operations to be performed on it.
pub struct Circuit
{
    /// The number of quantum bits in the system
    nr_qbits: usize,
    /// The number of classical bit in the system
    nr_cbits: usize,
    /// The quantum state of the system
    q_state: Option<crate::qustate::QuState>,
    /// The classial state of the system
    c_state: Option<ndarray::Array1<u64>>,
    /// The operations to perform on the state
    ops: Vec<CircuitOp>
}

impl Circuit
{
    /// Create a new circuit.
    ///
    /// Create a new (empty) quantum circuit, with `nr_qbits` quantum bits and
    /// `nr_cbits` classical bits.
    pub fn new(nr_qbits: usize, nr_cbits: usize) -> Self
    {
        Circuit
        {
            nr_qbits: nr_qbits,
            nr_cbits: nr_cbits,
            q_state: None,
            c_state: None,
            ops: vec![]
        }
    }

    /// The number of quantum bits in this circuit
    pub fn nr_qbits(&self) -> usize
    {
        self.nr_qbits
    }

    /// The number of classical bits in this circuit
    pub fn nr_cbits(&self) -> usize
    {
        self.nr_cbits
    }

    /// The classical register.
    ///
    /// Return a reference to the classical bit register, containing the results
    /// of any measurements made on the system. If no experiment has been run
    /// yet, `None` is returned.
    pub fn cstate(&self) -> Option<&ndarray::Array1<u64>>
    {
        self.c_state.as_ref()
    }

    /// Add a gate.
    ///
    /// Append a `n`-ary gate `gate`, operating on the `n` qubits in `bits`, to
    /// this circuit.
    pub fn add_gate<G: 'static>(&mut self, gate: G, bits: &[usize]) -> crate::error::Result<()>
    where G: CircuitGate
    {
        if let Some(&bit) = bits.iter().find(|&&b| b >= self.nr_qbits)
        {
            Err(crate::error::Error::InvalidQBit(bit))
        }
        else
        {
            self.ops.push(CircuitOp::Gate(Box::new(gate), bits.to_owned()));
            Ok(())
        }
    }

    /// Add a conditional gate.
    ///
    /// Append a `n`-ary gate `gate`, that will operate on the `n` qubits in
    /// `bits` to this circuit. The gate will only be applied only when the
    /// classical bits with indices from `control` form the target word `target`.
    /// The bit at the position of the first index in `control` is interpreted
    /// as the most significant bit to check.
    pub fn add_conditional_gate<G: 'static>(&mut self, control: &[usize],
        target: u64, gate: G, qbits: &[usize]) -> crate::error::Result<()>
    where G: CircuitGate
    {
        if let Some(&bit) = control.iter().find(|&&b| b >= self.nr_cbits)
        {
            Err(crate::error::Error::InvalidCBit(bit))
        }
        else if let Some(&bit) = qbits.iter().find(|&&b| b >= self.nr_qbits)
        {
            Err(crate::error::Error::InvalidQBit(bit))
        }
        else
        {
            self.ops.push(CircuitOp::ConditionalGate(control.to_owned(), target,
                Box::new(gate), qbits.to_owned()));
            Ok(())
        }
    }

    /// Add a measurement
    ///
    /// Add measurement of qubit `qbit` in basis `basis`, into classical bit
    /// `cbit`, to this circuit.
    pub fn measure_basis(&mut self, qbit: usize, cbit: usize, basis: Basis)
         -> crate::error::Result<()>
    {
        if qbit >= self.nr_qbits
        {
            Err(crate::error::Error::InvalidQBit(qbit))
        }
        else if cbit >= self.nr_cbits
        {
            Err(crate::error::Error::InvalidCBit(cbit))
        }
        else
        {
            self.ops.push(CircuitOp::Measure(qbit, cbit, basis));
            Ok(())
        }
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `X` basis, into classical
    /// bit `cbit` to this circuit.
    #[inline(always)]
    pub fn measure_x(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.measure_basis(qbit, cbit, Basis::X)
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Y` basis, into classical
    /// bit `cbit` to this circuit.
    #[inline(always)]
    pub fn measure_y(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.measure_basis(qbit, cbit, Basis::Y)
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Z` basis, into classical
    /// bit `cbit` to this circuit.
    #[inline(always)]
    pub fn measure_z(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.measure_basis(qbit, cbit, Basis::Z)
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` into classical bit `cbit` to this circuit.
    /// This is an alias for `measure_z()`.
    #[inline(always)]
    pub fn measure(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.measure_basis(qbit, cbit, Basis::Z)
    }

    /// Add a measurement.
    ///
    /// Add the measurement of all qubits in the quantum state into the classical
    /// bits `cbits`. Measurement is done in basis `basis`.
    pub fn measure_all_basis(&mut self, cbits: &[usize], basis: Basis)
         -> crate::error::Result<()>
    {
        if let Some(&bit) = cbits.iter().find(|&&b| b >= self.nr_cbits)
        {
            Err(crate::error::Error::InvalidCBit(bit))
        }
        else
        {
            self.ops.push(CircuitOp::MeasureAll(cbits.to_owned(), basis));
            Ok(())
        }
    }

    /// Add a measurement.
    ///
    /// Add the measurement of all qubits in the quantum state into the classical
    /// bits `cbits`. Measurement is done in the Pauli `Z` basis.
    #[inline(always)]
    pub fn measure_all(&mut self, cbits: &[usize])-> crate::error::Result<()>
    {
        self.measure_all_basis(cbits, Basis::Z)
    }

    /// Add a measurement.
    ///
    /// Add the measurement of qubit `qbit` in the quantum state into the
    /// classical bit `cbit`. Measurement is done in basis `basis`, without
    /// collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    pub fn peek_basis(&mut self, qbit: usize, cbit: usize, basis: Basis)
        -> crate::error::Result<()>
    {
        if qbit >= self.nr_qbits
        {
            Err(crate::error::Error::InvalidQBit(qbit))
        }
        else if cbit >= self.nr_cbits
        {
            Err(crate::error::Error::InvalidCBit(cbit))
        }
        else
        {
            self.ops.push(CircuitOp::Peek(qbit, cbit, basis));
            Ok(())
        }
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `X` basis, into classical
    /// bit `cbit` to this circuit, without collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    #[inline(always)]
    pub fn peek_x(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.peek_basis(qbit, cbit, Basis::X)
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Y` basis, into classical
    /// bit `cbit` to this circuit, without collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    #[inline(always)]
    pub fn peek_y(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.peek_basis(qbit, cbit, Basis::Y)
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Z` basis, into classical
    /// bit `cbit` to this circuit, without collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    #[inline(always)]
    pub fn peek_z(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.peek_basis(qbit, cbit, Basis::Z)
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Z` basis, into classical
    /// bit `cbit` to this circuit, without collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    #[inline(always)]
    pub fn peek(&mut self, qbit: usize, cbit: usize) -> crate::error::Result<()>
    {
        self.peek_basis(qbit, cbit, Basis::Z)
    }

    /// Add a measurement.
    ///
    /// Add the measurement of all qubits in the quantum state into the classical
    /// bits `cbits`. Measurement is done in basis `basis`, without
    /// collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    pub fn peek_all_basis(&mut self, cbits: &[usize], basis: Basis)
        -> crate::error::Result<()>
    {
        if let Some(&bit) = cbits.iter().find(|&&b| b >= self.nr_cbits)
        {
            Err(crate::error::Error::InvalidCBit(bit))
        }
        else
        {
            self.ops.push(CircuitOp::PeekAll(cbits.to_owned(), basis));
            Ok(())
        }
    }


    /// Add a measurement.
    ///
    /// Add the measurement of all qubits in the quantum state into the classical
    /// bits `cbits`. Measurement is done in the Pauli `Z` basis, without
    /// collapsing the quantum state.
    /// NOTE: this is not a physical process, and cannot be reproduced on a real
    /// quantum computer.
    #[inline(always)]
    pub fn peek_all(&mut self, cbits: &[usize]) -> crate::error::Result<()>
    {
        self.peek_all_basis(cbits, Basis::Z)
    }

    /// Reset a qubit
    ///
    /// Reset the qubit `qbit` to |0〉. This is done by measuring the bit, and
    /// flipping it if the result is `1`, so this is potentially an expensive
    /// operation.
    pub fn reset(&mut self, qbit: usize) -> crate::error::Result<()>
    {
        if qbit >= self.nr_qbits
        {
            Err(crate::error::Error::InvalidQBit(qbit))
        }
        else
        {
            self.ops.push(CircuitOp::Reset(qbit));
            Ok(())
        }
    }

    /// Reset all qubits
    ///
    /// Reset the entire quantum state of the circuit to |00...0〉. The classical
    /// register is not affected.
    pub fn reset_all(&mut self)
    {
        self.ops.push(CircuitOp::ResetAll);
    }

    /// Add a Hadamard gate.
    ///
    /// Add a Hadamard operating on qubit `qbit`, to this circuit.
    pub fn h(&mut self, qbit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::H::new(), &[qbit])
    }

    /// Add a Pauli X gate.
    ///
    /// Add a Pauli X gate operating on qubit `bit`, to this circuit.
    pub fn x(&mut self, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::X::new(), &[bit])
    }

    /// Add a Pauli Y gate.
    ///
    /// Add a Pauli Y gate operating on qubit `bit`, to this circuit.
    pub fn y(&mut self, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::Y::new(), &[bit])
    }

    /// Add a Pauli Z gate.
    ///
    /// Add a Pauli Z gate operating on qubit `bit`, to this circuit.
    pub fn z(&mut self, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::Z::new(), &[bit])
    }

    /// Add a phase gate
    ///
    /// Add an `S` phase gate operating on qubit `bit`, to this circuit.
    pub fn s(&mut self, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::S::new(), &[bit])
    }

    /// Add a phase gate
    ///
    /// Add an `S`<sup>`\dagger`</sup> phase gate operating on qubit `bit`, to this circuit.
    pub fn sdg(&mut self, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::Sdg::new(), &[bit])
    }

    /// Add a R<sub>X</sub> gate.
    ///
    /// Add a `R`<sub>`X`</sub>`(θ)` gate operating on qubit `bit`, to this circuit.
    pub fn rx(&mut self, theta: f64, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::RX::new(theta), &[bit])
    }

    /// Add a R<sub>Y</sub> gate.
    ///
    /// Add a `R`<sub>`Y`</sub>`(θ)` gate operating on qubit `bit`, to this circuit.
    pub fn ry(&mut self, theta: f64, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::RY::new(theta), &[bit])
    }

    /// Add a R<sub>Z</sub> gate.
    ///
    /// Add a `R`<sub>`Z`</sub>`(λ)` gate operating on qubit `bit`, to this circuit.
    pub fn rz(&mut self, lambda: f64, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::RZ::new(lambda), &[bit])
    }

    /// Add a U<sub>1</sub> gate.
    ///
    /// Add a `U`<sub>`1`</sub>`(λ)` gate operating on qubit `bit`, to this circuit.
    pub fn u1(&mut self, lambda: f64, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::U1::new(lambda), &[bit])
    }

    /// Add a U<sub>2</sub> gate.
    ///
    /// Add a `U`<sub>`2`</sub>`(ϕ, λ)` gate operating on qubit `bit`, to this circuit.
    pub fn u2(&mut self, phi: f64, lambda: f64, bit: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::U2::new(phi, lambda), &[bit])
    }

    /// Add a U<sub>3</sub> gate.
    ///
    /// Add a `U`<sub>`3`</sub>`(θ, ϕ, λ)` gate operating on qubit `bit`, to this circuit.
    pub fn u3(&mut self, theta: f64, phi: f64, lambda: f64, bit: usize)
         -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::U3::new(theta, phi, lambda), &[bit])
    }

    /// Add a C<sub>X</sub> gate.
    ///
    /// Add a `C`<sub>`X`</sub> gate, controlled by qubit `control` and
    /// operating on qubit `target`, to this circuit.
    pub fn cx(&mut self, control: usize, target: usize) -> crate::error::Result<()>
    {
        self.add_gate(crate::gates::CX::new(), &[control, target])
    }

    /// Add a barrier
    ///
    /// Add a barrier on the bits in `bits`. No transformations on these bits
    /// are allowed across this barrier.
    pub fn barrier(&mut self, qbits: &[usize]) -> crate::error::Result<()>
    {
        if let Some(&bit) = qbits.iter().find(|&&b| b >= self.nr_qbits)
        {
            Err(crate::error::Error::InvalidQBit(bit))
        }
        else
        {
            self.ops.push(CircuitOp::Barrier(qbits.to_vec()));
            Ok(())
        }
    }

    /// Execute this circuit
    ///
    /// Execute this circuit, performing its operations and measurements.
    /// Measurements are made over `nr_shots` executions of the circuit. This
    /// function clears any previous states of the system (quantum or classical).
    pub fn execute(&mut self, nr_shots: usize)
    {
        self.execute_with_rng(nr_shots, &mut rand::thread_rng());
    }

    /// Execute this circuit
    ///
    /// Execute this circuit, performing its operations and measurements.
    /// Measurements are made over `nr_shots` executions of the circuit, using
    /// random number generator `rng` for sampling. This function clears any
    /// previous states of the system (quantum or classical).
    pub fn execute_with_rng<R: rand::RngCore>(&mut self, nr_shots: usize, rng: &mut R)
    {
        self.q_state = Some(crate::qustate::QuState::new(self.nr_qbits, nr_shots));
        self.c_state = Some(ndarray::Array::zeros(nr_shots));
        // It shpuld not be possible to get an error from reexecute, so
        // ignore it for now.
        self.reexecute_with_rng(rng).unwrap();
    }

    /// Execute a circuit again.
    ///
    /// Run this circuit again, starting with the state from the previous
    /// execution. If this circuit has not been run before, a `NotExecuted`
    /// error is returned.
    pub fn reexecute(&mut self) -> crate::error::Result<()>
    {
        self.reexecute_with_rng(&mut rand::thread_rng())
    }

    /// Execute a circuit again.
    ///
    /// Run this circuit again, starting with the state from the previous
    /// execution, using random number generator `rng` for sampling. If this
    /// circuit has not been run before, a `NotExecuted` error is returned.
    pub fn reexecute_with_rng<R: rand::Rng>(&mut self, rng: &mut R)
        -> crate::error::Result<()>
    {
        if self.q_state.is_none() || self.c_state.is_none()
        {
            return Err(crate::error::Error::NotExecuted);
        }

        let q_state = self.q_state.as_mut().unwrap();
        let c_state = self.c_state.as_mut().unwrap();

        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    q_state.apply_gate(&**gate, bits.as_slice());
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    let mut cbits = vec![0; c_state.len()];
                    for (idst, &isrc) in control.iter().enumerate()
                    {
                        for (db, &sb) in cbits.iter_mut().zip(c_state.iter())
                        {
                            *db |= ((sb >> isrc) & 1) << idst;
                        }
                    }
                    let apply_gate: Vec<bool> = cbits.iter()
                        .map(|&b| b == target)
                        .collect();
                    q_state.apply_conditional_gate(&apply_gate, &**gate,
                        bits.as_slice());
                },
                CircuitOp::Measure(qbit, cbit, basis) => {
                    match basis
                    {
                        Basis::X => {
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                            q_state.measure_into(qbit, cbit, c_state, rng);
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                        },
                        Basis::Y => {
                            q_state.apply_gate(&crate::gates::Sdg::new(), &[qbit]);
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                            q_state.measure_into(qbit, cbit, c_state, rng);
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                            q_state.apply_gate(&crate::gates::S::new(), &[qbit]);
                        },
                        Basis::Z => {
                            q_state.measure_into(qbit, cbit, c_state, rng);
                        }
                    }
                }
                CircuitOp::MeasureAll(ref cbits, basis) => {
                    match basis
                    {
                        Basis::X => {
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                            q_state.measure_all_into(cbits, c_state, rng);
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                        },
                        Basis::Y => {
                            q_state.apply_unary_gate_all(&crate::gates::Sdg::new());
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                            q_state.measure_all_into(cbits, c_state, rng);
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                            q_state.apply_unary_gate_all(&crate::gates::S::new());
                        },
                        Basis::Z => {
                            q_state.measure_all_into(cbits, c_state, rng);
                        }
                    }
                },
                CircuitOp::Peek(qbit, cbit, basis) => {
                    match basis
                    {
                        Basis::X => {
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                            q_state.peek_into(qbit, cbit, c_state, rng);
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                        },
                        Basis::Y => {
                            q_state.apply_gate(&crate::gates::Sdg::new(), &[qbit]);
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                            q_state.peek_into(qbit, cbit, c_state, rng);
                            q_state.apply_gate(&crate::gates::H::new(), &[qbit]);
                            q_state.apply_gate(&crate::gates::S::new(), &[qbit]);
                        },
                        Basis::Z => {
                            q_state.peek_into(qbit, cbit, c_state, rng);
                        }
                    }
                },
                CircuitOp::PeekAll(ref cbits, basis) => {
                    match basis
                    {
                        Basis::X => {
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                            q_state.peek_all_into(cbits, c_state, rng);
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                        },
                        Basis::Y => {
                            q_state.apply_unary_gate_all(&crate::gates::Sdg::new());
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                            q_state.peek_all_into(cbits, c_state, rng);
                            q_state.apply_unary_gate_all(&crate::gates::H::new());
                            q_state.apply_unary_gate_all(&crate::gates::S::new());
                        },
                        Basis::Z => {
                            q_state.peek_all_into(cbits, c_state, rng);
                        }
                    }
                },
                CircuitOp::Reset(bit) => {
                    q_state.reset(bit, rng);
                },
                CircuitOp::ResetAll => {
                    q_state.reset_all();
                },
                CircuitOp::Barrier(_) => {
                    /* Nothing to be done */
                }
            }
        }

        Ok(())
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a single `u64` integer value. The
    /// first bit in the classical register (at index 0) corresponds to the
    /// least significant bit in the key; the last classical bit (at index `n-1`)
    /// to the most significant bit in the key. This function of course only works
    /// when there are at most 64 bits in the register. If there are more, use
    /// `histogram_string()`.
    pub fn histogram(&self) -> crate::error::Result<crate::idhash::U64HashMap<usize>>
    {
        if let Some(ref c_state) = self.c_state
        {
            let mut res = crate::idhash::new_u64_hash_map();
            for &key in c_state
            {
                let count = res.entry(key).or_insert(0);
                *count += 1;
            }
            Ok(res)
        }
        else
        {
            Err(crate::error::Error::NotExecuted)
        }
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a single `usize` integer value,
    /// which is used as an index in a vector. The bit order of the indices
    /// is the same as in the `histogram()` function. The vector is of length
    /// `2`<sub>`n`</sub>, so use this function only for reasonably small
    /// numbers of `n`. For sparse collections, using `histogram()` or
    /// `histogram_string` may be better.
    pub fn histogram_vec(&self) -> crate::error::Result<Vec<usize>>
    {
        if let Some(ref c_state) = self.c_state
        {
            let mut res = vec![0; 1 << self.nr_cbits];
            for &key in c_state
            {
                res[key as usize] += 1;
            }
            Ok(res)
        }
        else
        {
            Err(crate::error::Error::NotExecuted)
        }
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a string key, with the last character
    /// in the key corresponding to the first bit (at index 0) in the classical
    /// register and vice versa.
    pub fn histogram_string(&self)
        -> crate::error::Result<::std::collections::HashMap<String, usize>>
    {
        if let Some(ref c_state) = self.c_state
        {
            let mut res = ::std::collections::HashMap::new();
            for &key in c_state
            {
                let skey = format!("{:0width$b}", key, width=self.nr_cbits);
                let count = res.entry(skey).or_insert(0);
                *count += 1;
            }
            Ok(res)
        }
        else
        {
            Err(crate::error::Error::NotExecuted)
        }
    }

    fn is_full_register(&self, control: &[usize]) -> bool
    {
        let n = control.len();
        if n != self.nr_cbits
        {
            return false;
        }

        let mut scontrol = vec![0; n];
        scontrol.copy_from_slice(control);
        scontrol.sort();
        for i in 0..n
        {
            if scontrol[i] != i
            {
                return false;
            }
        }

        true
    }

    fn check_open_qasm_condition_bits(&self, control: &[usize])
        -> crate::error::ExportResult<()>
    {
        if !self.is_full_register(control)
        {
            Err(crate::error::ExportError::IncompleteConditionRegister)
        }
        else
        {
            Ok(())
        }
    }

    /// Export to OpenQasm
    ///
    /// Export this circuit to a program in OpenQasm format. On a successful
    /// conversion, the result is `Ok` with the program text. When the conversion
    /// to OpenQasm fails, `Err` with an error message is returned.
    pub fn open_qasm(&self) -> crate::error::Result<String>
    {
        let mut res = String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n");

        let mut qbit_names = vec![];
        if self.nr_qbits > 0
        {
            res += &format!("qreg q[{}];\n", self.nr_qbits);
            for i in 0..self.nr_qbits
            {
                qbit_names.push(format!("q[{}]", i));
            }
        }
        let mut cbit_names = vec![];
        if self.nr_cbits > 0
        {
            res += &format!("creg b[{}];\n", self.nr_cbits);
            for i in 0..self.nr_cbits
            {
                cbit_names.push(format!("b[{}]", i));
            }
        }

        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    res += &format!("{};\n", gate.open_qasm(&qbit_names, bits)?);
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    if control.is_empty()
                    {
                        res += &format!("{};\n", gate.open_qasm(&qbit_names, bits)?);
                    }
                    else
                    {
                        // We do require that the control bits span the entire classical
                        // register, but not necessarily in the order 0..#bits.
                        self.check_open_qasm_condition_bits(control)?;
                        let mut starget = 0;
                        for (tshift, sshift) in control.iter().enumerate()
                        {
                            starget |= ((target >> tshift) & 0x01) << sshift;
                        }
                        let condition = format!("b == {}", starget);
                        let gate_qasm = gate.conditional_open_qasm(&condition, &qbit_names, bits)?;
                        res += &format!("{};\n", gate_qasm);
                    }
                },
                CircuitOp::Measure(qbit, cbit, basis) => {
                    match basis
                    {
                        Basis::X => {
                            res += &format!("{};\n",
                                crate::gates::H::new().open_qasm(&qbit_names, &[qbit])?);
                        },
                        Basis::Y => {
                            res += &format!("{};\n",
                                crate::gates::Sdg::new().open_qasm(&qbit_names, &[qbit])?);
                            res += &format!("{};\n",
                                crate::gates::H::new().open_qasm(&qbit_names, &[qbit])?);
                        }
                        _ => {}
                    }
                    res += &format!("measure {} -> {};\n", qbit_names[qbit], cbit_names[cbit]);
                }
                CircuitOp::MeasureAll(ref cbits, basis) => {
                    match basis
                    {
                        Basis::X => {
                            let names = [String::from("q")];
                            res += &format!("{};\n",
                                crate::gates::H::new().open_qasm(&names, &[0])?);
                        },
                        Basis::Y => {
                            let names = [String::from("q")];
                            res += &format!("{};\n",
                                crate::gates::Sdg::new().open_qasm(&names, &[0])?);
                            res += &format!("{};\n",
                                crate::gates::H::new().open_qasm(&names, &[0])?);
                        }
                        _ => {}
                    }

                    if cbits.len() == self.nr_cbits
                        && cbits.iter().enumerate().all(|(i, &b)| i==b)
                    {
                        res += &format!("measure q -> b;\n");
                    }
                    else
                    {
                        for (qbit, &cbit) in cbits.iter().enumerate()
                        {
                            res += &format!("measure {} -> {};\n", qbit_names[qbit],
                                cbit_names[cbit]);
                        }
                    }
                },
                CircuitOp::Peek(_, _, _) => {
                    return Err(crate::error::Error::from(
                        crate::error::ExportError::ExportPeekInvalid("OpenQasm")
                    ));
                },
                CircuitOp::PeekAll(_, _) => {
                    return Err(crate::error::Error::from(
                        crate::error::ExportError::ExportPeekInvalid("OpenQasm")
                    ));
                },
                CircuitOp::Reset(qbit) => {
                    res += &format!("reset {};\n", qbit_names[qbit]);
                },
                CircuitOp::ResetAll => {
                    res += "reset q;\n";
                },
                CircuitOp::Barrier(ref qbits) => {
                    if qbits.len() == self.nr_qbits
                        && qbits.iter().enumerate().all(|(i, &b)| i==b)
                    {
                        res += "barrier q;\n";
                    }
                    else
                    {
                        res += &format!("barrier {};\n",
                            qbits.iter()
                            .map(|&b| qbit_names[b].as_str())
                            .collect::<Vec<&str>>()
                            .join(", "));
                    }
                }
            }
        }

        Ok(res)
    }

    fn check_c_qasm_measurement(qbit: usize, cbit: usize) -> crate::error::ExportResult<()>
    {
        if qbit != cbit
        {
            Err(crate::error::ExportError::NoClassicalRegister)
        }
        else
        {
            Ok(())
        }
    }

    /// Export to c-Qasm
    ///
    /// Export this circuit to a program in c-Qasm format. On a successful
    /// conversion, the result is `Ok` with the program text. When the conversion
    /// to OpenQasm fails, `Err` with an error message is returned.
    pub fn c_qasm(&self) -> crate::error::Result<String>
    {
        let mut res = String::from("version 1.0\n");

        let mut qbit_names = vec![];
        let mut cbit_names = vec![];
        if self.nr_qbits > 0
        {
            res += &format!("qubits {}\n", self.nr_qbits);
            for i in 0..self.nr_qbits
            {
                qbit_names.push(format!("q[{}]", i));
                cbit_names.push(format!("b[{}]", i));
            }
        }

        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    res += &format!("{}\n", gate.c_qasm(&qbit_names, bits)?);
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    if control.is_empty()
                    {
                        res += &format!("{}\n", gate.c_qasm(&qbit_names, bits)?);
                    }
                    else
                    {
                        let mut conditions = vec![];
                        for (shift, &idx) in control.iter().enumerate()
                        {
                            if target & (1 << shift) == 0
                            {
                                res += &format!("not {}\n", cbit_names[idx]);
                            }
                            conditions.push(cbit_names[idx].as_str());
                        }
                        let condition = conditions.join(", ");
                        let gate_qasm = gate.conditional_c_qasm(&condition,
                            &qbit_names, bits)?;
                        res += &format!("{}\n", gate_qasm);
                        for (shift, &idx) in control.iter().enumerate()
                        {
                            if target & (1 << shift) == 0
                            {
                                res += &format!("not {}\n", cbit_names[idx]);
                            }
                        }
                    }
                },
                CircuitOp::Measure(qbit, cbit, basis) => {
                    Self::check_c_qasm_measurement(qbit, cbit)?;
                    let op = match basis
                    {
                        Basis::X => "measure_x",
                        Basis::Y => "measure_y",
                        _        => "measure"
                    };
                    res += &format!("{} q[{}]\n", op, qbit);
                }
                CircuitOp::MeasureAll(ref cbits, basis) => {
                    for (qbit, &cbit) in cbits.iter().enumerate()
                    {
                        Self::check_c_qasm_measurement(qbit, cbit)?;
                    }
                    match basis
                    {
                        Basis::X => {
                            for bit in 0..self.nr_qbits
                            {
                                res += &format!("{}\n",
                                    crate::gates::H::new().c_qasm(&qbit_names, &[bit])?);
                            }
                        },
                        Basis::Y => {
                            for bit in 0..self.nr_qbits
                            {
                                res += &format!("{}\n",
                                    crate::gates::Sdg::new().c_qasm(&qbit_names, &[bit])?);
                                res += &format!("{}\n",
                                    crate::gates::H::new().c_qasm(&qbit_names, &[bit])?);
                            }
                        },
                        _ => {
                            /* do nothing */
                        }
                    }
                    res += &format!("measure_all\n");
                },
                CircuitOp::Peek(_, _, _) => {
                    return Err(crate::error::Error::from(
                        crate::error::ExportError::ExportPeekInvalid("c-Qasm")
                    ));
                },
                CircuitOp::PeekAll(_, _) => {
                    return Err(crate::error::Error::from(
                        crate::error::ExportError::ExportPeekInvalid("c-Qasm")
                    ));
                },
                CircuitOp::Reset(qbit) => {
                    res += &format!("prep_z {}\n", qbit_names[qbit]);
                },
                CircuitOp::ResetAll => {
                    for i in 0..self.nr_qbits
                    {
                        res += &format!("prep_z {}\n", qbit_names[i]);
                    }
                },
                CircuitOp::Barrier(_) => {
                    /* Not available */
                }
            }
        }

        Ok(res)
    }

    pub fn latex(&self) -> crate::error::Result<String>
    {
        let mut state = crate::export::LatexExportState::new(self.nr_qbits, self.nr_cbits);
        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    gate.latex(bits, &mut state)?;
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    state.start_range_op(bits, Some(control))?;
                    let controlled = state.set_controlled(true);
                    gate.latex(bits, &mut state)?;
                    state.set_controlled(controlled);
                    state.set_condition(control, target, bits)?;
                    state.end_range_op();
                },
                CircuitOp::Measure(qbit, cbit, basis) => {
                    let basis_lbl = match basis
                    {
                        Basis::X => Some("X"),
                        Basis::Y => Some("Y"),
                        _        => None
                    };
                    state.set_measurement(qbit, cbit, basis_lbl)?;
                }
                CircuitOp::MeasureAll(ref cbits, basis) => {
                    let basis_lbl = match basis
                    {
                        Basis::X => Some("X"),
                        Basis::Y => Some("Y"),
                        _        => None
                    };
                    for (qbit, &cbit) in cbits.iter().enumerate()
                    {
                        state.set_measurement(qbit, cbit, basis_lbl)?;
                    }
                },
                CircuitOp::Peek(_, _, _) => {
                    return Err(crate::error::Error::from(
                        crate::error::ExportError::NotImplemented("LaTeX",
                            String::from("peek")
                        )
                    ));
                },
                CircuitOp::PeekAll(_, _) => {
                    return Err(crate::error::Error::from(
                        crate::error::ExportError::NotImplemented("LaTeX",
                            String::from("peek all")
                        )
                    ));
                },
                CircuitOp::Reset(qbit) => {
                    state.set_reset(qbit)?;
                },
                CircuitOp::ResetAll => {
                    state.start_range_op(&[0, self.nr_qbits-1], None)?;
                    for qbit in 0..self.nr_qbits
                    {
                        state.set_reset(qbit)?;
                    }
                    state.end_range_op();
                },
                CircuitOp::Barrier(ref qbits) => {
                    state.set_barrier(qbits)?;
                }
            }
        }

        Ok(state.code())
    }
}

#[macro_export]
macro_rules! circuit_method_check
{
    ( add_conditional_gate $res:expr ) => { $res? };
    ( add_gate $res:expr ) => { $res? };
    ( barrier $res:expr ) => { $res? };
    ( cx $res:expr ) => { $res? };
    ( h $res:expr ) => { $res? };
    ( measure $res:expr ) => { $res? };
    ( measure_all $res:expr ) => { $res? };
    ( measure_all_basis $res:expr ) => { $res? };
    ( measure_x $res:expr ) => { $res? };
    ( measure_y $res:expr ) => { $res? };
    ( measure_z $res:expr ) => { $res? };
    ( peek $res:expr ) => { $res? };
    ( peek_x $res:expr ) => { $res? };
    ( peek_y $res:expr ) => { $res? };
    ( peek_z $res:expr ) => { $res? };
    ( peek_all $res:expr ) => { $res? };
    ( peek_all_basis $res:expr ) => { $res? };
    ( reset $res:expr ) => { $res? };
    ( s $res:expr ) => { $res? };
    ( sdg $res:expr ) => { $res? };
    ( x $res:expr ) => { $res? };
    ( y $res:expr ) => { $res? };
    ( z $res:expr ) => { $res? };
    ( $name:ident $res:expr ) => { $res };
}

#[macro_export]
macro_rules! circuit
{
    ($nr_qbits:expr, $nr_cbits:expr, { $( $method_name:ident ( $( $arg:expr ),* ) );* ; } ) => {
        {
            let generator = || {
                let mut circuit = $crate::circuit::Circuit::new($nr_qbits, $nr_cbits);
                $(
                    circuit_method_check!(
                        $method_name
                        circuit.$method_name($($arg),*)
                    );
                );*
                Ok(circuit) as $crate::error::Result<$crate::circuit::Circuit>
            };
            generator()
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::{Basis, Circuit, CircuitOp};
    use crate::gates::{CX, H, S, X};

    #[test]
    fn test_gate_methods()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;

        let mut circuit = Circuit::new(2, 0);
        assert_eq!(circuit.h(0), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), &array![[x, x], [x, -x]]);
                assert_eq!(bits, &vec![0]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not an H gate"),
            None => panic!("H gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.x(1), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[z, o], [o, z]]);
                assert_eq!(bits, &vec![1]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not an X gate"),
            None => panic!("X gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.y(0), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [i, z]]);
                assert_eq!(bits, &vec![0]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not a Y gate"),
            None => panic!("Y gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.z(1), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, -o]]);
                assert_eq!(bits, &vec![1]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not a Z gate"),
            None => panic!("Z gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.rx(::std::f64::consts::PI, 1), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [-i, z]]);
                assert_eq!(bits, &vec![1]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not an RX gate"),
            None => panic!("RX gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.ry(::std::f64::consts::PI, 0), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[z, -o], [o, z]]);
                assert_eq!(bits, &vec![0]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not an RY gate"),
            None => panic!("RY gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.rz(::std::f64::consts::PI, 1), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[-i, z], [z, i]]);
                assert_eq!(bits, &vec![1]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not an RZ gate"),
            None => panic!("RZ gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.u1(::std::f64::consts::FRAC_PI_4, 1), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, x*(o+i)]]);
                assert_eq!(bits, &vec![1]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not a U1 gate"),
            None => panic!("U1 gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.u2(::std::f64::consts::FRAC_PI_4,
            ::std::f64::consts::FRAC_PI_2, 0), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [x, -x*i],
                    [0.5*(o+i), 0.5*(-o+i)]
                ]);
                assert_eq!(bits, &vec![0]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not a U2 gate"),
            None => panic!("U2 gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.u3(::std::f64::consts::PI, ::std::f64::consts::FRAC_PI_4,
            ::std::f64::consts::FRAC_PI_2, 0), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [z, -i],
                    [x*(o+i), z]
                ]);
                assert_eq!(bits, &vec![0]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not a U3 gate"),
            None => panic!("U3 gate was not added")
            // LCOV_EXCL_STOP
        }

        assert_eq!(circuit.cx(1, 0), Ok(()));
        match circuit.ops.last()
        {
            Some(CircuitOp::Gate(gate, bits)) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, o, z, z],
                    [z, z, z, o],
                    [z, z, o, z]
                ]);
                assert_eq!(bits, &vec![1, 0]);
            },
            // LCOV_EXCL_START
            Some(_) => panic!("Value added was not a CX gate"),
            None => panic!("CX gate was not added")
            // LCOV_EXCL_STOP
        }
    }

    #[test]
    fn test_execute()
    {
        let nr_shots = 5;
        let mut circuit = circuit!(2, 2, {
            add_gate(X::new(), &[0]);
            add_gate(X::new(), &[1]);
            add_gate(CX::new(), &[0, 1]);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        assert_eq!(circuit.cstate(), Some(&array![0b01, 0b01, 0b01, 0b01, 0b01]));
    }

    #[test]
    fn test_measure()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            x(0);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            x(0);
            measure_x(0, 0);
            measure_x(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert!(hist.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));

        let mut circuit = circuit!(2, 2, {
            x(0);
            h(0);
            h(1);
            measure_x(0, 0);
            measure_x(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            x(0);
            measure_y(0, 0);
            measure_y(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert!(hist.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_peek()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(1, 3, {
            h(0);
            peek(0, 0);
            h(0);
            peek(0, 1);
            h(0);
            peek(0, 2);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram_vec().unwrap();
        // Results of first and third measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 0.
        let n00 = hist[0] + hist[2] + hist[4] + hist[6];
        assert!(crate::stats::measurement_ok(n00, nr_shots, 0.5, tol));
        let n10 = hist[0] + hist[1] + hist[4] + hist[5];
        assert!(n10 == nr_shots);
        let n20 = hist[0] + hist[1] + hist[2] + hist[3];
        assert!(crate::stats::measurement_ok(n20, nr_shots, 0.5, tol));

        let mut circuit = circuit!(2, 6, {
            h(0);
            h(1);
            peek(0, 0);
            h(0);
            peek(0, 1);
            h(0);
            peek(0, 2);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram().unwrap();
        // Results of first and third measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 0.
        let mut n0 = [0; 2];
        let mut n1 = [0; 2];
        let mut n2 = [0; 2];
        for (key, count) in hist
        {
            n0[key as usize & 1] += count;
            n1[(key as usize >> 1) & 1] += count;
            n2[(key as usize >> 2) & 1] += count;
        }
        assert!(n0.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.5, tol)
        ));
        assert_eq!(n1, [nr_shots, 0]);
        assert!(n2.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.5, tol)
        ));
    }

    #[test]
    fn test_peek_basis()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(1, 3, {
            peek_x(0, 0);
            h(0);
            peek_x(0, 1);
            h(0);
            peek_x(0, 2);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram_vec().unwrap();
        // Results of first and third measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 0.
        let n00 = hist[0] + hist[2] + hist[4] + hist[6];
        assert!(crate::stats::measurement_ok(n00, nr_shots, 0.5, tol));
        let n10 = hist[0] + hist[1] + hist[4] + hist[5];
        assert_eq!(n10, nr_shots);
        let n20 = hist[0] + hist[1] + hist[2] + hist[3];
        assert!(crate::stats::measurement_ok(n20, nr_shots, 0.5, tol));

        let mut circuit = circuit!(2, 6, {
            peek_y(0, 0);
            h(0);
            peek_y(0, 1);
            sdg(0);
            peek_y(0, 2);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram_vec().unwrap();
        // Results of first and second measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 1.
        let n00 = hist[0] + hist[2] + hist[4] + hist[6];
        assert!(crate::stats::measurement_ok(n00, nr_shots, 0.5, tol));
        let n10 = hist[0] + hist[1] + hist[4] + hist[5];
        assert!(crate::stats::measurement_ok(n10, nr_shots, 0.5, tol));
        let n20 = hist[0] + hist[1] + hist[2] + hist[3];
        assert_eq!(n20, 0);
    }

    #[test]
    fn test_conditional()
    {
        let mut circuit = circuit!(2, 2, {
            add_conditional_gate(&[0, 1], 1, X::new(), &[1]);
            measure_all(&[0, 1]);
        }).unwrap();
        circuit.execute(5);
        assert_eq!(circuit.c_state, Some(array![0b00, 0b00, 0b00, 0b00, 0b00]));

        let mut circuit = Circuit::new(2, 2);
        circuit.q_state = Some(crate::qustate::QuState::new(2, 5));
        circuit.c_state = Some(array![0b01, 0b10, 0b10, 0b11, 0b00]);
        circuit.add_conditional_gate(&[0, 1], 1, X::new(), &[1]).unwrap();
        circuit.measure_all(&[0, 1]).unwrap();
        circuit.reexecute().unwrap();
        assert_eq!(circuit.c_state, Some(array![0b10, 0b00, 0b00, 0b00, 0b00]));

        let mut circuit = Circuit::new(2, 2);
        circuit.q_state = Some(crate::qustate::QuState::new(2, 5));
        circuit.c_state = Some(array![0b01, 0b10, 0b10, 0b11, 0b00]);
        circuit.add_conditional_gate(&[0, 1], 2, X::new(), &[1]).unwrap();
        circuit.measure_all(&[0, 1]).unwrap();
        circuit.reexecute().unwrap();
        assert_eq!(circuit.c_state, Some(array![0b00, 0b10, 0b10, 0b00, 0b00]));

        let mut circuit = Circuit::new(2, 2);
        circuit.q_state = Some(crate::qustate::QuState::new(2, 5));
        circuit.c_state = Some(array![0b01, 0b10, 0b10, 0b11, 0b00]);
        circuit.add_conditional_gate(&[1], 1, X::new(), &[0]).unwrap();
        circuit.measure_all(&[0, 1]).unwrap();
        circuit.reexecute().unwrap();
        assert_eq!(circuit.c_state, Some(array![0b00, 0b01, 0b01, 0b01, 0b00]));
    }

    #[test]
    fn test_measure_all()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            x(0);
            measure_all(&[0, 1]);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            x(0);
            measure_all(&[1, 0]);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, 0, nr_shots, 0]);

        let mut circuit = circuit!(2, 2, {
            h(0);
            h(1);
            measure_all(&[0, 1]);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert!(hist.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_measure_all_basis()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            h(0);
            h(1);
            measure_all_basis(&[0, 1], Basis::X);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![nr_shots, 0, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            x(0);
            h(0);
            h(1);
            measure_all_basis(&[0, 1], Basis::X);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            x(0);
            h(0);
            h(1);
            add_gate(S::new(), &[0]);
            add_gate(S::new(), &[1]);
            measure_all_basis(&[0, 1], Basis::Y);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            measure_all_basis(&[0, 1], Basis::Y);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert!(hist.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_peek_all()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(1, 3, {
            h(0);
            peek_all(&[0]);
            h(0);
            peek_all(&[1]);
            h(0);
            peek_all(&[2]);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram_vec().unwrap();
        // Results of first and third measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 0.
        let n00 = hist[0] + hist[2] + hist[4] + hist[6];
        assert!(crate::stats::measurement_ok(n00, nr_shots, 0.5, tol));
        let n10 = hist[0] + hist[1] + hist[4] + hist[5];
        assert!(n10 == nr_shots);
        let n20 = hist[0] + hist[1] + hist[2] + hist[3];
        assert!(crate::stats::measurement_ok(n20, nr_shots, 0.5, tol));

        let mut circuit = circuit!(2, 6, {
            h(0);
            h(1);
            peek_all(&[0, 1]);
            h(0);
            peek_all(&[2, 3]);
            h(0);
            peek_all(&[4, 5]);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram().unwrap();
        // Results of first and third measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 0.
        let mut n0 = [0; 4];
        let mut n1 = [0; 4];
        let mut n2 = [0; 4];
        for (key, count) in hist
        {
            n0[key as usize & 0x03] += count;
            n1[(key as usize >> 2) & 0x03] += count;
            n2[(key as usize >> 4) & 0x03] += count;
        }
        assert!(n0.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
        assert_eq!(n1[1], 0);
        assert_eq!(n1[3], 0);
        assert!(crate::stats::measurement_ok(n1[0], nr_shots, 0.5, tol));
        assert!(crate::stats::measurement_ok(n1[2], nr_shots, 0.5, tol));
        assert!(n2.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_peek_all_basis()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(1, 3, {
            peek_all_basis(&[0], Basis::X);
            h(0);
            peek_all_basis(&[1], Basis::X);
            h(0);
            peek_all_basis(&[2], Basis::X);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram_vec().unwrap();
        // Results of first and third measurement should be approximately equally
        // distributed over 0 and 1, second should be pure 0.
        let n00 = hist[0] + hist[2] + hist[4] + hist[6];
        assert!(crate::stats::measurement_ok(n00, nr_shots, 0.5, tol));
        let n10 = hist[0] + hist[1] + hist[4] + hist[5];
        assert!(n10 == nr_shots);
        let n20 = hist[0] + hist[1] + hist[2] + hist[3];
        assert!(crate::stats::measurement_ok(n20, nr_shots, 0.5, tol));

        let mut circuit = circuit!(2, 6, {
            h(0);
            h(1);
            peek_all_basis(&[0, 1], Basis::Y);
            s(1);
            peek_all_basis(&[2, 3], Basis::Y);
            s(0);
            peek_all_basis(&[4, 5], Basis::Y);
        }).unwrap();
        circuit.execute(1024);
        let hist = circuit.histogram().unwrap();
        // Results of first measurement should be approximately equally
        // distributed over 0 and 1 for both qubits, second should be pure 0
        // for second qubit, third pure |00〉.
        let mut n0 = [0; 4];
        let mut n1 = [0; 4];
        let mut n2 = [0; 4];
        for (key, count) in hist
        {
            n0[key as usize & 0x03] += count;
            n1[(key as usize >> 2) & 0x03] += count;
            n2[(key as usize >> 4) & 0x03] += count;
        }
        assert!(n0.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
        assert_eq!(n1[2], 0);
        assert_eq!(n1[3], 0);
        assert!(crate::stats::measurement_ok(n1[0], nr_shots, 0.5, tol));
        assert!(crate::stats::measurement_ok(n1[1], nr_shots, 0.5, tol));
        assert_eq!(n2, [nr_shots, 0, 0, 0]);
    }

    #[test]
    fn test_histogram()
    {
        let nr_shots = 4096;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            add_gate(H::new(), &[0]);
            add_gate(H::new(), &[1]);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);

        let hist = circuit.histogram().unwrap();
        // With this many shots, we expect all keys to be present
        let mut keys: Vec<&u64> = hist.keys().collect();
        keys.sort();
        assert_eq!(keys, vec![&0, &1, &2, &3]);

        assert_eq!(hist.values().sum::<usize>(), nr_shots);
        assert!(hist.values().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_histogram_vec()
    {
        let nr_shots = 4096;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            add_gate(H::new(), &[0]);
            add_gate(H::new(), &[1]);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);

        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist.iter().sum::<usize>(), nr_shots);
        assert!(hist.iter().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_histogram_string()
    {
        let nr_shots = 4096;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            add_gate(H::new(), &[0]);
            add_gate(H::new(), &[1]);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);

        let hist = circuit.histogram_string().unwrap();
        // With this many shots, we expect all keys to be present
        let mut keys: Vec<&String> = hist.keys().collect();
        keys.sort();
        assert_eq!(keys, vec!["00", "01", "10", "11"]);

        assert_eq!(hist.values().sum::<usize>(), nr_shots);
        assert!(hist.values().all(
            |&count| crate::stats::measurement_ok(count, nr_shots, 0.25, tol)
        ));
    }

    #[test]
    fn test_reset()
    {
        let nr_shots = 1024;
        let tol = 1.0e-5;

        let mut circuit = circuit!(2, 2, {
            h(0);
            z(0);
            reset(0);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![nr_shots, 0, 0, 0]);

        let mut circuit = circuit!(2, 2, {
            h(0);
            z(0);
            x(1);
            reset(0);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist, vec![0, 0, nr_shots, 0]);

        let mut circuit = circuit!(2, 2, {
            h(0);
            z(0);
            h(1);
            reset(0);
            measure(0, 0);
            measure(1, 1);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert!(crate::stats::measurement_ok(hist[0], nr_shots, 0.5, tol));
        assert_eq!(hist[1], 0);
        assert!(crate::stats::measurement_ok(hist[2], nr_shots, 0.5, tol));
        assert_eq!(hist[3], 0);
    }

    #[test]
    fn test_reset_all()
    {
        let nr_shots = 1024;

        let mut circuit = circuit!(5, 5, {
            h(0);
            z(0);
            x(4);
            h(3);
            reset_all();
            measure_all(&[0, 1, 2, 3, 4]);
        }).unwrap();
        circuit.execute(nr_shots);
        let hist = circuit.histogram_vec().unwrap();
        assert_eq!(hist[0], nr_shots);
        assert!(hist[1..].iter().all(|&c| c == 0));
    }

    #[test]
    fn test_open_qasm()
    {
        let circuit = circuit!(2, 2, {
            x(0);
            cx(0, 1);
            barrier(&[0, 1]);
            cx(1, 0);
            barrier(&[1]);
            cx(0, 1);
            barrier(&[1, 0]);
            measure_x(0, 0);
            measure_y(1, 1);
        }).unwrap();
        assert_eq!(circuit.open_qasm(), Ok(String::from(
r#"OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg b[2];
x q[0];
cx q[0], q[1];
barrier q;
cx q[1], q[0];
barrier q[1];
cx q[0], q[1];
barrier q[1], q[0];
h q[0];
measure q[0] -> b[0];
sdg q[1];
h q[1];
measure q[1] -> b[1];
"#)));

        let circuit = circuit!(2, 2, {
            x(0);
            measure_all(&[0, 1]);
            measure_all(&[1, 0]);
            measure_all_basis(&[0, 1], Basis::X);
            measure_all_basis(&[0, 1], Basis::Y);
        }).unwrap();
        assert_eq!(circuit.open_qasm(), Ok(String::from(
r#"OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg b[2];
x q[0];
measure q -> b;
measure q[0] -> b[1];
measure q[1] -> b[0];
h q;
measure q -> b;
sdg q;
h q;
measure q -> b;
"#)));

        let circuit = circuit!(2, 0, {
            x(0);
            h(1);
            reset(0);
            x(0);
            reset_all();
        }).unwrap();
        assert_eq!(circuit.open_qasm(), Ok(String::from(
r#"OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
x q[0];
h q[1];
reset q[0];
x q[0];
reset q;
"#)));

        let circuit = circuit!(2, 2, {
            x(0);
            measure_all(&[0, 1]);
            add_conditional_gate(&[0, 1], 1, X::new(), &[0]);
            add_conditional_gate(&[], 1, X::new(), &[1]);
        }).unwrap();
        assert_eq!(circuit.open_qasm(), Ok(String::from(
r#"OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg b[2];
x q[0];
measure q -> b;
if (b == 1) x q[0];
x q[1];
"#)));

        let circuit = circuit!(2, 2, {
            add_conditional_gate(&[0], 1, X::new(), &[0]);
        }).unwrap();
        assert!(matches!(circuit.open_qasm(), Err(_)));
    }

    #[test]
    fn test_c_qasm()
    {
        let circuit = circuit!(3, 3, {
            x(0);
            cx(0, 1);
            cx(1, 0);
            cx(0, 1);
            measure(0, 0);
            measure_x(1, 1);
            measure_y(2, 2);
        }).unwrap();
        assert_eq!(circuit.c_qasm(), Ok(String::from(
r#"version 1.0
qubits 3
x q[0]
cnot q[0], q[1]
cnot q[1], q[0]
cnot q[0], q[1]
measure q[0]
measure_x q[1]
measure_y q[2]
"#)));

        let circuit = circuit!(2, 2, {
            x(0);
            h(1);
            measure_all(&[0, 1]);
            reset_all();
            measure_all_basis(&[0, 1], Basis::X);
            reset(1);
            measure_all_basis(&[0, 1], Basis::Y);
        }).unwrap();
        assert_eq!(circuit.c_qasm(), Ok(String::from(
r#"version 1.0
qubits 2
x q[0]
h q[1]
measure_all
prep_z q[0]
prep_z q[1]
h q[0]
h q[1]
measure_all
prep_z q[1]
sdag q[0]
h q[0]
sdag q[1]
h q[1]
measure_all
"#)));

        let circuit = circuit!(2, 2, {
            x(0);
            measure_all(&[0, 1]);
            add_conditional_gate(&[0, 1], 1, X::new(), &[0]);
            add_conditional_gate(&[], 1, X::new(), &[1]);
        }).unwrap();
        assert_eq!(circuit.c_qasm(), Ok(String::from(
r#"version 1.0
qubits 2
x q[0]
measure_all
not b[1]
c-x b[0], b[1], q[0]
not b[1]
x q[1]
"#)));

        let circuit = circuit!(2, 2, {
            measure(0, 1);
        }).unwrap();
        // c-Qasm only allows for measuring to the classical bit with the same index
        assert!(matches!(circuit.c_qasm(), Err(_)));
    }

    #[test]
    fn test_latex()
    {
        let circuit = circuit!(2, 2, {
            h(0);
            x(1);
            measure(0, 0);
            measure_x(1, 1);
            add_conditional_gate(&[0, 1], 2, X::new(), &[0]);
            reset_all();
            measure_all_basis(&[1, 0], Basis::Y);
            reset(0);
            measure_y(1, 0);
            barrier(&[1]);
        }).unwrap();

        assert_eq!(circuit.latex(), Ok(String::from(
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{H} & \meter & \qw & \targ & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \meterB{Y} & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \qw & \qw & \qw \\
    \lstick{\ket{0}} & \gate{X} & \qw & \meterB{X} & \qw & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \qw & \meterB{Y} & \meterB{Y} & \qw \barrier{0} & \qw \\
    \lstick{0} & \cw & \cw \cwx[-2] & \cw & \cctrlo{-2} & \cw & \cw & \cw \cwx[-1] & \cw \cwx[-1] & \cw & \cw \\
    \lstick{0} & \cw & \cw & \cw \cwx[-2] & \cctrl{-1} & \cw & \cw \cwx[-3] & \cw & \cw & \cw & \cw \\
}
"#)));
    }
}
