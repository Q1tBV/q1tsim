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


extern crate ndarray;

use export;
use gates;
use qustate;

use export::{CircuitGate, CQasm, OpenQasm};

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
    /// Prevent gate reordering on the associated bits across the barrier
    Barrier(Vec<usize>)
}

/// A quantum circuit
///
/// Struct Circuit represents a quantum circuit, holding a quantum state and the
/// operations to be performed on it.
pub struct Circuit
{
    /// The quantum state of the system
    q_state: qustate::QuState,
    /// The classial state of the system
    c_state: ndarray::Array2<u8>,
    /// The operations to perform on the state
    ops: Vec<CircuitOp>
}

impl Circuit
{
    /// Create a new circuit.
    ///
    /// Create a new (empty) quantum circuit, with `nr_qubits` quantum bits and
    /// `nr_cbits` classical bits, to be run `nr_shots` times.
    pub fn new(nr_qubits: usize, nr_cbits: usize, nr_shots: usize) -> Self
    {
        Circuit
        {
            q_state: qustate::QuState::new(nr_qubits, nr_shots),
            c_state: ndarray::Array::zeros((nr_cbits, nr_shots)),
            ops: vec![]
        }
    }

    /// The number of quantum bits in this circuit
    pub fn nr_qbits(&self) -> usize
    {
        self.q_state.nr_bits()
    }

    /// The number of classical bits in this circuit
    pub fn nr_cbits(&self) -> usize
    {
        self.c_state.rows()
    }

    /// The classical register.
    ///
    /// Return a reference to the classical bit register, containing the results
    /// of any measurements made on the system.
    pub fn cstate(&self) -> &ndarray::Array2<u8>
    {
        &self.c_state
    }

    /// Add a gate.
    ///
    /// Append a `n`-ary gate `gate`, operating on the `n` qubits in `bits`, to
    /// this circuit.
    pub fn add_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: CircuitGate
    {
        self.ops.push(CircuitOp::Gate(Box::new(gate), bits.to_owned()));
    }

    /// Add a conditional gate.
    ///
    /// Append a `n`-ary gate `gate`, that will operate on the `n` qubits in
    /// `bits` to this circuit. The gate will only be applied only when the
    /// classical bits with indices from `control` form the target word `target`.
    /// The bit at the position of the first index in `control` is interpreted
    /// as the most significant bit to check.
    pub fn add_conditional_gate<G: 'static>(&mut self, control: &[usize],
        target: u64, gate: G, bits: &[usize])
    where G: CircuitGate
    {
        self.ops.push(CircuitOp::ConditionalGate(control.to_owned(), target,
            Box::new(gate), bits.to_owned()));
    }

    /// Add a measurement
    ///
    /// Add measurement of qubit `qbit` in basis `basis`, into classical bit
    /// `cbit`, to this circuit.
    pub fn measure_basis(&mut self, qbit: usize, cbit: usize, basis: Basis)
    {
        self.ops.push(CircuitOp::Measure(qbit, cbit, basis));
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `X` basis, into classical
    /// bit `cbit` to this circuit.
    #[inline(always)]
    pub fn measure_x(&mut self, qbit: usize, cbit: usize)
    {
        self.measure_basis(qbit, cbit, Basis::X);
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Y` basis, into classical
    /// bit `cbit` to this circuit.
    #[inline(always)]
    pub fn measure_y(&mut self, qbit: usize, cbit: usize)
    {
        self.measure_basis(qbit, cbit, Basis::Y);
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Z` basis, into classical
    /// bit `cbit` to this circuit.
    #[inline(always)]
    pub fn measure_z(&mut self, qbit: usize, cbit: usize)
    {
        self.measure_basis(qbit, cbit, Basis::Z);
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` into classical bit `cbit` to this circuit.
    /// This is an alias for `measure_z()`.
    #[inline(always)]
    pub fn measure(&mut self, qbit: usize, cbit: usize)
    {
        self.measure_basis(qbit, cbit, Basis::Z);
    }

    /// Add a measurement.
    ///
    /// Add the measurement of all qubits in the quantum state into the classical
    /// bits `cbits`. Measurement is done in basis `basis`.
    pub fn measure_all_basis(&mut self, cbits: &[usize], basis: Basis)
    {
        self.ops.push(CircuitOp::MeasureAll(cbits.to_owned(), basis));
    }

    /// Add a measurement.
    ///
    /// Add the measurement of all qubits in the quantum state into the classical
    /// bits `cbits`. Measurement is done in the Pauli `Z` basis.
    #[inline(always)]
    pub fn measure_all(&mut self, cbits: &[usize])
    {
        self.measure_all_basis(cbits, Basis::Z);
    }

    /// Reset a qubit
    ///
    /// Reset the qubit `qbit` to |0〉. This is done by measuring the bit, and
    /// flipping it if the result is `1`, so this is potentially an expensive
    /// operation.
    pub fn reset(&mut self, qbit: usize)
    {
        self.ops.push(CircuitOp::Reset(qbit));
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
    /// Add a Hadamard operating on qubit `bit`, to this circuit.
    pub fn h(&mut self, bit: usize)
    {
        self.add_gate(gates::H::new(), &[bit]);
    }

    /// Add a Pauli X gate.
    ///
    /// Add a Pauli X gate operating on qubit `bit`, to this circuit.
    pub fn x(&mut self, bit: usize)
    {
        self.add_gate(gates::X::new(), &[bit]);
    }

    /// Add a Pauli Y gate.
    ///
    /// Add a Pauli Y gate operating on qubit `bit`, to this circuit.
    pub fn y(&mut self, bit: usize)
    {
        self.add_gate(gates::Y::new(), &[bit]);
    }

    /// Add a Pauli Z gate.
    ///
    /// Add a Pauli Z gate operating on qubit `bit`, to this circuit.
    pub fn z(&mut self, bit: usize)
    {
        self.add_gate(gates::Z::new(), &[bit]);
    }

    /// Add a R<sub>X</sub> gate.
    ///
    /// Add a `R`<sub>`X`</sub>`(θ)` gate operating on qubit `bit`, to this circuit.
    pub fn rx(&mut self, theta: f64, bit: usize)
    {
        self.add_gate(gates::RX::new(theta), &[bit]);
    }

    /// Add a R<sub>Y</sub> gate.
    ///
    /// Add a `R`<sub>`Y`</sub>`(θ)` gate operating on qubit `bit`, to this circuit.
    pub fn ry(&mut self, theta: f64, bit: usize)
    {
        self.add_gate(gates::RY::new(theta), &[bit]);
    }

    /// Add a R<sub>Z</sub> gate.
    ///
    /// Add a `R`<sub>`Z`</sub>`(λ)` gate operating on qubit `bit`, to this circuit.
    pub fn rz(&mut self, lambda: f64, bit: usize)
    {
        self.add_gate(gates::RZ::new(lambda), &[bit]);
    }

    /// Add a U<sub>1</sub> gate.
    ///
    /// Add a `U`<sub>`1`</sub>`(λ)` gate operating on qubit `bit`, to this circuit.
    pub fn u1(&mut self, lambda: f64, bit: usize)
    {
        self.add_gate(gates::U1::new(lambda), &[bit]);
    }

    /// Add a U<sub>2</sub> gate.
    ///
    /// Add a `U`<sub>`2`</sub>`(ϕ, λ)` gate operating on qubit `bit`, to this circuit.
    pub fn u2(&mut self, phi: f64, lambda: f64, bit: usize)
    {
        self.add_gate(gates::U2::new(phi, lambda), &[bit]);
    }

    /// Add a U<sub>3</sub> gate.
    ///
    /// Add a `U`<sub>`3`</sub>`(θ, ϕ, λ)` gate operating on qubit `bit`, to this circuit.
    pub fn u3(&mut self, theta: f64, phi: f64, lambda: f64, bit: usize)
    {
        self.add_gate(gates::U3::new(theta, phi, lambda), &[bit]);
    }

    /// Add a C<sub>X</sub> gate.
    ///
    /// Add a `C`<sub>`X`</sub> gate, controlled by qubit `control` and
    /// operating on qubit `target`, to this circuit.
    pub fn cx(&mut self, control: usize, target: usize)
    {
        self.add_gate(gates::CX::new(), &[control, target]);
    }

    /// Add a barrier
    ///
    /// Add a barrier on the bits in `bits`. No transformations on these bits
    /// are allowed across this barrier.
    pub fn barrier(&mut self, bits: &[usize])
    {
        self.ops.push(CircuitOp::Barrier(bits.to_vec()));
    }

    /// Execute this circuit
    ///
    /// Execute this circuit, performing its operations and measurements. Note
    /// that this does not reset the state before execution. In case multiple
    /// runs of the same circuit are to be done, call `reset()` between
    /// executions.
    pub fn execute(&mut self)
    {
        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    self.q_state.apply_gate(&**gate, bits.as_slice());
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    let nr_shots = self.c_state.cols();
                    let mut cbits = vec![0; nr_shots];
                    for (ibit, &irow) in control.iter().enumerate()
                    {
                        for (cb, &sb) in cbits.iter_mut().zip(self.c_state.row(irow).iter())
                        {
                            *cb |= (sb << ibit) as u64;
                        }
                    }
                    let apply_gate: Vec<bool> = cbits.iter().map(|&b| b == target).collect();
                    self.q_state.apply_conditional_gate(&apply_gate, &**gate,
                        bits.as_slice());
                },
                CircuitOp::Measure(qbit, cbit, basis) => {
                    match basis
                    {
                        Basis::X => {
                            self.q_state.apply_gate(&gates::H::new(), &[qbit]);
                        },
                        Basis::Y => {
                            self.q_state.apply_gate(&gates::Sdg::new(), &[qbit]);
                            self.q_state.apply_gate(&gates::H::new(), &[qbit]);
                        },
                        _ => {
                            /* do nothing */
                        }
                    }
                    self.q_state.measure_into(qbit, self.c_state.row_mut(cbit));
                }
                CircuitOp::MeasureAll(ref cbits, basis) => {
                    match basis
                    {
                        Basis::X => {
                            self.q_state.apply_unary_gate_all(&gates::H::new());
                        },
                        Basis::Y => {
                            self.q_state.apply_unary_gate_all(&gates::Sdg::new());
                            self.q_state.apply_unary_gate_all(&gates::H::new());
                        },
                        _ => {
                            /* do nothing */
                        }
                    }
                    let msr = self.q_state.measure_all();
                    for (&i, row) in cbits.iter().zip(msr.genrows())
                    {
                        self.c_state.row_mut(i).assign(&row);
                    }
                },
                CircuitOp::Reset(bit) => {
                    self.q_state.reset(bit);
                },
                CircuitOp::ResetAll => {
                    self.q_state.reset_all();
                },
                CircuitOp::Barrier(_) => {
                    /* Nothing to be done */
                }
            }
        }
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
    pub fn histogram(&self) -> ::std::collections::HashMap<u64, usize>
    {
        let mut res = ::std::collections::HashMap::new();
        for col in self.c_state.gencolumns()
        {
            let key = col.iter().rev().fold(0, |k, &b| (k << 1) | b as u64);
            let count = res.entry(key).or_insert(0);
            *count += 1;
        }
        res
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
    pub fn histogram_vec(&self) -> Vec<usize>
    {
        let mut res = vec![0; 1 << self.c_state.rows()];
        for col in self.c_state.gencolumns()
        {
            let key = col.iter().rev().fold(0, |k, &b| (k << 1) | b as usize);
            res[key] += 1;
        }
        res
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a string key, with the last character
    /// in the key corresponding to the first bit (at index 0) in the classical
    /// register and vice versa.
    pub fn histogram_string(&self) -> ::std::collections::HashMap<String, usize>
    {
        let mut res = ::std::collections::HashMap::new();
        for col in self.c_state.gencolumns()
        {
            let key = col.iter().rev()
                .map(|&b| ::std::char::from_digit(b as u32, 10).unwrap())
                .collect();
            let count = res.entry(key).or_insert(0);
            *count += 1;
        }
        res
    }

    fn is_full_register(&self, control: &[usize]) -> bool
    {
        let n = control.len();
        if n != self.c_state.rows()
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

    fn check_open_qasm_condition_bits(&self, control: &[usize]) -> Result<(), String>
    {
        if !self.is_full_register(control)
        {
            Err(String::from("OpenQasm can only perform conditional operations based on a complete classical register"))
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
    pub fn open_qasm(&self) -> Result<String, String>
    {
        let mut res = String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n");

        let mut qbit_names = vec![];
        let nr_qbits = self.nr_qbits();
        if nr_qbits > 0
        {
            res += &format!("qreg q[{}];\n", nr_qbits);
            for i in 0..nr_qbits
            {
                qbit_names.push(format!("q[{}]", i));
            }
        }
        let mut cbit_names = vec![];
        let nr_cbits = self.nr_cbits();
        if nr_cbits > 0
        {
            res += &format!("creg b[{}];\n", nr_cbits);
            for i in 0..nr_qbits
            {
                cbit_names.push(format!("b[{}]", i));
            }
        }

        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    res += &format!("{};\n", gate.open_qasm(&qbit_names, bits));
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    if control.is_empty()
                    {
                        res += &format!("{};\n", gate.open_qasm(&qbit_names, bits));
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
                                gates::H::new().open_qasm(&qbit_names, &[qbit]));
                        },
                        Basis::Y => {
                            res += &format!("{};\n",
                                gates::Sdg::new().open_qasm(&qbit_names, &[qbit]));
                            res += &format!("{};\n",
                                gates::H::new().open_qasm(&qbit_names, &[qbit]));
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
                                gates::H::new().open_qasm(&names, &[0]));
                        },
                        Basis::Y => {
                            let names = [String::from("q")];
                            res += &format!("{};\n",
                                gates::Sdg::new().open_qasm(&names, &[0]));
                            res += &format!("{};\n",
                                gates::H::new().open_qasm(&names, &[0]));
                        }
                        _ => {}
                    }

                    if cbits.len() == self.nr_cbits()
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
                CircuitOp::Reset(qbit) => {
                    res += &format!("reset {};\n", qbit_names[qbit]);
                },
                CircuitOp::ResetAll => {
                    res += "reset q;\n";
                },
                CircuitOp::Barrier(ref qbits) => {
                    if qbits.len() == self.nr_qbits()
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

    fn check_c_qasm_measurement(qbit: usize, cbit: usize) -> Result<(), String>
    {
        if qbit != cbit
        {
            Err(String::from("In cQasm, no classical registers can be specified. Measurements must be made to a classical bit with the same index as the qubit"))
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
    pub fn c_qasm(&self) -> Result<String, String>
    {
        let mut res = String::from("version 1.0\n");

        let mut qbit_names = vec![];
        let mut cbit_names = vec![];
        let nr_qbits = self.nr_qbits();
        if nr_qbits > 0
        {
            res += &format!("qubits {}\n", nr_qbits);
            for i in 0..nr_qbits
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
                    res += &format!("{}\n", gate.c_qasm(&qbit_names, bits));
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    if control.is_empty()
                    {
                        res += &format!("{}\n", gate.c_qasm(&qbit_names, bits));
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
                        let gate_qasm = gate.conditional_c_qasm(&condition, &qbit_names, bits)?;
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
                            for bit in 0..self.nr_qbits()
                            {
                                res += &format!("{}\n",
                                    gates::H::new().c_qasm(&qbit_names, &[bit]));
                            }
                        },
                        Basis::Y => {
                            for bit in 0..self.nr_qbits()
                            {
                                res += &format!("{}\n",
                                    gates::Sdg::new().c_qasm(&qbit_names, &[bit]));
                                res += &format!("{}\n",
                                    gates::H::new().c_qasm(&qbit_names, &[bit]));
                            }
                        },
                        _ => {
                            /* do nothing */
                        }
                    }
                    res += &format!("measure_all\n");
                },
                CircuitOp::Reset(qbit) => {
                    res += &format!("prep_z {}\n", qbit_names[qbit]);
                },
                CircuitOp::ResetAll => {
                    for i in 0..nr_qbits
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

    pub fn latex(&self) -> String
    {
        let nr_qbits = self.nr_qbits();
        let nr_cbits = self.nr_cbits();
        let mut state = export::LatexExportState::new(nr_qbits, nr_cbits);
        for op in self.ops.iter()
        {
            match *op
            {
                CircuitOp::Gate(ref gate, ref bits) => {
                    gate.latex_checked(bits, &mut state);
                },
                CircuitOp::ConditionalGate(ref control, target, ref gate, ref bits) => {
                    state.reserve_range(bits, Some(control));
                    let controlled = state.set_controlled(true);
                    gate.latex(bits, &mut state);
                    state.set_controlled(controlled);
                    state.set_condition(control, target, bits);
                },
                CircuitOp::Measure(qbit, cbit, basis) => {
                    let basis_lbl = match basis
                    {
                        Basis::X => Some("X"),
                        Basis::Y => Some("Y"),
                        _        => None
                    };
                    state.set_measurement(qbit, cbit, basis_lbl);
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
                        state.set_measurement(qbit, cbit, basis_lbl);
                    }
                },
                CircuitOp::Reset(qbit) => {
                    state.set_reset(qbit);
                },
                CircuitOp::ResetAll => {
                    state.reserve_range(&[0, nr_qbits-1], None);
                    for qbit in 0..nr_qbits
                    {
                        state.set_reset(qbit);
                    }
                },
                CircuitOp::Barrier(ref qbits) => {
                    state.set_barrier(qbits);
                }
            }
        }

        state.code()
    }
}

#[cfg(test)]
mod tests
{
    use super::{Basis, Circuit, CircuitOp};
    use cmatrix;
    use gates;

    #[test]
    fn test_gate_methods()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;

        let mut circuit = Circuit::new(2, 0, 1);
        circuit.h(0);
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

        circuit.x(1);
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

        circuit.y(0);
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

        circuit.z(1);
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

        circuit.rx(::std::f64::consts::PI, 1);
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

        circuit.ry(::std::f64::consts::PI, 0);
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

        circuit.rz(::std::f64::consts::PI, 1);
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

        circuit.u1(::std::f64::consts::FRAC_PI_4, 1);
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

        circuit.u2(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::FRAC_PI_2, 0);
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

        circuit.u3(::std::f64::consts::PI, ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::FRAC_PI_2, 0);
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

        circuit.cx(1, 0);
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
        let mut circuit = Circuit::new(2, 2, 5);
        circuit.add_gate(gates::X::new(), &[0]);
        circuit.add_gate(gates::X::new(), &[1]);
        circuit.add_gate(gates::CX::new(), &[0, 1]);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();
        assert_eq!(circuit.cstate(), &array![[1, 1, 1, 1, 1], [0, 0, 0, 0, 0]]);
    }

    #[test]
    fn test_measure()
    {
        let nr_shots = 1024;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 196;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.measure_x(0, 0);
        circuit.measure_x(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert!(*hist.iter().min().unwrap() >= min_count);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.h(0);
        circuit.h(1);
        circuit.measure_x(0, 0);
        circuit.measure_x(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.measure_y(0, 0);
        circuit.measure_y(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert!(*hist.iter().min().unwrap() >= min_count);
    }

    #[test]
    fn test_conditional()
    {
        let mut circuit = Circuit::new(2, 2, 5);
        circuit.add_conditional_gate(&[0, 1], 1, gates::X::new(), &[1]);
        circuit.measure_all(&[0, 1]);
        circuit.execute();
        assert_eq!(circuit.c_state, array![
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0]
        ]);

        let mut circuit = Circuit::new(2, 2, 5);
        circuit.c_state.assign(&array![
            [1, 0, 0, 1, 0],
            [0, 1, 1, 1, 0]
        ]);
        circuit.add_conditional_gate(&[0, 1], 1, gates::X::new(), &[1]);
        circuit.measure_all(&[0, 1]);
        circuit.execute();
        assert_eq!(circuit.c_state, array![
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0]
        ]);

        let mut circuit = Circuit::new(2, 2, 5);
        circuit.c_state.assign(&array![
            [1, 0, 0, 1, 0],
            [0, 1, 1, 1, 0]
        ]);
        circuit.add_conditional_gate(&[0, 1], 2, gates::X::new(), &[1]);
        circuit.measure_all(&[0, 1]);
        circuit.execute();
        assert_eq!(circuit.c_state, array![
            [0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0]
        ]);

        let mut circuit = Circuit::new(2, 2, 5);
        circuit.c_state.assign(&array![
            [1, 0, 0, 1, 0],
            [0, 1, 1, 1, 0]
        ]);
        circuit.add_conditional_gate(&[1], 1, gates::X::new(), &[0]);
        circuit.measure_all(&[0, 1]);
        circuit.execute();
        assert_eq!(circuit.c_state, array![
            [0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0]
        ]);
    }

    #[test]
    fn test_measure_all()
    {
        let nr_shots = 1024;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 196;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.measure_all(&[0, 1]);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.measure_all(&[1, 0]);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, 0, nr_shots, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.h(0);
        circuit.h(1);
        circuit.measure_all(&[0, 1]);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert!(*hist.iter().min().unwrap() >= min_count);
    }

    #[test]
    fn test_measure_all_basis()
    {
        let nr_shots = 1024;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 196;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.h(0);
        circuit.h(1);
        circuit.measure_all_basis(&[0, 1], Basis::X);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![nr_shots, 0, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.h(0);
        circuit.h(1);
        circuit.measure_all_basis(&[0, 1], Basis::X);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.h(0);
        circuit.h(1);
        circuit.add_gate(gates::S::new(), &[0]);
        circuit.add_gate(gates::S::new(), &[1]);
        circuit.measure_all_basis(&[0, 1], Basis::Y);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, nr_shots, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.measure_all_basis(&[0, 1], Basis::Y);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert!(*hist.iter().min().unwrap() >= min_count);
    }

    #[test]
    fn test_histogram()
    {
        let nr_shots = 4096;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 906;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.add_gate(gates::H::new(), &[0]);
        circuit.add_gate(gates::H::new(), &[1]);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();

        let hist = circuit.histogram();
        // With this many shots, we expect all keys to be present
        let mut keys: Vec<&u64> = hist.keys().collect();
        keys.sort();
        assert_eq!(keys, vec![&0, &1, &2, &3]);

        assert_eq!(hist.values().sum::<usize>(), nr_shots);
        assert!(*hist.values().min().unwrap() >= min_count);
    }

    #[test]
    fn test_histogram_vec()
    {
        let nr_shots = 4096;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 906;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.add_gate(gates::H::new(), &[0]);
        circuit.add_gate(gates::H::new(), &[1]);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();

        let hist = circuit.histogram_vec();
        assert_eq!(hist.iter().sum::<usize>(), nr_shots);
        assert!(*hist.iter().min().unwrap() >= min_count);
    }

    #[test]
    fn test_histogram_string()
    {
        let nr_shots = 4096;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 906;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.add_gate(gates::H::new(), &[0]);
        circuit.add_gate(gates::H::new(), &[1]);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();

        let hist = circuit.histogram_string();
        // With this many shots, we expect all keys to be present
        let mut keys: Vec<&String> = hist.keys().collect();
        keys.sort();
        assert_eq!(keys, vec!["00", "01", "10", "11"]);

        assert_eq!(hist.values().sum::<usize>(), nr_shots);
        assert!(*hist.values().min().unwrap() >= min_count);
    }

    #[test]
    fn test_reset()
    {
        let nr_shots = 1024;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 443;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.h(0);
        circuit.z(0);
        circuit.reset(0);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![nr_shots, 0, 0, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.h(0);
        circuit.z(0);
        circuit.x(1);
        circuit.reset(0);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist, vec![0, 0, nr_shots, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.h(0);
        circuit.z(0);
        circuit.h(1);
        circuit.reset(0);
        circuit.measure(0, 0);
        circuit.measure(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert!(hist[0] > min_count);
        assert!(hist[2] > min_count);
        assert_eq!(hist[1], 0);
        assert_eq!(hist[3], 0);
    }

    #[test]
    fn test_reset_all()
    {
        let nr_shots = 1024;

        let mut circuit = Circuit::new(5, 5, nr_shots);
        circuit.h(0);
        circuit.z(0);
        circuit.x(4);
        circuit.h(3);
        circuit.reset_all();
        circuit.measure_all(&[0, 1, 2, 3, 4]);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert_eq!(hist[0], nr_shots);
        assert!(hist[1..].iter().all(|&c| c == 0));
    }

    #[test]
    fn test_open_qasm()
    {
        let mut circuit = Circuit::new(2, 2, 1);
        circuit.x(0);
        circuit.cx(0, 1);
        circuit.barrier(&[0, 1]);
        circuit.cx(1, 0);
        circuit.barrier(&[1]);
        circuit.cx(0, 1);
        circuit.barrier(&[1, 0]);
        circuit.measure_x(0, 0);
        circuit.measure_y(1, 1);
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

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.x(0);
        circuit.measure_all(&[0, 1]);
        circuit.measure_all(&[1, 0]);
        circuit.measure_all_basis(&[0, 1], Basis::X);
        circuit.measure_all_basis(&[0, 1], Basis::Y);
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

        let mut circuit = Circuit::new(2, 0, 1);
        circuit.x(0);
        circuit.h(1);
        circuit.reset(0);
        circuit.x(0);
        circuit.reset_all();
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

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.x(0);
        circuit.measure_all(&[0, 1]);
        circuit.add_conditional_gate(&[0, 1], 1, gates::X::new(), &[0]);
        circuit.add_conditional_gate(&[], 1, gates::X::new(), &[1]);
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

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.add_conditional_gate(&[0], 1, gates::X::new(), &[0]);
        assert!(matches!(circuit.open_qasm(), Err(_)));

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.add_conditional_gate(&[1, 2], 1, gates::X::new(), &[0]);
        assert!(matches!(circuit.open_qasm(), Err(_)));
    }

    #[test]
    fn test_c_qasm()
    {
        let mut circuit = Circuit::new(3, 3, 10);
        circuit.x(0);
        circuit.cx(0, 1);
        circuit.cx(1, 0);
        circuit.cx(0, 1);
        circuit.measure(0, 0);
        circuit.measure_x(1, 1);
        circuit.measure_y(2, 2);
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

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.x(0);
        circuit.h(1);
        circuit.measure_all(&[0, 1]);
        circuit.reset_all();
        circuit.measure_all_basis(&[0, 1], Basis::X);
        circuit.reset(1);
        circuit.measure_all_basis(&[0, 1], Basis::Y);
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

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.x(0);
        circuit.measure_all(&[0, 1]);
        circuit.add_conditional_gate(&[0, 1], 1, gates::X::new(), &[0]);
        circuit.add_conditional_gate(&[], 1, gates::X::new(), &[1]);
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

        let mut circuit = Circuit::new(2, 2, 1);
        circuit.measure(0, 1);
        // c-Qasm only allows for measuring to the classical bit with the same index
        assert!(matches!(circuit.c_qasm(), Err(_)));
    }

    #[test]
    fn test_latex()
    {
        let mut circuit = Circuit::new(2, 2, 1024);

        circuit.h(0);
        circuit.x(1);
        circuit.measure(0, 0);
        circuit.measure_x(1, 1);
        circuit.add_conditional_gate(&[0, 1], 2, gates::X::new(), &[0]);
        circuit.reset_all();
        circuit.measure_all_basis(&[1, 0], Basis::Y);
        circuit.reset(0);
        circuit.measure_y(1, 0);
        circuit.barrier(&[1]);

        assert_eq!(circuit.latex(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{H} & \meter & \qw & \targ & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \meterB{Y} & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \qw & \qw & \qw \\
    \lstick{\ket{0}} & \gate{X} & \qw & \meterB{X} & \qw & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \qw & \meterB{Y} & \meterB{Y} & \qw \barrier{0} & \qw \\
    \lstick{0} & \cw & \cw \cwx[-2] & \cw & \cctrlo{-2} & \cw & \cw & \cw \cwx[-1] & \cw \cwx[-1] & \cw & \cw \\
    \lstick{0} & \cw & \cw & \cw \cwx[-2] & \cctrl{-1} & \cw & \cw \cwx[-3] & \cw & \cw & \cw & \cw \\
}
"#);
    }
}
