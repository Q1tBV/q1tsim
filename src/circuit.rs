extern crate ndarray;

use gates;
use qustate;

use gates::Gate;

/// A single operation in a circuit
enum CircuitOp
{
    /// Apply a gate to the state
    Gate(Box<gates::Gate>, Vec<usize>),
    /// Conditionally apply a gate, depending on classical bits
    ConditionalGate(Vec<usize>, u64, Box<gates::Gate>, Vec<usize>),
    /// Measure a qubit in the Pauli X basis
    MeasureX(usize, usize),
    /// Measure a qubit in the Pauli Y basis
    MeasureY(usize, usize),
    /// Measure a qubit in the Pauli Z basis
    MeasureZ(usize, usize)
}

/// A quantum circuit
///
/// Struct Circuit represents a quantum circuit, holding a quantum state and the
/// operations to be performed on it.
pub struct Circuit
{
    /// The quantum state of the system
    pub q_state: qustate::QuState,
    /// The classial state of the system
    pub c_state: ndarray::Array2<u8>,
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
    where G: gates::Gate
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
    where G: gates::Gate
    {
        self.ops.push(CircuitOp::ConditionalGate(control.to_owned(), target,
            Box::new(gate), bits.to_owned()));
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `X` basis, into classical
    /// bit `cbit` to this circuit.
    pub fn measure_x(&mut self, qbit: usize, cbit: usize)
    {
        self.ops.push(CircuitOp::MeasureX(qbit, cbit));
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Y` basis, into classical
    /// bit `cbit` to this circuit.
    pub fn measure_y(&mut self, qbit: usize, cbit: usize)
    {
        self.ops.push(CircuitOp::MeasureY(qbit, cbit));
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` in the Pauli `Z` basis, into classical
    /// bit `cbit` to this circuit.
    pub fn measure_z(&mut self, qbit: usize, cbit: usize)
    {
        self.ops.push(CircuitOp::MeasureZ(qbit, cbit));
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` into classical bit `cbit` to this circuit.
    /// This is an alias for `measure_z()`.
    #[inline(always)]
    pub fn measure(&mut self, qbit: usize, cbit: usize)
    {
        self.measure_z(qbit, cbit);
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
        self.add_gate(gates::U3::new(phi, theta, lambda), &[bit]);
    }

    /// Add a C<sub>X</sub> gate.
    ///
    /// Add a `C`<sub>`X`</sub> gate, controlled by qubit `control` and
    /// operating on qubit `target`, to this circuit.
    pub fn cx(&mut self, control: usize, target: usize)
    {
        self.add_gate(gates::CX::new(), &[control, target]);
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
                    for &i in control
                    {
                        for (cb, &sb) in cbits.iter_mut().zip(self.c_state.row(i).iter())
                        {
                            *cb = (*cb << 1) | sb as u64;
                        }
                    }
                    let apply_gate: Vec<bool> = cbits.iter().map(|&b| b == target).collect();
                    self.q_state.apply_conditional_gate(&apply_gate, &**gate,
                        bits.as_slice());
                },
                CircuitOp::MeasureX(qbit, cbit)     => {
                    self.q_state.apply_gate(&gates::H::new(), &[qbit]);
                    let msr = self.q_state.measure(qbit);
                    self.c_state.row_mut(cbit).assign(&msr);
                }
                CircuitOp::MeasureY(qbit, cbit)     => {
                    self.q_state.apply_gate(&gates::Sdg::new(), &[qbit]);
                    self.q_state.apply_gate(&gates::H::new(), &[qbit]);
                    let msr = self.q_state.measure(qbit);
                    self.c_state.row_mut(cbit).assign(&msr);
                }
                CircuitOp::MeasureZ(qbit, cbit)     => {
                    let msr = self.q_state.measure(qbit);
                    self.c_state.row_mut(cbit).assign(&msr);
                }
            }
        }
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a single `u64` integer value. The
    /// most significant bit `n-1` in the histogram key corresponds to the first
    /// bit in the classical register, and the least significant bit in the key
    /// to the last bit in the register. This function of course only works
    /// when there are at most 64 bits in the register. If there are more, use
    /// `histogram_string()`.
    pub fn histogram(&self) -> ::std::collections::HashMap<u64, usize>
    {
        let mut res = ::std::collections::HashMap::new();
        for col in self.c_state.gencolumns()
        {
            let key = col.iter().fold(0, |k, &b| (k << 1) | b as u64);
            let count = res.entry(key).or_insert(0);
            *count += 1;
        }
        res
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a single `usize` integer value,
    /// which is used as an index in a vector. The vector is of length
    /// `2`<sub>`n`</sub>, so use this function only for reasonably small
    /// numbers of `n`. For sparse collections, using `histogram()` or
    /// `histogram_string` may be better.
    pub fn histogram_vec(&self) -> Vec<usize>
    {
        let mut res = vec![0; 1 << self.c_state.rows()];
        for col in self.c_state.gencolumns()
        {
            let key = col.iter().fold(0, |k, &b| (k << 1) | b as usize);
            res[key] += 1;
        }
        res
    }

    /// Create a histogram of measurements.
    ///
    /// Create a histogram of the measured classical bits. The `n` bits in the
    /// classical register are collected in a string key, with the first character
    /// in the key corresponding to the first bit in the classical register.
    pub fn histogram_string(&self) -> ::std::collections::HashMap<String, usize>
    {
        let mut res = ::std::collections::HashMap::new();
        for col in self.c_state.gencolumns()
        {
            let key = col.iter().map(|&b| ::std::char::from_digit(b as u32, 10).unwrap()).collect();
            let count = res.entry(key).or_insert(0);
            *count += 1;
        }
        res
    }

    pub fn open_qasm(&self) -> String
    {
        let mut res = String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n");

        let mut qbit_names = vec![];
        let nr_qbits = self.q_state.nr_bits();
        if nr_qbits > 0
        {
            res += &format!("qreg q[{}];\n", nr_qbits);
            for i in 0..nr_qbits
            {
                qbit_names.push(format!("q[{}]", i));
            }
        }
        let mut cbit_names = vec![];
        let nr_cbits = self.c_state.rows();
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
                    let instrs = gate.open_qasm(&qbit_names, bits);
                    for instr in instrs.split(';')
                    {
                    // XXX
//                         let s = instr.strip();
//                         res += &format!("if ({} == 1) {};\n", cbit_names[control], s);
                    }
                },
                CircuitOp::MeasureX(qbit, cbit)   => {
                    res += &format!("{};\n", gates::H::new().open_qasm(&qbit_names, &[qbit]));
                    res += &format!("measure q[{}] -> b[{}];\n", qbit, cbit);
                }
                CircuitOp::MeasureY(qbit, cbit)   => {
                    res += &format!("{};\n", gates::Sdg::new().open_qasm(&qbit_names, &[qbit]));
                    res += &format!("{};\n", gates::H::new().open_qasm(&qbit_names, &[qbit]));
                    res += &format!("measure q[{}] -> b[{}];\n", qbit, cbit);
                }
                CircuitOp::MeasureZ(qbit, cbit)   => {
                    res += &format!("measure q[{}] -> b[{}];\n", qbit, cbit);
                }
            }
        }

        res
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

    fn check_c_asm_binary_controlled_gate<G>(gate: &G) -> Result<(), String>
    where G: gates::Gate
    {
        if gate.nr_affected_bits() != 1
        {
            Err(String::from("Binary control can only be used with single-bit quantum operations"))
        }
        else
        {
            Ok(())
        }
    }

    pub fn c_qasm(&self) -> Result<String, String>
    {
        let mut res = String::from("version 1.0\n");

        let mut qbit_names = vec![];
        let mut cbit_names = vec![];
        let nr_qbits = self.q_state.nr_bits();
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
                // XXX
//                     Self::check_c_asm_binary_controlled_gate(gate)?;
//                     let instrs = gate.open_qasm(&qbit_names, bits);
//                     for instr in instrs.split('\n')
//                     {
//                         res += &format!("c-{}, {};\n", instr.strip(), cbit_names[control]);
//                     }
                },
                CircuitOp::MeasureX(qbit, cbit)     => {
                    Self::check_c_qasm_measurement(qbit, cbit)?;
                    res += &format!("measure_x q[{}]\n", qbit);
                }
                CircuitOp::MeasureY(qbit, cbit)     => {
                    Self::check_c_qasm_measurement(qbit, cbit)?;
                    res += &format!("measure_y q[{}]\n", qbit);
                }
                CircuitOp::MeasureZ(qbit, cbit)     => {
                    Self::check_c_qasm_measurement(qbit, cbit)?;
                    res += &format!("measure q[{}]\n", qbit);
                }
            }
        }

        Ok(res)
    }
}

#[cfg(test)]
mod tests
{
    use circuit::Circuit;
    use gates;

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
        assert_eq!(hist, vec![0, 0, 1024, 0]);

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
        assert_eq!(hist, vec![0, 0, 1024, 0]);

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.x(0);
        circuit.measure_y(0, 0);
        circuit.measure_y(1, 1);
        circuit.execute();
        let hist = circuit.histogram_vec();
        assert!(*hist.iter().min().unwrap() >= min_count);
    }

    #[test]
    fn test_open_qasm()
    {
        let mut circuit = Circuit::new(2, 2, 10);
        circuit.x(0);
        circuit.cx(0, 1);
        circuit.cx(1, 0);
        circuit.cx(0, 1);

        circuit.measure(0, 0);
        circuit.measure(1, 1);

        let asm = circuit.open_qasm();
        assert_eq!(asm, concat!(
            "OPENQASM 2.0;\n",
            "include \"qelib1.inc\";\n",
            "qreg q[2];\n",
            "creg b[2];\n",
            "x q[0];\n",
            "cx q[0], q[1];\n",
            "cx q[1], q[0];\n",
            "cx q[0], q[1];\n",
            "measure q[0] -> b[0];\n",
            "measure q[1] -> b[1];\n"
        ));
    }

    #[test]
    fn test_c_qasm()
    {
        let mut circuit = Circuit::new(2, 2, 10);
        circuit.x(0);
        circuit.cx(0, 1);
        circuit.cx(1, 0);
        circuit.cx(0, 1);

        circuit.measure(0, 0);
        circuit.measure(1, 1);

        let asm = circuit.c_qasm();
        assert_eq!(asm, Ok(String::from(concat!(
            "version 1.0\n",
            "qubits 2\n",
            "x q[0]\n",
            "cnot q[0], q[1]\n",
            "cnot q[1], q[0]\n",
            "cnot q[0], q[1]\n",
            "measure q[0]\n",
            "measure q[1]\n"
        ))));
    }
}
