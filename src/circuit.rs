extern crate rulinalg;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

use gates;
use qustate;

/// A single operation in a circuit
enum CircuitOp
{
    /// Apply a unary gate to a single qubit
    Unary(Box<gates::UnaryGate>, usize),
    /// Apply a binary gate on two qubits
    Binary(Box<gates::BinaryGate>, usize, usize),
    /// Apply a multi-bit gate
    Nary(Box<gates::NaryGate>, Vec<usize>),
    /// Measure a qubit
    Measure(usize, usize)
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
    c_state: rulinalg::matrix::Matrix<u8>,
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
            c_state: rulinalg::matrix::Matrix::zeros(nr_cbits, nr_shots),
            ops: vec![]
        }
    }

    /// The classical register.
    ///
    /// Return a reference to the classical bit register, containing the results
    /// of any measurements made on the system.
    pub fn cstate(&self) -> &rulinalg::matrix::Matrix<u8>
    {
        &self.c_state
    }

    /// Add a unary gate.
    ///
    /// Append a unary gate `gate` operating on qubit `bit` to this circuit.
    pub fn add_unary_gate<G: 'static>(&mut self, gate: G, bit: usize)
    where G: gates::UnaryGate
    {
        self.ops.push(CircuitOp::Unary(Box::new(gate), bit));
    }

    /// Add a binary gate.
    ///
    /// Append a binary gate `gate` operating on qubits `bit0` and `bit1` to
    /// this circuit.
    pub fn add_binary_gate<G: 'static>(&mut self, gate: G, bit0: usize, bit1: usize)
    where G: gates::BinaryGate
    {
        self.ops.push(CircuitOp::Binary(Box::new(gate), bit0, bit1));
    }

    /// Add a multi-bit gate.
    ///
    /// Append a `n`-ary gate `gate`, operating on the `n` qubits in `bits`, to
    /// this circuit.
    pub fn add_n_ary_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: gates::NaryGate
    {
        self.ops.push(CircuitOp::Nary(Box::new(gate), bits.to_owned()));
    }

    /// Add a measurement.
    ///
    /// Add measurement of qubit `qbit` into classical bit `cbit` to this circuit.
    pub fn add_measurement(&mut self, qbit: usize, cbit: usize)
    {
        self.ops.push(CircuitOp::Measure(qbit, cbit));
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
                CircuitOp::Unary(ref gate, bit)         => {
                    self.q_state.apply_unary_gate(&**gate, bit);
                },
                CircuitOp::Binary(ref gate, bit0, bit1) => {
                    self.q_state.apply_binary_gate(&**gate, bit0, bit1);
                },
                CircuitOp::Nary(ref gate, ref bits) => {
                    self.q_state.apply_n_ary_gate(&**gate, bits.as_slice());
                },
                CircuitOp::Measure(qbit, cbit)          => {
                    let msr = self.q_state.measure(qbit);
                    self.c_state.row_mut(cbit).raw_slice_mut().copy_from_slice(msr.data());
                }
            }
        }
    }

    pub fn histogram(&self) -> ::std::collections::HashMap<String, usize>
    {
        let mut res = ::std::collections::HashMap::new();
        for col in self.c_state.col_iter()
        {
            let key = col.iter().map(|&b| ::std::char::from_digit(b as u32, 10).unwrap()).collect();
            let count = res.entry(key).or_insert(0);
            *count += 1;
        }
        res
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
        circuit.add_unary_gate(gates::X::new(), 0);
        circuit.add_unary_gate(gates::X::new(), 1);
        circuit.add_binary_gate(gates::CX::new(), 0, 1);
        circuit.add_measurement(0, 0);
        circuit.add_measurement(1, 1);
        circuit.execute();
        assert_matrix_eq!(circuit.cstate(), matrix![1, 1, 1, 1, 1; 0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_histogram()
    {
        let nr_shots = 4096;
        // chance of individual count being less than min_count is less than 10^-5
        // (assuming normal distribution)
        let min_count = 906;

        let mut circuit = Circuit::new(2, 2, nr_shots);
        circuit.add_unary_gate(gates::Hadamard::new(), 0);
        circuit.add_unary_gate(gates::Hadamard::new(), 1);
        circuit.add_measurement(0, 0);
        circuit.add_measurement(1, 1);
        circuit.execute();

        let hist = circuit.histogram();
        // With this many shots, we expect all keys to be present
        let mut keys: Vec<&String> = hist.keys().collect();
        keys.sort();
        assert_eq!(keys, vec!["00", "01", "10", "11"]);

        assert_eq!(hist.values().sum::<usize>(), nr_shots);
        assert!(*hist.values().min().unwrap() >= min_count);
    }
}