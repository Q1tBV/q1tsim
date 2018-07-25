extern crate rulinalg;

use cmatrix;
use gates;

/// Operation in a Custom gate.
enum CustomOp
{
    /// Unary gate operating on a single qubit.
    Unary(Box<gates::UnaryGate>, usize),
    /// Binary gate operating on two qubits.
    Binary(Box<gates::BinaryGate>, usize, usize),
    /// Gate operating on multiple qubits.
    Nary(Box<gates::NaryGate>, Vec<usize>)
}

/// Custom gate.
///
/// Struct Custom provides for user-defined gates that are made out of a
/// sequence of more primitive gates.
pub struct Custom
{
    // The name of the gate
    name: String,
    // The number of gates on which this gate operates
    nr_bits: usize,
    // The operations making up the gate
    ops: Vec<CustomOp>
}

impl Custom
{
    /// Create a new custom gate.
    ///
    /// Initialize a new custom gate with name `name` for operating on `nr_bits`
    /// qubits at a time. The gates making up the operation should be added
    /// using the `add_unary_gate()`, `add_binary_gate()` and `add_n_ary_gate()` functions.
    pub fn new(name: &str, nr_bits: usize) -> Self
    {
        Custom
        {
            name: name.to_owned(),
            nr_bits: nr_bits,
            ops: vec![]
        }
    }

    /// Add a unary gate.
    ///
    /// Append a unary gate `gate` operating on qubit `bit` to this gate.
    pub fn add_unary_gate<G: 'static>(&mut self, gate: G, bit: usize)
    where G: gates::UnaryGate
    {
        assert!(bit < self.nr_bits, "Invalid bit index {} for {}-bit gate", bit, self.nr_bits);
        self.ops.push(CustomOp::Unary(Box::new(gate), bit));
    }

    /// Add a binary gate.
    ///
    /// Append a binary gate `gate` operating on qubits `bit0` and `bit1` to
    /// this gate.
    pub fn add_binary_gate<G: 'static>(&mut self, gate: G, bit0: usize, bit1: usize)
    where G: gates::BinaryGate
    {
        self.ops.push(CustomOp::Binary(Box::new(gate), bit0, bit1));
    }

    /// Add a multi-bit gate.
    ///
    /// Append a `n`-ary gate `gate`, operating on the `n` qubits in `bits`, to
    /// this gate.
    pub fn add_n_ary_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: gates::NaryGate
    {
        self.ops.push(CustomOp::Nary(Box::new(gate), bits.to_owned()));
    }
}

impl gates::Gate for Custom
{
    fn description(&self) -> &str
    {
        &self.name
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let mut res = cmatrix::CMatrix::eye_sq(1 << self.nr_bits);
        for op in self.ops.iter()
        {
            match *op
            {
                CustomOp::Unary(ref gate, bit)         => {
                    res = gate.expanded_matrix(bit, self.nr_bits) * res;
                },
                CustomOp::Binary(ref gate, bit0, bit1) => {
                    res = gate.expanded_matrix(bit0, bit1, self.nr_bits) * res;
                },
                CustomOp::Nary(ref gate, ref bits)     => {
                    res = gate.expanded_matrix(bits, self.nr_bits) * res;
                }
            }
        }

        res
    }
}

#[cfg(test)]
mod tests
{
    use cmatrix;
    use gates::Gate;
    use gates::{Custom, CCX, CX, Hadamard, X};
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let gate = Custom::new("G", 3);
        assert_eq!(gate.description(), "G");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let mut gate = Custom::new("CZ", 2);
        gate.add_unary_gate(Hadamard::new(), 1);
        gate.add_binary_gate(CX::new(), 0, 1);
        gate.add_unary_gate(Hadamard::new(), 1);
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![
            o, z, z,  z;
            z, o, z,  z;
            z, z, o,  z;
            z, z, z, -o
        ]);

        let mut gate = Custom::new("Inc", 2);
        gate.add_unary_gate(Hadamard::new(), 0);
        gate.add_unary_gate(Hadamard::new(), 1);
        gate.add_binary_gate(CX::new(), 0, 1);
        gate.add_unary_gate(Hadamard::new(), 1);
        gate.add_unary_gate(Hadamard::new(), 0);
        gate.add_unary_gate(X::new(), 1);
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![
            z, z, z, o;
            o, z, z, z;
            z, o, z, z;
            z, z, o, z
        ]);

        let mut gate = Custom::new("Inc", 3);
        gate.add_unary_gate(Hadamard::new(), 0);
        gate.add_unary_gate(Hadamard::new(), 2);
        gate.add_n_ary_gate(CCX::new(), &[0, 1, 2]);
        gate.add_unary_gate(Hadamard::new(), 1);
        gate.add_binary_gate(CX::new(), 1, 2);
        gate.add_unary_gate(Hadamard::new(), 2);
        gate.add_unary_gate(Hadamard::new(), 0);
        gate.add_unary_gate(Hadamard::new(), 1);
        gate.add_unary_gate(X::new(), 2);
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![
            z, z, z, z, z, z, z, o;
            o, z, z, z, z, z, z, z;
            z, o, z, z, z, z, z, z;
            z, z, o, z, z, z, z, z;
            z, z, z, o, z, z, z, z;
            z, z, z, z, o, z, z, z;
            z, z, z, z, z, o, z, z;
            z, z, z, z, z, z, o, z
        ]);
    }
}
