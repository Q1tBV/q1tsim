use cmatrix;
use gates;

/// Operation in a Custom gate.
struct CustomOp
{
    /// The gate
    gate: Box<gates::Gate>,
    /// The bits on which the gate acts
    bits: Vec<usize>
}

impl CustomOp
{
    /// Create a new Custom gate operation
    fn new<G>(gate: G, bits: &[usize]) -> Self
    where G: 'static + gates::Gate
    {
        CustomOp
        {
            gate: Box::new(gate),
            bits: bits.to_owned()
        }
    }
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
    /// using the `add_gate()`, `add_gate()` and `add_gate()` functions.
    pub fn new(name: &str, nr_bits: usize) -> Self
    {
        Custom
        {
            name: name.to_owned(),
            nr_bits: nr_bits,
            ops: vec![]
        }
    }

    /// Add a gate.
    ///
    /// Append a `n`-ary gate `gate`, operating on the `n` qubits in `bits`, to
    /// this gate.
    pub fn add_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: gates::Gate
    {
        self.ops.push(CustomOp::new(gate, bits));
    }
}

impl gates::Gate for Custom
{
    fn cost(&self) -> f64
    {
        self.ops.iter().map(|op| op.gate.cost()).sum()
    }

    fn description(&self) -> &str
    {
        &self.name
    }

    fn nr_affected_bits(&self) -> usize
    {
        self.nr_bits
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let mut res = cmatrix::CMatrix::eye(1 << self.nr_bits);
        for op in self.ops.iter()
        {
            res = op.gate.expanded_matrix(&op.bits, self.nr_bits).dot(&res);
        }

        res
    }
}

#[cfg(test)]
mod tests
{
    use cmatrix;
    use gates::{Gate, Custom, CCX, CX, H, X};

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
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(H::new(), &[1]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z,  z],
            [z, o, z,  z],
            [z, z, o,  z],
            [z, z, z, -o]
        ]);

        let mut gate = Custom::new("Inc", 2);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(X::new(), &[1]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [z, z, z, o],
            [o, z, z, z],
            [z, o, z, z],
            [z, z, o, z]
        ]);

        let mut gate = Custom::new("Inc", 3);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(CCX::new(), &[0, 1, 2]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(CX::new(), &[1, 2]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(X::new(), &[2]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [z, z, z, z, z, z, z, o],
            [o, z, z, z, z, z, z, z],
            [z, o, z, z, z, z, z, z],
            [z, z, o, z, z, z, z, z],
            [z, z, z, o, z, z, z, z],
            [z, z, z, z, o, z, z, z],
            [z, z, z, z, z, o, z, z],
            [z, z, z, z, z, z, o, z]
        ]);
    }
}
