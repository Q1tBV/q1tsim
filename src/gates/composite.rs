use cmatrix;
use gates;

/// Operation in a composite gate.
struct SubGate
{
    /// The gate
    gate: Box<gates::Gate>,
    /// The bits on which the gate acts
    bits: Vec<usize>
}

impl SubGate
{
    /// Create a new composite gate operation
    fn new<G>(gate: G, bits: &[usize]) -> Self
    where G: 'static + gates::Gate
    {
        SubGate
        {
            gate: Box::new(gate),
            bits: bits.to_owned()
        }
    }
}

/// Compsite gate.
///
/// Struct Composite provides for user-defined gates that are made out of a
/// sequence of more primitive gates.
pub struct Composite
{
    // The name of the gate
    name: String,
    // The number of gates on which this gate operates
    nr_bits: usize,
    // The operations making up the gate
    ops: Vec<SubGate>
}

impl Composite
{
    /// Create a new composite gate.
    ///
    /// Initialize a new composite gate with name `name` for operating on `nr_bits`
    /// qubits at a time. The gates making up the operation should be added
    /// using the `add_gate()` function.
    pub fn new(name: &str, nr_bits: usize) -> Self
    {
        Composite
        {
            name: name.to_owned(),
            nr_bits: nr_bits,
            ops: vec![]
        }
    }

    /// Add a gate.
    ///
    /// Append a `n`-ary subgate `gate`, operating on the `n` qubits in `bits`,
    /// to this composite gate.
    pub fn add_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: gates::Gate
    {
        self.ops.push(SubGate::new(gate, bits));
    }
}

impl gates::Gate for Composite
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
    use super::Composite;
    use gates::{Gate, CCX, CX, H, X};

    #[test]
    fn test_description()
    {
        let gate = Composite::new("G", 3);
        assert_eq!(gate.description(), "G");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let mut gate = Composite::new("CZ", 2);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(H::new(), &[1]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z,  z],
            [z, o, z,  z],
            [z, z, o,  z],
            [z, z, z, -o]
        ]);

        let mut gate = Composite::new("Inc", 2);
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

        let mut gate = Composite::new("Inc", 3);
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
