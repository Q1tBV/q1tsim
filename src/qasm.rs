use gates;

use gates::Gate;

/// Trait for gates that can be represented in OpenQasm.
pub trait OpenQasm
{
    /// OpenQasm representation
    ///
    /// Return an OpenQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits.
    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String;
}

/// Trait for gates that can be represented in c-Qasm.
pub trait CQasm
{
    /// cQasm representation
    ///
    /// Return an cQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits.
    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String;
}


/// Trait combining the traits necessary for a gate in a quantum circuit
pub trait CircuitGate: gates::Gate + OpenQasm + CQasm {}

impl<G: Gate + OpenQasm + CQasm> CircuitGate for G {}
