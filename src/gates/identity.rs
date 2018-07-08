use cmatrix;
use gates;

/// The identity gate
///
/// The identity gate leaves the qubits on which it acts unchanged.
pub struct Identity
{
}

impl Identity
{
    /// Create a new identity gate.
    pub fn new() -> Self
    {
        Identity { }
    }
}

impl gates::Gate for Identity
{
    fn description(&self) -> &str
    {
        return "I";
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::CMatrix::eye(2)
    }
}
