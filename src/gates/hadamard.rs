use cmatrix;
use gates;

/// The Hadamard gate.
///
/// The Hadamard gate maps the zero state |0> to the symmetric combination of
/// |0> and |1>, and the |1> state to the anti-symmetric combination.
pub struct Hadamard
{
}

impl Hadamard
{
    /// Create a new Hadamard gate.
    pub fn new() -> Self
    {
        Hadamard { }
    }
}

impl gates::Gate for Hadamard
{
    fn description(&self) -> &str
    {
        return "H";
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let s = ::std::f64::consts::FRAC_1_SQRT_2;
        cmatrix::CMatrix::new_real(matrix![s, s; s, -s])
    }
}
