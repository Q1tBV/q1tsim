use cmatrix;
use gates;

/// The Hadamard gate.
///
/// The Hadamard gate maps the zero state |0&rang; to the symmetric combination
/// of |0&rang; and |1&rang;, and the |1&rang; state to the anti-symmetric
/// combination.
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
        "H"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let s = ::std::f64::consts::FRAC_1_SQRT_2;
        cmatrix::CMatrix::new_real(2, 2, vec![s, s, s, -s])
    }
}

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::Hadamard;

    #[test]
    fn test_description()
    {
        let h = Hadamard::new();
        assert_eq!(h.description(), "H");
    }

    #[test]
    fn test_matrix()
    {
        let h = Hadamard::new();
        let s = ::std::f64::consts::FRAC_1_SQRT_2;
        assert_eq!(h.matrix().real().data(), &vec![s, s, s, -s]);
        assert_eq!(h.matrix().imag().data(), &vec![0.0; 4]);
    }
}
