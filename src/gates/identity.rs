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
        "I"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::CMatrix::eye(2, 2)
    }
}

impl gates::UnaryGate for Identity
{
    fn apply_unary(&self, _s: &mut cmatrix::CMatrix)
    {
        // Identity, leave state unchanged, so do nothing
    }
}

#[cfg(test)]
mod tests
{
    use gates::{Gate, UnaryGate};
    use gates::Identity;
    use cmatrix;

    #[test]
    fn test_description()
    {
        let i = Identity::new();
        assert_eq!(i.description(), "I");
    }

    #[test]
    fn test_matrix()
    {
        let i = Identity::new();
        assert_matrix_eq!(i.matrix().real(), matrix![1.0, 0.0; 0.0, 1.0], comp=float);
        assert_matrix_eq!(i.matrix().imag(), matrix![0.0, 0.0; 0.0, 0.0], comp=float);
    }

    #[test]
    fn test_apply_unary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        let xc = o*x;
        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, xc, xc, z, o, xc, -xc]);

        Identity::new().apply_unary(&mut state);
        assert_matrix_eq!(state.real(), matrix![1.0, 0.0, x, x; 0.0, 1.0, x, -x], comp=float);
        assert_matrix_eq!(state.imag(), matrix![0.0, 0.0, 0.0, 0.0; 0.0, 0.0, 0.0, 0.0], comp=float);
    }

}
