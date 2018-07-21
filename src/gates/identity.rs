extern crate num_complex;
extern crate rulinalg;

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
    fn apply_unary<T>(&self, _s: &mut T)
    where T: rulinalg::matrix::BaseMatrixMut<num_complex::Complex64>
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
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let i = Identity::new();
        assert_eq!(i.description(), "I");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = Identity::new();
        assert_complex_matrix_eq!(i.matrix().as_ref(), matrix![o, z; z, o]);
    }

    #[test]
    fn test_apply_unary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, x, x, z, o, x, -x]);

        Identity::new().apply_unary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![o, z, x, x; z, o, x, -x]);
    }

}
