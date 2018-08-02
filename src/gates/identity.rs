extern crate num_complex;

use cmatrix;
use gates;

/// The identity gate
///
/// The identity gate leaves the qubits on which it acts unchanged.
pub struct I
{
}

impl I
{
    /// Create a new identity gate.
    pub fn new() -> Self
    {
        I { }
    }
}

impl gates::Gate for I
{
    fn cost(&self) -> f64
    {
        0.0
    }

    fn description(&self) -> &str
    {
        "I"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::CMatrix::eye(2)
    }

    fn apply_slice(&self, _state: &mut cmatrix::CVecSliceMut)
    {
        // Identity, leave state unchanged, so do nothing
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, I};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let i = I::new();
        assert_eq!(i.description(), "I");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = I::new();
        assert_complex_matrix_eq!(i.matrix(), array![[o, z], [z, o]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [z, o, x, -x]];
        gate_test(I::new(), &mut state, &result);
    }

}
