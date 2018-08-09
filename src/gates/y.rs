extern crate num_complex;

use cmatrix;
use gates;

/// The Pauli Y gate.
///
/// The Y gate rotates the state over π radians around the `y` axis of the Bloch
/// sphere.
pub struct Y
{
}

impl Y
{
    /// Create a new Pauli Y gate.
    pub fn new() -> Self
    {
        Y { }
    }
}

impl gates::Gate for Y
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "Y"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let i = cmatrix::COMPLEX_I;
        array![[z, -i], [i, z]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        gates::X::transform(state);
        let n = state.len() / 2;
        {
            let mut slice = state.slice_mut(s![..n]);
            slice *= -cmatrix::COMPLEX_I;
        }
        {
            let mut slice = state.slice_mut(s![n..]);
            slice *=  cmatrix::COMPLEX_I;
        }
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        gates::X::transform_mat(state);
        let n = state.rows() / 2;
        {
            let mut slice = state.slice_mut(s![..n, ..]);
            slice *= -cmatrix::COMPLEX_I;
        }
        {
            let mut slice = state.slice_mut(s![n.., ..]);
            slice *=  cmatrix::COMPLEX_I;
        }
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, Y};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let y = Y::new();
        assert_eq!(y.description(), "Y");
    }

    #[test]
    fn test_matrix()
    {
        let y = Y::new();
        let z = cmatrix::COMPLEX_ZERO;
        let i = cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(y.matrix(), array![[z, -i], [i, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[z, -i, -i*x, i*x], [i, z, i*x, i*x]];
        gate_test(Y::new(), &mut state, &result);
    }
}
