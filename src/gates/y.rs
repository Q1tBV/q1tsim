extern crate num_complex;

use cmatrix;
use gates;
use qasm;

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

impl qasm::OpenQasm for Y
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("y {}", bit_names[bits[0]])
    }
}

impl qasm::CQasm for Y
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("y {}", bit_names[bits[0]])
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, Y};
    use qasm::{OpenQasm, CQasm};
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

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Y::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "y qb");
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Y::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "y qb");
    }
}
