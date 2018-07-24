extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

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
        let s = cmatrix::COMPLEX_HSQRT2;
        cmatrix::CMatrix::new(2, 2, vec![s, s, s, -s])
    }
}

impl gates::UnaryGate for Hadamard
{
    fn apply_unary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let m = state.cols();
        let s0 = state.sub_slice([0, 0], n, m).into_matrix() * cmatrix::COMPLEX_HSQRT2;
        let s1 = state.sub_slice([n, 0], n, m).into_matrix() * cmatrix::COMPLEX_HSQRT2;

        state.sub_slice_mut([0, 0], n, m).set_to(&s0 + &s1);
        state.sub_slice_mut([n, 0], n, m).set_to(&s0 - &s1);
    }
}


#[cfg(test)]
mod tests
{
    use gates::{Gate, UnaryGate};
    use gates::Hadamard;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

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
        let s = cmatrix::COMPLEX_HSQRT2;
        assert_complex_matrix_eq!(h.matrix().as_ref(), matrix![s, s; s, -s]);
    }

    #[test]
    fn test_apply_unary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, x, x, z, o, x, -x]);

        Hadamard::new().apply_unary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![x, x, o, z; x, -x, z, o]);
    }
}
