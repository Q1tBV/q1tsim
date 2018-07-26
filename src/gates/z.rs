extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// The Pauli Z gate.
///
/// The Z gate rotates the state over π radians around the `z` axis of
/// the Bloch sphere, i.e. it flips the sign of the |1〉 components of the qubit.
pub struct Z
{
}

impl Z
{
    /// Create a new Pauli Z gate.
    pub fn new() -> Self
    {
        Z { }
    }
}

impl gates::Gate for Z
{
    fn description(&self) -> &str
    {
        "Z"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        cmatrix::CMatrix::new(2, 2, vec![o, z, z, -o])
    }

    fn apply_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let m = state.cols();
        state.sub_slice_mut([n, 0], n, m).apply(&|c: num_complex::Complex64| -c);
    }
}

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::Z;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let z = Z::new();
        assert_eq!(z.description(), "Z");
    }

    #[test]
    fn test_matrix()
    {
        let gate = Z::new();
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![o, z; z, -o]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, x, x, z, o, x, -x]);

        Z::new().apply(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![o, z, x, x; z, -o, -x, x]);
    }
}
