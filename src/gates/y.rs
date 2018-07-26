extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

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
        cmatrix::CMatrix::new(2, 2, vec![z, -i, i, z])
    }

    fn apply_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let m = state.cols();
        for i in 0..n
        {
            state.swap_rows(i, i+n);
        }
        {
            let mut slice = state.sub_slice_mut([0, 0], n, m);
            slice *= -cmatrix::COMPLEX_I;
        }
        {
            let mut slice = state.sub_slice_mut([n, 0], n, m);
            slice *=  cmatrix::COMPLEX_I;
        }
    }
}

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::Y;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

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
        assert_complex_matrix_eq!(y.matrix().as_ref(), matrix![z, -i; i, z]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, x, x, z, o, x, -x]);

        Y::new().apply(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![z, -i, -i*x, i*x; i, z, i*x, i*x]);
    }
}
