extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// The C<sub>H</sub> gate.
///
/// The C<sub>H</sub> gate is a controlled Hadamard gate: when the control bit
// is zero, it leaves the target unchanged; when the control bit is one, it
// applies a Hadamard transform to the target bit.
pub struct CH
{
}

impl CH
{
    /// Create a new C<sub>H</sub> gate.
    pub fn new() -> Self
    {
        CH { }
    }
}

impl gates::Gate for CH
{
    fn description(&self) -> &str
    {
        "CH"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        cmatrix::CMatrix::new(4, 4, vec![
            o, z, z,  z,
            z, o, z,  z,
            z, z, x,  x,
            z, z, x, -x
        ])
    }
}

impl gates::BinaryGate for CH
{
    fn apply_binary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 4 == 0, "Number of rows is not a multiple of four.");

        let n = state.rows() / 2;
        let m = state.cols();
        gates::H::transform(&mut state.sub_slice_mut([n, 0], n, m));
    }
}


#[cfg(test)]
mod tests
{
    use gates::{Gate, BinaryGate};
    use gates::CH;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let gate = CH::new();
        assert_eq!(gate.description(), "CH");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let gate = CH::new();
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![
            o, z, z,  z;
            z, o, z,  z;
            z, z, x,  x;
            z, z, x, -x
        ]);
    }

    #[test]
    fn test_apply_binary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let h = o * 0.5;

        let mut state = cmatrix::CMatrix::new(4, 4, vec![
            o, z,  h,  z,
            z, z, -h,  z,
            z, o,  h,  x,
            z, z, -h, -x
        ]);

        CH::new().apply_binary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            o, z,  h,  z;
            z, z, -h,  z;
            z, x,  z,  z;
            z, x,  x,  o
        ]);
    }
}
