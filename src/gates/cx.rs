extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// The C<sub>X</sub> gate.
///
/// The C<sub>X</sub> or CNOT gate flips a qubit based on a control bit: when
/// the control bit is zero, it leaves the target unchanged; when the control
/// bit is one, it flips the |0〉 and |1〉 components of the target bit.
pub struct CX
{
}

impl CX
{
    /// Create a new C<sub>X</sub> gate.
    pub fn new() -> Self
    {
        CX { }
    }
}

impl gates::Gate for CX
{
    fn description(&self) -> &str
    {
        "CX"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        cmatrix::CMatrix::new(4, 4, vec![
            o, z, z, z,
            z, o, z, z,
            z, z, z, o,
            z, z, o, z
        ])
    }
}

impl gates::BinaryGate for CX
{
    fn apply_binary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 4 == 0, "Number of rows is not a multiple of four.");

        let n = state.rows() / 4;
        for i in 2*n..3*n
        {
            state.swap_rows(i, i+n);
        }
    }
}


#[cfg(test)]
mod tests
{
    use gates::{Gate, BinaryGate};
    use gates::CX;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let h = CX::new();
        assert_eq!(h.description(), "CX");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let cx = CX::new();
        assert_complex_matrix_eq!(cx.matrix().as_ref(), matrix![
            o, z, z, z;
            z, o, z, z;
            z, z, z, o;
            z, z, o, z
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

        CX::new().apply_binary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            o, z,  h,  z;
            z, z, -h,  z;
            z, z, -h, -x;
            z, o,  h,  x
        ]);
    }
}
