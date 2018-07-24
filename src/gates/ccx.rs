extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// The CC<sub>X</sub> gate.
///
/// The CC<sub>X</sub> or Toffoli gate flips a qubit based on two control bits:
/// when either control bit is zero, it leaves the target unchanged; when both
/// control bits are one, it flips the |0〉 and |1〉 components of the target bit.
pub struct CCX
{
}

impl CCX
{
    /// Create a new CC<sub>X</sub> gate.
    pub fn new() -> Self
    {
        CCX { }
    }
}

impl gates::Gate for CCX
{
    fn description(&self) -> &str
    {
        "CCX"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        cmatrix::CMatrix::new(8, 8, vec![
            o, z, z, z, z, z, z, z,
            z, o, z, z, z, z, z, z,
            z, z, o, z, z, z, z, z,
            z, z, z, o, z, z, z, z,
            z, z, z, z, o, z, z, z,
            z, z, z, z, z, o, z, z,
            z, z, z, z, z, z, z, o,
            z, z, z, z, z, z, o, z
        ])
    }
}

impl gates::NaryGate for CCX
{
    fn nr_affected_bits(&self) -> usize
    {
        3
    }

    fn apply_n_ary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 8 == 0, "Number of rows is not a multiple of eight.");

        let n = state.rows() / 8;
        for i in 6*n..7*n
        {
            state.swap_rows(i, i+n);
        }
    }
}


#[cfg(test)]
mod tests
{
    use gates::{Gate, NaryGate};
    use gates::CCX;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let h = CCX::new();
        assert_eq!(h.description(), "CCX");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let cx = CCX::new();
        assert_complex_matrix_eq!(cx.matrix().as_ref(), matrix![
            o, z, z, z, z, z, z, z;
            z, o, z, z, z, z, z, z;
            z, z, o, z, z, z, z, z;
            z, z, z, o, z, z, z, z;
            z, z, z, z, o, z, z, z;
            z, z, z, z, z, o, z, z;
            z, z, z, z, z, z, z, o;
            z, z, z, z, z, z, o, z
        ]);
    }

    #[test]
    fn test_nr_affected_bits()
    {
        assert_eq!(CCX::new().nr_affected_bits(), 3);
    }

    #[test]
    fn test_apply_n_ary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let h = o * 0.5;

        let mut state = cmatrix::CMatrix::new(8, 4, vec![
            o, z,  h,  z,
            z, z, -h,  z,
            z, o,  h,  x,
            z, z, -h, -x,
            o, z,  h,  z,
            z, z, -h,  z,
            z, o,  h,  x,
            z, z, -h, -x
        ]);

        CCX::new().apply_n_ary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            o, z,  h,  z;
            z, z, -h,  z;
            z, o,  h,  x;
            z, z, -h, -x;
            o, z,  h,  z;
            z, z, -h,  z;
            z, z, -h, -x;
            z, o,  h,  x
        ]);
    }
}
