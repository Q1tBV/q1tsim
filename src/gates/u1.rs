extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// Phase gate.
///
/// The `U`<sub>`1`</sub>`(λ)` gate shifts the phase of the |1〉 component of a
/// qubit over an angle `λ`.
pub struct U1
{
    lambda: f64,
    desc: String
}

impl U1
{
    /// Create a new `U`<sub>`1`</sub> gate.
    pub fn new(lambda: f64) -> Self
    {
        U1 { lambda: lambda, desc: format!("U1({:.4})", lambda) }
    }
}

impl gates::Gate for U1
{
    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        cmatrix::CMatrix::new(2, 2, vec![o, z, z, num_complex::Complex::from_polar(&1.0, &self.lambda)])
    }

    fn apply_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let m = state.cols();
        let mut slice = state.sub_slice_mut([n, 0], n, m);
        slice *= num_complex::Complex::from_polar(&1.0, &self.lambda);
    }
}

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::U1;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        assert_eq!(gate.description(), "U1(0.7854)");
    }

    #[test]
    fn test_matrix()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_2);
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![o, z; z, i]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, x, x, z, o, x, -x]);

        let gate = U1::new(::std::f64::consts::FRAC_PI_4);

        gate.apply(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![o, z, x, x; z, x*(o+i), 0.5*(o+i), -0.5*(o+i)]);
    }
}
