extern crate num_complex;

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
        let s = num_complex::Complex::new(::std::f64::consts::FRAC_1_SQRT_2, 0.0);
        cmatrix::CMatrix::new(2, 2, vec![s, s, s, -s])
    }
}

impl gates::UnaryGate for Hadamard
{
    fn apply_unary(&self, state: &mut cmatrix::CMatrix)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let m = state.cols();
        let s0 = state.sub_slice([0, 0], n, m).into_matrix();
        let s1 = state.sub_slice([n, 0], n, m).into_matrix();

        state.sub_slice_mut([0, 0], n, m).set_to(&s0 + &s1);
        state.sub_slice_mut([n, 0], n, m).set_to(&s0 - &s1);
        *state *= ::std::f64::consts::FRAC_1_SQRT_2;
    }
}


#[cfg(test)]
mod tests
{
    use gates::{Gate, UnaryGate};
    use gates::Hadamard;
    use cmatrix;

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
        let s = ::std::f64::consts::FRAC_1_SQRT_2;
        assert_matrix_eq!(h.matrix().real(), matrix![s, s; s, -s]);
        assert_matrix_eq!(h.matrix().imag(), matrix![0.0, 0.0; 0.0, 0.0]);
    }

    #[test]
    fn test_apply_unary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        let xc = cmatrix::COMPLEX_HSQRT2;
        let mut state = cmatrix::CMatrix::new(2, 4, vec![o, z, xc, xc, z, o, xc, -xc]);

        Hadamard::new().apply_unary(&mut state);
        assert_matrix_eq!(state.real(), matrix![x, x, 1.0, 0.0; x, -x, 0.0, 1.0], comp=float);
        assert_matrix_eq!(state.imag(), matrix![0.0, 0.0, 0.0, 0.0; 0.0, 0.0, 0.0, 0.0], comp=float);
    }
}
