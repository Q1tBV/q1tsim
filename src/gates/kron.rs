extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// Gate describing the Kronecker product of two other gates operating on
/// different qubits.
pub struct Kron<G0, G1>
{
    g0: G0,
    g1: G1,
    desc: String
}

impl<G0, G1> Kron<G0, G1>
where G0: gates::Gate, G1: gates::Gate
{
    /// Create a new Kronecker product gate `g1` ⊗ `g2`.
    pub fn new(g0: G0, g1: G1) -> Self
    {
        let desc = format!("{}⊗{}", g0.description(), g1.description());
        Kron { g0: g0, g1: g1, desc: desc }
    }
}

impl<G0, G1> gates::Gate for Kron<G0, G1>
where G0: gates::Gate, G1: gates::Gate
{
    fn description(&self) -> &str
    {
        &self.desc
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        self.g0.matrix().kron(&self.g1.matrix())
    }
}

impl<G0, G1> gates::BinaryGate for Kron<G0, G1>
where G0: gates::UnaryGate, G1: gates::UnaryGate
{
    fn apply_binary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 4 == 0, "Number of rows is not a multiple of four.");

        let n = state.rows() / 2;
        let m = state.cols();

        self.g0.apply_unary_slice(state);
        self.g1.apply_unary_slice(&mut state.sub_slice_mut([0, 0], n, m));
        self.g1.apply_unary_slice(&mut state.sub_slice_mut([n, 0], n, m));
    }
}

#[cfg(test)]
mod tests
{
    use gates::{Gate, BinaryGate};
    use gates::H;
    use gates::Identity;
    use gates::Kron;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let ih = Kron::new(Identity::new(), H::new());
        assert_eq!(ih.description(), "I⊗H");
        let hh = Kron::new(H::new(), H::new());
        assert_eq!(hh.description(), "H⊗H");
        let hih = Kron::new(H::new(), Kron::new(Identity::new(), H::new()));
        assert_eq!(hih.description(), "H⊗I⊗H");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let s = cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * cmatrix::COMPLEX_ONE;

        let ih = Kron::new(Identity::new(), H::new());
        assert_complex_matrix_eq!(ih.matrix().as_ref(), matrix![
            s,  s, z,  z;
            s, -s, z,  z;
            z,  z, s,  s;
            z,  z, s, -s
        ]);

        let hh = Kron::new(H::new(), H::new());
        assert_complex_matrix_eq!(hh.matrix().as_ref(), matrix![
            h,  h,  h,  h;
            h, -h,  h, -h;
            h,  h, -h, -h;
            h, -h, -h,  h
        ]);

        let hih = Kron::new(H::new(), Kron::new(Identity::new(), H::new()));
        assert_complex_matrix_eq!(hih.matrix().as_ref(), matrix![
            h,  h,  z,  z,  h,  h,  z,  z;
            h, -h,  z,  z,  h, -h,  z,  z;
            z,  z,  h,  h,  z,  z,  h,  h;
            z,  z,  h, -h,  z,  z,  h, -h;
            h,  h,  z,  z, -h, -h,  z,  z;
            h, -h,  z,  z, -h,  h,  z,  z;
            z,  z,  h,  h,  z,  z, -h, -h;
            z,  z,  h, -h,  z,  z, -h,  h
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
        let ih = Kron::new(Identity::new(), H::new());
        ih.apply_binary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            x, z, z, z;
            x, z, x, z;
            z, x, z, z;
            z, x, x, o
        ]);

        let mut state = cmatrix::CMatrix::new(4, 4, vec![
            o, z,  h,  z,
            z, z, -h,  z,
            z, o,  h,  x,
            z, z, -h, -x
        ]);
        let hi = Kron::new(H::new(), Identity::new());
        hi.apply_binary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            x,  x,  x,  h;
            z,  z, -x, -h;
            x, -x,  z, -h;
            z,  z,  z,  h
        ]);

        let mut state = cmatrix::CMatrix::new(4, 4, vec![
            o, z,  h,  z,
            z, z, -h,  z,
            z, o,  h,  x,
            z, z, -h, -x
        ]);
        let hh = Kron::new(H::new(), H::new());
        hh.apply_binary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            h,  h, z,  z;
            h,  h, o,  x;
            h, -h, z,  z;
            h, -h, z, -x
        ]);
    }
}
