extern crate num_complex;
extern crate rulinalg;

pub use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// Controlled gates.
///
/// A controlled gate is a gate which acts upon a control but: when
/// the control bit is zero, it leaves the target unchanged; when the control
/// bit is one, the gate is applied.
pub struct C<G>
where G: gates::Gate
{
    gate: G,
    desc: String
}

impl<G> C<G>
where G: gates::Gate
{
    /// Create a new controlled gate for `gate`.
    pub fn new(gate: G) -> Self
    {
        let desc = format!("C{}", gate.description());
        C { gate: gate, desc: desc }
    }
}

impl<G> gates::Gate for C<G>
where G: gates::Gate
{
    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        1 + self.gate.nr_affected_bits()
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let gm = self.gate.matrix();
        let gsize = gm.as_ref().rows();

        let mut res = cmatrix::CMatrix::eye_sq(2*gsize);
        res.as_mut().sub_slice_mut([gsize, gsize], gsize, gsize).set_to(gm.as_matrix());

        res
    }

    fn apply_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        let n = state.rows() / 2;
        let m = state.cols();
        self.gate.apply_slice(&mut state.sub_slice_mut([n, 0], n, m));
    }
}

#[macro_export]
macro_rules! declare_controlled
{
    ($name:ident, $gate_type:ty) => {
        pub struct $name(gates::C<$gate_type>);
        impl $name
        {
            pub fn new() -> Self
            {
                type T = $gate_type;
                $name(gates::C::new(T::new()))
            }
        }
        impl gates::Gate for $name
        {
            fn description(&self) -> &str { self.0.description() }
            fn nr_affected_bits(&self) -> usize { self.0.nr_affected_bits() }
            fn matrix(&self) -> $crate::cmatrix::CMatrix { self.0.matrix() }
            fn apply_slice(&self,
                state: &mut $crate::cmatrix::RLMatrixSliceMut)
            {
                self.0.apply_slice(state);
            }
        }
    }
}

declare_controlled!(CX, gates::X);
declare_controlled!(CZ, gates::Z);
declare_controlled!(CH, gates::H);
declare_controlled!(CCX, gates::CX);
declare_controlled!(CCZ, gates::CZ);

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::{C, CH, X, H};
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let gate = C::new(X::new());
        assert_eq!(gate.description(), "CX");
        let gate = CH::new();
        assert_eq!(gate.description(), "CH");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let gate = C::new(X::new());
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![
            o, z, z, z;
            z, o, z, z;
            z, z, z, o;
            z, z, o, z
        ]);

        let gate = C::new(H::new());
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
        C::new(X::new()).apply(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            o, z,  h,  z;
            z, z, -h,  z;
            z, z, -h, -x;
            z, o,  h,  x
        ]);

        let mut state = cmatrix::CMatrix::new(4, 4, vec![
            o, z,  h,  z,
            z, z, -h,  z,
            z, o,  h,  x,
            z, z, -h, -x
        ]);
        C::new(H::new()).apply(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            o, z,  h,  z;
            z, z, -h,  z;
            z, x,  z,  z;
            z, x,  x,  o
        ]);
    }

    #[test]
    fn test_apply_n_ary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = cmatrix::CMatrix::from_matrix(matrix![
            o,  z,  z ;
            z,  z,  z ;
            z,  z,  z ;
            z,  x,  z ;
            z,  z,  z ;
            z,  z,  z ;
            z,  z, -x ;
            z, -x,  x
        ]);
        C::new(C::new(X::new())).apply(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
            o,  z,  z ;
            z,  z,  z ;
            z,  z,  z ;
            z,  x,  z ;
            z,  z,  z ;
            z,  z,  z ;
            z, -x,  x ;
            z,  z, -x
        ]);
    }
}
