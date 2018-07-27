extern crate num_complex;

pub use cmatrix;
use gates;

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
        let gsize = gm.rows();

        let mut res = cmatrix::CMatrix::eye(2*gsize);
        res.slice_mut(s![gsize.., gsize..]).assign(&gm);

        res
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        let n = state.len() / 2;
        self.gate.apply_slice(&mut state.slice_mut(s![n..]));
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
            fn apply_slice(&self, state: &mut $crate::cmatrix::CVecSliceMut)
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
    use gates::{gate_test, Gate, C, CH, X, H};
    use cmatrix;

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
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z, z],
            [z, o, z, z],
            [z, z, z, o],
            [z, z, o, z]
        ]);

        let gate = C::new(H::new());
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z,  z],
            [z, o, z,  z],
            [z, z, x,  x],
            [z, z, x, -x]
        ]);
    }

    #[test]
    fn test_apply_binary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let h = o * 0.5;

        let mut state = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, o,  h,  x],
            [z, z, -h, -x]
        ];
        let result = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, z, -h, -x],
            [z, o,  h,  x]
        ];
        gate_test(C::new(X::new()), &mut state, &result);

        let mut state = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, o,  h,  x],
            [z, z, -h, -x]
        ];
        let result = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, x,  z,  z],
            [z, x,  x,  o]
        ];
        gate_test(C::new(H::new()), &mut state, &result);
    }

    #[test]
    fn test_apply_n_ary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = array![
            [o,  z,  z],
            [z,  z,  z],
            [z,  z,  z],
            [z,  x,  z],
            [z,  z,  z],
            [z,  z,  z],
            [z,  z, -x],
            [z, -x,  x]
        ];
        let result = array![
            [o,  z,  z],
            [z,  z,  z],
            [z,  z,  z],
            [z,  x,  z],
            [z,  z,  z],
            [z,  z,  z],
            [z, -x,  x],
            [z,  z, -x]
        ];
        gate_test(C::new(C::new(X::new())), &mut state, &result);
    }
}