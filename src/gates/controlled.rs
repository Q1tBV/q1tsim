extern crate num_complex;

pub use cmatrix;
use gates;

/// Controlled gates.
///
/// A controlled gate is a gate which acts upon a control bit: when
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
    fn cost(&self) -> f64
    {
        // Wild guess, probably wildly wrong
        2.0 * self.gate.cost()
    }

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

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        let n = state.len() / 2;
        self.gate.apply_mat_slice(&mut state.slice_mut(s![n.., ..]));
    }
}

#[macro_export]
macro_rules! declare_controlled_type
{
    ($(#[$attr:meta])* $name:ident, $gate_type:ty) => {
        $(#[$attr])*
        pub struct $name(gates::C<$gate_type>);
    }
}

#[macro_export]
macro_rules! declare_controlled_impl
{
    ($name:ident, $gate_type:ty $(, $arg:ident)*) => {
        impl $name
        {
            pub fn new($($arg: f64, )*) -> Self
            {
                $name(gates::C::new(<$gate_type>::new($($arg, )*)))
            }
        }
    };
    ($name:ident, $gate_type:ty, cost=$cost:expr $(, $arg:ident)*) => {
        impl $name
        {
            pub fn new($($arg: f64, )*) -> Self
            {
                $name(gates::C::new(<$gate_type>::new($($arg, )*)))
            }
            pub fn cost() -> f64
            {
                $cost
            }
        }
    };
}

#[macro_export]
macro_rules! declare_controlled_impl_gate
{
    ($name:ident, $gate_type:ty) => {
        impl gates::Gate for $name
        {
            fn cost(&self) -> f64 { self.0.cost() }
            fn description(&self) -> &str { self.0.description() }
            fn nr_affected_bits(&self) -> usize { self.0.nr_affected_bits() }
            fn matrix(&self) -> $crate::cmatrix::CMatrix { self.0.matrix() }
            fn apply_slice(&self, state: &mut $crate::cmatrix::CVecSliceMut)
            {
                self.0.apply_slice(state);
            }
        }
    };
    ($name:ident, $gate_type:ty, cost=$cost:expr) => {
        impl gates::Gate for $name
        {
            fn cost(&self) -> f64 { $cost }
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

#[macro_export]
macro_rules! declare_controlled
{
    ($(#[$attr:meta])* $name:ident, $gate_type:ty) => {
        declare_controlled_type!($(#[$attr])* $name, $gate_type);
        declare_controlled_impl!($name, $gate_type);
        declare_controlled_impl_gate!($name, $gate_type);
    };
    ($(#[$attr:meta])* $name:ident, $gate_type:ty, cost=$cost:expr $(, $arg:ident)*) => {
        declare_controlled_type!($(#[$attr])* $name, $gate_type);
        declare_controlled_impl!($name, $gate_type, cost=$cost $(, $arg)*);
        declare_controlled_impl_gate!($name, $gate_type, cost=Self::cost());
    };
}

declare_controlled!(
    /// Controlled Hadamard gate.
    CH, gates::H,
    cost=2.0*CX::cost() + 5.0*gates::U1::cost() + 3.0*gates::U2::cost() + gates::U3::cost());

declare_controlled!(
    /// Controlled `R`<sub>`X`</sub> gate.
    CRX, gates::RX,
    cost=2.0*CX::cost() + gates::U1::cost() + 2.0*gates::U3::cost(),
    theta);
declare_controlled!(
    /// Controlled `R`<sub>`Y`</sub> gate.
    CRY, gates::RY,
    cost=2.0*CX::cost() + gates::U1::cost() + 2.0*gates::U3::cost(),
    theta);
declare_controlled!(
    /// Controlled `R`<sub>`Z`</sub> gate.
    CRZ, gates::RZ,
    cost=2.0*CX::cost() + 2.0*gates::U1::cost(),
    lambda);

declare_controlled!(
    /// Controlled `S` gate.
    CS, gates::S, cost=2.0*CX::cost() + 3.0*gates::U1::cost());
declare_controlled!(
    /// Controlled `S`<sup>`†`</sup> gate.
    CSdg, gates::Sdg, cost=2.0*CX::cost() + 3.0*gates::U1::cost());

declare_controlled!(
    /// Controlled `U`<sub>`1`</sub> gate.
    CU1, gates::U1,
    cost=2.0*CX::cost() + 3.0*gates::U1::cost(),
    lambda);
declare_controlled!(
    /// Controlled `U`<sub>`2`</sub> gate.
    CU2, gates::U2,
    cost=2.0*CX::cost() + 2.0*gates::U1::cost() + gates::U2::cost(),
    phi, lambda);
declare_controlled!(
    /// Controlled `U`<sub>`3`</sub> gate.
    CU3, gates::U3,
    cost=2.0*CX::cost() + gates::U1::cost() + 2.0*gates::U3::cost(),
    theta, phi, lambda);

declare_controlled!(
    /// Controlled `V` gate.
    CV, gates::V,
    cost=2.0*CX::cost() + gates::U1::cost() + 2.0*gates::U3::cost());
declare_controlled!(
    /// Controlled `V`<sup>`†`</sup> gate.
    CVdg, gates::Vdg,
    cost=2.0*CX::cost() + gates::U1::cost() + 2.0*gates::U3::cost());

declare_controlled!(
    /// Controlled `X` gate.
    CX, gates::X, cost=1001.0);
declare_controlled!(
    /// Controlled `Y` gate.
    CY, gates::Y, cost=CX::cost() + 2.0*gates::U1::cost());
declare_controlled!(
    /// Controlled `Z` gate.
    CZ, gates::Z, cost=CX::cost() + 2.0*gates::U2::cost());

declare_controlled!(
    /// Doubly controlled `X` gate.
    CCX, gates::CX,
    cost=6.0*CX::cost() + 7.0*gates::U1::cost() + 2.0*gates::U2::cost());
declare_controlled!(
    /// Doubly controlled `Z` gate.
    CCZ, gates::CZ);

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, H, X};
    use super::{C, CCX, CCZ, CH, CRX, CRY, CRZ, CX, CY, CZ};
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
        gate_test(CX::new(), &mut state, &result);

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
        gate_test(CH::new(), &mut state, &result);
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
        gate_test(CCX::new(), &mut state, &result);
    }

    #[test]
    fn test_cost()
    {
        assert_eq!(CX::new().cost(), 1001.0);
        assert_eq!(CY::new().cost(), 1015.0);
        assert_eq!(CZ::new().cost(), 1209.0);
        assert_eq!(CH::new().cost(), 2550.0);
        assert_eq!(CRX::cost(), 2411.0);
        assert_eq!(CRY::cost(), 2411.0);
        assert_eq!(CRZ::cost(), 2016.0);
        assert_eq!(CCX::new().cost(), 6263.0);
        assert_eq!(CCZ::new().cost(), 2418.0);
    }
}
