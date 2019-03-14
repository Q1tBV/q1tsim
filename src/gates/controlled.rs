// Copyright 2019 Q1t BV
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


extern crate num_complex;

pub use cmatrix;
use gates;
use error;
use export;

use gates::Gate;

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

impl<G> export::Latex for C<G>
where G: gates::Gate + export::Latex
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;

        // We can only reasonably draw this if the controlled gate
        // operates on bits all above or all below the control bit
        let control = bits[0];
        let min = *bits[1..].iter().min().unwrap();
        let max = *bits[1..].iter().max().unwrap();
        if min > control && max > control
        {
            state.set_field(control, format!(r"\ctrl{{{}}}", min - control));
        }
        else if min < control && max < control
        {
            state.set_field(control, format!(r"\ctrl{{{}}}", max as isize - control as isize));
        }
        else
        {
            // XXX FIXME: figure out how to do this, and if really not possible,
            // return an error.
            panic!("Unable to draw controlled gate with control in the middle");
        }

        let controlled = state.set_controlled(true);
        self.gate.latex(&bits[1..], state)?;
        state.set_controlled(controlled);

        Ok(())
    }

    fn latex_checked(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        state.reserve_range(bits, None)?;
        self.latex(bits, state)?;
        state.claim_range(bits, None)
    }
}

#[macro_export]
macro_rules! declare_controlled_type
{
    ($(#[$attr:meta])* $name:ident, $gate_type:ty $(, $arg:ident)*) => {
        $(#[$attr])*
        pub struct $name
        {
            $( #[allow(dead_code)] $arg: f64, )*
            cgate: gates::C<$gate_type>
        }
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
                $name
                {
                    $( $arg: $arg, )*
                    cgate: gates::C::new(<$gate_type>::new($($arg, )*))
                }
            }
        }
    };
    ($name:ident, $gate_type:ty, cost=$cost:expr $(, $arg:ident)*) => {
        impl $name
        {
            pub fn new($($arg: f64, )*) -> Self
            {
                $name
                {
                    $( $arg: $arg, )*
                    cgate: gates::C::new(<$gate_type>::new($($arg, )*))
                }
            }
            pub fn cost() -> f64
            {
                $cost
            }
        }
    };
}

#[macro_export]
macro_rules! declare_controlled_cost
{
    ($cost:expr) => {
        fn cost(&self) -> f64 { $cost }
    };
    () => {
        fn cost(&self) -> f64 { self.cgate.cost() }
    };
}

#[macro_export]
macro_rules! declare_controlled_qasm
{
    ($trait_name:ident, $gate_name:ident, $method_name:ident $(, arg=$arg:ident)*) => {
        impl export::$trait_name for $gate_name
        {
            fn $method_name(&self, bit_names: &[String], bits: &[usize])
                -> error::ExportResult<String>
            {
                let mut res = stringify!($gate_name).to_lowercase();
                if bits.len() > 0
                {
                    res += &format!(" {}", bit_names[bits[0]]);
                    for &bit in bits[1..].iter()
                    {
                        res += &format!(", {}", &bit_names[bit]);
                    }
                }
                Ok(res)
            }
        }
    };
    ($trait_name:ident, $gate_name:ident, $method_name: ident, qasm=$qasm:expr $(, arg=$arg:ident)*) => {
        impl export::$trait_name for $gate_name
        {
            fn $method_name(&self, bit_names: &[String], bits: &[usize])
                -> error::ExportResult<String>
            {
                let mut res = String::from($qasm);
                for (i, &bit) in bits.iter().enumerate()
                {
                    let pattern = format!("{{{}}}", i);
                    res = res.replace(&pattern, &bit_names[bit]);
                }
                $(
                    let pattern = concat!("{", stringify!($arg), "}");
                    res = res.replace(pattern, &self.$arg.to_string());
                    let pattern = concat!("{half_", stringify!($arg), "}");
                    res = res.replace(pattern, &(0.5 * self.$arg).to_string());
                    let pattern = concat!("{quarter_", stringify!($arg), "}");
                    res = res.replace(pattern, &(0.25 * self.$arg).to_string());
                )*

                Ok(res)
            }
        }
    };
}

#[macro_export]
macro_rules! declare_controlled_latex
{
    ($gate_name:ident) => {
        impl export::Latex for $gate_name
        {
            fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
                -> error::Result<()>
            {
                self.cgate.latex(bits, state)
            }

            fn latex_checked(&self, bits: &[usize], state: &mut export::LatexExportState)
                -> error::Result<()>
            {
                self.cgate.latex_checked(bits, state)
            }
        }
    }
}

#[macro_export]
macro_rules! declare_controlled_impl_gate
{
    ($name:ident, $gate_type:ty $(, cost=$cost:expr)*) => {
        impl gates::Gate for $name
        {
            declare_controlled_cost!($($cost)*);
            fn description(&self) -> &str { self.cgate.description() }
            fn nr_affected_bits(&self) -> usize { self.cgate.nr_affected_bits() }
            fn matrix(&self) -> $crate::cmatrix::CMatrix { self.cgate.matrix() }
            fn apply_slice(&self, state: &mut $crate::cmatrix::CVecSliceMut)
            {
                self.cgate.apply_slice(state);
            }
        }
    };
}

#[macro_export]
macro_rules! declare_controlled
{
    ($(#[$attr:meta])* $name:ident, $gate_type:ty) => {
        declare_controlled_type!($(#[$attr])* $name, $gate_type);
        declare_controlled_impl!($name, $gate_type);
        declare_controlled_impl_gate!($name, $gate_type);
        declare_controlled_qasm!(OpenQasm, $name, open_qasm);
        declare_controlled_qasm!(CQasm, $name, c_qasm);
        declare_controlled_latex!($name);
    };
    ($(#[$attr:meta])* $name:ident, $gate_type:ty, cost=$cost:expr $(, arg=$arg:ident)* $(, open_qasm=$open_qasm:expr)* $(, c_qasm=$c_qasm:expr)*) => {
        declare_controlled_type!($(#[$attr])* $name, $gate_type $(, $arg)*);
        declare_controlled_impl!($name, $gate_type, cost=$cost $(, $arg)*);
        declare_controlled_impl_gate!($name, $gate_type, cost=Self::cost());
        declare_controlled_qasm!(OpenQasm, $name, open_qasm $(, qasm=$open_qasm)* $(, arg=$arg)*);
        declare_controlled_qasm!(CQasm, $name, c_qasm $(, qasm=$c_qasm)* $(, arg=$arg)*);
        declare_controlled_latex!($name);
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
    arg=theta);
declare_controlled!(
    /// Controlled `R`<sub>`Y`</sub> gate.
    CRY, gates::RY,
    cost=2.0*CX::cost() + 2.0*gates::U3::cost(),
    arg=theta,
    open_qasm="cx {0}, {1}; u3(-{theta}/2, 0, 0) {1}; cx {0}, {1}; u3({theta}/2, 0, 0) {1}",
    c_qasm="cnot {0}, {1}\nry {1}, -{half_theta}\ncnot {0}, {1}\nry {1}, {half_theta}");
declare_controlled!(
    /// Controlled `R`<sub>`Z`</sub> gate.
    CRZ, gates::RZ,
    cost=2.0*CX::cost() + 2.0*gates::U1::cost(),
    arg=lambda);

declare_controlled!(
    /// Controlled `S` gate.
    CS, gates::S, cost=2.0*CX::cost() + 3.0*gates::U1::cost());
declare_controlled!(
    /// Controlled `S`<sup>`†`</sup> gate.
    CSdg, gates::Sdg, cost=2.0*CX::cost() + 3.0*gates::U1::cost());

declare_controlled!(
    /// Controlled `T` gate.
    CT, gates::T, cost=2.0*CX::cost() + 3.0*gates::U1::cost());
declare_controlled!(
    /// Controlled `T`<sup>`†`</sup> gate.
    CTdg, gates::Tdg, cost=2.0*CX::cost() + 3.0*gates::U1::cost());

declare_controlled!(
    /// Controlled `U`<sub>`1`</sub> gate.
    CU1, gates::U1,
    cost=2.0*CX::cost() + 3.0*gates::U1::cost(),
    arg=lambda);
declare_controlled!(
    /// Controlled `U`<sub>`2`</sub> gate.
    CU2, gates::U2,
    cost=2.0*CX::cost() + 2.0*gates::U1::cost() + gates::U2::cost(),
    arg=phi, arg=lambda);
declare_controlled!(
    /// Controlled `U`<sub>`3`</sub> gate.
    CU3, gates::U3,
    cost=2.0*CX::cost() + gates::U1::cost() + 2.0*gates::U3::cost(),
    arg=theta, arg=phi, arg=lambda);

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
    CX, gates::X, cost=1001.0,
    c_qasm="cnot {0}, {1}");
declare_controlled!(
    /// Controlled `Y` gate.
    CY, gates::Y, cost=CX::cost() + 2.0*gates::U1::cost());
declare_controlled!(
    /// Controlled `Z` gate.
    CZ, gates::Z, cost=CX::cost() + 2.0*gates::U2::cost());

declare_controlled!(
    /// Doubly controlled `R`<sub>`Y`</sub> gate.
    CCRY, gates::CRY,
    cost=6.0 * gates::U3::cost() + 8.0*CX::cost(),
    arg=theta,
    open_qasm="cx {1}, {2}; u3(-{theta}/4, 0, 0) {2}; cx {1}, {2}; u3({theta}/4, 0, 0) {2}; cx {0}, {1}; cx {1}, {2}; u3({theta}/4, 0, 0) {2}; cx {1}, {2}; u3(-{theta}/4, 0, 0) {2}; cx {0}, {1}; cx {0}, {2}; u3(-{theta}/4, 0, 0) {2}; cx {0}, {2}; u3({theta}/4, 0, 0) {2}",
    c_qasm=r#"cnot {1}, {2}
ry {2}, -{quarter_theta}
cnot {1}, {2}
ry {2}, {quarter_theta}
cnot {0}, {1}
cnot {1}, {2}
ry {2}, {quarter_theta}
cnot {1}, {2}
ry {2}, -{quarter_theta}
cnot {0}, {1}
cnot {0}, {2}
ry {2}, -{quarter_theta}
cnot {0}, {2}
ry {2}, {quarter_theta}"#);

declare_controlled!(
    /// Doubly controlled `X` gate.
    CCX, gates::CX,
    cost=6.0*CX::cost() + 7.0*gates::U1::cost() + 2.0*gates::U2::cost(),
    c_qasm="toffoli {0}, {1}, {2}");
declare_controlled!(
    /// Doubly controlled `Z` gate.
    CCZ, gates::CZ,
    cost=CCX::cost() + 2.0*gates::H::cost(),
    open_qasm="h {2}; ccx {0}, {1}, {2}; h {2}",
    c_qasm="h {2}\ntoffoli {0}, {1}, {2}\nh {2}");

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, H, X};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use super::{C, CCRY, CCX, CCZ, CH, CRX, CRY, CRZ, CS, CV, CX, CY, CZ};
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
        assert_eq!(CRY::cost(), 2404.0);
        assert_eq!(CRZ::cost(), 2016.0);
        assert_eq!(CCX::new().cost(), 6263.0);
        assert_eq!(CCZ::new().cost(), 6471.0);
        assert_eq!(CCRY::new(1.6).cost(), 9214.0);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let open_qasm = CX::new().open_qasm(&bit_names, &[0, 1]);
        assert_eq!(open_qasm, Ok(String::from("cx qb0, qb1")));

        let bit_names = [String::from("qb0"), String::from("qb1"), String::from("qb2")];
        let open_qasm = CCZ::new().open_qasm(&bit_names, &[0, 1, 2]);
        assert_eq!(open_qasm, Ok(String::from("h qb2; ccx qb0, qb1, qb2; h qb2")));

        let bit_names = [String::from("qb0"), String::from("qb1"), String::from("qb2")];
        let open_qasm = CCRY::new(1.6).open_qasm(&bit_names, &[1, 2, 0]);
        assert_eq!(open_qasm, Ok(String::from("cx qb2, qb0; u3(-1.6/4, 0, 0) qb0; cx qb2, qb0; u3(1.6/4, 0, 0) qb0; cx qb1, qb2; cx qb2, qb0; u3(1.6/4, 0, 0) qb0; cx qb2, qb0; u3(-1.6/4, 0, 0) qb0; cx qb1, qb2; cx qb1, qb0; u3(-1.6/4, 0, 0) qb0; cx qb1, qb0; u3(1.6/4, 0, 0) qb0")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let c_qasm = CX::new().c_qasm(&bit_names, &[0, 1]);
        assert_eq!(c_qasm, Ok(String::from("cnot qb0, qb1")));

        let bit_names = [String::from("qb0"), String::from("qb1"), String::from("qb2")];
        let c_qasm = CCZ::new().c_qasm(&bit_names, &[0, 1, 2]);
        assert_eq!(c_qasm, Ok(String::from("h qb2\ntoffoli qb0, qb1, qb2\nh qb2")));

        let bit_names = [String::from("qb0"), String::from("qb1"), String::from("qb2")];
        let c_qasm = CCRY::new(1.6).c_qasm(&bit_names, &[1, 2, 0]);
        assert_eq!(c_qasm, Ok(String::from(
r#"cnot qb2, qb0
ry qb0, -0.4
cnot qb2, qb0
ry qb0, 0.4
cnot qb1, qb2
cnot qb2, qb0
ry qb0, 0.4
cnot qb2, qb0
ry qb0, -0.4
cnot qb1, qb2
cnot qb1, qb0
ry qb0, -0.4
cnot qb1, qb0
ry qb0, 0.4"#)));
    }

    #[test]
    fn test_latex()
    {
        let gate = CS::new();
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex_checked(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \ctrl{1} & \qw \\
    \lstick{\ket{0}} & \gate{S} & \qw \\
}
"#);

        let gate = CV::new();
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex_checked(&[1, 0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{V} & \qw \\
    \lstick{\ket{0}} & \ctrl{-1} & \qw \\
}
"#);
    }

    #[test]
    #[should_panic]
    fn test_latex_error()
    {
        let gate = CCX::new();
        let mut state = LatexExportState::new(3, 0);
        let _ltx = gate.latex_checked(&[1, 2, 0], &mut state);
    }
}
