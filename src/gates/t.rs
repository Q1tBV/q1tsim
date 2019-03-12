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

use cmatrix;
use gates;
use error;
use export;

/// The `T` gate
///
/// The `T` gate rotates the state over π/4 radians around the `z` axis of
/// the Bloch sphere. It is the square root of the `S` gate.
pub struct T
{
}

impl T
{
    /// Create a new `T` gate.
    pub fn new() -> Self
    {
        T { }
    }
}

impl gates::Gate for T
{
    fn cost(&self) -> f64
    {
        gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "T"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        array![[o, z], [z, x+x*i]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);
    }
}

impl export::OpenQasm for T
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("t {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for T
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("t {}", bit_names[bits[0]]))
    }
}

impl export::Latex for T
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
    {
        assert!(bits.len() == 1, "T gate operates on a single bit");
        state.set_field(bits[0], String::from(r"\gate{T}"));
    }
}

/// Conjugate of `T` gate
///
/// The `T`<sup>`†`</sup> gate rotates the state over -π/4 radians around the
/// `z` axis of the Bloch sphere. It is the conjugate of the `T` gate.
pub struct Tdg
{
}

impl Tdg
{
    /// Create a new `T`<sup>`†`</sup> gate.
    pub fn new() -> Self
    {
        Tdg { }
    }
}

impl gates::Gate for Tdg
{
    fn cost(&self) -> f64
    {
        gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "T†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        array![[o, z], [z, x-x*i]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= num_complex::Complex::from_polar(&1.0, &-::std::f64::consts::FRAC_PI_4);
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= num_complex::Complex::from_polar(&1.0, &-::std::f64::consts::FRAC_PI_4);
    }
}

impl export::OpenQasm for Tdg
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("tdg {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for Tdg
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("tdag {}", bit_names[bits[0]]))
    }
}

impl export::Latex for Tdg
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
    {
        assert!(bits.len() == 1, "Tdg gate operates on a single bit");
        state.set_field(bits[0], String::from(r"\gate{T^\dagger}"));
    }
}

#[cfg(test)]
mod tests
{
    extern crate num_complex;

    use super::{T, Tdg};
    use gates::{gate_test, Gate};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = T::new();
        assert_eq!(gate.description(), "T");
        let gate = Tdg::new();
        assert_eq!(gate.description(), "T†");
    }

    #[test]
    fn test_cost()
    {
        let gate = T::new();
        assert_eq!(gate.cost(), 7.0);
        let gate = Tdg::new();
        assert_eq!(gate.cost(), 7.0);
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let t = num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);

        let gate = T::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, t]]);

        let gate = Tdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, t.conj()]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;
        let x = cmatrix::COMPLEX_HSQRT2;
        let t = num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);
        let td = t.conj();

        let mut state = array![
            [o, z, x,  h, z],
            [z, o, z, -h, z],
            [z, z, x,  h, z],
            [z, z, z, -h, o]
        ];
        let result = array![
            [o, z,   x,    h, z],
            [z, o,   z,   -h, z],
            [z, z, t*x,  t*h, z],
            [z, z,   z, -t*h, t]
        ];
        gate_test(T::new(), &mut state, &result);

        let mut state = array![
            [o, z, x,  h, z],
            [z, o, z, -h, z],
            [z, z, x,  h, z],
            [z, z, z, -h, o]
        ];
        let result = array![
            [o, z,    x,     h,  z],
            [z, o,    z,    -h,  z],
            [z, z, td*x,  td*h,  z],
            [z, z,    z, -td*h, td]
        ];
        gate_test(Tdg::new(), &mut state, &result);
    }

    #[test]
    fn test_apply_mat()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;
        let x = cmatrix::COMPLEX_HSQRT2;
        let t = num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);
        let td = t.conj();

        let mut state = array![
            [o, z, x,  h, z],
            [z, o, z, -h, z],
            [z, z, x,  h, z],
            [z, z, z, -h, o]
        ];
        T::new().apply_mat(&mut state);
        assert_complex_matrix_eq!(&state, &array![
            [o, z,   x,    h, z],
            [z, o,   z,   -h, z],
            [z, z, t*x,  t*h, z],
            [z, z,   z, -t*h, t]
        ]);

        let mut state = array![
            [o, z, x,  h, z],
            [z, o, z, -h, z],
            [z, z, x,  h, z],
            [z, z, z, -h, o]
        ];
        Tdg::new().apply_mat(&mut state);
        assert_complex_matrix_eq!(&state, &array![
            [o, z,    x,     h,  z],
            [z, o,    z,    -h,  z],
            [z, z, td*x,  td*h,  z],
            [z, z,    z, -td*h, td]
        ]);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = T::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("t qb")));
        let qasm = Tdg::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("tdg qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = T::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("t qb")));
        let qasm = Tdg::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("tdag qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = T::new();
        let mut state = LatexExportState::new(1, 0);
        gate.latex_checked(&[0], &mut state);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{T} & \qw \\
}
"#);

        let gate = Tdg::new();
        let mut state = LatexExportState::new(1, 0);
        gate.latex_checked(&[0], &mut state);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{T^\dagger} & \qw \\
}
"#);
    }
}
