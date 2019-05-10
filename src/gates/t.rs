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

use crate::gates::Gate;

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

impl crate::gates::Gate for T
{
    fn cost(&self) -> f64
    {
        crate::gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "T"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;
        array![[o, z], [z, x+x*i]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);
    }
}

impl crate::export::OpenQasm for T
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("t {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for T
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("t {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for T
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, "T")
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

impl crate::gates::Gate for Tdg
{
    fn cost(&self) -> f64
    {
        crate::gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "T†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;
        array![[o, z], [z, x-x*i]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= num_complex::Complex::from_polar(&1.0, &-::std::f64::consts::FRAC_PI_4);
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= num_complex::Complex::from_polar(&1.0, &-::std::f64::consts::FRAC_PI_4);
    }
}

impl crate::export::OpenQasm for Tdg
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("tdg {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for Tdg
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("tdag {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for Tdg
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, r"T^\dagger")
    }
}

#[cfg(test)]
mod tests
{
    use super::{T, Tdg};
    use crate::gates::{gate_test, Gate};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};

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
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let t = num_complex::Complex::from_polar(&1.0, &::std::f64::consts::FRAC_PI_4);

        let gate = T::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, t]]);

        let gate = Tdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, t.conj()]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
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
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
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
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{T} & \qw \\
}
"#);

        let gate = Tdg::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{T^\dagger} & \qw \\
}
"#);
    }
}
