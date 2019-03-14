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

use gates::Gate;

/// The Clifford `S` gate
///
/// The `S` gate rotates the state over π/2 radians around the `z` axis of
/// the Bloch sphere. It is the square root of the `Z` gate. The associated
/// matrix is
/// ```text
/// ┌     ┐
/// │ 1 0 │
/// │ 0 i │
/// └     ┘
/// ```
pub struct S
{
}

impl S
{
    /// Create a new `S` gate.
    pub fn new() -> Self
    {
        S { }
    }
}

impl gates::Gate for S
{
    fn cost(&self) -> f64
    {
        gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "S"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        array![[o, z], [z, i]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= cmatrix::COMPLEX_I;
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= cmatrix::COMPLEX_I;
    }
}

impl export::OpenQasm for S
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("s {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for S
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("s {}", bit_names[bits[0]]))
    }
}

impl export::Latex for S
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;

        state.set_field(bits[0], String::from(r"\gate{S}"));
        Ok(())
    }
}

/// Conjugate of Clifford `S` gate
///
/// The `S`<sup>`†`</sup> gate rotates the state over -π/2 radians around the
/// `z` axis of the Bloch sphere. It is the conjugate of the `S` gate. The
/// associated matrix is
/// ```text
/// ┌      ┐
/// │ 1  0 │
/// │ 0 -i │
/// └      ┘
/// ```
pub struct Sdg
{
}

impl Sdg
{
    /// Create a new `S`<sup>`†`</sup> gate.
    pub fn new() -> Self
    {
        Sdg { }
    }
}

impl gates::Gate for Sdg
{
    fn cost(&self) -> f64
    {
        gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "S†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        array![[o, z], [z, -i]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= -cmatrix::COMPLEX_I;
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= -cmatrix::COMPLEX_I;
    }
}

impl export::OpenQasm for Sdg
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("sdg {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for Sdg
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("sdag {}", bit_names[bits[0]]))
    }
}

impl export::Latex for Sdg
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;

        state.set_field(bits[0], String::from(r"\gate{S^\dagger}"));
        Ok(())
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, S, Sdg};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = S::new();
        assert_eq!(gate.description(), "S");
        let gate = Sdg::new();
        assert_eq!(gate.description(), "S†");
    }

    #[test]
    fn test_cost()
    {
        let gate = S::new();
        assert_eq!(gate.cost(), 7.0);
        let gate = Sdg::new();
        assert_eq!(gate.cost(), 7.0);
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;

        let gate = S::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, i]]);

        let gate = Sdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, -i]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [o, z,   x,    x],
            [z, i, x*i, -x*i]
        ];
        gate_test(S::new(), &mut state, &result);

        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [o,  z,    x,   x],
            [z, -i, -x*i, x*i]
        ];
        gate_test(Sdg::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = S::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("s qb")));
        let qasm = Sdg::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("sdg qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = S::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("s qb")));
        let qasm = Sdg::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("sdag qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = S::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{S} & \qw \\
}
"#);

        let gate = Sdg::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{S^\dagger} & \qw \\
}
"#);
    }
}
