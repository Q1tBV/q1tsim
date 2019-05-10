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
use crate::stabilizer::PauliOp;

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

impl crate::gates::Gate for S
{
    fn cost(&self) -> f64
    {
        crate::gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "S"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;
        array![[o, z], [z, i]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= crate::cmatrix::COMPLEX_I;
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= crate::cmatrix::COMPLEX_I;
    }

    fn conjugate(&self, ops: &mut [PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        let (op, phase) = match ops[0]
        {
            PauliOp::I => (PauliOp::I, false),
            PauliOp::Z => (PauliOp::Z, false),
            PauliOp::X => (PauliOp::Y, false),
            PauliOp::Y => (PauliOp::X, true)
        };
        ops[0] = op;
        Ok(phase)
    }
}

impl crate::export::OpenQasm for S
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("s {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for S
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("s {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for S
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, "S")
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

impl crate::gates::Gate for Sdg
{
    fn cost(&self) -> f64
    {
        crate::gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "S†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;
        array![[o, z], [z, -i]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= -crate::cmatrix::COMPLEX_I;
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let mut slice = state.slice_mut(s![n.., ..]);
        slice *= -crate::cmatrix::COMPLEX_I;
    }

    fn conjugate(&self, ops: &mut [PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        let (op, phase) = match ops[0]
        {
            PauliOp::I => (PauliOp::I, false),
            PauliOp::Z => (PauliOp::Z, false),
            PauliOp::X => (PauliOp::Y, true),
            PauliOp::Y => (PauliOp::X, false)
        };
        ops[0] = op;
        Ok(phase)
    }
}

impl crate::export::OpenQasm for Sdg
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("sdg {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for Sdg
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("sdag {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for Sdg
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, r"S^\dagger")
    }
}

#[cfg(test)]
mod tests
{
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use crate::gates::{gate_test, Gate, S, Sdg};
    use crate::stabilizer::PauliOp;

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
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;

        let gate = S::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, i]]);

        let gate = Sdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, -i]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;
        let x = crate::cmatrix::COMPLEX_HSQRT2;

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
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{S} & \qw \\
}
"#);

        let gate = Sdg::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{S^\dagger} & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let mut op = [PauliOp::I];
        assert_eq!(S::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::I]);

        let mut op = [PauliOp::Z];
        assert_eq!(S::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::Z]);

        let mut op = [PauliOp::X];
        assert_eq!(S::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::Y]);

        let mut op = [PauliOp::Y];
        assert_eq!(S::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::X]);

        let mut op = [PauliOp::I];
        assert_eq!(Sdg::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::I]);

        let mut op = [PauliOp::Z];
        assert_eq!(Sdg::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::Z]);

        let mut op = [PauliOp::X];
        assert_eq!(Sdg::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::Y]);

        let mut op = [PauliOp::Y];
        assert_eq!(Sdg::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::X]);
    }
}
