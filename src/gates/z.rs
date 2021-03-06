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

/// The Pauli Z gate.
///
/// The Z gate rotates the state over π radians around the `z` axis of
/// the Bloch sphere, i.e. it flips the sign of the |1⟩ components of the qubit.
#[derive(Clone)]
pub struct Z
{
}

impl Z
{
    /// Create a new Pauli Z gate.
    pub fn new() -> Self
    {
        Z { }
    }
}

impl crate::gates::Gate for Z
{
    fn cost(&self) -> f64
    {
        crate::gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "Z"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        array![[o, z], [z, -o]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        state.slice_mut(s![n..]).mapv_inplace(|c| -c);
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        state.slice_mut(s![n.., ..]).mapv_inplace(|c| -c);
    }

    fn is_stabilizer(&self) -> bool
    {
        true
    }

    fn conjugate(&self, ops: &mut [PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        Ok(ops[0] == PauliOp::X || ops[0] == PauliOp::Y)
    }
}

impl crate::export::OpenQasm for Z
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("z {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for Z
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("z {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for Z
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;

        let symbol = if state.is_controlled() { r"\control \qw" } else { r"\gate{Z}" };
        state.set_field(bits[0], String::from(symbol))
    }
}

impl crate::arithmetic::Square for Z
{
    type SqType = crate::gates::I;

    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        Ok(crate::gates::I::new())
    }
}

#[cfg(test)]
mod tests
{
    use super::Z;
    use crate::gates::{gate_test, Gate};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use crate::arithmetic::Square;
    use crate::stabilizer::PauliOp;

    #[test]
    fn test_description()
    {
        let gate = Z::new();
        assert_eq!(gate.description(), "Z");
    }

    #[test]
    fn test_cost()
    {
        let gate = Z::new();
        assert_eq!(gate.cost(), 7.0);
    }

    #[test]
    fn test_matrix()
    {
        let gate = Z::new();
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, -o]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [ z, -o, -x, x]];
        gate_test(Z::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Z::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("z qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Z::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("z qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = Z::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{Z} & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let mut op = [PauliOp::I];
        assert_eq!(Z::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::I]);

        let mut op = [PauliOp::Z];
        assert_eq!(Z::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::Z]);

        let mut op = [PauliOp::X];
        assert_eq!(Z::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::X]);

        let mut op = [PauliOp::Y];
        assert_eq!(Z::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::Y]);
    }

    #[test]
    fn test_square()
    {
        let gate = Z::new();
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);
    }
}
