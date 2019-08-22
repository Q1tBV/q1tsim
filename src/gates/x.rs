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

/// The Pauli X gate.
///
/// The X, or NOT, gate rotates the state over π radians around the `x` axis of
/// the Bloch sphere, i.e. it swaps the |0⟩ and |1⟩ components of the qubit.
#[derive(Clone)]
pub struct X
{
}

impl X
{
    /// Create a new Pauli X gate.
    pub fn new() -> Self
    {
        X { }
    }

    pub fn transform(mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        for i in 0..n
        {
            state.swap(i, i+n);
        }
    }

    pub fn transform_mat(mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let m = state.cols();
        for i in 0..n
        {
            for j in 0..m
            {
                state.swap((i, j), (i+n, j));
            }
        }
    }
}

impl crate::gates::Gate for X
{
    fn cost(&self) -> f64
    {
        crate::gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "X"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        array![[z, o], [o, z]]
    }

    fn apply_slice(&self, state: crate::cmatrix::CVecSliceMut)
    {
        Self::transform(state);
    }

    fn apply_mat_slice(&self, state: crate::cmatrix::CMatSliceMut)
    {
        Self::transform_mat(state);
    }

    fn is_stabilizer(&self) -> bool
    {
        true
    }

    fn conjugate(&self, ops: &mut [PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        Ok(ops[0] == PauliOp::Z || ops[0] == PauliOp::Y)
    }
}

impl crate::export::OpenQasm for X
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("x {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for X
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("x {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for X
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;

        let symbol = if state.is_controlled() { r"\targ" } else { r"\gate{X}" };
        state.set_field(bits[0], String::from(symbol))
    }
}

impl crate::arithmetic::Square for X
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
    use crate::export::{LatexExportState, Latex, OpenQasm, CQasm};
    use crate::gates::{gate_test, Gate, X};
    use crate::arithmetic::Square;
    use crate::stabilizer::PauliOp;

    #[test]
    fn test_description()
    {
        let x = X::new();
        assert_eq!(x.description(), "X");
    }

    #[test]
    fn test_matrix()
    {
        let x = X::new();
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(x.matrix(), array![[z, o], [o, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[z, o, x, -x],[o, z, x, x]];
        gate_test(X::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = X::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("x qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = X::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("x qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = X::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{X} & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let mut op = [PauliOp::I];
        assert_eq!(X::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::I]);

        let mut op = [PauliOp::Z];
        assert_eq!(X::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::Z]);

        let mut op = [PauliOp::X];
        assert_eq!(X::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::X]);

        let mut op = [PauliOp::Y];
        assert_eq!(X::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::Y]);
    }

    #[test]
    fn test_square()
    {
        let gate = X::new();
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);
    }
}
