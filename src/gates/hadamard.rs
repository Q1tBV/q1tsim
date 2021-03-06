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

/// The Hadamard gate.
///
/// The Hadamard gate maps the zero state |0&rang; to the symmetric combination
/// of |0&rang; and |1&rang;, and the |1&rang; state to the anti-symmetric
/// combination.
#[derive(Clone)]
pub struct H
{
}

impl H
{
    /// Create a new Hadamard gate.
    pub fn new() -> Self
    {
        H { }
    }

    pub fn transform(mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let s1_copy = state.slice(s![n..]).to_owned();

        let (mut s0, mut s1) = state.view_mut().split_at(ndarray::Axis(0), n);

        s1 -= &s0;
        s1 *= -crate::cmatrix::COMPLEX_HSQRT2;
        s0 += &s1_copy;
        s0 *= crate::cmatrix::COMPLEX_HSQRT2;
    }

    pub fn transform_mat(mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let s1_copy = state.slice(s![n.., ..]).to_owned();

        let (mut s0, mut s1) = state.view_mut().split_at(ndarray::Axis(0), n);

        s1 -= &s0;
        s1 *= -crate::cmatrix::COMPLEX_HSQRT2;
        s0 += &s1_copy;
        s0 *= crate::cmatrix::COMPLEX_HSQRT2;
    }

    pub fn cost() -> f64
    {
        crate::gates::U2::cost()
    }
}

impl crate::gates::Gate for H
{
    fn cost(&self) -> f64
    {
        Self::cost()
    }

    fn description(&self) -> &str
    {
        "H"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        array![[x, x], [x, -x]]
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
        let (op, phase) = match ops[0]
            {
                PauliOp::I => (PauliOp::I, false),
                PauliOp::Z => (PauliOp::X, false),
                PauliOp::X => (PauliOp::Z, false),
                PauliOp::Y => (PauliOp::Y, true)
            };
        ops[0] = op;
        Ok(phase)
    }
}

impl crate::export::OpenQasm for H
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("h {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for H
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("h {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for H
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, "H")
    }
}

impl crate::arithmetic::Square for H
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
    use crate::gates::{gate_test, Gate, H};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use crate::arithmetic::Square;
    use crate::stabilizer::PauliOp;

    #[test]
    fn test_description()
    {
        let gate = H::new();
        assert_eq!(gate.description(), "H");
    }

    #[test]
    fn test_cost()
    {
        let gate = H::new();
        assert_eq!(gate.cost(), 104.0);
    }

    #[test]
    fn test_matrix()
    {
        let h = H::new();
        let s = crate::cmatrix::COMPLEX_HSQRT2;
        assert_complex_matrix_eq!(h.matrix(), array![[s, s], [s, -s]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[x, x, o, z], [x, -x, z, o]];
        gate_test(H::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = H::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("h qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = H::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("h qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = H::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{H} & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let gate = H::new();

        let mut ops = [PauliOp::I];
        assert_eq!(gate.conjugate(&mut ops), Ok(false));
        assert_eq!(ops, [PauliOp::I]);

        let mut ops = [PauliOp::Z];
        assert_eq!(gate.conjugate(&mut ops), Ok(false));
        assert_eq!(ops, [PauliOp::X]);

        let mut ops = [PauliOp::X];
        assert_eq!(gate.conjugate(&mut ops), Ok(false));
        assert_eq!(ops, [PauliOp::Z]);

        let mut ops = [PauliOp::Y];
        assert_eq!(gate.conjugate(&mut ops), Ok(true));
        assert_eq!(ops, [PauliOp::Y]);
    }

    #[test]
    fn test_square()
    {
        let gate = H::new();
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);
    }
}
