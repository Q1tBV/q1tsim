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

/// Controlled `X` gate.
pub struct CX
{
    cgate: crate::gates::C<crate::gates::X>
}

impl CX
{
    pub fn new() -> Self
    {
        CX
        {
            cgate: crate::gates::C::new(crate::gates::X::new())
        }
    }

    pub fn cost() -> f64
    {
        1001.0
    }
}

impl crate::gates::Gate for CX
{
    fn cost(&self) -> f64 { Self::cost() }
    fn description(&self) -> &str { self.cgate.description() }
    fn nr_affected_bits(&self) -> usize { self.cgate.nr_affected_bits() }
    fn matrix(&self) -> crate::cmatrix::CMatrix { self.cgate.matrix() }
    fn apply_slice(&self, state: crate::cmatrix::CVecSliceMut)
    {
        self.cgate.apply_slice(state);
    }
    fn apply_mat_slice(&self, state: crate::cmatrix::CMatSliceMut)
    {
        self.cgate.apply_mat_slice(state);
    }
    fn is_stabilizer(&self) -> bool
    {
        true
    }
    fn conjugate(&self, ops: &mut [crate::stabilizer::PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        let (phase, op0, op1) = match (ops[0], ops[1])
        {
            (PauliOp::I, PauliOp::I) => (false, PauliOp::I, PauliOp::I),
            (PauliOp::I, PauliOp::Z) => (false, PauliOp::Z, PauliOp::Z),
            (PauliOp::I, PauliOp::X) => (false, PauliOp::I, PauliOp::X),
            (PauliOp::I, PauliOp::Y) => (false, PauliOp::Z, PauliOp::Y),
            (PauliOp::Z, PauliOp::I) => (false, PauliOp::Z, PauliOp::I),
            (PauliOp::Z, PauliOp::Z) => (false, PauliOp::I, PauliOp::Z),
            (PauliOp::Z, PauliOp::X) => (false, PauliOp::Z, PauliOp::X),
            (PauliOp::Z, PauliOp::Y) => (false, PauliOp::I, PauliOp::Y),
            (PauliOp::X, PauliOp::I) => (false, PauliOp::X, PauliOp::X),
            (PauliOp::X, PauliOp::Z) => (true,  PauliOp::Y, PauliOp::Y),
            (PauliOp::X, PauliOp::X) => (false, PauliOp::X, PauliOp::I),
            (PauliOp::X, PauliOp::Y) => (false, PauliOp::Y, PauliOp::Z),
            (PauliOp::Y, PauliOp::I) => (false, PauliOp::Y, PauliOp::X),
            (PauliOp::Y, PauliOp::Z) => (false, PauliOp::X, PauliOp::Y),
            (PauliOp::Y, PauliOp::X) => (false, PauliOp::Y, PauliOp::I),
            (PauliOp::Y, PauliOp::Y) => (true,  PauliOp::X, PauliOp::Z),
        };
        ops[0] = op0;
        ops[1] = op1;
        Ok(phase)
    }
}

impl crate::export::OpenQasm for CX
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        self.check_nr_bits(bits.len())?;
        Ok(format!("cx {}, {}", bit_names[bits[0]], bit_names[bits[1]]))
    }
}

impl crate::export::CQasm for CX
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        self.cgate.check_nr_bits(bits.len())?;
        Ok(format!("cnot {}, {}", bit_names[bits[0]], bit_names[bits[1]]))
    }
}

impl crate::export::Latex for CX
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.cgate.latex(bits, state)
    }
}

#[cfg(test)]
mod tests
{
    use super::CX;
    use crate::cmatrix;
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use crate::gates::{gate_test, Gate};
    use crate::stabilizer::PauliOp;

    #[test]
    fn test_description()
    {
        let gate = CX::new();
        assert_eq!(gate.description(), "CX");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let gate = CX::new();
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z, z],
            [z, o, z, z],
            [z, z, z, o],
            [z, z, o, z]
        ]);
    }

    #[test]
    fn test_apply()
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
    }

    #[test]
    fn test_cost()
    {
        assert_eq!(CX::new().cost(), 1001.0);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let open_qasm = CX::new().open_qasm(&bit_names, &[0, 1]);
        assert_eq!(open_qasm, Ok(String::from("cx qb0, qb1")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let c_qasm = CX::new().c_qasm(&bit_names, &[0, 1]);
        assert_eq!(c_qasm, Ok(String::from("cnot qb0, qb1")));
    }

    #[test]
    fn test_latex()
    {
        let gate = CX::new();
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \ctrl{1} & \qw \\
    \lstick{\ket{0}} & \targ & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        const OPS_X: [(PauliOp, PauliOp, PauliOp, PauliOp, bool); 16] = [
            (PauliOp::I, PauliOp::I, PauliOp::I, PauliOp::I, false),
            (PauliOp::I, PauliOp::Z, PauliOp::Z, PauliOp::Z, false),
            (PauliOp::I, PauliOp::X, PauliOp::I, PauliOp::X, false),
            (PauliOp::I, PauliOp::Y, PauliOp::Z, PauliOp::Y, false),
            (PauliOp::Z, PauliOp::I, PauliOp::Z, PauliOp::I, false),
            (PauliOp::Z, PauliOp::Z, PauliOp::I, PauliOp::Z, false),
            (PauliOp::Z, PauliOp::X, PauliOp::Z, PauliOp::X, false),
            (PauliOp::Z, PauliOp::Y, PauliOp::I, PauliOp::Y, false),
            (PauliOp::X, PauliOp::I, PauliOp::X, PauliOp::X, false),
            (PauliOp::X, PauliOp::Z, PauliOp::Y, PauliOp::Y, true ),
            (PauliOp::X, PauliOp::X, PauliOp::X, PauliOp::I, false),
            (PauliOp::X, PauliOp::Y, PauliOp::Y, PauliOp::Z, false),
            (PauliOp::Y, PauliOp::I, PauliOp::Y, PauliOp::X, false),
            (PauliOp::Y, PauliOp::Z, PauliOp::X, PauliOp::Y, false),
            (PauliOp::Y, PauliOp::X, PauliOp::Y, PauliOp::I, false),
            (PauliOp::Y, PauliOp::Y, PauliOp::X, PauliOp::Z, true )
        ];

        let mut op = [PauliOp::I; 2];
        for &(op0, op1, r0, r1, p) in &OPS_X
        {
            op[0] = op0;
            op[1] = op1;
            let phase = CX::new().conjugate(&mut op).unwrap();
            assert_eq!(op[0], r0);
            assert_eq!(op[1], r1);
            assert_eq!(phase, p);
        }
    }
}
