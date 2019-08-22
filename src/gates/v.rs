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

/// The `V` gate
///
/// The `V` gate is the square root of the `X` gate. The associated matrix is
/// ```text
///     ┌         ┐
/// 1/2 │ 1+i 1-i │
///     │ 1-i 1+i │
///     └         ┘
/// ```
#[derive(Clone)]
pub struct V
{
}

impl V
{
    /// Create a new `V` gate.
    pub fn new() -> Self
    {
        V { }
    }
}

impl crate::gates::Gate for V
{
    fn cost(&self) -> f64
    {
        crate::gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "V"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let h = 0.5 * crate::cmatrix::COMPLEX_ONE;
        let hi = 0.5 * crate::cmatrix::COMPLEX_I;
        array![[h+hi, h-hi], [h-hi, h+hi]]
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
            PauliOp::Z => (PauliOp::Y, true),
            PauliOp::X => (PauliOp::X, false),
            PauliOp::Y => (PauliOp::Z, false)
        };
        ops[0] = op;
        Ok(phase)
    }
}

impl crate::export::OpenQasm for V
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("u3(pi/2, -pi/2, pi/2) {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for V
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("x90 {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for V
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, "V")
    }
}

impl crate::arithmetic::Square for V
{
    type SqType = crate::gates::X;

    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        Ok(crate::gates::X::new())
    }
}

/// Conjugate of `V` gate.
///
/// The `V`<sup>`†`</sup> is the conjugate of the `V` gate.
#[derive(Clone)]
pub struct Vdg
{
}

impl Vdg
{
    /// Create a new `V`<sup>`†`</sup> gate.
    pub fn new() -> Self
    {
        Vdg { }
    }
}

impl crate::gates::Gate for Vdg
{
    fn cost(&self) -> f64
    {
        crate::gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "V†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let h = 0.5 * crate::cmatrix::COMPLEX_ONE;
        let hi = 0.5 * crate::cmatrix::COMPLEX_I;
        array![[h-hi, h+hi], [h+hi, h-hi]]
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
            PauliOp::Z => (PauliOp::Y, false),
            PauliOp::X => (PauliOp::X, false),
            PauliOp::Y => (PauliOp::Z, true)
        };
        ops[0] = op;
        Ok(phase)
    }
}

impl crate::export::OpenQasm for Vdg
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("u3(pi/2, pi/2, -pi/2) {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for Vdg
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("mx90 {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for Vdg
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, r"V^\dagger")
    }
}

impl crate::arithmetic::Square for Vdg
{
    type SqType = crate::gates::X;

    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        Ok(crate::gates::X::new())
    }
}

#[cfg(test)]
mod tests
{
    use super::{V, Vdg};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use crate::gates::{gate_test, Gate};
    use crate::arithmetic::Square;
    use crate::stabilizer::PauliOp;

    #[test]
    fn test_description()
    {
        let gate = V::new();
        assert_eq!(gate.description(), "V");
        let gate = Vdg::new();
        assert_eq!(gate.description(), "V†");
    }

    #[test]
    fn test_cost()
    {
        let gate = V::new();
        assert_eq!(gate.cost(), 201.0);
        let gate = Vdg::new();
        assert_eq!(gate.cost(), 201.0);
    }

    #[test]
    fn test_matrix()
    {
        let h = 0.5*crate::cmatrix::COMPLEX_ONE;
        let hi = 0.5*crate::cmatrix::COMPLEX_I;

        let gate = V::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[h+hi, h-hi], [h-hi, h+hi]]);

        let gate = Vdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[h-hi, h+hi], [h+hi, h-hi]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * o;
        let hi = 0.5 * i;

        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [h+hi, h-hi, x,  x*i],
            [h-hi, h+hi, x, -x*i]
        ];
        gate_test(V::new(), &mut state, &result);

        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [h-hi, h+hi, x, -x*i],
            [h+hi, h-hi, x,  x*i]
        ];
        gate_test(Vdg::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = V::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("u3(pi/2, -pi/2, pi/2) qb")));
        let qasm = Vdg::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("u3(pi/2, pi/2, -pi/2) qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = V::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("x90 qb")));
        let qasm = Vdg::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("mx90 qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = V::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{V} & \qw \\
}
"#);

        let gate = Vdg::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{V^\dagger} & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let mut op = [PauliOp::I];
        assert_eq!(V::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::I]);

        let mut op = [PauliOp::Z];
        assert_eq!(V::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::Y]);

        let mut op = [PauliOp::X];
        assert_eq!(V::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::X]);

        let mut op = [PauliOp::Y];
        assert_eq!(V::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::Z]);

        let mut op = [PauliOp::I];
        assert_eq!(Vdg::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::I]);

        let mut op = [PauliOp::Z];
        assert_eq!(Vdg::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::Y]);

        let mut op = [PauliOp::X];
        assert_eq!(Vdg::new().conjugate(&mut op), Ok(false));
        assert_eq!(op, [PauliOp::X]);

        let mut op = [PauliOp::Y];
        assert_eq!(Vdg::new().conjugate(&mut op), Ok(true));
        assert_eq!(op, [PauliOp::Z]);
    }

    #[test]
    fn test_square()
    {
        let gate = V::new();
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);

        let gate = Vdg::new();
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);
    }
}
