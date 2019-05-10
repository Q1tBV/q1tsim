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

/// The `Swap` gate
///
/// The `Swap` gate swap two qubits.
pub struct Swap
{
}

impl Swap
{
    /// Create a `Swap` gate
    pub fn new() -> Self
    {
        Swap { }
    }

    pub fn transform(mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 4 == 0, "Number of rows is not a mutiple of 4.");

        let n = state.len() / 4;
        for i in n..2*n
        {
            state.swap(i, i+n);
        }
    }

    pub fn transform_mat(mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.len() % 4 == 0, "Number of rows is not a multiple of 4.");

        let n = state.rows() / 4;
        let m = state.cols();
        for i in n..2*n
        {
            for j in 0..m
            {
                state.swap((i, j), (i+n, j));
            }
        }
    }
}

impl crate::gates::Gate for Swap
{
    fn cost(&self) -> f64
    {
        3.0 * crate::gates::CX::cost()
    }

    fn description(&self) -> &str
    {
        "Swap"
    }

    fn nr_affected_bits(&self) -> usize
    {
        2
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        array![
            [o, z, z, z],
            [z, z, o, z],
            [z, o, z, z],
            [z, z, z, o]
        ]
    }

    fn apply_slice(&self, state: crate::cmatrix::CVecSliceMut)
    {
        Self::transform(state);
    }

    fn apply_mat_slice(&self, state: crate::cmatrix::CMatSliceMut)
    {
        Self::transform_mat(state);
    }

    fn conjugate(&self, ops: &mut [PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        ops.swap(0, 1);
        Ok(false)
    }
}

impl crate::export::OpenQasm for Swap
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        let b0 = &bit_names[bits[0]];
        let b1 = &bit_names[bits[1]];
        Ok(format!("cx {}, {}; cx {}, {}; cx {}, {}", b0, b1, b1, b0, b0, b1))
    }
}

impl crate::export::CQasm for Swap
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("swap {}, {}", bit_names[bits[0]], bit_names[bits[1]]))
    }
}

impl crate::export::Latex for Swap
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;

        let (mut b0, mut b1) = (bits[0], bits[1]);
        if b1 < b0
        {
            ::std::mem::swap(&mut b0, &mut b1);
        }

        state.start_range_op(bits, None)?;
        state.set_field(b0, format!(r"\qswap \qwx[{}]", b1-b0))?;
        state.set_field(b1, String::from(r"\qswap"))?;
        state.end_range_op();

        Ok(())
    }
}


#[cfg(test)]
mod tests
{
    use super::Swap;
    use crate::export::{LatexExportState, Latex, OpenQasm, CQasm};
    use crate::gates::{gate_test, Gate};
    use crate::stabilizer::PauliOp;

    #[test]
    fn test_description()
    {
        let gate = Swap::new();
        assert_eq!(gate.description(), "Swap");
    }

    #[test]
    fn test_cost()
    {
        let gate = Swap::new();
        assert_eq!(gate.cost(), 3.0 * 1001.0);
    }

    #[test]
    fn test_matrix()
    {
        let gate = Swap::new();
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z, z],
            [z, z, o, z],
            [z, o, z, z],
            [z, z, z, o]
        ]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * o;
        let mut state = array![
            [o, z, x,  x, x,  h, z],
            [z, o, x, -x, z, -h, z],
            [z, z, z,  z, x,  h, o],
            [z, z, z,  z, z, -h, z]
        ];
        let result = array![
            [o, z, x,  x, x,  h, z],
            [z, z, z,  z, x,  h, o],
            [z, o, x, -x, z, -h, z],
            [z, z, z,  z, z, -h, z]
        ];
        gate_test(Swap::new(), &mut state, &result);
    }

    #[test]
    fn test_apply_mat()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * o;
        let mut state = array![
            [o, z, x,  x, x,  h, z],
            [z, o, x, -x, z, -h, z],
            [z, z, z,  z, x,  h, o],
            [z, z, z,  z, z, -h, z]
        ];
        Swap::new().apply_mat(&mut state);
        let result = array![
            [o, z, x,  x, x,  h, z],
            [z, z, z,  z, x,  h, o],
            [z, o, x, -x, z, -h, z],
            [z, z, z,  z, z, -h, z]
        ];
        assert_complex_matrix_eq!(&state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = Swap::new().open_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("cx qb0, qb1; cx qb1, qb0; cx qb0, qb1")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = Swap::new().c_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("swap qb0, qb1")));
    }

    #[test]
    fn test_latex()
    {
        let gate = Swap::new();
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qswap \qwx[1] & \qw \\
    \lstick{\ket{0}} & \qswap & \qw \\
}
"#);

        let gate = Swap::new();
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[1, 0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qswap \qwx[1] & \qw \\
    \lstick{\ket{0}} & \qswap & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let mut ops = [PauliOp::I, PauliOp::X];
        assert_eq!(Swap::new().conjugate(&mut ops), Ok(false));
        assert_eq!(ops, [PauliOp::X, PauliOp::I]);
    }
}
