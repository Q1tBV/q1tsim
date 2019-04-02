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

/// The identity gate
///
/// The identity gate leaves the qubits on which it acts unchanged.
pub struct I
{
}

impl I
{
    /// Create a new identity gate.
    pub fn new() -> Self
    {
        I { }
    }
}

impl crate::gates::Gate for I
{
    fn cost(&self) -> f64
    {
        0.0
    }

    fn description(&self) -> &str
    {
        "I"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        crate::cmatrix::CMatrix::eye(2)
    }

    fn apply_slice(&self, _state: crate::cmatrix::CVecSliceMut)
    {
        // Identity, leave state unchanged, so do nothing
    }
}

impl crate::export::OpenQasm for I
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("id {}", bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for I
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("i {}", bit_names[bits[0]]))
    }
}

impl crate::export::Latex for I
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits)?;
        state.set_field(bits[0], String::from(r"\qw"))
    }
}

#[cfg(test)]
mod tests
{
    use crate::gates::{gate_test, Gate, I};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};

    #[test]
    fn test_description()
    {
        let gate = I::new();
        assert_eq!(gate.description(), "I");
    }

    #[test]
    fn test_nr_affected_bits()
    {
        let gate = I::new();
        assert_eq!(gate.nr_affected_bits(), 1);
    }

    #[test]
    fn test_cost()
    {
        let gate = I::new();
        assert_eq!(gate.cost(), 0.0);
    }

    #[test]
    fn test_matrix()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = I::new();
        assert_complex_matrix_eq!(i.matrix(), array![[o, z], [z, o]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;

        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [z, o, x, -x]];
        gate_test(I::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = I::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("id qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = I::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("i qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = I::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qw & \qw \\
}
"#);
    }
}
