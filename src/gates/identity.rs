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
use export;

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

impl gates::Gate for I
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

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::CMatrix::eye(2)
    }

    fn apply_slice(&self, _state: &mut cmatrix::CVecSliceMut)
    {
        // Identity, leave state unchanged, so do nothing
    }
}

impl export::OpenQasm for I
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("id {}", bit_names[bits[0]])
    }
}

impl export::CQasm for I
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("i {}", bit_names[bits[0]])
    }
}

impl export::Latex for I
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
    {
        assert!(bits.len() == 1, "I gate operates on a single bit");
        state.set_field(bits[0], String::from(r"\qw"));
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, I};
    use export::{OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let i = I::new();
        assert_eq!(i.description(), "I");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = I::new();
        assert_complex_matrix_eq!(i.matrix(), array![[o, z], [z, o]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [z, o, x, -x]];
        gate_test(I::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = I::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "id qb");
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = I::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "i qb");
    }
}
