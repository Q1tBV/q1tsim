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

/// The Pauli Z gate.
///
/// The Z gate rotates the state over π radians around the `z` axis of
/// the Bloch sphere, i.e. it flips the sign of the |1〉 components of the qubit.
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

impl gates::Gate for Z
{
    fn cost(&self) -> f64
    {
        gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "Z"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        array![[o, z], [z, -o]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        state.slice_mut(s![n..]).mapv_inplace(|c| -c);
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        state.slice_mut(s![n.., ..]).mapv_inplace(|c| -c);
    }
}

impl export::OpenQasm for Z
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("z {}", bit_names[bits[0]])
    }
}

impl export::CQasm for Z
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("z {}", bit_names[bits[0]])
    }
}

impl export::Latex for Z
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
    {
        assert!(bits.len() == 1, "Z gate operates on a single bit");
        let symbol = if state.is_controlled() { r"\control \qw" } else { r"\gate{Z}" };
        state.set_field(bits[0], String::from(symbol));
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, Z};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

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
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, -o]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [ z, -o, -x, x]];
        gate_test(Z::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Z::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "z qb");
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Z::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "z qb");
    }

    #[test]
    fn test_latex()
    {
        let gate = Z::new();
        let mut state = LatexExportState::new(1, 0);
        gate.latex_checked(&[0], &mut state);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{Z} & \qw \\
}
"#);
    }
}
