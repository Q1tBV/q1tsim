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
use error;
use export;

/// The Pauli X gate.
///
/// The X, or NOT, gate rotates the state over π radians around the `x` axis of
/// the Bloch sphere, i.e. it swaps the |0〉 and |1〉 components of the qubit.
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

    pub fn transform(state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        for i in 0..n
        {
            state.swap(i, i+n);
        }
    }

    pub fn transform_mat(state: &mut cmatrix::CMatSliceMut)
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

impl gates::Gate for X
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "X"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        array![[z, o], [o, z]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        Self::transform(state);
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        Self::transform_mat(state);
    }
}

impl export::OpenQasm for X
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("x {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for X
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("x {}", bit_names[bits[0]]))
    }
}

impl export::Latex for X
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
    {
        assert!(bits.len() == 1, "X gate operates on a single bit");
        let symbol = if state.is_controlled() { r"\targ" } else { r"\gate{X}" };
        state.set_field(bits[0], String::from(symbol));
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, X};
    use export::{OpenQasm, CQasm};
    use cmatrix;

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
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(x.matrix(), array![[z, o], [o, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
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
}
