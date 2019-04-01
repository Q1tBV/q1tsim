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

use gates::Gate;

/// The Pauli Y gate.
///
/// The Y gate rotates the state over π radians around the `y` axis of the Bloch
/// sphere.
pub struct Y
{
}

impl Y
{
    /// Create a new Pauli Y gate.
    pub fn new() -> Self
    {
        Y { }
    }
}

impl gates::Gate for Y
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "Y"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let i = cmatrix::COMPLEX_I;
        array![[z, -i], [i, z]]
    }

    fn apply_slice(&self, mut state: cmatrix::CVecSliceMut)
    {
        gates::X::transform(state.view_mut());
        let n = state.len() / 2;
        {
            let mut slice = state.slice_mut(s![..n]);
            slice *= -cmatrix::COMPLEX_I;
        }
        {
            let mut slice = state.slice_mut(s![n..]);
            slice *=  cmatrix::COMPLEX_I;
        }
    }

    fn apply_mat_slice(&self, mut state: cmatrix::CMatSliceMut)
    {
        gates::X::transform_mat(state.view_mut());
        let n = state.rows() / 2;
        {
            let mut slice = state.slice_mut(s![..n, ..]);
            slice *= -cmatrix::COMPLEX_I;
        }
        {
            let mut slice = state.slice_mut(s![n.., ..]);
            slice *=  cmatrix::COMPLEX_I;
        }
    }
}

impl export::OpenQasm for Y
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("y {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for Y
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("y {}", bit_names[bits[0]]))
    }
}

impl export::Latex for Y
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;
        state.add_block_gate(bits, "Y")
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, Y};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = Y::new();
        assert_eq!(gate.description(), "Y");
    }

    #[test]
    fn test_cost()
    {
        let gate = Y::new();
        assert_eq!(gate.cost(), 201.0);
    }

    #[test]
    fn test_matrix()
    {
        let y = Y::new();
        let z = cmatrix::COMPLEX_ZERO;
        let i = cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(y.matrix(), array![[z, -i], [i, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[z, -i, -i*x, i*x], [i, z, i*x, i*x]];
        gate_test(Y::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Y::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("y qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = Y::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("y qb")));
    }

    #[test]
    fn test_latex()
    {
        let gate = Y::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{Y} & \qw \\
}
"#);
    }
}
