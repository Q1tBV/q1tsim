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

/// Phase gate.
///
/// The `U`<sub>`1`</sub>`(λ)` gate shifts the phase of the |1〉 component of a
/// qubit over an angle `λ`. The associated matrix is
/// ```text
/// ┌              ┐
/// │ 1          0 │
/// │              │
/// │ 0    exp(iλ) │
/// └              ┘
/// ```

pub struct U1
{
    lambda: f64,
    desc: String
}

impl U1
{
    /// Create a new `U`<sub>`1`</sub> gate.
    pub fn new(lambda: f64) -> Self
    {
        U1 { lambda: lambda, desc: format!("U1({:.4})", lambda) }
    }

    pub fn cost() -> f64
    {
        7.0
    }
}

impl gates::Gate for U1
{
    fn cost(&self) -> f64
    {
        Self::cost()
    }

    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let p = num_complex::Complex::from_polar(&1.0, &self.lambda);
        array![[o, z], [z, p]]
    }

    fn apply_slice(&self, mut state: cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= num_complex::Complex::from_polar(&1.0, &self.lambda);
    }
}

impl export::OpenQasm for U1
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("u1({}) {}", self.lambda, bit_names[bits[0]]))
    }
}

impl export::CQasm for U1
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        // U1 is R_Z up to a phase
        Ok(format!("rz {}, {}", bit_names[bits[0]], self.lambda))
    }
}

impl export::Latex for U1
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;
        let contents = format!("U_1({:.4})", self.lambda);
        state.add_block_gate(bits, &contents)
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, U1};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        assert_eq!(gate.description(), "U1(0.7854)");
    }

    #[test]
    fn test_cost()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        assert_eq!(gate.cost(), 7.0);
    }

    #[test]
    fn test_matrix()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_2);
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, i]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [z, x*(o+i), 0.5*(o+i), -0.5*(o+i)]];
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        gate_test(gate, &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = U1::new(::std::f64::consts::PI).open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("u1(3.141592653589793) qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = U1::new(::std::f64::consts::PI).c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("rz qb, 3.141592653589793")));
    }

    #[test]
    fn test_latex()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{U_1(0.7854)} & \qw \\
}
"#);

        let gate = U1::new(-1.2);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{U_1(-1.2000)} & \qw \\
}
"#);
    }
}
