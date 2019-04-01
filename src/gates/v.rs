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

/// The `V` gate
///
/// The `V` gate is the square root of the `X` gate. The associated matrix is
/// ```text
///     ┌         ┐
/// 1/2 │ 1+i 1-i │
///     │ 1-i 1+i │
///     └         ┘
/// ```
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

impl gates::Gate for V
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "V"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let h = 0.5 * cmatrix::COMPLEX_ONE;
        let hi = 0.5 * cmatrix::COMPLEX_I;
        array![[h+hi, h-hi], [h-hi, h+hi]]
    }
}

impl export::OpenQasm for V
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("u3(pi/2, -pi/2, pi/2) {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for V
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("x90 {}", bit_names[bits[0]]))
    }
}

impl export::Latex for V
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;
        state.add_block_gate(bits, "V")
    }
}

/// Conjugate of `V` gate.
///
/// The `V`<sup>`†`</sup> is the conjugate of the `V` gate.
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

impl gates::Gate for Vdg
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "V†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let h = 0.5 * cmatrix::COMPLEX_ONE;
        let hi = 0.5 * cmatrix::COMPLEX_I;
        array![[h-hi, h+hi], [h+hi, h-hi]]
    }
}

impl export::OpenQasm for Vdg
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("u3(pi/2, pi/2, -pi/2) {}", bit_names[bits[0]]))
    }
}

impl export::CQasm for Vdg
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        Ok(format!("mx90 {}", bit_names[bits[0]]))
    }
}

impl export::Latex for Vdg
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;
        state.add_block_gate(bits, r"V^\dagger")
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, V, Vdg};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

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
        let h = 0.5*cmatrix::COMPLEX_ONE;
        let hi = 0.5*cmatrix::COMPLEX_I;

        let gate = V::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[h+hi, h-hi], [h-hi, h+hi]]);

        let gate = Vdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[h-hi, h+hi], [h+hi, h-hi]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        let x = cmatrix::COMPLEX_HSQRT2;
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
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{V} & \qw \\
}
"#);

        let gate = Vdg::new();
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{V^\dagger} & \qw \\
}
"#);
    }
}
