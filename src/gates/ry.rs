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

/// Rotation around `y` axis.
///
/// The `R`<sub>`Y`</sub>`(θ)` gate rotates the qubit around the `y` axis of the
/// Bloch sphere over an angle `θ`. The associated matrix is
/// ```text
/// ┌                    ┐
/// │ cos(θ/2) -sin(θ/2) │
/// │                    │
/// │ sin(θ/2)  cos(θ/2) │
/// └                    ┘
/// ```
pub struct RY
{
    theta: f64,
    desc: String
}

impl RY
{
    /// Create a new `R`<sub>`Y`</sub> gate.
    pub fn new(theta: f64) -> Self
    {
        RY { theta: theta, desc: format!("RY({:.4})", theta) }
    }
}

impl gates::Gate for RY
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
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
        let c = num_complex::Complex::new((0.5 * self.theta).cos(), 0.0);
        let s = num_complex::Complex::new((0.5 * self.theta).sin(), 0.0);
        array![[c, -s], [s, c]]
    }

    fn apply_slice(&self, mut state: cmatrix::CVecSliceMut)
    {
        let cos_t = num_complex::Complex::new((0.5 * self.theta).cos(), 0.0);
        let sin_t = num_complex::Complex::new((0.5 * self.theta).sin(), 0.0);

        let mut s = state.to_owned();
        s *= sin_t;
        state *= cos_t;

        let n = state.len() / 2;
        {
            let mut slice = state.slice_mut(s![..n]);
            slice -= &s.slice(s![n..]);
        }
        {
            let mut slice = state.slice_mut(s![n..]);
            slice += &s.slice(s![..n]);
        }
    }

    fn apply_mat_slice(&self, mut state: cmatrix::CMatSliceMut)
    {
        let cos_t = num_complex::Complex::new((0.5 * self.theta).cos(), 0.0);
        let sin_t = num_complex::Complex::new((0.5 * self.theta).sin(), 0.0);

        let mut s = state.to_owned();
        s *= sin_t;
        state *= cos_t;

        let n = state.rows() / 2;
        {
            let mut slice = state.slice_mut(s![..n, ..]);
            slice -= &s.slice(s![n.., ..]);
        }
        {
            let mut slice = state.slice_mut(s![n.., ..]);
            slice += &s.slice(s![..n, ..]);
        }
    }
}

impl export::OpenQasm for RY
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        // For some reason, the web interface on QX claims RY is not a defined
        // gate, even though it is defined in the specification. Replace by U3.
        //format!("ry({}) {}", self.theta, bit_names[bits[0]])
        Ok(format!("u3({}, 0, 0) {}", self.theta, bit_names[bits[0]]))
    }
}

impl export::CQasm for RY
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
    {
        Ok(format!("ry {}, {}", bit_names[bits[0]], self.theta))
    }
}

impl export::Latex for RY
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;

        let contents = format!(r"\gate{{R_y({:.4})}}", self.theta);
        state.set_field(bits[0], contents);
        Ok(())
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, RY};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = RY::new(0.21675627161);
        assert_eq!(gate.description(), "RY(0.2168)");
    }

    #[test]
    fn test_cost()
    {
        let gate = RY::new(0.21675627161);
        assert_eq!(gate.cost(), 201.0);
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let gate = RY::new(::std::f64::consts::FRAC_PI_2);
        assert_complex_matrix_eq!(gate.matrix(), array![[x, -x], [x, x]]);

        let gate = RY::new(::std::f64::consts::PI);
        assert_complex_matrix_eq!(gate.matrix(), array![[z, -o], [o, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [x, -x, z, o],
            [x,  x, o, z]
        ];
        let gate = RY::new(::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = RY::new(2.25).open_qasm(&bit_names, &[0]);
        //assert_eq!(qasm, "ry(2.25) qb");
        assert_eq!(qasm, Ok(String::from("u3(2.25, 0, 0) qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = RY::new(2.25).c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("ry qb, 2.25")));
    }

    #[test]
    fn test_latex()
    {
        let gate = RY::new(::std::f64::consts::FRAC_PI_2);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{R_y(1.5708)} & \qw \\
}
"#);

        let gate = RY::new(-24.0);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex_checked(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{R_y(-24.0000)} & \qw \\
}
"#);
    }
}
