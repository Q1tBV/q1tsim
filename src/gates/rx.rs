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

/// Rotation around `x` axis.
///
/// The `R`<sub>`X`</sub>`(θ)` gate rotates the qubit around the `x` axis of the
/// Bloch sphere over an angle `θ`. The associated matrix is
/// ```text
/// ┌                     ┐
/// │  cos(θ/2) -isin(θ/2)│
/// │                     │
/// │-isin(θ/2)   cos(θ/2)│
/// └                     ┘
/// ```
#[derive(Clone)]
pub struct RX
{
    theta: crate::gates::Parameter,
    desc: String
}

impl RX
{
    /// Create a new `R`<sub>`X`</sub> gate.
    pub fn new<T>(theta: T) -> Self
    where crate::gates::Parameter: From<T>
    {
        let param = crate::gates::Parameter::from(theta);
        let desc = format!("RX({:.4})", param);
        RX { theta: param, desc: desc }
    }
}

impl crate::gates::Gate for RX
{
    fn cost(&self) -> f64
    {
        crate::gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let htheta = 0.5 * self.theta.value();
        let c  = num_complex::Complex::new((htheta).cos(), 0.0);
        let si = num_complex::Complex::new(0.0, (htheta).sin());
        array![[c, -si], [-si, c]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        let htheta = 0.5 * self.theta.value();
        let c  = num_complex::Complex::new((htheta).cos(), 0.0);
        let si = num_complex::Complex::new(0.0, (htheta).sin());

        let mut s = state.to_owned();
        s *= si;
        state *= c;

        let n = state.len() / 2;
        {
            let mut slice = state.slice_mut(s![..n]);
            slice -= &s.slice(s![n..]);
        }
        {
            let mut slice = state.slice_mut(s![n..]);
            slice -= &s.slice(s![..n]);
        }
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        let htheta = 0.5 * self.theta.value();
        let cos_t   = num_complex::Complex::new((htheta).cos(), 0.0);
        let sin_t_i = num_complex::Complex::new(0.0, (htheta).sin());

        let mut s = state.to_owned();
        s *= sin_t_i;
        state *= cos_t;

        let n = state.rows() / 2;
        {
            let mut slice = state.slice_mut(s![..n, ..]);
            slice -= &s.slice(s![n.., ..]);
        }
        {
            let mut slice = state.slice_mut(s![n.., ..]);
            slice -= &s.slice(s![..n, ..]);
        }
    }
}

impl crate::export::OpenQasm for RX
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("rx({}) {}", self.theta, bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for RX
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("rx {}, {}", bit_names[bits[0]], self.theta))
    }
}

impl crate::export::Latex for RX
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        let contents = format!("R_x({:.4})", self.theta);
        state.add_block_gate(bits, &contents)
    }
}

impl crate::arithmetic::Square for RX
{
    type SqType = Self;

    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        match self.theta
        {
            crate::gates::Parameter::Direct(x) => Ok(Self::new(2.0 * x)),
            _                                  => Err(crate::error::Error::ReferenceArithmetic)
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::RX;
    use crate::arithmetic::Square;
    use crate::gates::{gate_test, Gate};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};

    #[test]
    fn test_description()
    {
        let gate = RX::new(::std::f64::consts::PI);
        assert_eq!(gate.description(), "RX(3.1416)");
    }

    #[test]
    fn test_cost()
    {
        let gate = RX::new(::std::f64::consts::PI);
        assert_eq!(gate.cost(), 201.0);
    }

    #[test]
    fn test_matrix()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;

        let gate = RX::new(::std::f64::consts::FRAC_PI_2);
        assert_complex_matrix_eq!(gate.matrix(), array![[x, -x*i], [-x*i, x]]);

        let gate = RX::new(::std::f64::consts::PI);
        assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [-i, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;
        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [   x, -x*i, 0.5*(o-i),  0.5*(o+i)],
            [-x*i,    x, 0.5*(o-i), -0.5*(o+i)]
        ];
        let gate = RX::new(::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = RX::new(2.25).open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("rx(2.25) qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = RX::new(2.25).c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("rx qb, 2.25")));
    }

    #[test]
    fn test_latex()
    {
        let gate = RX::new(::std::f64::consts::FRAC_PI_2);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{R_x(1.5708)} & \qw \\
}
"#);

        let gate = RX::new(-24.0);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{R_x(-24.0000)} & \qw \\
}
"#);
    }

    #[test]
    fn test_square()
    {
        let gate = RX::new(0.0);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);

        let gate = RX::new(1.3);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);

        let gate = RX::new(-2.5);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);
    }
}
