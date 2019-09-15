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

/// U<sub>2</sub> gate.
///
/// The `U`<sub>`2`</sub>`(ϕ, λ)` gate transforms a qubit by the matrix
/// ```text
///           ┌                      ┐
///           │      1       -exp(iλ)│
/// 1/sqrt(2) │                      │
///           │exp(iϕ)    exp(i(λ+ϕ))│
///           └                      ┘
/// ```
#[derive(Clone)]
pub struct U2
{
    phi: crate::gates::Parameter,
    lambda: crate::gates::Parameter,
    desc: String
}

impl U2
{
    /// Create a new `U`<sub>`2`</sub> gate.
    pub fn new<Tp, Tl>(phi: Tp, lambda: Tl) -> Self
    where crate::gates::Parameter: From<Tp> + From<Tl>
    {
        let pphi = crate::gates::Parameter::from(phi);
        let plambda = crate::gates::Parameter::from(lambda);
        let desc = format!("U2({:.4}, {:.4})", pphi, plambda);
        U2 { phi: pphi, lambda: plambda, desc: desc }
    }

    pub fn cost() -> f64
    {
        104.0
    }
}

impl crate::gates::Gate for U2
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

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let lambda = self.lambda.value();
        let phi = self.phi.value();
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        array![[ num_complex::Complex::new(x, 0.0),
                -num_complex::Complex::from_polar(&x, &lambda)],
               [ num_complex::Complex::from_polar(&x, &phi),
                 num_complex::Complex::from_polar(&x, &(phi+lambda))]]
    }
}

impl crate::export::OpenQasm for U2
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("u2({}, {}) {}", self.phi, self.lambda, bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for U2
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        let name = &bit_names[bits[0]];
        Ok(format!("rz {}, {}\nh {}\nrz {} {}", name,
            self.lambda.value() + ::std::f64::consts::PI, name, name, self.phi))
    }
}

impl crate::export::Latex for U2
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        let contents = format!("U_2({:.4}, {:.4})", self.phi, self.lambda);
        state.add_block_gate(bits, &contents)
    }
}

impl crate::arithmetic::Square for U2
{
    type SqType = crate::gates::U3;

    // NOTE: the U3 given below as the square of U2 is exact up to a global
    // phase exp(i(lambda+phi-pi)/2).
    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        match (&self.phi, &self.lambda)
        {
            (&crate::gates::Parameter::Direct(p), &crate::gates::Parameter::Direct(l)) => {
                let pi = ::std::f64::consts::PI;
                Ok(crate::gates::U3::new(l+p-pi, p-0.5*pi, l-0.5*pi))
            },
            _ => Err(crate::error::Error::ReferenceArithmetic)
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::U2;
    use crate::arithmetic::Square;
    use crate::gates::{gate_test, Gate};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use num_complex::Complex;

    #[test]
    fn test_description()
    {
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_eq!(gate.description(), "U2(0.7854, 0.6931)");
    }

    #[test]
    fn test_cost()
    {
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_eq!(gate.cost(), 104.0);
    }

    #[test]
    fn test_matrix()
    {
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [Complex::new(0.7071067811865476, 0.0), Complex::new(-0.5439340435069544, -0.4518138513969824)],
            [Complex::new(               0.5, 0.5), Complex::new(0.06513881252516862,  0.7041000888388035)]
        ]);
    }

    #[test]
    fn test_apply()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;
        let mut state = array![[x, x], [z, -x], [x, z], [z, z]];
        let result = array![
            [    o-i,        o],
            [      z,       -o],
            [2.0*i*x,  (o+i)*x],
            [      z, -(o+i)*x]
        ] * (0.5 * o);
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = U2::new(1.0, 2.25).open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("u2(1, 2.25) qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = U2::new(1.0, 2.25).c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("rz qb, 5.391592653589793\nh qb\nrz qb 1")));
    }

    #[test]
    fn test_latex()
    {
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{U_2(0.7854, 0.6931)} & \qw \\
}
"#);

        let gate = U2::new(-1.2, ::std::f64::consts::LN_2);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{U_2(-1.2000, 0.6931)} & \qw \\
}
"#);

        let gate = U2::new(::std::f64::consts::FRAC_PI_4, -3.14);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{U_2(0.7854, -3.1400)} & \qw \\
}
"#);
    }

    #[test]
    fn test_square()
    {
        let phi = 1.2;
        let lambda = -2.3;
        let phase = num_complex::Complex::from_polar(&1.0,
            &(0.5*(lambda+phi-::std::f64::consts::PI)));
        let gate = U2::new(phi, lambda);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(phase * gate.square().unwrap().matrix(), &sq_mat);

        let phi = 2.4;
        let lambda = 0.0;
        let phase = num_complex::Complex::from_polar(&1.0,
            &(0.5*(lambda+phi-::std::f64::consts::PI)));
        let gate = U2::new(phi, lambda);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(phase * gate.square().unwrap().matrix(), &sq_mat);

        let phi = 0.0;
        let lambda = -2.3;
        let phase = num_complex::Complex::from_polar(&1.0,
            &(0.5*(lambda+phi-::std::f64::consts::PI)));
        let gate = U2::new(phi, lambda);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(phase * gate.square().unwrap().matrix(), &sq_mat);

        let phi = 2.28;
        let lambda = 17.12;
        let phase = num_complex::Complex::from_polar(&1.0,
            &(0.5*(lambda+phi-::std::f64::consts::PI)));
        let gate = U2::new(phi, lambda);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(phase * gate.square().unwrap().matrix(), &sq_mat);
    }
}
