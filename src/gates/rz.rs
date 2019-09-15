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

/// Rotation around `z` axis.
///
/// The `R`<sub>`Z`</sub>`(λ)` gate rotates the qubit around the `z` axis of the
/// Bloch sphere over an angle `λ`. It is equivalent to the `U`<sub>`1`</sub>
/// gate, up to an overall phase. The associated matrix is
/// ```text
/// ┌                     ┐
/// │ exp(-iλ/2)        0 │
/// │                     │
/// │          0 exp(iλ/2)│
/// └                     ┘
/// ```
#[derive(Clone)]
pub struct RZ
{
    lambda: crate::gates::Parameter,
    desc: String
}

impl RZ
{
    /// Create a new `R`<sub>`Z`</sub> gate with fixed angle `lambda`
    pub fn new<T>(lambda: T) -> Self
    where crate::gates::Parameter: From<T>
    {
        let param = crate::gates::Parameter::from(lambda);
        let desc = format!("RZ({:.4})", param);
        RZ { lambda: param, desc: desc }
    }
}

impl crate::gates::Gate for RZ
{
    fn cost(&self) -> f64
    {
        crate::gates::U1::cost()
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
        let z = crate::cmatrix::COMPLEX_ZERO;
        let p = num_complex::Complex::from_polar(&1.0, &(0.5 * self.lambda.value()));
        array![[p.conj(), z], [z, p]]
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let hlambda = 0.5 * self.lambda.value();
        {
            let mut slice = state.slice_mut(s![..n]);
            slice *= num_complex::Complex::from_polar(&1.0, &(-hlambda));
        }
        {
            let mut slice = state.slice_mut(s![n..]);
            slice *= num_complex::Complex::from_polar(&1.0, &( hlambda));
        }
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        let hlambda = 0.5 * self.lambda.value();
        {
            let mut slice = state.slice_mut(s![..n, ..]);
            slice *= num_complex::Complex::from_polar(&1.0, &(-hlambda));
        }
        {
            let mut slice = state.slice_mut(s![n.., ..]);
            slice *= num_complex::Complex::from_polar(&1.0, &( hlambda));
        }
    }
}

impl crate::export::OpenQasm for RZ
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("rz({}) {}", self.lambda, bit_names[bits[0]]))
    }
}

impl crate::export::CQasm for RZ
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        Ok(format!("rz {}, {}", bit_names[bits[0]], self.lambda))
    }
}

impl crate::export::Latex for RZ
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        let contents = format!("R_z({:.4})", self.lambda);
        state.add_block_gate(bits, &contents)
    }
}

impl crate::arithmetic::Square for RZ
{
    type SqType = Self;

    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        match self.lambda
        {
            crate::gates::Parameter::Direct(x) => Ok(Self::new(2.0 * x)),
            _                                  => Err(crate::error::Error::ReferenceArithmetic)
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::RZ;
    use crate::arithmetic::Square;
    use crate::gates::{gate_test, Gate};
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};

    #[test]
    fn test_description()
    {
        let gate = RZ::new(::std::f64::consts::FRAC_PI_4);
        assert_eq!(gate.description(), "RZ(0.7854)");
    }

    #[test]
    fn test_cost()
    {
        let gate = RZ::new(0.21675627161);
        assert_eq!(gate.cost(), 7.0);
    }

    #[test]
    fn test_matrix()
    {
        let gate = RZ::new(::std::f64::consts::PI);
        let z = crate::cmatrix::COMPLEX_ZERO;
        let i = crate::cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(gate.matrix(), array![[-i, z], [z, i]]);
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
            [x*(o-i), z,       0.5*(o-i),  0.5*(o-i)],
            [z      , x*(o+i), 0.5*(o+i), -0.5*(o+i)]
        ];
        let gate = RZ::new(::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = RZ::new(2.25).open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("rz(2.25) qb")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = RZ::new(2.25).c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::from("rz qb, 2.25")));
    }

    #[test]
    fn test_latex()
    {
        let gate = RZ::new(::std::f64::consts::FRAC_PI_2);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{R_z(1.5708)} & \qw \\
}
"#);

        let gate = RZ::new(-24.0);
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(gate.latex(&[0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{R_z(-24.0000)} & \qw \\
}
"#);
    }

    #[test]
    fn test_square()
    {
        let gate = RZ::new(0.0);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);

        let gate = RZ::new(1.3);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);

        let gate = RZ::new(-2.5);
        let mat = gate.matrix();
        let sq_mat = mat.dot(&mat);
        assert_complex_matrix_eq!(gate.square().unwrap().matrix(), &sq_mat);
    }
}
