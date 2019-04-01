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
use error;
use gates;
use export;

use gates::Gate;

/// Gate describing the Kronecker product of two other gates operating on
/// different qubits.
pub struct Kron<G0, G1>
{
    g0: G0,
    g1: G1,
    desc: String
}

impl<G0, G1> Kron<G0, G1>
where G0: gates::Gate, G1: gates::Gate
{
    /// Create a new Kronecker product gate `g1` ⊗ `g2`.
    pub fn new(g0: G0, g1: G1) -> Self
    {
        let desc = format!("{}⊗{}", g0.description(), g1.description());
        Kron { g0: g0, g1: g1, desc: desc }
    }
}

impl<G0, G1> gates::Gate for Kron<G0, G1>
where G0: gates::Gate, G1: gates::Gate
{
    fn cost(&self) -> f64
    {
        self.g0.cost() + self.g1.cost()
    }

    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        self.g0.nr_affected_bits() + self.g1.nr_affected_bits()
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::kron_mat(&self.g0.matrix(), &self.g1.matrix())
    }

    fn apply_slice(&self, mut state: cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 4 == 0, "Number of rows is not a multiple of four.");

        let n = state.len() / 2;

        self.g0.apply_slice(state.view_mut());
        self.g1.apply_slice(state.slice_mut(s![..n]));
        self.g1.apply_slice(state.slice_mut(s![n..]));
    }
}

impl<G0, G1> export::OpenQasm for Kron<G0, G1>
where G0: export::OpenQasm, G1: export::OpenQasm
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        let n0 = self.g0.nr_affected_bits();
        let op0 = self.g0.open_qasm(bit_names, &bits[..n0])?;
        let op1 = self.g1.open_qasm(bit_names, &bits[n0..])?;
        Ok(format!("{}; {}", op0, op1))
    }

    fn conditional_open_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> error::Result<String>
    {
        let n0 = self.g0.nr_affected_bits();
        let op0 = self.g0.conditional_open_qasm(condition, bit_names, &bits[..n0])?;
        let op1 = self.g1.conditional_open_qasm(condition, bit_names, &bits[n0..])?;
        Ok(op0 + "; " + &op1)
    }
}

impl<G0, G1> export::CQasm for Kron<G0, G1>
where G0: export::CQasm, G1: export::CQasm
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        let n0 = self.g0.nr_affected_bits();
        let op0 = self.g0.c_qasm(bit_names, &bits[..n0])?;
        let op1 = self.g1.c_qasm(bit_names, &bits[n0..])?;
        Ok(format!("{{ {} | {} }}", op0, op1))
    }

    fn conditional_c_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> error::Result<String>
    {
        let n0 = self.g0.nr_affected_bits();
        let op0 = self.g0.conditional_c_qasm(condition, bit_names, &bits[..n0])?;
        let op1 = self.g1.conditional_c_qasm(condition, bit_names, &bits[n0..])?;
        Ok(op0 + "\n" + &op1)
    }
}

impl<G0, G1> export::Latex for Kron<G0, G1>
where G0: export::Latex, G1: export::Latex
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;

        let n0 = self.g0.nr_affected_bits();
        self.g0.latex_checked(&bits[..n0], state)?;
        self.g1.latex_checked(&bits[n0..], state)
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, CX, Gate, H, I, Kron, X};
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use cmatrix;

    #[test]
    fn test_cost()
    {
        let gate = Kron::new(X::new(), H::new());
        assert_eq!(gate.cost(), 201.0 + 104.0);

        let gate = Kron::new(Kron::new(X::new(), H::new()), CX::new());
        assert_eq!(gate.cost(), 201.0 + 104.0 + 1001.0);
    }

    #[test]
    fn test_description()
    {
        let ih = Kron::new(I::new(), H::new());
        assert_eq!(ih.description(), "I⊗H");
        let hh = Kron::new(H::new(), H::new());
        assert_eq!(hh.description(), "H⊗H");
        let hih = Kron::new(H::new(), Kron::new(I::new(), H::new()));
        assert_eq!(hih.description(), "H⊗I⊗H");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let s = cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * cmatrix::COMPLEX_ONE;

        let ih = Kron::new(I::new(), H::new());
        assert_complex_matrix_eq!(ih.matrix(), array![
            [s,  s, z,  z],
            [s, -s, z,  z],
            [z,  z, s,  s],
            [z,  z, s, -s]
        ]);

        let hh = Kron::new(H::new(), H::new());
        assert_complex_matrix_eq!(hh.matrix(), array![
            [h,  h,  h,  h],
            [h, -h,  h, -h],
            [h,  h, -h, -h],
            [h, -h, -h,  h]
        ]);

        let hih = Kron::new(H::new(), Kron::new(I::new(), H::new()));
        assert_complex_matrix_eq!(hih.matrix(), array![
            [h,  h,  z,  z,  h,  h,  z,  z],
            [h, -h,  z,  z,  h, -h,  z,  z],
            [z,  z,  h,  h,  z,  z,  h,  h],
            [z,  z,  h, -h,  z,  z,  h, -h],
            [h,  h,  z,  z, -h, -h,  z,  z],
            [h, -h,  z,  z, -h,  h,  z,  z],
            [z,  z,  h,  h,  z,  z, -h, -h],
            [z,  z,  h, -h,  z,  z, -h,  h]
        ]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let h = o * 0.5;

        let mut state = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, o,  h,  x],
            [z, z, -h, -x]
        ];
        let result = array![
            [x, z, z, z],
            [x, z, x, z],
            [z, x, z, z],
            [z, x, x, o]
        ];
        gate_test(Kron::new(I::new(), H::new()), &mut state, &result);

        let mut state = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, o,  h,  x],
            [z, z, -h, -x]
        ];
        let result = array![
            [x,  x,  x,  h],
            [z,  z, -x, -h],
            [x, -x,  z, -h],
            [z,  z,  z,  h]
        ];
        gate_test(Kron::new(H::new(), I::new()), &mut state, &result);

        let mut state = array![
            [o, z,  h,  z],
            [z, z, -h,  z],
            [z, o,  h,  x],
            [z, z, -h, -x]
        ];
        let result = array![
            [h,  h, z,  z],
            [h,  h, o,  x],
            [h, -h, z,  z],
            [h, -h, z, -x]
        ];
        gate_test(Kron::new(H::new(), H::new()), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = Kron::new(H::new(), I::new()).open_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("h qb0; id qb1")));
    }

    #[test]
    fn test_conditional_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = Kron::new(H::new(), X::new())
            .conditional_open_qasm("b == 1", &bit_names, &[1, 0]);
        assert_eq!(qasm, Ok(String::from("if (b == 1) h qb1; if (b == 1) x qb0")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = Kron::new(H::new(), I::new()).c_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("{ h qb0 | i qb1 }")));
    }

    #[test]
    fn test_conditional_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = Kron::new(H::new(), X::new())
            .conditional_c_qasm("b == 1", &bit_names, &[1, 0]);
        assert_eq!(qasm, Ok(String::from(
r#"c-h b == 1, qb1
c-x b == 1, qb0"#)));
    }

    #[test]
    fn test_latex()
    {
        let gate = Kron::new(H::new(), X::new());
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{H} & \qw \\
    \lstick{\ket{0}} & \gate{X} & \qw \\
}
"#);

        let gate = Kron::new(H::new(), X::new());
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[1, 0], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{X} & \qw \\
    \lstick{\ket{0}} & \gate{H} & \qw \\
}
"#);

        let gate = Kron::new(CX::new(), X::new());
        let mut state = LatexExportState::new(3, 0);
        assert_eq!(gate.latex(&[0, 2, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \ctrl{2} & \qw & \qw \\
    \lstick{\ket{0}} & \qw & \gate{X} & \qw \\
    \lstick{\ket{0}} & \targ & \qw & \qw \\
}
"#);
    }
}
