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


use cmatrix;
use error;
use gates;
use export;

use gates::Gate;

/// Static loop gate
///
/// The `Loop` gate represents a static loop, i.e. a set of instructions that
/// is executed a fixed number of times.
pub struct Loop
{
    /// The loop label, naming the loop
    label: String,
    /// The number of times to execute the loop body
    nr_iterations: usize,
    /// The instructions to loop
    body: gates::Composite,
    /// A descriptions string, describing the loop
    desc: String
}

impl Loop
{
    /// Create a new static loop.
    ///
    /// Initialize a new static loop executing the instructions in `body`,
    /// `nr_iterations` times.
    pub fn new(label: &str, nr_iterations: usize, body: gates::Composite) -> Self
    {
        let desc = format!("{}({})", nr_iterations, body.description());
        Loop
        {
            label: String::from(label),
            nr_iterations: nr_iterations,
            body: body,
            desc: desc
        }
    }
}

impl gates::Gate for Loop
{
    fn cost(&self) -> f64
    {
        self.nr_iterations as f64 * self.body.cost()
    }

    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        self.body.nr_affected_bits()
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let mut res = cmatrix::CMatrix::eye(1 << self.nr_affected_bits());
        self.apply_mat_slice(res.view_mut());
        res
    }

    fn apply_slice(&self, mut state: cmatrix::CVecSliceMut)
    {
        for _ in 0..self.nr_iterations
        {
            self.body.apply_slice(state.view_mut());
        }
    }

    fn apply_mat_slice(&self, mut state: cmatrix::CMatSliceMut)
    {
        for _ in 0..self.nr_iterations
        {
            self.body.apply_mat_slice(state.view_mut());
        }
    }
}

impl export::OpenQasm for Loop
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        if self.nr_iterations == 0
        {
            Ok(String::new())
        }
        else
        {
            let qasm_body = self.body.open_qasm(bit_names, bits)?;
            let mut res = qasm_body.clone();
            for _ in 1..self.nr_iterations
            {
                res += ";\n";
                res += &qasm_body;
            }
            Ok(res)
        }
    }

    fn conditional_open_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> error::Result<String>
    {
        if self.nr_iterations == 0
        {
            Ok(String::new())
        }
        else
        {
            let qasm_body = self.body.conditional_open_qasm(condition,
                bit_names, bits)?;
            let mut res = qasm_body.clone();
            for _ in 1..self.nr_iterations
            {
                res += ";\n";
                res += &qasm_body;
            }
            Ok(res)
        }
    }
}

impl export::CQasm for Loop
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::Result<String>
    {
        let body_qasm = self.body.c_qasm(bit_names, bits)?;
        Ok(format!(".{}({})\n{}\n.end", self.label, self.nr_iterations,
            body_qasm))
    }

    fn conditional_c_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> error::Result<String>
    {
        if self.nr_iterations == 0
        {
            Ok(String::new())
        }
        else
        {
            // We can't create a classically controlled loop in c-Qasm. We can
            // repeat the classically controlled loop body, however.
            let qasm_body = self.body.conditional_c_qasm(condition,
                bit_names, bits)?;
            let mut res = qasm_body.clone();
            for _ in 1..self.nr_iterations
            {
                res += "\n";
                res += &qasm_body;
            }
            Ok(res)
        }
    }
}

impl export::Latex for Loop
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
        -> error::Result<()>
    {
        self.check_nr_bits(bits)?;

        if self.nr_iterations == 0
        {
            Ok(())
        }
        else if self.nr_iterations == 1
        {
            self.body.latex(bits, state)
        }
        else if self.nr_iterations == 2
        {
            self.body.latex(bits, state)?;
            self.body.latex_checked(bits, state)
        }
        else
        {
            let min = *bits.iter().min().unwrap();
            let max = *bits.iter().max().unwrap();

            state.start_loop(self.nr_iterations);
            self.body.latex(bits, state)?;
            state.add_cds(min, max - min, r"\cdots");
            self.body.latex_checked(bits, state)?;
            state.end_loop()
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::Loop;
    use gates::{gate_test, Composite, Gate};
    use export::{CQasm, OpenQasm, Latex, LatexExportState};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let body = Composite::from_string("body", "RX(1.0471975511965976) 0").unwrap();
        let gate = Loop::new("myloop", 3, body);
        assert_eq!(gate.description(), "3(body)");
    }

    #[test]
    fn test_cost()
    {
        let body = Composite::from_string("body", "H 0").unwrap();
        let count = 15;
        let gate = Loop::new("myloop", count, body);
        assert_eq!(gate.cost(), count as f64 * 104.0);

        let body = Composite::from_string("body", "H 1; CX 0 1; H 1").unwrap();
        let count = 9;
        let gate = Loop::new("myloop", count, body);
        assert_eq!(gate.cost(), count as f64 * 1209.0);
    }

    #[test]
    fn test_matrix()
    {
        let body = Composite::from_string("body", "RX(1.0471975511965976) 0").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let z = cmatrix::COMPLEX_ZERO;
        let i = cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [-i, z]]);
    }

    #[test]
    fn test_apply()
    {
        let body = Composite::from_string("body", "RX(1.0471975511965976) 0; RX(-1.0471975511965976) 1").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = array![
            [o, z, x, x],
            [z, o, z, x],
            [z, z, z, z],
            [z, z, x, z],
        ];
        let result = array![
            [ z, z, x, z ],
            [ z, z, z, z ],
            [ z, o, z, x ],
            [ o, z, x, x ]
        ];
        gate_test(gate, &mut state, &result);
    }

    #[test]
    fn test_apply_mat()
    {
        let body = Composite::from_string("body", "Sdg 0; I 1").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut state = array![
            [o, z, x, z],
            [z, z, z, z],
            [z, o, z, x],
            [z, z, x, x],
        ];
        let result = array![
            [ o, z,   x,   z ],
            [ z, z,   z,   z ],
            [ z, i,   z, i*x ],
            [ z, z, i*x, i*x ]
        ];
        gate.apply_mat(&mut state);
        assert_complex_matrix_eq!(&state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let body = Composite::from_string("body", "H 0").unwrap();
        let gate = Loop::new("myloop", 0, body);
        let bit_names = [String::from("qb0")];
        let qasm = gate.open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, Ok(String::new()));

        let body = Composite::from_string("body", "H 0; H 1; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = gate.open_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from(
r#"h qb0; h qb1; cx qb0, qb1;
h qb0; h qb1; cx qb0, qb1;
h qb0; h qb1; cx qb0, qb1"#)));
    }

    #[test]
    fn test_conditional_open_qasm()
    {
        let body = Composite::from_string("body", "H 0; H 1; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let res = gate.conditional_open_qasm("c == 3", &bit_names, &[0, 1]);
        let expected = String::from(
r#"if (c == 3) h qb0; if (c == 3) h qb1; if (c == 3) cx qb0, qb1;
if (c == 3) h qb0; if (c == 3) h qb1; if (c == 3) cx qb0, qb1;
if (c == 3) h qb0; if (c == 3) h qb1; if (c == 3) cx qb0, qb1"#);
        assert_eq!(res, Ok(expected));
    }

    #[test]
    fn test_c_qasm()
    {
        let body = Composite::from_string("body", "H 0; H 1; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let qasm = gate.c_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm,
            Ok(String::from(".myloop(3)\nh qb0\nh qb1\ncnot qb0, qb1\n.end")));
    }

    #[test]
    fn test_conditional_c_qasm()
    {
        let body = Composite::from_string("body", "H 0; H 1; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 3, body);
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let res = gate.conditional_c_qasm("c == 3", &bit_names, &[0, 1]);
        let expected = String::from(
r#"c-h c == 3, qb0
c-h c == 3, qb1
c-cnot c == 3, qb0, qb1
c-h c == 3, qb0
c-h c == 3, qb1
c-cnot c == 3, qb0, qb1
c-h c == 3, qb0
c-h c == 3, qb1
c-cnot c == 3, qb0, qb1"#
        );
        assert_eq!(res, Ok(expected));
    }

    #[test]
    fn test_latex()
    {
        let body = Composite::from_string("body", "H 0; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 0, body);
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qw \\
    \lstick{\ket{0}} & \qw \\
}
"#);

        let body = Composite::from_string("body", "H 0; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 1, body);
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{H} & \ctrl{1} & \qw \\
    \lstick{\ket{0}} & \qw & \targ & \qw \\
}
"#);

        let body = Composite::from_string("body", "H 0; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 2, body);
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{H} & \ctrl{1} & \gate{H} & \ctrl{1} & \qw \\
    \lstick{\ket{0}} & \qw & \targ & \qw & \targ & \qw \\
}
"#);

        let body = Composite::from_string("body", "H 0; CX 0 1").unwrap();
        let gate = Loop::new("myloop", 15, body);
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    & \mbox{} \POS"2,2"."2,2"."2,6"."2,6"!C*+<.7em>\frm{^\}},+U*++!D{15\times}\\
    & & & & & & \\
    \lstick{\ket{0}} & \gate{H} & \ctrl{1} & \cds{1}{\cdots} & \gate{H} & \ctrl{1} & \qw \\
    \lstick{\ket{0}} & \qw & \targ & \qw & \qw & \targ & \qw \\
}
"#);
    }
}
