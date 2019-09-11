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

use crate::export::CircuitGate;
use crate::stabilizer::PauliOp;
use super::*;

/// Structure for a description of a subgate.
#[derive(Debug)]
struct SubGateDesc
{
    /// Name of the gate.
    name: String,
    /// Parameters to the gate.
    args: Vec<f64>,
    /// Bits this gate will operate on.
    bits: Vec<usize>
}

impl SubGateDesc
{
    /// Create a new subgate description.
    fn new(name: &str, args: Vec<f64>, bits: Vec<usize>) -> Self
    {
        SubGateDesc
        {
            name: String::from(name),
            args: args,
            bits: bits
        }
    }
}

/// Operation in a composite gate.
#[derive(Clone)]
struct SubGate
{
    /// The gate
    gate: Box<CircuitGate>,
    /// The bits on which the gate acts
    bits: Vec<usize>
}

impl SubGate
{
    /// Create a new composite gate operation
    fn new<G>(gate: G, bits: &[usize]) -> Self
    where G: 'static + CircuitGate
    {
        SubGate
        {
            gate: Box::new(gate),
            bits: bits.to_owned()
        }
    }
}

/// Composite gate.
///
/// Struct Composite provides for user-defined gates that are made out of a
/// sequence of more primitive gates.
#[derive(Clone)]
pub struct Composite
{
    // The name of the gate
    name: String,
    // The number of gates on which this gate operates
    nr_bits: usize,
    // The operations making up the gate
    ops: Vec<SubGate>
}

impl Composite
{
    /// Create a new composite gate.
    ///
    /// Initialize a new composite gate with name `name` for operating on `nr_bits`
    /// qubits at a time. The gates making up the operation should be added
    /// using the `add_gate()` function.
    pub fn new(name: &str, nr_bits: usize) -> Self
    {
        Composite
        {
            name: name.to_owned(),
            nr_bits: nr_bits,
            ops: vec![]
        }
    }

    /// Parse the subgate name.
    ///
    /// Try to retrieve the name of the subgate from `desc`. On success,
    /// return the name, and the remainder of the subgate description to be
    /// parsed. On failure, return ParseError::NoGateName.
    fn parse_gate_name(desc: &str) -> crate::error::ParseResult<(&str, &str)>
    {
        let re = regex::Regex::new(r"(?i)^\s*([a-z][a-z0-9]*)").unwrap();
        if let Some(captures) = re.captures(desc)
        {
            let m = captures.get(1).unwrap();
            let rest = &desc[m.end()..];
            Ok((m.as_str(), rest))
        }
        else
        {
            Err(crate::error::ParseError::NoGateName(String::from(desc)))
        }
    }

    /// Parse arguments to a subgate.
    ///
    /// Parse arguments to a subgate, if any, from description string `desc`.
    /// If no parenthesized argument list is found, an emmpty argument vector
    /// is returned. If there is an argument list, then if it can be parsed
    /// successfully, the arguments are returned,¸together with the rest of the
    /// description string that needs to be parsed for bit numbers. On failure,
    /// a ParseError is returned.
    fn parse_gate_args(desc: &str) -> crate::error::ParseResult<(Vec<f64>, &str)>
    {
        let open_args = regex::Regex::new(r"^\s*\(").unwrap();
        let sep_args = regex::Regex::new(r"^\s*,").unwrap();
        let close_args = regex::Regex::new(r"^\s*\)").unwrap();
        if let Some(m) = open_args.find(desc)
        {
            let (arg, mut rest) = crate::expression::Expression::parse(&desc[m.end()..])?;
            let mut args = vec![];
            match arg.eval()
            {
                Ok(x) => { args.push(x); },
                _     => {
                    return Err(
                        crate::error::ParseError::InvalidArgument(String::from(m.as_str()))
                    );
                }
            }

            while let Some(m) = sep_args.find(rest)
            {
                let (arg, new_rest) = crate::expression::Expression::parse(&rest[m.end()..])?;
                match arg.eval()
                {
                    Ok(x) => { args.push(x); },
                    _     => {
                        return Err(
                            crate::error::ParseError::InvalidArgument(String::from(m.as_str()))
                        );
                    }
                }
                rest = new_rest;
            }

            if let Some(m) = close_args.find(rest)
            {
                Ok((args, &rest[m.end()..]))
            }
            else
            {
                Err(crate::error::ParseError::UnclosedParentheses(String::from(desc)))
            }
        }
        else
        {
            Ok((vec![], desc))
        }
    }

    /// Parse the bit numbers for a subgate.
    ///
    /// Parse the bit numbers on which the subgate operates from description
    /// string `desc`. Return the bits and the unparsed remainder of the
    /// description string on success, or a ParseError on failure.
    fn parse_gate_bits<'a>(desc: &'a str, name: &str)
        -> crate::error::ParseResult<(Vec<usize>, &'a str)>
    {
        let re = regex::Regex::new(r"^\s*(\d+)").unwrap();
        let mut rest = desc;
        let mut bits = vec![];
        while let Some(captures) = re.captures(rest)
        {
            let m = captures.get(0).unwrap();
            rest = &rest[m.end()..];

            let bit_txt = captures[1].trim();
            if let Ok(bit) = bit_txt.parse()
            {
                bits.push(bit);
            }
            else
            {
                return Err(crate::error::ParseError::InvalidBit(String::from(bit_txt)));
            }
        }

        if bits.is_empty()
        {
            Err(crate::error::ParseError::NoBits(String::from(name)))
        }
        else
        {
            Ok((bits, rest))
        }
    }

    /// Parse a gate description.
    ///
    /// Parse the subgate description string `desc`. Returns the subgate
    /// description on success, or a ParseError on failure.
    fn parse_gate_desc(desc: &str) -> crate::error::ParseResult<SubGateDesc>
    {
        let (name, rest) = Self::parse_gate_name(desc)?;
        let (args, rest) = Self::parse_gate_args(rest)?;
        let (bits, rest) = Self::parse_gate_bits(rest, name)?;

        let rest = rest.trim();
        if !rest.is_empty()
        {
            Err(crate::error::ParseError::TrailingText(String::from(rest)))
        }
        else
        {
            Ok(SubGateDesc::new(name, args, bits))
        }
    }

    /// Ensure correct number of arguments and bits.
    ///
    /// Ensure that the number of arguments in the subgate description `desc`
    /// matches `nr_args`, and that the number of bits in `desc` is equal to
    //// `nr_bits`. Return a ParseError on failure.
    fn assert_nr_args_bits(nr_args: usize, nr_bits: usize, desc: &SubGateDesc)
        -> crate::error::ParseResult<()>
    {
        if nr_args != desc.args.len()
        {
            Err(crate::error::ParseError::InvalidNrArguments(desc.args.len(), nr_args, desc.name.clone()))
        }
        else if nr_bits != desc.bits.len()
        {
            Err(crate::error::ParseError::InvalidNrBits(desc.bits.len(), nr_bits, desc.name.clone()))
        }
        else
        {
            Ok(())
        }
    }

    /// Create a new composite gate from a description string.
    ///
    /// Create a new composite gate with name `name`, based on the description
    /// in `desc`. The format of the description is as follows:
    /// * One or more subgate descriptions, separated by semicolons.
    /// * A subgate description consists of the name of the gate; optionally
    ///   followed by comma-separated parameter list in parentheses; followed by
    ///   one or more bit numbers on which the sub gate operates, separated by
    ///   white space."Failed to parse argument \"{}\"", text
    /// * Currently, only real numbers are allowed for parameters.
    /// Examples:
    /// ```text
    /// H 1; CX 0 1; H 1
    /// RY(4.7124) 1; CX 1 0; RY(1.5708) 1; X1
    /// ```
    pub fn from_string(name: &str, desc: &str) -> crate::error::ParseResult<Self>
    {
        let mut gates = vec![];
        let mut max_bit = 0;
        for part in desc.split(';')
        {
            let gate = Self::parse_gate_desc(part)?;
            max_bit = ::std::cmp::max(max_bit, *gate.bits.iter().max().unwrap());
            gates.push(gate);
        }

        let mut composite = Self::new(name, max_bit+1);
        for gate in gates
        {
            match gate.name.to_lowercase().as_str()
            {
                "ccrx" => {
                    Self::assert_nr_args_bits(1, 3, &gate)?;
                    composite.add_gate(CCRX::new(gate.args[0]), &gate.bits);
                },
                "ccry" => {
                    Self::assert_nr_args_bits(1, 3, &gate)?;
                    composite.add_gate(CCRY::new(gate.args[0]), &gate.bits);
                },
                "ccrz" => {
                    Self::assert_nr_args_bits(1, 3, &gate)?;
                    composite.add_gate(CCRZ::new(gate.args[0]), &gate.bits);
                },
                "ccx" => {
                    Self::assert_nr_args_bits(0, 3, &gate)?;
                    composite.add_gate(CCX::new(), &gate.bits);
                },
                "ccz" => {
                    Self::assert_nr_args_bits(0, 3, &gate)?;
                    composite.add_gate(CCZ::new(), &gate.bits);
                },
                "ch" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CH::new(), &gate.bits);
                },
                "crx" => {
                    Self::assert_nr_args_bits(1, 2, &gate)?;
                    composite.add_gate(CRX::new(gate.args[0]), &gate.bits);
                },
                "cry" => {
                    Self::assert_nr_args_bits(1, 2, &gate)?;
                    composite.add_gate(CRY::new(gate.args[0]), &gate.bits);
                },
                "crz" => {
                    Self::assert_nr_args_bits(1, 2, &gate)?;
                    composite.add_gate(CRZ::new(gate.args[0]), &gate.bits);
                },
                "cs" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CS::new(), &gate.bits);
                },
                "csdg" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CSdg::new(), &gate.bits);
                },
                "ct" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CT::new(), &gate.bits);
                },
                "ctdg" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CTdg::new(), &gate.bits);
                },
                "cu1" => {
                    Self::assert_nr_args_bits(1, 2, &gate)?;
                    composite.add_gate(CU1::new(gate.args[0]), &gate.bits);
                },
                "cu2" => {
                    Self::assert_nr_args_bits(2, 2, &gate)?;
                    composite.add_gate(CU2::new(gate.args[0], gate.args[1]), &gate.bits);
                },
                "cu3" => {
                    Self::assert_nr_args_bits(3, 2, &gate)?;
                    composite.add_gate(CU3::new(gate.args[0], gate.args[1], gate.args[2]), &gate.bits);
                },
                "cv" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CV::new(), &gate.bits);
                },
                "cvdg" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CVdg::new(), &gate.bits);
                },
                "cx" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CX::new(), &gate.bits);
                },
                "cy" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CY::new(), &gate.bits);
                },
                "cz" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(CZ::new(), &gate.bits);
                },
                "h" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(H::new(), &gate.bits);
                },
                "i" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(I::new(), &gate.bits);
                },
                "rx" => {
                    Self::assert_nr_args_bits(1, 1, &gate)?;
                    composite.add_gate(RX::new(gate.args[0]), &gate.bits);
                },
                "ry" => {
                    Self::assert_nr_args_bits(1, 1, &gate)?;
                    composite.add_gate(RY::new(gate.args[0]), &gate.bits);
                },
                "rz" => {
                    Self::assert_nr_args_bits(1, 1, &gate)?;
                    composite.add_gate(RZ::new(gate.args[0]), &gate.bits);
                },
                "s" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(S::new(), &gate.bits);
                },
                "sdg" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(Sdg::new(), &gate.bits);
                },
                "t" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(T::new(), &gate.bits);
                },
                "tdg" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(Tdg::new(), &gate.bits);
                },
                "swap" => {
                    Self::assert_nr_args_bits(0, 2, &gate)?;
                    composite.add_gate(Swap::new(), &gate.bits);
                },
                "u1" => {
                    Self::assert_nr_args_bits(1, 1, &gate)?;
                    composite.add_gate(U1::new(gate.args[0]), &gate.bits);
                },
                "u2" => {
                    Self::assert_nr_args_bits(2, 1, &gate)?;
                    composite.add_gate(U2::new(gate.args[0], gate.args[1]), &gate.bits);
                },
                "u3" => {
                    Self::assert_nr_args_bits(3, 1, &gate)?;
                    composite.add_gate(U3::new(gate.args[0], gate.args[1], gate.args[2]), &gate.bits);
                },
                "v" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(V::new(), &gate.bits);
                },
                "vdg" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(Vdg::new(), &gate.bits);
                },
                "x" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(X::new(), &gate.bits);
                },
                "y" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(Y::new(), &gate.bits);
                },
                "z" => {
                    Self::assert_nr_args_bits(0, 1, &gate)?;
                    composite.add_gate(Z::new(), &gate.bits);
                },
                _ => { return Err(crate::error::ParseError::UnknownGate(gate.name)); }
            }
        }

        Ok(composite)
    }

    /// Add a gate.
    ///
    /// Append a `n`-ary subgate `gate`, operating on the `n` qubits in `bits`,
    /// to this composite gate.
    pub fn add_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: CircuitGate
    {
        self.ops.push(SubGate::new(gate, bits));
    }
}

impl crate::gates::Gate for Composite
{
    fn cost(&self) -> f64
    {
        self.ops.iter().map(|op| op.gate.cost()).sum()
    }

    fn description(&self) -> &str
    {
        &self.name
    }

    fn nr_affected_bits(&self) -> usize
    {
        self.nr_bits
    }

    fn matrix(&self) -> crate::cmatrix::CMatrix
    {
        let mut res = crate::cmatrix::CMatrix::eye(1 << self.nr_bits);
        self.apply_mat_slice(res.view_mut());
        res
    }

    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        for op in self.ops.iter()
        {
            apply_gate_slice(state.view_mut(), &*op.gate, &op.bits, self.nr_bits);
        }
    }

    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        for op in self.ops.iter()
        {
            apply_gate_mat_slice(state.view_mut(), &*op.gate, &op.bits, self.nr_bits);
        }
    }

    fn is_stabilizer(&self) -> bool
    {
        self.ops.iter().all(|op| op.gate.is_stabilizer())
    }

    fn conjugate(&self, ops: &mut [PauliOp]) -> crate::error::Result<bool>
    {
        self.check_nr_bits(ops.len())?;
        let mut flip_sign = false;
        for op in self.ops.iter()
        {
            let mut gate_ops: Vec<PauliOp> = op.bits.iter().map(|&b| ops[b]).collect();
            flip_sign ^= op.gate.conjugate(&mut gate_ops)?;
            for (&i, gop) in op.bits.iter().zip(gate_ops)
            {
                ops[i] = gop;
            }
        }
        Ok(flip_sign)
    }
}

impl crate::export::OpenQasm for Composite
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        let mut res = String::new();
        if self.ops.len() > 0
        {
            let gate_bits: Vec<usize> = self.ops[0].bits.iter().map(|&b| bits[b]).collect();
            res = self.ops[0].gate.open_qasm(bit_names, &gate_bits)?;
            for op in self.ops[1..].iter()
            {
                let gate_bits: Vec<usize> = op.bits.iter().map(|&b| bits[b]).collect();
                let qasm = op.gate.open_qasm(bit_names, &gate_bits)?;
                res += &format!("; {}", qasm);
            }
        }
        Ok(res)
    }

    fn conditional_open_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> crate::error::Result<String>
    {
        let mut res = String::new();
        if self.ops.len() > 0
        {
            let gate_bits: Vec<usize> = self.ops[0].bits.iter().map(|&b| bits[b]).collect();
            res = self.ops[0].gate.conditional_open_qasm(condition, bit_names,
                &gate_bits)?;
            for op in self.ops[1..].iter()
            {
                let gate_bits: Vec<usize> = op.bits.iter().map(|&b| bits[b]).collect();
                let qasm = op.gate.conditional_open_qasm(condition, bit_names,
                    &gate_bits)?;
                res += "; ";
                res += &qasm;
            }
        }
        Ok(res)
    }
}

impl crate::export::CQasm for Composite
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> crate::error::Result<String>
    {
        let mut res = String::new();
        if self.ops.len() > 0
        {
            let gate_bits: Vec<usize> = self.ops[0].bits.iter()
                .map(|&b| bits[b]).collect();
            res = self.ops[0].gate.c_qasm(bit_names, &gate_bits)?;
            for op in self.ops[1..].iter()
            {
                let gate_bits: Vec<usize> = op.bits.iter().map(|&b| bits[b]).collect();
                let qasm = op.gate.c_qasm(bit_names, &gate_bits)?;
                res += &format!("\n{}", qasm);
            }
        }
        Ok(res)
    }

    fn conditional_c_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> crate::error::Result<String>
    {
        let mut res = String::new();
        if self.ops.len() > 0
        {
            let gate_bits: Vec<usize> = self.ops[0].bits.iter()
                .map(|&b| bits[b]).collect();
            res = self.ops[0].gate.conditional_c_qasm(condition, bit_names,
                &gate_bits)?;
            for op in self.ops[1..].iter()
            {
                let gate_bits: Vec<usize> = op.bits.iter().map(|&b| bits[b]).collect();
                let qasm = op.gate.conditional_c_qasm(condition, bit_names,
                    &gate_bits)?;
                res += &format!("\n{}", qasm);
            }
        }
        Ok(res)
    }
}

impl crate::export::Latex for Composite
{
    fn latex(&self, bits: &[usize], state: &mut crate::export::LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;

        if state.expand_composite()
        {
            for op in self.ops.iter()
            {
                let gate_bits: Vec<usize> = op.bits.iter().map(|&b| bits[b]).collect();
                op.gate.latex(&gate_bits, state)?;
            }
            Ok(())
        }
        else
        {
            state.add_block_gate(bits, self.description())
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::Composite;
    use crate::export::{Latex, LatexExportState, OpenQasm, CQasm};
    use crate::gates::{Gate, CCX, CX, H, X};
    use crate::stabilizer::PauliOp;
    use num_complex::Complex;

    #[test]
    fn test_description()
    {
        let gate = Composite::new("G", 3);
        assert_eq!(gate.description(), "G");
    }

    #[test]
    fn test_cost()
    {
        let mut gate = Composite::new("Inc2", 2);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(X::new(), &[1]);
        assert_eq!(gate.cost(), 1001.0 + 201.0);

        let mut gate = Composite::new("H3", 3);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(H::new(), &[2]);
        assert_eq!(gate.cost(), 3.0 * 104.0);
    }

    #[test]
    fn test_matrix()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;

        let mut gate = Composite::new("CZ", 2);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(H::new(), &[1]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z,  z],
            [z, o, z,  z],
            [z, z, o,  z],
            [z, z, z, -o]
        ]);

        let mut gate = Composite::new("Inc", 2);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(X::new(), &[1]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [z, z, z, o],
            [o, z, z, z],
            [z, o, z, z],
            [z, z, o, z]
        ]);

        let mut gate = Composite::new("Inc", 3);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(CCX::new(), &[0, 1, 2]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(CX::new(), &[1, 2]);
        gate.add_gate(H::new(), &[2]);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(H::new(), &[1]);
        gate.add_gate(X::new(), &[2]);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [z, z, z, z, z, z, z, o],
            [o, z, z, z, z, z, z, z],
            [z, o, z, z, z, z, z, z],
            [z, z, o, z, z, z, z, z],
            [z, z, z, o, z, z, z, z],
            [z, z, z, z, o, z, z, z],
            [z, z, z, z, z, o, z, z],
            [z, z, z, z, z, z, o, z]
        ]);
    }

    #[test]
    fn test_from_string()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;

        // Test composition
        match Composite::from_string("Inc3", "CCX 2 1 0; CX 2 1; X 2")
        {
            Ok(gate) => {
                assert_eq!(gate.description(), "Inc3");
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [z, z, z, z, z, z, z, o],
                    [o, z, z, z, z, z, z, z],
                    [z, o, z, z, z, z, z, z],
                    [z, z, o, z, z, z, z, z],
                    [z, z, z, o, z, z, z, z],
                    [z, z, z, z, o, z, z, z],
                    [z, z, z, z, z, o, z, z],
                    [z, z, z, z, z, z, o, z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!(err); }
            // LCOV_EXCL_STOP
        }

        // Test arguments
        match Composite::from_string("Y", "U3(3.141592653589793,1.570796326794897,1.570796326794897) 0")
        {
            Ok(gate) => {
                assert_eq!(gate.description(), "Y");
                assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [i, z]]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!(err); }
            // LCOV_EXCL_STOP
        }
    }

    #[test]
    fn test_from_string_gates()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;

        match Composite::from_string("G", "CCRX(3.141592653589793) 0 1 2")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z, z, z,  z,  z],
                    [z, o, z, z, z, z,  z,  z],
                    [z, z, o, z, z, z,  z,  z],
                    [z, z, z, o, z, z,  z,  z],
                    [z, z, z, z, o, z,  z,  z],
                    [z, z, z, z, z, o,  z,  z],
                    [z, z, z, z, z, z,  z, -i],
                    [z, z, z, z, z, z, -i,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CCRY(3.141592653589793) 0 1 2")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z, z, z, z,  z],
                    [z, o, z, z, z, z, z,  z],
                    [z, z, o, z, z, z, z,  z],
                    [z, z, z, o, z, z, z,  z],
                    [z, z, z, z, o, z, z,  z],
                    [z, z, z, z, z, o, z,  z],
                    [z, z, z, z, z, z, z, -o],
                    [z, z, z, z, z, z, o,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CCRZ(3.141592653589793) 0 1 2")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z, z, z,  z, z],
                    [z, o, z, z, z, z,  z, z],
                    [z, z, o, z, z, z,  z, z],
                    [z, z, z, o, z, z,  z, z],
                    [z, z, z, z, o, z,  z, z],
                    [z, z, z, z, z, o,  z, z],
                    [z, z, z, z, z, z, -i, z],
                    [z, z, z, z, z, z,  z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CCX 0 1 2")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z, z, z, z, z],
                    [z, o, z, z, z, z, z, z],
                    [z, z, o, z, z, z, z, z],
                    [z, z, z, o, z, z, z, z],
                    [z, z, z, z, o, z, z, z],
                    [z, z, z, z, z, o, z, z],
                    [z, z, z, z, z, z, z, o],
                    [z, z, z, z, z, z, o, z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CCZ 0 1 2")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z, z, z, z,  z],
                    [z, o, z, z, z, z, z,  z],
                    [z, z, o, z, z, z, z,  z],
                    [z, z, z, o, z, z, z,  z],
                    [z, z, z, z, o, z, z,  z],
                    [z, z, z, z, z, o, z,  z],
                    [z, z, z, z, z, z, o,  z],
                    [z, z, z, z, z, z, z, -o]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CH 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,  z],
                    [z, o, z,  z],
                    [z, z, x,  x],
                    [z, z, x, -x]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CRX(3.141592653589793) 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z,  z,  z],
                    [z, o,  z,  z],
                    [z, z,  z, -i],
                    [z, z, -i,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CRY(3.141592653589793) 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,  z],
                    [z, o, z,  z],
                    [z, z, z, -o],
                    [z, z, o,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CRZ(3.141592653589793) 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z,  z, z],
                    [z, o,  z, z],
                    [z, z, -i, z],
                    [z, z,  z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CS 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, o, z, z],
                    [z, z, o, z],
                    [z, z, z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CSdg 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,  z],
                    [z, o, z,  z],
                    [z, z, o,  z],
                    [z, z, z, -i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CT 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,       z],
                    [z, o, z,       z],
                    [z, z, o,       z],
                    [z, z, z, x*(o+i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CTdg 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,       z],
                    [z, o, z,       z],
                    [z, z, o,       z],
                    [z, z, z, x*(o-i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CU1(1.570796326794897) 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, o, z, z],
                    [z, z, o, z],
                    [z, z, z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CU2(0.7853981633974483, 0.6931471805599453) 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, o, z, z],
                    [z, z, Complex::new(0.7071067811865476, 0.0), Complex::new(-0.5439340435069544, -0.4518138513969824)],
                    [z, z, Complex::new(               0.5, 0.5), Complex::new(0.06513881252516862,  0.7041000888388035)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CU3(0.32, 0.7853981633974483, 0.6931471805599453) 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, o, z, z],
                    [z, z, Complex::new(0.9872272833756269,                0.0), Complex::new( -0.1225537622232209, -0.1017981646382380)],
                    [z, z, Complex::new(0.1126549842634128, 0.1126549842634128), Complex::new(0.09094356700076842,  0.9830294892130130)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CV 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z,         z,         z],
                    [z, o,         z,         z],
                    [z, z, 0.5*(o+i), 0.5*(o-i)],
                    [z, z, 0.5*(o-i), 0.5*(o+i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!(err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CVdg 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z,         z,         z],
                    [z, o,         z,         z],
                    [z, z, 0.5*(o-i), 0.5*(o+i)],
                    [z, z, 0.5*(o+i), 0.5*(o-i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!(err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CX 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, o, z, z],
                    [z, z, z, o],
                    [z, z, o, z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CY 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,  z],
                    [z, o, z,  z],
                    [z, z, z, -i],
                    [z, z, i,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "CZ 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z,  z],
                    [z, o, z,  z],
                    [z, z, o,  z],
                    [z, z, z, -o]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "H 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [x,  x],
                    [x, -x]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "I 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z],
                    [z, o]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "RX(3.141592653589793) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [ z, -i],
                    [-i,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "RY(3.141592653589793) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [z, -o],
                    [o,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "RZ(3.141592653589793) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [-i, z],
                    [ z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "S 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z],
                    [z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "Sdg 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,  z],
                    [z, -i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "T 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,       z],
                    [z, x*(o+i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "Tdg 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,       z],
                    [z, x*(o-i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "Swap 0 1")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z, z, z],
                    [z, z, o, z],
                    [z, o, z, z],
                    [z, z, z, o]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(1.570796326794897) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o, z],
                    [z, i]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U2(0.7853981633974483, 0.6931471805599453) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [Complex::new(0.7071067811865476, 0.0), Complex::new(-0.5439340435069544, -0.4518138513969824)],
                    [Complex::new(               0.5, 0.5), Complex::new(0.06513881252516862,  0.7041000888388035)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U3(0.32, 0.7853981633974483, 0.6931471805599453) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [Complex::new(0.9872272833756269,                0.0), Complex::new( -0.1225537622232209, -0.1017981646382380)],
                    [Complex::new(0.1126549842634128, 0.1126549842634128), Complex::new(0.09094356700076842,  0.9830294892130130)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "V 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [0.5*(o+i), 0.5*(o-i)],
                    [0.5*(o-i), 0.5*(o+i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "Vdg 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [0.5*(o-i), 0.5*(o+i)],
                    [0.5*(o+i), 0.5*(o-i)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "X 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [z, o],
                    [o, z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "Y 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [z, -i],
                    [i,  z]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "Z 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,  z],
                    [z, -o]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }
    }

    #[test]
    fn test_from_string_arguments()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;

        match Composite::from_string("G", "U1(0.23) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                z],
                    [z, Complex::from_polar(&1.0, &0.23)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(0.23+0.16) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                z],
                    [z, Complex::from_polar(&1.0, &0.39)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(1.23-0.16) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                z],
                    [z, Complex::from_polar(&1.0, &1.07)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(0.8*0.6) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                z],
                    [z, Complex::from_polar(&1.0, &0.48)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(1.38/2) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                z],
                    [z, Complex::from_polar(&1.0, &0.69)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(-0.23) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                 z],
                    [z, Complex::from_polar(&1.0, &-0.23)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(1.03^5) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                        z],
                    [z, Complex::from_polar(&1.0, &1.1592740743)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(sin(1.0)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                              z],
                    [z, Complex::from_polar(&1.0, &0.8414709848078965)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(cos(1.0)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                              z],
                    [z, Complex::from_polar(&1.0, &0.5403023058681398)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(tan(1.0)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                              z],
                    [z, Complex::from_polar(&1.0, &1.5574077246549023)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(exp(1.0)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                             z],
                    [z, Complex::from_polar(&1.0, &2.718281828459045)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(ln(0.8)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                               z],
                    [z, Complex::from_polar(&1.0, &-0.2231435513142097)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(sqrt(0.8)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                              z],
                    [z, Complex::from_polar(&1.0, &0.8944271909999159)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(0.5*(1.0-0.26)) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                                z],
                    [z, Complex::from_polar(&1.0, &0.37)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(3) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,                               z],
                    [z, Complex::from_polar(&1.0, &3.0)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U1(pi) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [o,  z],
                    [z, -o]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }

        match Composite::from_string("G", "U3(pi/2, 1.03^2^(1.05-0.23), -1.78) 0")
        {
            Ok(gate) => {
                assert_complex_matrix_eq!(gate.matrix(), array![
                    [Complex::new(0.7071067811865476, 0.0),
                        Complex::new(0.1468526445611853, 0.6916894540076393)],
                    [Complex::new(0.3496446456390944, 0.6146125786020914),
                        Complex::new(0.5285973886766332, -0.4696645618782032)]
                ]);
            },
            // LCOV_EXCL_START
            Err(err) => { panic!("{}", err); }
            // LCOV_EXCL_STOP
        }
    }

    #[test]
    fn test_from_string_errors()
    {
        // Invalid gate name
        let res = Composite::from_string("XXX", "XYZ 0");
        assert!(matches!(res, Err(crate::error::ParseError::UnknownGate(_))));

        // Missing gate name
        let res = Composite::from_string("XXX", "X 1; 0");
        assert!(matches!(res, Err(crate::error::ParseError::NoGateName(_))));

        // Invalid nr of arguments
        let res = Composite::from_string("XXX", "RX(1.2, 3.4) 1");
        assert!(matches!(res, Err(crate::error::ParseError::InvalidNrArguments(2, 1, _))));

        // Invalid nr of bits to operate on
        let res = Composite::from_string("XXX", "H 0 1");
        assert!(matches!(res, Err(crate::error::ParseError::InvalidNrBits(2, 1, _))));

        // Invalid arguments
        let res = Composite::from_string("XXX", "RX(abc) 1");
        assert!(matches!(res, Err(crate::error::ParseError::InvalidArgument(_))));

        let res = Composite::from_string("G", "U1(12897231928172918729136192817936) 0");
        assert!(matches!(res, Err(crate::error::ParseError::InvalidArgument(_))));

        // Missing bit number
        let res = Composite::from_string("XXX", "H 0; X");
        assert!(matches!(res, Err(crate::error::ParseError::NoBits(_))));

        // Invalid bit number
        let res = Composite::from_string("XXX", "H 117356715625188271521875");
        assert!(matches!(res, Err(crate::error::ParseError::InvalidBit(_))));

        // Trailing junk
        let res = Composite::from_string("XXX", "H 0 and something");
        assert!(matches!(res, Err(crate::error::ParseError::TrailingText(_))));

        // Unclosed arguments
        let res = Composite::from_string("XXX", "RX(1.2a) 1");
        assert!(matches!(res, Err(crate::error::ParseError::UnclosedParentheses(_))));

        let res = Composite::from_string("XXX", "RX(1.2*(1+2 1");
        assert!(matches!(res, Err(crate::error::ParseError::UnclosedParentheses(_))));

        let res = Composite::from_string("XXX", "RX(sin(1.2 1");
        assert!(matches!(res, Err(crate::error::ParseError::UnclosedParentheses(_))));
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let mut gate = Composite::new("Inc2", 2);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(X::new(), &[1]);
        let qasm = gate.open_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("cx qb0, qb1; x qb1")));
    }

    #[test]
    fn test_conditional_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];

        let mut gate = Composite::new("XXX", 2);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(X::new(), &[1]);
        let qasm = gate.conditional_open_qasm("b == 3", &bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("if (b == 3) h qb0; if (b == 3) x qb1")));
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];
        let mut gate = Composite::new("Inc2", 2);
        gate.add_gate(CX::new(), &[0, 1]);
        gate.add_gate(X::new(), &[1]);
        let qasm = gate.c_qasm(&bit_names, &[0, 1]);
        assert_eq!(qasm, Ok(String::from("cnot qb0, qb1\nx qb1")));
    }

    #[test]
    fn test_conditional_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];

        let mut gate = Composite::new("XXX", 2);
        gate.add_gate(H::new(), &[0]);
        gate.add_gate(X::new(), &[1]);
        let qasm = gate.conditional_c_qasm("b == 3", &bit_names, &[0, 1]);
        let expected = String::from(
r#"c-h b == 3, qb0
c-x b == 3, qb1"#);
        assert_eq!(qasm, Ok(expected));
    }

    #[test]
    fn test_latex()
    {
        let gate = Composite::from_string("CZ", "H 1; CX 0 1; H 1").unwrap();
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qw & \ctrl{1} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{H} & \targ & \gate{H} & \qw \\
}
"#);

        let gate = Composite::from_string("CZ", "H 1; CX 0 1; H 1").unwrap();
        let mut state = LatexExportState::new(2, 0);
        state.set_expand_composite(false);
        assert_eq!(gate.latex(&[0, 1], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \multigate{1}{CZ} & \qw \\
    \lstick{\ket{0}} & \ghost{CZ} & \qw \\
}
"#);

        let gate = Composite::from_string("CZ", "H 1; CX 0 1; H 1").unwrap();
        let mut state = LatexExportState::new(3, 0);
        state.set_expand_composite(false);
        assert_eq!(gate.latex(&[0, 2], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{CZ} & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{CZ} \qwx[-2] & \qw \\
}
"#);

        let gate = Composite::from_string("G", "H 0; CX 1 2; H 0").unwrap();
        let mut state = LatexExportState::new(4, 0);
        state.set_expand_composite(false);
        assert_eq!(gate.latex(&[0, 2, 3], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{G} & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
    \lstick{\ket{0}} & \multigate{1}{G} \qwx[-2] & \qw \\
    \lstick{\ket{0}} & \ghost{G} & \qw \\
}
"#);

        let gate = Composite::from_string("G", "H 0; CX 1 2; H 0").unwrap();
        let mut state = LatexExportState::new(4, 0);
        state.set_expand_composite(false);
        assert_eq!(gate.latex(&[0, 1, 3], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \multigate{1}{G} & \qw \\
    \lstick{\ket{0}} & \ghost{G} & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{G} \qwx[-2] & \qw \\
}
"#);
    }

    #[test]
    fn test_conjugate()
    {
        let gate = Composite::from_string("G", "CX 0 1; X 0").unwrap();
        let mut ops = [PauliOp::X, PauliOp::Y];
        assert_eq!(gate.conjugate(&mut ops), Ok(true));
        assert_eq!(ops, [PauliOp::Y, PauliOp::Z]);

        let gate = Composite::from_string("G", "CX 0 1; T 0").unwrap();
        let mut ops = [PauliOp::X, PauliOp::Y];
        assert!(matches!(gate.conjugate(&mut ops), Err(crate::error::Error::NotAStabilizer(_))));
    }
}
