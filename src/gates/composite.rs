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


extern crate ndarray;
extern crate regex;

use cmatrix;
use error;
use gates;
use export;
use support;

use export::CircuitGate;
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
    fn parse_gate_name(desc: &str) -> error::ParseResult<(&str, &str)>
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
            Err(error::ParseError::NoGateName(String::from(desc)))
        }
    }

    /// Parse a real number
    ///
    /// Parse a real number. This can be either in the form of a real literal,
    /// an integer literal (which is converted to a real number), or the string
    /// "pi", which obviously indicates the number π. On success, the parsed
    /// number is returned, together with the remainder of the string to be
    /// parsed. On failure, a ParseError::InvalidArgument error is returned.
    fn parse_real_literal(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let real = regex::Regex::new(
            r"^\s*((?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)"
        ).unwrap();
        let integer = regex::Regex::new(r"^\s*([1-9][0-9]*|0)").unwrap();
        let pi = regex::Regex::new(r"^\s*pi").unwrap();

        if let Some(captures) = real.captures(expr)
        {
            let m = captures.get(0).unwrap();
            if let Ok(nr) = captures[1].parse::<f64>()
            {
                Ok((nr, &expr[m.end()..]))
            }
            else
            {
                // It's probably impossible to get here with the above regular
                // expression. Too large numbers are mapped to Inf.
                Err(error::ParseError::InvalidArgument(String::from(expr)))
            }
        }
        else if let Some(captures) = integer.captures(expr)
        {
            let m = captures.get(0).unwrap();
            if let Ok(nr) = captures[1].parse::<u64>()
            {
                Ok((nr as f64, &expr[m.end()..]))
            }
            else
            {
                Err(error::ParseError::InvalidArgument(String::from(expr)))
            }
        }
        else if let Some(m) = pi.find(expr)
        {
            Ok((::std::f64::consts::PI, &expr[m.end()..]))
        }
        else
        {
            Err(error::ParseError::InvalidArgument(String::from(expr)))
        }
    }

    /// Parse a possibly parenthesized expression
    ///
    /// Parse an argument to a gate. This function expects either an expression
    /// in parentheses, or a literal real number. On success, the parsed
    /// number is returned, together with the remainder of the string to be
    /// parsed. On failure, a ParseError is returned.
    fn parse_parenthesized_expression(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let fun_open = regex::Regex::new(r"^\s*\(").unwrap();
        let fun_close = regex::Regex::new(r"^\s*\)").unwrap();
        if let Some(m) = fun_open.find(expr)
        {
            let (result, rest) = Self::parse_sum_expression(&expr[m.end()..])?;
            if let Some(m) = fun_close.find(rest)
            {
                Ok((result, &rest[m.end()..]))
            }
            else
            {
                Err(error::ParseError::UnclosedParentheses(String::from(expr)))
            }
        }
        else
        {
            Self::parse_real_literal(expr)
        }
    }

    /// Parse a function call
    ///
    /// Parse an argument to a gate. If a function call is found, it is parsed
    /// in this function, otherwise control is passed on to
    /// `parse_parenthesized_expression()`. The functions that are recognised
    /// are `sin`, `cos`, `tan`, `exp`, `ln`, and `sqrt`. On success, the parsed
    /// number is returned, together with the remainder of the string to be
    /// parsed. On failure, a ParseError is returned.
    fn parse_function_expression(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let fun_open = regex::Regex::new(r"^\s*(sin|cos|tan|exp|ln|sqrt)\s*\(").unwrap();
        let fun_close = regex::Regex::new(r"^\s*\)").unwrap();
        if let Some(captures) = fun_open.captures(expr)
        {
            let m = captures.get(0).unwrap();
            let (arg, new_rest) = Self::parse_sum_expression(&expr[m.end()..])?;
            if let Some(m) = fun_close.find(new_rest)
            {
                let result = match &captures[1]
                {
                    "sin"  => arg.sin(),
                    "cos"  => arg.cos(),
                    "tan"  => arg.tan(),
                    "exp"  => arg.exp(),
                    "ln"   => arg.ln(),
                    "sqrt" => arg.sqrt(),
                    // LCOV_EXCL_START
                    _      => { unreachable!() }
                    // LCOV_EXCL_STOP
                };
                Ok((result, &new_rest[m.end()..]))
            }
            else
            {
                Err(error::ParseError::UnclosedParentheses(String::from(expr)))
            }
        }
        else
        {
            Self::parse_parenthesized_expression(expr)
        }
    }

    /// Parse a power-raising expression
    ///
    /// Parse an argument to a gate, in the form of a number possibly raised
    /// to another number.  On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_power_expression(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let op = regex::Regex::new(r"^\s*\^").unwrap();
        let (left, rest) = Self::parse_function_expression(expr)?;
        if let Some(m) = op.find(rest)
        {
            let (right, new_rest) = Self::parse_power_expression(&rest[m.end()..])?;
            Ok((left.powf(right), new_rest))
        }
        else
        {
            Ok((left, rest))
        }
    }

    /// Parse a negated expression
    ///
    /// Parse an argument to a gate, in the form of a number that is zero or
    /// more times negated.  On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_negative_expression(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let op = regex::Regex::new(r"^\s*\-").unwrap();
        let mut rest = expr;
        let mut flip_sign = false;
        while let Some(m) = op.find(rest)
        {
            rest = &rest[m.end()..];
            flip_sign = !flip_sign;
        }

        let (mut result, new_rest) = Self::parse_power_expression(rest)?;
        if flip_sign
        {
            result = -result;
        }

        Ok((result, new_rest))
    }

    /// Parse a product expression
    ///
    /// Parse an argument to a gate, in the form of a product or quotient of one
    /// or more expressions. On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_product_expression(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let op = regex::Regex::new(r"^\s*([*/])").unwrap();
        let (mut left, mut rest) = Self::parse_negative_expression(expr)?;
        while let Some(captures) = op.captures(rest)
        {
            let m = captures.get(0).unwrap();
            let (right, new_rest) = Self::parse_negative_expression(&rest[m.end()..])?;
            if &captures[1] == "*"
            {
                left *= right;
            }
            else
            {
                left /= right;
            }
            rest = new_rest;
        }

        Ok((left, rest))
    }

    /// Parse a sum expression
    ///
    /// Parse an argument to a gate, in the form of a sum or difference of one
    /// or more expressions. On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_sum_expression(expr: &str) -> error::ParseResult<(f64, &str)>
    {
        let op = regex::Regex::new(r"^\s*([-+])").unwrap();
        let (mut left, mut rest) = Self::parse_product_expression(expr)?;
        while let Some(captures) = op.captures(rest)
        {
            let m = captures.get(0).unwrap();
            let (right, new_rest) = Self::parse_product_expression(&rest[m.end()..])?;
            if &captures[1] == "+"
            {
                left += right;
            }
            else
            {
                left -= right;
            }
            rest = new_rest;
        }

        Ok((left, rest))
    }

    /// Parse arguments to a subgate.
    ///
    /// Parse arguments to a subgate, if any, from description string `desc`.
    /// If no parenthesized argument list is found, an emmpty argument vector
    /// is returned. If there is an argument list, then if it can be parsed
    /// successfully, the arguments are returned,¸together with the rest of the
    /// description string that needs to be parsed for bit numbers. On failure,
    /// ParseError::InvalidArgument is returned.
    fn parse_gate_args(desc: &str) -> error::ParseResult<(Vec<f64>, &str)>
    {
        let open_args = regex::Regex::new(r"^\s*\(").unwrap();
        let sep_args = regex::Regex::new(r"^\s*,").unwrap();
        let close_args = regex::Regex::new(r"^\s*\)").unwrap();
        if let Some(m) = open_args.find(desc)
        {
            let (arg, mut rest) = Self::parse_sum_expression(&desc[m.end()..])?;
            let mut args = vec![arg];
            while let Some(m) = sep_args.find(rest)
            {
                let (arg, new_rest) = Self::parse_sum_expression(&rest[m.end()..])?;
                args.push(arg);
                rest = new_rest;
            }

            if let Some(m) = close_args.find(rest)
            {
                Ok((args, &rest[m.end()..]))
            }
            else
            {
                Err(error::ParseError::UnclosedParentheses(String::from(desc)))
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
        -> error::ParseResult<(Vec<usize>, &'a str)>
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
                return Err(error::ParseError::InvalidBit(String::from(bit_txt)));
            }
        }

        if bits.is_empty()
        {
            Err(error::ParseError::NoBits(String::from(name)))
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
    fn parse_gate_desc(desc: &str) -> error::ParseResult<SubGateDesc>
    {
        let (name, rest) = Self::parse_gate_name(desc)?;
        let (args, rest) = Self::parse_gate_args(rest)?;
        let (bits, rest) = Self::parse_gate_bits(rest, name)?;

        let rest = rest.trim();
        if !rest.is_empty()
        {
            Err(error::ParseError::TrailingText(String::from(rest)))
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
        -> error::ParseResult<()>
    {
        if nr_args != desc.args.len()
        {
            Err(error::ParseError::InvalidNrArguments(desc.name.clone()))
        }
        else if nr_bits != desc.bits.len()
        {
            Err(error::ParseError::InvalidNrBits(desc.name.clone()))
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
    pub fn from_string(name: &str, desc: &str) -> error::ParseResult<Self>
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
                "ccry" => {
                    Self::assert_nr_args_bits(1, 3, &gate)?;
                    composite.add_gate(CCRY::new(gate.args[0]), &gate.bits);
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
                _ => { return Err(error::ParseError::UnknownGate(gate.name)); }
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

impl gates::Gate for Composite
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

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let mut res = cmatrix::CMatrix::eye(1 << self.nr_bits);
        for op in self.ops.iter()
        {
            apply_gate(&mut res, &*op.gate, &op.bits, self.nr_bits);
        }

        res
    }
}

impl export::OpenQasm for Composite
{
    fn open_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
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
        bits: &[usize]) -> error::ExportResult<String>
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

impl export::CQasm for Composite
{
    fn c_qasm(&self, bit_names: &[String], bits: &[usize])
        -> error::ExportResult<String>
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
        bits: &[usize]) -> error::ExportResult<String>
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

impl export::Latex for Composite
{
    fn latex(&self, bits: &[usize], state: &mut export::LatexExportState)
    {
        if state.expand_composite()
        {
            for op in self.ops.iter()
            {
                let gate_bits: Vec<usize> = op.bits.iter().map(|&b| bits[b]).collect();
                // The composite gate may take up more than one column, so for
                // each gate in the composition, we need to ensure there's
                // place available. Therefore, use latex_checked() instead of
                // the latex() method on the sub-gate.
                op.gate.latex_checked(&gate_bits, state);
            }
        }
        else
        {
            let ranges = support::get_ranges(bits);

            let (first, last) = ranges[0];
            if last == first
            {
                state.set_field(first, format!(r"\gate{{{}}}", self.description()));
            }
            else
            {
                state.set_field(first,
                    format!(r"\multigate{{{}}}{{{}}}", last-first, self.description()));
                for bit in first+1..last+1
                {
                    state.set_field(bit, format!(r"\ghost{{{}}}", self.description()));
                }
            }

            let mut prev_last = last;
            for &(first, last) in ranges[1..].iter()
            {
                if last == first
                {
                    state.set_field(first,
                        format!(r"\gate{{{}}} \qwx[{}]", self.description(),
                            prev_last as isize - first as isize));
                }
                else
                {
                    state.set_field(first,
                        format!(r"\multigate{{{}}}{{{}}} \qwx[{}]",
                            last - first, self.description(),
                            prev_last as isize - first as isize));
                    for bit in first+1..last+1
                    {
                        state.set_field(bit, format!(r"\ghost{{{}}}", self.description()));
                    }
                }

                prev_last = last;
            }
        }
    }
}

#[cfg(test)]
mod tests
{
    extern crate num_complex;

    use super::Composite;
    use cmatrix;
    use error;
    use export::{Latex, LatexExportState, OpenQasm, CQasm};
    use gates::{Gate, CCX, CX, H, X};
    use self::num_complex::Complex;

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
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

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
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;

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
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;

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
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

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
        assert!(matches!(res, Err(error::ParseError::UnknownGate(_))));

        // Missing gate name
        let res = Composite::from_string("XXX", "X 1; 0");
        assert!(matches!(res, Err(error::ParseError::NoGateName(_))));

        // Invalid nr of arguments
        let res = Composite::from_string("XXX", "RX(1.2, 3.4) 1");
        assert!(matches!(res, Err(error::ParseError::InvalidNrArguments(_))));

        // Invalid nr of bits to operate on
        let res = Composite::from_string("XXX", "H 0 1");
        assert!(matches!(res, Err(error::ParseError::InvalidNrBits(_))));

        // Invalid arguments
        let res = Composite::from_string("XXX", "RX(abc) 1");
        assert!(matches!(res, Err(error::ParseError::InvalidArgument(_))));

        let res = Composite::from_string("G", "U1(12897231928172918729136192817936) 0");
        assert!(matches!(res, Err(error::ParseError::InvalidArgument(_))));

        // Missing bit number
        let res = Composite::from_string("XXX", "H 0; X");
        assert!(matches!(res, Err(error::ParseError::NoBits(_))));

        // Invalid bit number
        let res = Composite::from_string("XXX", "H 117356715625188271521875");
        assert!(matches!(res, Err(error::ParseError::InvalidBit(_))));

        // Trailing junk
        let res = Composite::from_string("XXX", "H 0 and something");
        assert!(matches!(res, Err(error::ParseError::TrailingText(_))));

        // Unclosed arguments
        let res = Composite::from_string("XXX", "RX(1.2a) 1");
        assert!(matches!(res, Err(error::ParseError::UnclosedParentheses(_))));

        let res = Composite::from_string("XXX", "RX(1.2*(1+2 1");
        assert!(matches!(res, Err(error::ParseError::UnclosedParentheses(_))));

        let res = Composite::from_string("XXX", "RX(sin(1.2 1");
        assert!(matches!(res, Err(error::ParseError::UnclosedParentheses(_))));
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
        gate.latex_checked(&[0, 1], &mut state);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qw & \ctrl{1} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{H} & \targ & \gate{H} & \qw \\
}
"#);

        let gate = Composite::from_string("CZ", "H 1; CX 0 1; H 1").unwrap();
        let mut state = LatexExportState::new(2, 0);
        state.set_expand_composite(false);
        gate.latex_checked(&[0, 1], &mut state);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \multigate{1}{CZ} & \qw \\
    \lstick{\ket{0}} & \ghost{CZ} & \qw \\
}
"#);

        let gate = Composite::from_string("CZ", "H 1; CX 0 1; H 1").unwrap();
        let mut state = LatexExportState::new(3, 0);
        state.set_expand_composite(false);
        gate.latex_checked(&[0, 2], &mut state);
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
        gate.latex_checked(&[0, 2, 3], &mut state);
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
        gate.latex_checked(&[0, 1, 3], &mut state);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \multigate{1}{G} & \qw \\
    \lstick{\ket{0}} & \ghost{G} & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{G} \qwx[-2] & \qw \\
}
"#);
    }
}
