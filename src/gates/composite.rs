use cmatrix;
use gates;
extern crate regex;

use super::*;

/// Strucutre for a description of a subgate.
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
    gate: Box<gates::Gate>,
    /// The bits on which the gate acts
    bits: Vec<usize>
}

impl SubGate
{
    /// Create a new composite gate operation
    fn new<G>(gate: G, bits: &[usize]) -> Self
    where G: 'static + gates::Gate
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
    /// parsed. On failure, return an error string.
    fn parse_gate_name(desc: &str) -> Result<(&str, &str), String>
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
            Err(format!("Failed to find gate name in \"{}\"", desc))
        }
    }

    /// Parse arguments to a subgate.
    ///
    /// Parse arguments to a subgate, if any, from description string `desc`.
    /// If no parenthesized argument list is found, an emmpty argument vector
    /// is returned. If there is an argument list, then if it can be parsed
    /// successfully, the arguments are returned,¸together with the rest of the
    /// description string that needs to be parsed for bit numbers. On failure,
    /// an error string is returned.
    fn parse_gate_args(desc: &str) -> Result<(Vec<f64>, &str), String>
    {
        let re = regex::Regex::new(r"^\s*\(\s*([^\)]*)\s*\)").unwrap();
        if let Some(captures) = re.captures(desc)
        {
            let m = captures.get(0).unwrap();
            let rest = &desc[m.end()..];
            let mut args = vec![];

            for arg_txt in captures[1].split(',')
            {
                if let Ok(arg) = arg_txt.trim().parse()
                {
                    args.push(arg);
                }
                else
                {
                    return Err(format!("Failed to parse argument \"{}\"", arg_txt));
                }
            }

            Ok((args, rest))
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
    /// description string on success, or an error string on failure.
    fn parse_gate_bits<'a>(desc: &'a str, name: &str) -> Result<(Vec<usize>, &'a str), String>
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
                return Err(format!("Failed to parse bit number in \"{}\"", bit_txt));
            }
        }

        if bits.is_empty()
        {
            Err(format!("Unable to find the bits gate {} operates on", name))
        }
        else
        {
            Ok((bits, rest))
        }
    }

    /// Parse a gate description.
    ///
    /// Parse the subgate description string `desc`. Returns the subgate
    /// description on success, or an error string on failure.
    fn parse_gate_desc(desc: &str) -> Result<SubGateDesc, String>
    {
        let (name, rest) = Self::parse_gate_name(desc)?;
        let (args, rest) = Self::parse_gate_args(rest)?;
        let (bits, rest) = Self::parse_gate_bits(rest, name)?;

        let rest = rest.trim();
        if !rest.is_empty()
        {
            Err(format!("Trailing text after gate description: \"{}\"", rest))
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
    //// `nr_bits`. Return an error string on failure.
    fn assert_nr_args_bits(nr_args: usize, nr_bits: usize, desc: &SubGateDesc)
        -> Result<(), String>
    {
        if nr_args != desc.args.len()
        {
            Err(format!("Invalid number of arguments for \"{}\" gate", desc.name))
        }
        else if nr_bits != desc.bits.len()
        {
            Err(format!("Invalid number of bits for \"{}\" gate", desc.name))
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
    ///   white space.
    /// * Currently, only real numbers are allowed for parameters.
    /// Examples:
    /// ```text
    /// H 1; CX 0 1; H 1
    /// RY(4.7124) 1; CX 1 0; RY(1.5708) 1; X1
    /// ```
    pub fn from_string(name: &str, desc: &str) -> Result<Self, String>
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
                _ => { return Err(format!("Unknown gate \"{}\"", gate.name)); }
            }
        }

        Ok(composite)
    }

    /// Add a gate.
    ///
    /// Append a `n`-ary subgate `gate`, operating on the `n` qubits in `bits`,
    /// to this composite gate.
    pub fn add_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: gates::Gate
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
            res = op.gate.expanded_matrix(&op.bits, self.nr_bits).dot(&res);
        }

        res
    }
}

#[cfg(test)]
mod tests
{
    use cmatrix;
    use super::Composite;
    use gates::{Gate, CCX, CX, H, X};

    #[test]
    fn test_description()
    {
        let gate = Composite::new("G", 3);
        assert_eq!(gate.description(), "G");
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
            Err(err) => { panic!(err); }
        }

        match Composite::from_string("Y", "U3(3.141592653589793,1.570796326794897,1.570796326794897) 0")
        {
            Ok(gate) => {
                assert_eq!(gate.description(), "Y");
                assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [i, z]]);
            },
            Err(err) => { panic!(err); }
        }
    }
}
