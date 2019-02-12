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


use gates;

use gates::Gate;

/// Trait for gates that can be represented in OpenQasm.
pub trait OpenQasm
{
    /// OpenQasm representation
    ///
    /// Return an OpenQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits.
    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String;

    /// OpenQasm representation of conditional gate.
    ///
    /// Return the OpenQasm representation of a gate that is only executed when
    /// the condition `condition` on the classical bits of the program state
    /// holds. The default implementation only works for a single gate,
    /// composite gates (like `Composite` or `Kron`) should overwrite this
    /// default. On success, returns `Ok` with the instruction string. On error,
    /// returns `Err` with an error message.
    fn conditional_open_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> Result<String, String>
    {
        Ok(format!("if ({}) {}", condition, self.open_qasm(bit_names, bits)))
    }
}

/// Trait for gates that can be represented in c-Qasm.
pub trait CQasm
{
    /// cQasm representation
    ///
    /// Return an cQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits.
    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String;

    /// cQasm representation of conditional gate.
    ///
    /// Return the cQasm representation of a gate that is only executed when
    /// the condition `condition` on the classical bits of the program state
    /// holds. The default implementation only works for a single gate,
    /// composite gates (like `Composite` or `Kron`) should overwrite this
    /// default. On success, returns `Ok` with the instruction string. On error,
    /// returns `Err` with an error message.
    fn conditional_c_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> Result<String, String>
    {
        let unc_qasm = self.c_qasm(bit_names, bits);
        let parts: Vec<&str> = unc_qasm.splitn(2, " ").collect();
        if parts.len() != 2
        {
            Err(format!("Unable to find gate name or argument in \"{}\"", unc_qasm))
        }
        else
        {
            Ok(format!("c-{} {}, {}", parts[0], condition, parts[1]))
        }
    }
}

pub struct LatexExportState
{
    matrix: Vec<Vec<Option<String>>>,
    in_use: Vec<bool>,
    nr_qbits: usize,
    nr_cbits: usize,
    // Settings for output
    add_init: bool,
    expand_composite: bool,
    // Runtime variables
    controlled: bool,
    loops: Vec<(usize, usize, usize)>,
    open_loops: Vec<(usize, usize)>
}

impl LatexExportState
{
    pub fn new(nr_qbits: usize, nr_cbits: usize) -> Self
    {
        LatexExportState
        {
            matrix: vec![],
            in_use: vec![true; nr_qbits + nr_cbits],
            nr_qbits: nr_qbits,
            nr_cbits: nr_cbits,
            add_init: true,
            expand_composite: true,
            controlled: false,
            loops: vec![],
            open_loops: vec![]
        }
    }

    fn total_nr_bits(&self) -> usize
    {
        self.nr_qbits + self.nr_cbits
    }

    fn add_column(&mut self)
    {
        let nr_bits = self.total_nr_bits();
        self.matrix.push(vec![None; nr_bits]);
        self.in_use.clear();
        self.in_use.resize(nr_bits, false);
    }

    pub fn reserve(&mut self, qbits: &[usize], cbits: Option<&[usize]>)
    {
        let mut bits = qbits.to_vec();
        if let Some(cbs) = cbits
        {
            bits.extend(cbs.iter().map(|&b| self.nr_qbits + b));
        }

        if bits.iter().any(|&b| self.in_use[b])
        {
            self.add_column();
        }
    }

    pub fn reserve_range(&mut self, qbits: &[usize], cbits: Option<&[usize]>)
    {
        let mut bits = qbits.to_vec();
        if let Some(cbs) = cbits
        {
            bits.extend(cbs.iter().map(|&b| self.nr_qbits + b));
        }

        if let Some(&start) = bits.iter().min()
        {
            let end = *bits.iter().max().unwrap();
            if self.in_use[start..end+1].contains(&true)
            {
                self.add_column();
            }
        }
    }

    pub fn claim_range(&mut self, qbits: &[usize], cbits: Option<&[usize]>)
    {
        let mut bits = qbits.to_vec();
        if let Some(cbs) = cbits
        {
            bits.extend(cbs.iter().map(|&b| self.nr_qbits + b));
        }

        if let Some(&start) = bits.iter().min()
        {
            let end = *bits.iter().max().unwrap();
            for bit in start..end
            {
                self.in_use[bit] = true;
            }
        }
    }

    pub fn set_field(&mut self, bit: usize, contents: String)
    {
        let col = self.matrix.last_mut().unwrap();
        col[bit] = Some(contents);
        self.in_use[bit] = true;
    }

    pub fn set_measurement(&mut self, qbit: usize, cbit: usize, basis: Option<&str>)
    {
        let cbit_idx = self.nr_qbits + cbit;
        self.reserve_range(&[qbit], Some(&[cbit]));
        let meter = if let Some(b) = basis
            {
                format!(r"\meterB{{{}}}", b)
            }
            else
            {
                String::from(r"\meter")
            };
        self.set_field(qbit, meter);
        self.set_field(cbit_idx, format!(r"\cw \cwx[{}]", qbit as isize - cbit_idx as isize));
        self.claim_range(&[qbit], Some(&[cbit]));
    }

    pub fn set_reset(&mut self, qbit: usize)
    {
        self.reserve(&[qbit], None);
        self.set_field(qbit, String::from(r"\push{~\ket{0}~} \ar @{|-{}} [0,-1]"));
    }

    pub fn set_condition(&mut self, control: &[usize], target: u64, qbits: &[usize])
    {
        if qbits.is_empty()
        {
            return;
        }

        let mut pbit = *qbits.iter().max().unwrap();
        let mut bp: Vec<(usize, usize)> = control.iter().enumerate()
            .map(|(pos, &idx)| (self.nr_qbits + idx, pos))
            .collect();
        bp.sort();
        for (bit, pos) in bp
        {
            let ctrl = if (target & (1 << pos)) == 0 { r"\cctrlo" } else { r"\cctrl" };
            self.set_field(bit, format!("{}{{{}}}", ctrl, pbit as isize - bit as isize));
            pbit = bit;
        }
    }

    pub fn start_loop(&mut self, count: usize)
    {
        let last_cbit = self.nr_cbits - 1;
        self.reserve_range(&[0], Some(&[last_cbit]));
        self.open_loops.push((self.matrix.len() - 1, count));
    }

    pub fn end_loop(&mut self)
    {
        if let Some((start, count)) = self.open_loops.pop()
        {
            let end = self.matrix.len() - 1;
            self.loops.push((start, end, count));
            self.add_column();
        }
        else
        {
            panic!("Unable to close loop, because no loop is currently open");
        }
    }

    pub fn add_cds(&mut self, bit: usize, count: usize, label: &str)
    {
        let last_cbit = self.nr_cbits - 1;
        self.reserve_range(&[0], Some(&[last_cbit]));
        self.set_field(bit, format!(r"\cds{{{}}}{{{}}}", count, label));
        self.add_column();
    }

    pub fn code(&self) -> String
    {
        let mut res = String::from("\\Qcircuit @C=1em @R=.7em {\n");

        if !self.loops.is_empty()
        {
            let mut prev_idx = 0;
            res += r"    & ";
            for (start, end, count) in self.loops.iter()
            {
                res += r"& ".repeat(start - prev_idx).as_str();
                res += format!("\\mbox{{}} \\POS\"{},{}\".\"{},{}\".\"{},{}\".\"{},{}\"!C*+<.7em>\\frm{{^\\}}}},+U*++!D{{{}\\times}}",
                    2, start+2, 2, start+2, /*self.total_nr_bits()+*/2, end+2,
                    /*self.total_nr_bits()+*/2, end+2, count).as_str();
                prev_idx = *start;
            }
            res += "\\\\\n";

            res += r"    ";
            res += r"& ".repeat(self.matrix.len()).as_str();
            res += "\\\\\n";
        }

        for i in 0..self.total_nr_bits()
        {
            if self.add_init
            {
                if i < self.nr_qbits
                {
                    res += r"    \lstick{\ket{0}}";
                }
                else
                {
                    res += r"    \lstick{0}";
                }
            }
            else
            {
                res += r"    ";
            }
            for row in self.matrix.iter()
            {
                res += " & ";
                if let Some(ref s) = row[i]
                {
                    res += s.as_str();
                }
                else if i < self.nr_qbits
                {
                    res += r"\qw";
                }
                else
                {
                    res += r"\cw";
                }
            }

            res += r" & ";
            res += if i < self.nr_qbits { r"\qw " } else { r"\cw " };
            res += "\\\\\n";
        }
        res += "}\n";

        res
    }

    pub fn set_controlled(&mut self, controlled: bool) -> bool
    {
        let res = self.controlled;
        self.controlled = controlled;
        res
    }

    pub fn is_controlled(&self) -> bool
    {
        self.controlled
    }

    pub fn expand_composite(&self) -> bool
    {
        self.expand_composite
    }
}

// Trait for gates that can be drawn in LaTeX
pub trait Latex
{
    fn latex(&self, bits: &[usize], state: &mut LatexExportState);

    fn latex_checked(&self, bits: &[usize], state: &mut LatexExportState)
    {
        state.reserve(bits, None);
        self.latex(bits, state);
    }
}

/// Trait combining the traits necessary for a gate in a quantum circuit
pub trait CircuitGate: gates::Gate + OpenQasm + CQasm + Latex {}

impl<G: Gate + OpenQasm + CQasm + Latex> CircuitGate for G {}


#[cfg(test)]
mod tests
{
    use gates::H;

    use super::{CQasm, OpenQasm};

    #[test]
    fn test_conditional_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];

        let res = H::new().conditional_open_qasm("b == 0", &bit_names, &[1]);
        assert_eq!(res, Ok(String::from("if (b == 0) h qb1")));
    }

    #[test]
    fn test_conditional_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];

        let res = H::new().conditional_c_qasm("b[0]", &bit_names, &[1]);
        assert_eq!(res, Ok(String::from("c-h b[0], qb1")));
    }
}
