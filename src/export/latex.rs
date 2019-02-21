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

use support;

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

    fn reserve_all(&mut self)
    {
        if self.in_use.contains(&true)
        {
            self.add_column();
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
        self.reserve_all();
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
        self.reserve_all();
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

    pub fn set_barrier(&mut self, bits: &[usize])
    {
        let ranges = support::get_ranges(bits);

        self.add_column();
        for (first, last) in ranges
        {
            self.set_field(first, format!(r"\qw \barrier{{{}}}", last - first))
        }
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

/// Trait for gates that can be drawn in LaTeX
pub trait Latex
{
    fn latex(&self, bits: &[usize], state: &mut LatexExportState);

    fn latex_checked(&self, bits: &[usize], state: &mut LatexExportState)
    {
        state.reserve(bits, None);
        self.latex(bits, state);
    }
}

#[cfg(test)]
mod tests
{
    use gates::H;

    use super::Latex;
}
