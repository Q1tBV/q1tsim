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

/// Structure to build up contents of LaTeX export
///
/// Struct `LatexExportState` is used to build up the matrix containing the
/// Qcircuit code for the export of a `Circuit` to LaTeX.
pub struct LatexExportState
{
    // Variables relating to the circuit

    /// The number of quantum bits in the circuit.
    nr_qbits: usize,
    /// The number of classical bits in the circuit.
    nr_cbits: usize,

    // Settings for output

    /// If `true` (the default), add initialization of the bits to the circuit.
    add_init: bool,
    /// If `true` (the default), composite gates are expanded into primitive
    /// gates in the export.
    expand_composite: bool,

    // Runtime variables

    /// Matrix containing LaTex code for each individual gate. Every row in
    /// the matrix corresponds to a column in the exported circuit.
    matrix: Vec<Vec<Option<String>>>,
    /// Vector containing which fields in the last row are currently occupied.
    /// Unoccupied fields can be used, if a gate operates on an occupied field,
    /// a new row must be added.
    in_use: Vec<bool>,
    /// Whether the gate currently being exported is being controlled or not.
    /// This has an impact on the representation of certain gates, e.g. the `X`
    /// gate which is normally represented as a boxed X character, but is
    /// represented by a circled plus (⊕) when it is controlled.
    controlled: bool,
    /// start and end row, and nr of iterations, of static loops.
    loops: Vec<(usize, usize, usize)>,
    /// Start index and nr of iterations of currently unfinished static loops.
    /// Vector because noops may be nested.
    open_loops: Vec<(usize, usize)>
}

impl LatexExportState
{
    /// Create a new LatexExportState
    ///
    /// Create a new `LatexExportState`, for a circuit with `nr_qbits` quantum
    /// bits and `nr_cbits` classical bits.
    pub fn new(nr_qbits: usize, nr_cbits: usize) -> Self
    {
        LatexExportState
        {
            nr_qbits: nr_qbits,
            nr_cbits: nr_cbits,
            add_init: true,
            expand_composite: true,
            matrix: vec![],
            in_use: vec![true; nr_qbits + nr_cbits],
            controlled: false,
            loops: vec![],
            open_loops: vec![]
        }
    }

    /// The total number of bits (quantum or classical) in the circuit.
    fn total_nr_bits(&self) -> usize
    {
        self.nr_qbits + self.nr_cbits
    }

    /// Add a new column.
    ///
    /// Add a new column to the export. Used when a new gate operates on a bit
    /// that is already in use.
    fn add_column(&mut self)
    {
        let nr_bits = self.total_nr_bits();
        self.matrix.push(vec![None; nr_bits]);
        self.in_use.clear();
        self.in_use.resize(nr_bits, false);
    }

    /// Ensure that fields are free.
    ///
    /// Ensure that the fields for the bits in `qbits` and (optionally) `cbits`
    /// are currently unoccupied. If not, add a new column to the export.
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

    /// Ensure that fields are free.
    ///
    /// Ensure that the fields for the bits in `qbits` and (optionally) `cbits`,
    /// as weel as all field in the range between the minimum and maximum bit,
    /// are currently unoccupied. If not, add a new column to the export.
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

    /// Ensure all fields are free.
    ///
    /// Ensure that the last column is currently fully empty. If not, add a new
    /// column to the export.
    fn reserve_all(&mut self)
    {
        if self.in_use.contains(&true)
        {
            self.add_column();
        }
    }

    /// Mark fields as in use.
    ///
    /// Mark the fields corresponding to the quantum bits in `qbits` and
    /// optionally the classical bits in `cbits`m as well as all other bits
    /// between them, as being currently in use. This is usually done for
    /// operations like controlled gates, which connect the control bit
    /// with a controlled operation bit, and for which no operation should be
    /// drawn between them.
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

    /// Set the contents of a field
    ///
    /// Set the contents of the field corresponding to bit `bit` to the LaTeX
    /// code in `contents`.
    pub fn set_field(&mut self, bit: usize, contents: String)
    {
        // Don't crash when user forgets to rerserve space
        if self.matrix.is_empty()
        {
            self.add_column();
        }

        let col = self.matrix.last_mut().unwrap();
        col[bit] = Some(contents);
        self.in_use[bit] = true;
    }

    /// Add a measurement
    ///
    /// Add a measurement of quantum bit `qbit` to classical bit `cbit` in basis
    /// `basis` to the export. If `basis` is `None`, no basis string is drawn in
    /// the measurement.
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

    /// Add a reset
    ///
    /// Add the reset of quantum bit `qbit` to the export.
    pub fn set_reset(&mut self, qbit: usize)
    {
        self.reserve(&[qbit], None);
        self.set_field(qbit, String::from(r"\push{~\ket{0}~} \ar @{|-{}} [0,-1]"));
    }

    /// Add classical control
    ///
    /// Add the control of an operation on quantum bits `qbits` by classical
    /// bits control to the export state. This function only adds the control
    /// part, the actual quantum operation should be drawn elsewhere. The bits
    /// in `control` make up a register, whose value should match `target`.
    /// The first bit in `control` corresponds to the least significant bit of
    /// `target`, the last bit in `control` to the most significant bit.
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

    /// Open a loop
    ///
    /// Open a loop of `count` ieterations at the current row in the export
    /// state. This loop should later be closed by a call to `end_loop()`, at
    /// which point the loop will be added to the export state.
    pub fn start_loop(&mut self, count: usize)
    {
        self.reserve_all();
        self.open_loops.push((self.matrix.len() - 1, count));
    }

    /// Close a loop
    ///
    /// Close the loop opened last by a call to `start_loop()`.
    pub fn end_loop(&mut self)
    {
        if let Some((start, count)) = self.open_loops.pop()
        {
            let end = self.matrix.len() - 1;
            self.loops.push((start, end, count));
            self.reserve_all();
        }
        else
        {
            panic!("Unable to close loop, because no loop is currently open");
        }
    }

    /// Add dots
    ///
    /// This function adds the string in `label` in the middle of the range of
    /// qbits starting at `bit` and going `count` bits down. This is usually
    /// used to add the dots used to indicate a repeated subcircuit in loops.
    pub fn add_cds(&mut self, bit: usize, count: usize, label: &str)
    {
        self.reserve_all();
        self.set_field(bit, format!(r"\cds{{{}}}{{{}}}", count, label));
        self.reserve_all();
    }

    /// Add a barrier
    ///
    /// Add a barrier for the quantum bits in `qbits`. Note that the placement
    /// of barriers may sometimes be off because the spacing between elements
    /// is not constant. It may therefore need some manual adjustment.
    pub fn set_barrier(&mut self, qbits: &[usize])
    {
        let ranges = support::get_ranges(qbits);

        self.add_column();
        for (first, last) in ranges
        {
            self.set_field(first, format!(r"\qw \barrier{{{}}}", last - first))
        }
    }

    /// Export to LaTeX
    ///
    /// This code exports the matrix that was built up in this state to LaTeX
    /// code. It uses the qcircuit package to do so.
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

        let last_col_used = self.in_use.contains(&true);
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

            if last_col_used
            {
                res += r" & ";
                res += if i < self.nr_qbits { r"\qw" } else { r"\cw" };
            }
            res += " \\\\\n";
        }
        res += "}\n";

        res
    }

    /// Set whether gates are controlled
    ///
    /// This sets the option to draw gates in their normal layout
    /// (`controlled = false`) or in their controlled layout (when
    /// `controlled = true`). This can make a difference for e.g. the X gate,
    /// which is drawn as a boxed X character normally, but as an exclusive
    /// or symbol (⊕) when used in a `CX` gate.
    pub fn set_controlled(&mut self, controlled: bool) -> bool
    {
        let res = self.controlled;
        self.controlled = controlled;
        res
    }

    /// Return whether gates should currently be drawn normally, or in their
    /// controlled form.
    pub fn is_controlled(&self) -> bool
    {
        self.controlled
    }

    /// Set whether to expand composite gates.
    ///
    /// Set whether composite gates should be drawn as individual components.
    /// If `expand` is `true`, composite gates are drawn by drawing their
    /// components. If `expand` is false, composite gates are drawn as a single
    /// block gate.
    pub fn set_expand_composite(&mut self, expand: bool)
    {
        self.expand_composite = expand;
    }

    /// Whether to expand composite gates.
    ///
    /// Return whether composite gates should be drawn as individual components
    /// (in which case `true` is returned), or as a single, possibly multi-bit,
    /// operation (when the result is `false`).
    pub fn expand_composite(&self) -> bool
    {
        self.expand_composite
    }
}

/// Trait for gates that can be drawn in LaTeX
pub trait Latex
{
    /// Add this gate to the export state.
    ///
    /// Add the execution of this gate on the bits in `bits`, to the export
    /// state `state`.
    fn latex(&self, bits: &[usize], state: &mut LatexExportState);

    /// Checked add to the export state.
    ///
    /// This function should first check if the fields needed for drawing this
    /// gate are free, and if not, add a new row in the export state `state`.
    /// The default implementation merely check if the fields corresponding to
    /// the bits in `bits` are free. Gates that need other fields free as well
    /// (e.g. controlled gates in which all fields between the control and the
    /// operation are occupied as well), should provide their own implementation
    /// of this function.
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
