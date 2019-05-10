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
    /// When the state is currently in the process of drawing an operation
    /// operating on a range of bits with possible gaps, this holds the minimum
    /// and maximum bit index of said range.
    reserved_ranges: Vec<(usize, usize)>,
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
            reserved_ranges: vec![],
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

    fn get_bit_indices(&self, qbits: &[usize], cbits: Option<&[usize]>)
        -> crate::error::Result<Vec<usize>>
    {
        if let Some(&bit) = qbits.iter().find(|&&b| b >= self.nr_qbits)
        {
            return Err(crate::error::Error::InvalidQBit(bit));
        }

        let mut bits = qbits.to_vec();
        if let Some(cbs) = cbits
        {
            if let Some(&bit) = cbs.iter().find(|&&b| b >= self.nr_cbits)
            {
                return Err(crate::error::Error::InvalidCBit(bit));
            }
            bits.extend(cbs.iter().map(|&b| self.nr_qbits + b));
        }

        Ok(bits)
    }

    /// Ensure that fields are free.
    ///
    /// Ensure that the fields for the bits in `qbits` and (optionally) `cbits`
    /// are currently unoccupied. If not, add a new column to the export.
    pub fn reserve(&mut self, qbits: &[usize], cbits: Option<&[usize]>)
        -> crate::error::Result<()>
    {
        let bits = self.get_bit_indices(qbits, cbits)?;
        if bits.iter().any(|&b| self.in_use[b])
        {
            self.add_column();
        }

        Ok(())
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

    /// Ensure that fields are free.
    ///
    /// Ensure that the fields for the bits in `qbits` and (optionally) `cbits`,
    /// as well as all field in the range between the minimum and maximum bit,
    /// are currently unoccupied. If not, add a new column to the export.
    pub fn start_range_op(&mut self, qbits: &[usize], cbits: Option<&[usize]>)
        -> crate::error::Result<()>
    {
        let bits = self.get_bit_indices(qbits, cbits)?;
        if let Some(&first) = bits.iter().min()
        {
            let last = *bits.iter().max().unwrap();

            if self.reserved_ranges.is_empty()
            {
                if self.in_use[first..=last].contains(&true)
                {
                    self.add_column();
                }
                self.reserved_ranges.push((first, last));
            }
            else
            {
                let &(outer_first, outer_last) = self.reserved_ranges.last().unwrap();
                if outer_first <= first && outer_last >= last
                {
                    self.reserved_ranges.push((first, last));
                }
                else
                {
                    return Err(crate::error::Error::from(crate::error::ExportError::RangeAlreadyOpen));
                }
            }

        }

        Ok(())
    }

    /// Mark fields as in use.
    ///
    /// Mark the fields in the fields in the range last opened by
    /// `start_range_op()` as in use, and close the range.
    pub fn end_range_op(&mut self)
    {
        if let Some((first, last)) = self.reserved_ranges.pop()
        {
            for bit in first..=last
            {
                self.in_use[bit] = true;
            }
        }
    }

    /// Set the contents of a field
    ///
    /// Set the contents of the field corresponding to bit `bit` to the LaTeX
    /// code in `contents`.
    pub fn set_field(&mut self, bit: usize, contents: String) -> crate::error::Result<()>
    {
        if self.reserved_ranges.is_empty()
        {
            self.reserve(&[bit], None)?;
        }

        let col = self.matrix.last_mut().unwrap();
        col[bit] = Some(contents);
        self.in_use[bit] = true;
        Ok(())
    }

    /// Add a measurement
    ///
    /// Add a measurement of quantum bit `qbit` to classical bit `cbit` in basis
    /// `basis` to the export. If `basis` is `None`, no basis string is drawn in
    /// the measurement.
    pub fn set_measurement(&mut self, qbit: usize, cbit: usize, basis: Option<&str>)
        -> crate::error::Result<()>
    {
        let cbit_idx = self.nr_qbits + cbit;
        self.start_range_op(&[qbit], Some(&[cbit]))?;
        let meter = if let Some(b) = basis
            {
                format!(r"\meterB{{{}}}", b)
            }
            else
            {
                String::from(r"\meter")
            };
        self.set_field(qbit, meter)?;
        self.set_field(cbit_idx, format!(r"\cw \cwx[{}]", qbit as isize - cbit_idx as isize))?;
        self.end_range_op();

        Ok(())
    }

    /// Add a reset
    ///
    /// Add the reset of quantum bit `qbit` to the export.
    pub fn set_reset(&mut self, qbit: usize) -> crate::error::Result<()>
    {
        self.set_field(qbit, String::from(r"\push{~\ket{0}~} \ar @{|-{}} [0,-1]"))
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
        -> crate::error::Result<()>
    {
        if let Some(&bit) = qbits.iter().find(|&&b| b >= self.nr_qbits)
        {
            return Err(crate::error::Error::InvalidQBit(bit));
        }
        if let Some(&bit) = control.iter().find(|&&b| b >= self.nr_cbits)
        {
            return Err(crate::error::Error::InvalidCBit(bit));
        }

        if qbits.is_empty()
        {
            return Ok(());
        }

        let mut pbit = *qbits.iter().max().unwrap();
        let mut bp: Vec<(usize, usize)> = control.iter().enumerate()
            .map(|(pos, &idx)| (self.nr_qbits + idx, pos))
            .collect();
        bp.sort();
        for (bit, pos) in bp
        {
            let ctrl = if (target & (1 << pos)) == 0 { r"\cctrlo" } else { r"\cctrl" };
            self.set_field(bit, format!("{}{{{}}}", ctrl, pbit as isize - bit as isize))?;
            pbit = bit;
        }

        Ok(())
    }

    /// Add a simple block gate
    ///
    /// Add a gate with description `desc` operating on qubits `qbits` as a
    /// simple box containing the description. If the gate opertes on mutiple
    /// disjoint ranges of qbits, it is drawn using mutiple boxes connected by
    /// wires.
    pub fn add_block_gate(&mut self, qbits: &[usize], desc: &str) -> crate::error::Result<()>
    {
        let ranges = crate::support::get_ranges(qbits);
        if !ranges.is_empty()
        {
            self.start_range_op(qbits, None)?;

            let (first, last) = ranges[0];
            if last == first
            {
                self.set_field(first, format!(r"\gate{{{}}}", desc))?;
            }
            else
            {
                self.set_field(first,
                    format!(r"\multigate{{{}}}{{{}}}", last-first, desc))?;
                for bit in first+1..last+1
                {
                    self.set_field(bit, format!(r"\ghost{{{}}}", desc))?;
                }
            }

            let mut prev_last = last;
            for &(first, last) in ranges[1..].iter()
            {
                if last == first
                {
                    self.set_field(first,
                        format!(r"\gate{{{}}} \qwx[{}]", desc, prev_last as isize - first as isize))?;
                }
                else
                {
                    self.set_field(first,
                        format!(r"\multigate{{{}}}{{{}}} \qwx[{}]",
                            last - first, desc, prev_last as isize - first as isize))?;
                    for bit in first+1..last+1
                    {
                        self.set_field(bit, format!(r"\ghost{{{}}}", desc))?;
                    }
                }

                prev_last = last;
            }

            self.end_range_op();
        }

        Ok(())
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
    pub fn end_loop(&mut self) -> crate::error::Result<()>
    {
        if let Some((start, count)) = self.open_loops.pop()
        {
            let end = self.matrix.len() - 1;
            self.loops.push((start, end, count));
            self.reserve_all();
            Ok(())
        }
        else
        {
            Err(crate::error::Error::from(crate::error::ExportError::CantCloseLoop))
        }
    }

    /// Add dots
    ///
    /// This function adds the string in `label` in the middle of the range of
    /// qbits starting at `bit` and going `count` bits down. This is usually
    /// used to add the dots used to indicate a repeated subcircuit in loops.
    pub fn add_cds(&mut self, bit: usize, count: usize, label: &str) -> crate::error::Result<()>
    {
        self.reserve_all();
        self.set_field(bit, format!(r"\cds{{{}}}{{{}}}", count, label))?;
        self.reserve_all();
        Ok(())
    }

    /// Add a barrier
    ///
    /// Add a barrier for the quantum bits in `qbits`. Note that the placement
    /// of barriers may sometimes be off because the spacing between elements
    /// is not constant. It may therefore need some manual adjustment.
    pub fn set_barrier(&mut self, qbits: &[usize]) -> crate::error::Result<()>
    {
        if let Some(&bit) = qbits.iter().find(|&&b| b >= self.nr_qbits)
        {
            return Err(crate::error::Error::InvalidQBit(bit));
        }

        let ranges = crate::support::get_ranges(qbits);

        self.add_column();
        for (first, last) in ranges
        {
            self.set_field(first, format!(r"\qw \barrier{{{}}}", last - first))?;
        }

        Ok(())
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

    /// Set whether to add initialization strings
    ///
    /// When `add_init` is `true`, initialization strings are added to the
    /// bits, when `false` they are omitted.
    pub fn set_add_init(&mut self, add_init: bool)
    {
        self.add_init = add_init;
    }
}

/// Trait for gates that can be drawn in LaTeX
pub trait Latex: crate::gates::Gate
{
    /// Add this gate to the export state.
    ///
    /// Add the execution of this gate on the bits in `bits`, to the export
    /// state `state`. The default implementation returns a NotImplemented error.
    fn latex(&self, bits: &[usize], state: &mut LatexExportState)
        -> crate::error::Result<()>
    {
        self.check_nr_bits(bits.len())?;
        state.add_block_gate(bits, self.description())
    }
}

#[cfg(test)]
mod tests
{
    use super::{LatexExportState, Latex};

    struct NoLatexGate;

    impl crate::gates::Gate for NoLatexGate
    {
        fn cost(&self) -> f64 { 0.0 }
        fn description(&self) -> &str { "NLG" }
        fn nr_affected_bits(&self) -> usize { 3 }
        fn matrix(&self) -> crate::cmatrix::CMatrix { unimplemented!() }
    }
    impl Latex for NoLatexGate {}

    #[test]
    fn test_new()
    {
        let nr_qbits = 5;
        let nr_cbits = 2;
        let state = LatexExportState::new(nr_qbits, nr_cbits);
        assert_eq!(state.nr_qbits, nr_qbits);
        assert_eq!(state.nr_cbits, nr_cbits);
        assert_eq!(state.add_init, true);
        assert_eq!(state.expand_composite, true);
        assert_eq!(state.matrix, Vec::<Vec<Option<String>>>::new());
        assert_eq!(state.in_use, vec![true; nr_qbits+nr_cbits]);
        assert_eq!(state.controlled, false);
        assert_eq!(state.loops, vec![]);
        assert_eq!(state.open_loops, vec![]);
    }

    #[test]
    fn test_total_nr_bits()
    {
        let state = LatexExportState::new(5, 2);
        assert_eq!(state.total_nr_bits(), 7);

        let state = LatexExportState::new(3, 0);
        assert_eq!(state.total_nr_bits(), 3);

        let state = LatexExportState::new(2, 8);
        assert_eq!(state.total_nr_bits(), 10);
    }

    #[test]
    fn test_add_column()
    {
        let mut state = LatexExportState::new(3, 1);
        state.add_column();
        assert_eq!(state.matrix, vec![vec![None, None, None, None]]);

        state.matrix[0][1] = Some(String::from("x"));
        state.add_column();
        assert_eq!(state.matrix, vec![
            vec![None, Some(String::from("x")), None, None],
            vec![None, None, None, None]
        ]);
    }

    #[test]
    fn test_reserve()
    {
        let mut state = LatexExportState::new(2, 2);
        assert_eq!(state.reserve(&[0], None), Ok(()));
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]]);

        state.in_use[0] = true;
        assert_eq!(state.reserve(&[1], None), Ok(()));
        assert_eq!(state.in_use, vec![true, false, false, false]);
        assert_eq!(state.matrix, vec![vec![None, None, None, None]]);

        assert_eq!(state.reserve(&[0], None), Ok(()));
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 2]);

        state.in_use[3] = true;
        assert_eq!(state.reserve(&[0, 1], None), Ok(()));
        assert_eq!(state.in_use, vec![false, false, false, true]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 2]);

        assert_eq!(state.reserve(&[0, 1], Some(&[0])), Ok(()));
        assert_eq!(state.in_use, vec![false, false, false, true]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 2]);

        assert_eq!(state.reserve(&[0, 1], Some(&[1])), Ok(()));
        assert_eq!(state.in_use, vec![false, false, false, false]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 3]);
    }

    #[test]
    fn test_reserve_all()
    {
        let mut state = LatexExportState::new(2, 2);
        state.reserve_all();
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]]);

        state.reserve_all();
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]]);

        state.in_use[0] = true;
        state.reserve_all();
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 2]);
    }

    #[test]
    fn test_range_op()
    {
        let mut state = LatexExportState::new(2, 2);
        assert_eq!(state.start_range_op(&[0], None), Ok(()));
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]]);
        state.end_range_op();
        assert_eq!(state.in_use, vec![true, false, false, false]);
        assert_eq!(state.matrix, vec![vec![None; 4]]);

        assert_eq!(state.start_range_op(&[0, 1], None), Ok(()));
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 2]);
        state.end_range_op();
        assert_eq!(state.in_use, vec![true, true, false, false]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 2]);

        assert_eq!(state.start_range_op(&[0], None), Ok(()));
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 3]);
        state.end_range_op();
        assert_eq!(state.in_use, vec![true, false, false, false]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 3]);

        assert_eq!(state.start_range_op(&[1], None), Ok(()));
        assert_eq!(state.in_use, vec![true, false, false, false]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 3]);
        state.end_range_op();
        assert_eq!(state.in_use, vec![true, true, false, false]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 3]);

        assert_eq!(state.start_range_op(&[0], Some(&[1])), Ok(()));
        assert_eq!(state.in_use, vec![false; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 4]);
        state.end_range_op();
        assert_eq!(state.in_use, vec![true; 4]);
        assert_eq!(state.matrix, vec![vec![None; 4]; 4]);
    }

    #[test]
    fn test_set_field()
    {
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(state.set_field(0, String::from("hello")), Ok(()));
        assert_eq!(state.in_use, vec![true, false]);
        assert_eq!(state.matrix, vec![
            vec![Some(String::from("hello")), None]
        ]);

        assert_eq!(state.set_field(1, String::from("world")), Ok(()));
        assert_eq!(state.in_use, vec![true, true]);
        assert_eq!(state.matrix, vec![
            vec![Some(String::from("hello")), Some(String::from("world"))]
        ]);

        assert_eq!(state.set_field(0, String::from("hi there")), Ok(()));
        assert_eq!(state.in_use, vec![true, false]);
        assert_eq!(state.matrix, vec![
            vec![Some(String::from("hello")), Some(String::from("world"))],
            vec![Some(String::from("hi there")), None]
        ]);

        state.add_column();
        assert_eq!(state.set_field(1, String::from("planet Mars")), Ok(()));
        assert_eq!(state.in_use, vec![false, true]);
        assert_eq!(state.matrix, vec![
            vec![Some(String::from("hello")), Some(String::from("world"))],
            vec![Some(String::from("hi there")), None],
            vec![None, Some(String::from("planet Mars"))]
        ]);
    }

    #[test]
    fn test_set_measurement()
    {
        let mut state = LatexExportState::new(2, 2);
        assert_eq!(state.set_measurement(0, 1, None), Ok(()));
        assert_eq!(state.set_measurement(1, 0, Some("X")), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \meter & \qw & \qw \\
    \lstick{\ket{0}} & \qw & \meterB{X} & \qw \\
    \lstick{0} & \cw & \cw \cwx[-1] & \cw \\
    \lstick{0} & \cw \cwx[-3] & \cw & \cw \\
}
"#);
    }

    #[test]
    fn test_set_reset()
    {
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(state.set_reset(0), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \push{~\ket{0}~} \ar @{|-{}} [0,-1] & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
}
"#);
    }

    #[test]
    fn test_set_condition()
    {
        let mut state = LatexExportState::new(2, 2);
        assert_eq!(state.start_range_op(&[], None), Ok(()));
        assert_eq!(state.set_condition(&[0, 1], 2, &[]), Ok(()));
        state.end_range_op();

        assert_eq!(state.start_range_op(&[0], Some(&[0, 1])), Ok(()));
        assert_eq!(state.set_field(0, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_condition(&[0, 1], 2, &[0]), Ok(()));
        state.end_range_op();

        assert_eq!(state.start_range_op(&[1], Some(&[0, 1])), Ok(()));
        assert_eq!(state.set_field(1, String::from(r"\gate{H}")), Ok(()));
        assert_eq!(state.set_condition(&[0, 1], 1, &[1]), Ok(()));
        state.end_range_op();

        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{X} & \qw & \qw \\
    \lstick{\ket{0}} & \qw & \gate{H} & \qw \\
    \lstick{0} & \cctrlo{-2} & \cctrl{-1} & \cw \\
    \lstick{0} & \cctrl{-1} & \cctrlo{-1} & \cw \\
}
"#);
    }

    #[test]
    fn test_add_block_gate()
    {
        let mut state = LatexExportState::new(1, 0);
        assert_eq!(state.add_block_gate(&[0], "G"), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{G} & \qw \\
}
"#);

        let mut state = LatexExportState::new(3, 0);
        assert_eq!(state.add_block_gate(&[1, 0, 2], "G"), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \multigate{2}{G} & \qw \\
    \lstick{\ket{0}} & \ghost{G} & \qw \\
    \lstick{\ket{0}} & \ghost{G} & \qw \\
}
"#);

        let mut state = LatexExportState::new(5, 0);
        assert_eq!(state.add_block_gate(&[1, 0, 3], "G"), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \multigate{1}{G} & \qw \\
    \lstick{\ket{0}} & \ghost{G} & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{G} \qwx[-2] & \qw \\
    \lstick{\ket{0}} & \qw & \qw \\
}
"#);
    }

    #[test]
    fn test_loop()
    {
        let mut state = LatexExportState::new(2, 0);
        state.start_loop(23);
        assert_eq!(state.reserve(&[0, 1], None), Ok(()));
        assert_eq!(state.set_field(0, String::from(r"\gate{H}")), Ok(()));
        assert_eq!(state.set_field(1, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.add_cds(0, 1, r"\leftrightarrow"), Ok(()));
        assert_eq!(state.reserve(&[0, 1], None), Ok(()));
        assert_eq!(state.set_field(0, String::from(r"\gate{H}")), Ok(()));
        assert_eq!(state.set_field(1, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.end_loop(), Ok(()));

        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    & \mbox{} \POS"2,2"."2,2"."2,4"."2,4"!C*+<.7em>\frm{^\}},+U*++!D{23\times}\\
    & & & & \\
    \lstick{\ket{0}} & \gate{H} & \cds{1}{\leftrightarrow} & \gate{H} & \qw \\
    \lstick{\ket{0}} & \gate{X} & \qw & \gate{X} & \qw \\
}
"#);
    }

    #[test]
    fn test_loop_close_error()
    {
        let mut state = LatexExportState::new(2, 0);
        assert_eq!(state.end_loop(),
            Err(crate::error::Error::ExportError(crate::error::ExportError::CantCloseLoop))
        );
    }

    #[test]
    fn test_set_barrier()
    {
        let mut state = LatexExportState::new(3, 0);
        assert_eq!(state.set_field(0, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_field(1, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_field(2, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_barrier(&[0]), Ok(()));

        assert_eq!(state.set_field(0, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_field(1, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_field(2, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_barrier(&[0, 2]), Ok(()));

        assert_eq!(state.set_field(0, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_field(1, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_field(2, String::from(r"\gate{X}")), Ok(()));
        assert_eq!(state.set_barrier(&[0, 1, 2]), Ok(()));

        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \gate{X} & \qw \barrier{0} & \gate{X} & \qw \barrier{0} & \gate{X} & \qw \barrier{2} & \qw \\
    \lstick{\ket{0}} & \gate{X} & \qw & \gate{X} & \qw & \gate{X} & \qw & \qw \\
    \lstick{\ket{0}} & \gate{X} & \qw & \gate{X} & \qw \barrier{0} & \gate{X} & \qw & \qw \\
}
"#);
    }

    #[test]
    fn test_no_init()
    {
        let mut state = LatexExportState::new(1, 1);
        assert_eq!(state.set_measurement(0, 0, None), Ok(()));

        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \meter & \qw \\
    \lstick{0} & \cw \cwx[-1] & \cw \\
}
"#);

        state.set_add_init(false);
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
     & \meter & \qw \\
     & \cw \cwx[-1] & \cw \\
}
"#);
    }

    #[test]
    fn test_default_trait_impl()
    {
        let gate = NoLatexGate {};
        let mut state = LatexExportState::new(5, 0);
        assert_eq!(gate.latex(&[1, 2, 4], &mut state), Ok(()));
        assert_eq!(gate.latex(&[0, 2, 3], &mut state), Ok(()));
        assert_eq!(state.code(),
r#"\Qcircuit @C=1em @R=.7em {
    \lstick{\ket{0}} & \qw & \gate{NLG} & \qw \\
    \lstick{\ket{0}} & \multigate{1}{NLG} & \qw & \qw \\
    \lstick{\ket{0}} & \ghost{NLG} & \multigate{1}{NLG} \qwx[-2] & \qw \\
    \lstick{\ket{0}} & \qw & \ghost{NLG} & \qw \\
    \lstick{\ket{0}} & \gate{NLG} \qwx[-2] & \qw & \qw \\
}
"#);
    }
}
