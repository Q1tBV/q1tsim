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

use crate::stabilizer::{MeasurementInfo, StabilizerTableau};

pub struct StabilizerState
{
    /// The number of qubits in this state
    nr_bits: usize,
    /// The number of separate runs for evolving this state
    nr_shots: usize,
    /// Run counts for each quantum state
    counts: Vec<usize>,
    /// The quantum states themselves
    tableaus: Vec<StabilizerTableau>
}

impl StabilizerState
{
    /// Create a new quantum state of `nr_bits` qubits, all initialized to |0⟩,
    /// which will be measured `nr_shots` times.
    pub fn new(nr_bits: usize, nr_shots: usize) -> Self
    {
        StabilizerState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            counts: vec![nr_shots],
            tableaus: vec![StabilizerTableau::new(nr_bits)]
        }
    }
}

impl crate::qustate::QuState for StabilizerState
{
    fn apply_gate<G>(&mut self, gate: &G, bits: &[usize]) -> crate::error::Result<()>
    where G: crate::gates::Gate + ?Sized
    {
        for t in self.tableaus.iter_mut()
        {
            t.apply_gate(gate, bits)?;
        }
        Ok(())
    }

    fn apply_unary_gate_all<G>(&mut self, gate: &G) -> crate::error::Result<()>
    where G: crate::gates::Gate
    {
        // XXX FIXME: this can possibly be done smarter
        for bit in 0..self.nr_bits
        {
            self.apply_gate(gate, &[bit])?;
        }
        Ok(())
    }

    fn apply_conditional_gate<G>(&mut self, control: &[bool], gate: &G,
        bits: &[usize]) -> crate::error::Result<()>
    where G: crate::gates::Gate + ?Sized
    {
        if control.len() != self.nr_shots
        {
            return Err(crate::error::Error::InvalidNrControlBits(control.len(),
                self.nr_shots, String::from(gate.description())));
        }

        // XXX FIXME: computing the ranges is the same in vectorstate. Take it
        // out and make a common function.
        let mut ranges = vec![];
        let mut off = 0;
        for (icol, &count) in self.counts.iter().enumerate()
        {
            let mut begin = off;
            let mut prev = control[off];
            for ibit in off+1..off+count
            {
                if control[ibit] != prev
                {
                    ranges.push((icol, ibit-begin, prev));
                    begin = ibit;
                    prev = !prev;
                }
            }
            if begin < off+count
            {
                ranges.push((icol, off+count-begin, prev));
            }

            off += count;
        }

        let mut new_tableaus = Vec::with_capacity(ranges.len());
        for &(icol, _, apply) in ranges.iter()
        {
            let mut tableau = self.tableaus[icol].clone();
            if apply
            {
                tableau.apply_gate(gate, bits)?;
            }
            new_tableaus.push(tableau);
        }

        self.tableaus = new_tableaus;
        self.counts = ranges.iter().map(|t| t.1).collect();

        Ok(())
    }

    fn measure<R: rand::Rng>(&mut self, qbit: usize, rng: &mut R)
        -> crate::error::Result<ndarray::Array1<u64>>
    {
        let mut res = ndarray::Array1::zeros(self.nr_shots);
        self.measure_into(qbit, 0, &mut res, rng)?;
        Ok(res)
    }

    fn measure_into<R: rand::Rng>(&mut self, qbit: usize, cbit: usize,
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>
    {
        if qbit >= self.nr_bits
        {
            return Err(crate::error::Error::InvalidQBit(qbit));
        }
        if res.len() < self.nr_shots
        {
            return Err(crate::error::Error::NotEnoughSpace(res.len(), self.nr_shots));
        }

        let one_mask = 1 << cbit;
        let zero_mask = !one_mask;

        let mut new_tableaus = vec![];
        let mut new_counts = vec![];
        let mut res_start = 0;
        for (mut tableau, &count) in self.tableaus.drain(..).zip(self.counts.iter())
        {
            match tableau.measure(qbit)
            {
                MeasurementInfo::Deterministic(value) => {
                    new_tableaus.push(tableau);
                    new_counts.push(count);

                    // Store the result.
                    if value
                    {
                        res.slice_mut(s![res_start..res_start+count]).map_inplace(
                            |b| *b |= one_mask
                        );
                    }
                    else
                    {
                        res.slice_mut(s![res_start..res_start+count]).map_inplace(
                            |b| *b &= zero_mask
                        );
                    }
                },
                MeasurementInfo::Random(i) => {
                    let distribution = rand::distributions::Binomial::new(count as u64, 0.5);
                    let n0 = rng.sample(distribution) as usize;

                    if n0 == 0
                    {
                        tableau.collapse(i, qbit, true);
                        new_tableaus.push(tableau);
                        new_counts.push(count);
                    }
                    else if n0 == count
                    {
                        tableau.collapse(i, qbit, false);
                        new_tableaus.push(tableau);
                        new_counts.push(count);
                    }
                    else
                    {
                        let mut t0 = tableau.clone();
                        t0.collapse(i, qbit, false);
                        tableau.collapse(i, qbit, true);
                        new_tableaus.push(t0);
                        new_tableaus.push(tableau);
                        new_counts.push(n0);
                        new_counts.push(count - n0);
                    }

                    // Store the result.
                    res.slice_mut(s![res_start..res_start+n0]).map_inplace(
                        |b| *b &= zero_mask
                    );
                    res.slice_mut(s![res_start+n0..res_start+count]).map_inplace(
                        |b| *b |= one_mask
                    );
                }
            }

            res_start += count;
        }

        self.counts = new_counts;
        self.tableaus = new_tableaus;

        Ok(())
    }

    fn measure_all<R: rand::Rng>(&mut self, rng: &mut R)
        -> crate::error::Result<ndarray::Array1<u64>>
    {
        let mut res = ndarray::Array1::zeros(self.nr_shots);
        let cbits: Vec<usize> = (0..self.nr_bits).collect();
        self.measure_all_into(&cbits, &mut res, rng)?;
        Ok(res)
    }

    fn measure_all_into<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>
    {
        // There does not seem to be a simpler way to measure all qubits from
        // a tableau, as there is for a vector state. So simply loop over all
        // bits
        for (qbit, &cbit) in cbits.iter().enumerate()
        {
            self.measure_into(qbit, cbit, res, rng)?;
        }
        Ok(())
    }

    fn peek_into<R: rand::Rng>(&self, qbit: usize, cbit: usize,
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>
    {
        if qbit >= self.nr_bits
        {
            return Err(crate::error::Error::InvalidQBit(qbit));
        }
        if res.len() < self.nr_shots
        {
            return Err(crate::error::Error::NotEnoughSpace(res.len(), self.nr_shots));
        }

        let one_mask = 1 << cbit;
        let zero_mask = !one_mask;

        let mut res_start = 0;
        for (tableau, &count) in self.tableaus.iter().zip(self.counts.iter())
        {
            let n0 = match tableau.measure(qbit)
                {
                    MeasurementInfo::Deterministic(false) => count,
                    MeasurementInfo::Deterministic(true) => 0,
                    MeasurementInfo::Random(_) => {
                        let distribution = rand::distributions::Binomial::new(count as u64, 0.5);
                        rng.sample(distribution) as usize
                    }
                };

            if n0 > 0
            {
                res.slice_mut(s![res_start..res_start+n0]).map_inplace(
                    |b| *b &= zero_mask
                );
            }
            if n0 < count
            {
                res.slice_mut(s![res_start+n0..res_start+count]).map_inplace(
                    |b| *b |= one_mask
                );
            }

            res_start += count;
        }

        Ok(())
    }

    fn peek_all_into<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>
    {
        // There does not seem to be a simpler way to measure all qubits from
        // a tableau, as there is for a vector state. So simply loop over all
        // bits
        for (qbit, &cbit) in cbits.iter().enumerate()
        {
            self.peek_into(qbit, cbit, res, rng)?;
        }
        Ok(())
    }

    fn reset<R: rand::Rng>(&mut self, bit: usize, _rng: &mut R) -> crate::error::Result<()>
    {
        for tableau in self.tableaus.iter_mut()
        {
            tableau.reset(bit);
        }
        Ok(())
    }

    fn reset_all(&mut self)
    {
        self.tableaus = vec![StabilizerTableau::new(self.nr_bits)];
        self.counts = vec![self.nr_shots];
    }
}

#[cfg(test)]
mod tests
{
    use super::StabilizerState;
    use crate::gates::{CX, H, X};
    use crate::qustate::QuState;

    #[test]
    fn test_new()
    {
        let state = StabilizerState::new(4, 1024);
        assert_eq!(state.nr_bits, 4);
        assert_eq!(state.nr_shots, 1024);
        assert_eq!(state.counts, vec![1024]);
        assert_eq!(state.tableaus.len(), 1);
        assert_eq!(format!("{}", state.tableaus[0]), String::from(
r#"+ZIII
+IZII
+IIZI
+IIIZ"#));
    }

    #[test]
    fn test_apply_conditional_gate()
    {
        let mut s = StabilizerState::new(2, 5);
        assert_eq!(s.apply_conditional_gate(&[false, false, true, true, false], &X::new(), &[1]), Ok(()));
        assert_eq!(s.counts, vec![2, 2, 1]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+ZI
+IZ

+ZI
-IZ

+ZI
+IZ"#));

        let mut s = StabilizerState::new(2, 5);
        assert_eq!(s.apply_conditional_gate(&[false, false, true, true, true], &X::new(), &[0]), Ok(()));
        assert_eq!(s.counts, vec![2, 3]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+ZI
+IZ

-ZI
+IZ"#));

        let mut s = StabilizerState::new(2, 5);
        assert_eq!(s.apply_conditional_gate(&[true, false, true, true, false], &H::new(), &[1]), Ok(()));
        assert_eq!(s.counts, vec![1, 1, 2, 1]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+IX
+ZI

+ZI
+IZ

+IX
+ZI

+ZI
+IZ"#));

        let mut s = StabilizerState::new(2, 5);
        assert_eq!(s.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(s.apply_conditional_gate(&[true, false, true, true, false], &CX::new(), &[1, 0]), Ok(()));
        assert_eq!(s.counts, vec![1, 1, 2, 1]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+XX
+ZZ

+IX
+ZI

+XX
+ZZ

+IX
+ZI"#));


        let mut s = StabilizerState::new(2, 5);
        assert_eq!(s.apply_conditional_gate(&[true, true, true, false, false], &H::new(), &[0]), Ok(()));
        assert_eq!(s.counts, vec![3, 2]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+XI
+IZ

+ZI
+IZ"#));
        assert_eq!(s.apply_conditional_gate(&[false, false, true, true, true], &H::new(), &[0]), Ok(()));
        assert_eq!(s.counts, vec![2, 1, 2]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+XI
+IZ

+ZI
+IZ

+XI
+IZ"#));
    }

    #[test]
    fn test_measure()
    {
        let mut rng = rand::thread_rng();

        // |0⟩
        let mut s = StabilizerState::new(1, 3);
        let m = s.measure(0, &mut rng).unwrap();
        assert_eq!(m, array![0, 0, 0]);
        assert_eq!(s.counts, vec![3]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+Z"));

        // |0⟩⊗|0⟩
        let mut s = StabilizerState::new(2, 3);
        let m = s.measure(1, &mut rng).unwrap();
        assert_eq!(m, array![0, 0, 0]);
        assert_eq!(s.counts, vec![3]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+ZI\n+IZ"));
        let m = s.measure(0, &mut rng).unwrap();
        assert_eq!(m, array![0, 0, 0]);
        assert_eq!(s.counts, vec![3]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+ZI\n+IZ"));

        // (H|0⟩)⊗(H|0⟩), unnormalized
        let mut s = StabilizerState::new(2, 1024);
        assert_eq!(s.apply_unary_gate_all(&H::new()), Ok(()));
        let m0 = s.measure(0, &mut rng).unwrap();
        let mut prev_b = m0[0];
        let mut sc_idx = 0;
        for &b in m0.iter()
        {
            if b != prev_b
            {
                sc_idx += 1;
                prev_b = b;
            }
            match b
            {
                0 => assert_eq!(format!("{}", s.tableaus[sc_idx]), String::from("+IX\n+ZI")),
                1 => assert_eq!(format!("{}", s.tableaus[sc_idx]), String::from("+IX\n-ZI")),
                // LCOV_EXCL_START
                _ => panic!("Invalid value {} for bit", b)
                // LCOV_EXCL_STOP
            }
        }

        // After collapse, a new measurement should yield the same result
        let m0b = s.measure(0, &mut rng).unwrap();
        assert_eq!(m0b, m0);

        // Measure second bit
        let m1 = s.measure(1, &mut rng).unwrap();
        let mut prev_b0 = m0[0];
        let mut prev_b1 = m1[0];
        let mut sc_idx = 0;
        for j in 0..s.nr_shots
        {
            let b0 = m0[j];
            let b1 = m1[j];
            if b0 != prev_b0 || b1 != prev_b1
            {
                sc_idx += 1;
                prev_b0 = b0;
                prev_b1 = b1;
            }
            match (b0, b1)
            {
                (0, 0) => assert_eq!(format!("{}", s.tableaus[sc_idx]), String::from("+ZI\n+IZ")),
                (0, 1) => assert_eq!(format!("{}", s.tableaus[sc_idx]), String::from("+ZI\n-IZ")),
                (1, 0) => assert_eq!(format!("{}", s.tableaus[sc_idx]), String::from("-ZI\n+IZ")),
                (1, 1) => assert_eq!(format!("{}", s.tableaus[sc_idx]), String::from("-ZI\n-IZ")),
                // LCOV_EXCL_START
                _      => panic!("Invalid value {:?} for bits", (b0, b1))
                // LCOV_EXCL_STOP
            }
        }
    }

    #[test]
    fn test_peek_into()
    {
        let nr_shots = 1024;
        let mut measurements = ndarray::Array1::zeros(nr_shots);

        let mut rng = rand::thread_rng();

        // |0⟩
        let s = StabilizerState::new(1, nr_shots);
        assert_eq!(s.peek_into(0, 0, &mut measurements, &mut rng), Ok(()));
        assert!(measurements.iter().all(|&bits| bits == 0));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+Z"));

        // H|0⟩
        let mut s = StabilizerState::new(1, nr_shots);
        assert_eq!(s.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(s.peek_into(0, 0, &mut measurements, &mut rng), Ok(()));
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+X"));

        // H|0⟩⊗ H|0⟩
        let mut s = StabilizerState::new(2, nr_shots);
        assert_eq!(s.apply_unary_gate_all(&H::new()), Ok(()));
        assert_eq!(s.peek_into(0, 0, &mut measurements, &mut rng), Ok(()));
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        measurements.fill(0);
        assert_eq!(s.peek_into(1, 0, &mut measurements, &mut rng), Ok(()));
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+XI\n+IX"));

        // H|0⟩⊗ |1⟩
        let mut s = StabilizerState::new(2, nr_shots);
        assert_eq!(s.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(s.apply_gate(&X::new(), &[1]), Ok(()));
        assert_eq!(s.peek_into(0, 0, &mut measurements, &mut rng), Ok(()));
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        measurements.fill(0);
        assert_eq!(s.peek_into(1, 0, &mut measurements, &mut rng), Ok(()));
        assert_eq!(measurements.sum() as usize, nr_shots);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+XI\n-IZ"));
    }

    #[test]
    fn test_measure_all()
    {
        let nr_bits = 3;
        let nr_shots = 5;

        let mut rng = rand::thread_rng();

        let mut s = StabilizerState::new(nr_bits, nr_shots);
        assert_eq!(s.apply_unary_gate_all(&X::new()), Ok(()));
        let result = s.measure_all(&mut rng).unwrap();
        assert_eq!(result.shape(), [nr_shots]);
        assert!(result.iter().all(|&b| b == 0b111));

        let mut s = StabilizerState::new(nr_bits, nr_shots);
        assert_eq!(s.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(s.apply_gate(&X::new(), &[1]), Ok(()));
        let result = s.measure_all(&mut rng).unwrap();
        assert_eq!(result.shape(), [nr_shots]);
        assert!(result.iter().all(|&b| b == 0b011));

        let mut s = StabilizerState::new(nr_bits, nr_shots);
        assert_eq!(s.apply_gate(&H::new(), &[2]), Ok(()));
        let result = s.measure_all(&mut rng).unwrap();
        assert_eq!(result.shape(), [nr_shots]);
        assert!(result.iter().all(|&b| (b & 0b011) == 0));
    }

    #[test]
    fn test_peek_all()
    {
        let tol = 1.0e-5;
        let nr_bits = 3;
        let nr_shots = 1024;

        let mut rng = rand::thread_rng();

        let mut s = StabilizerState::new(nr_bits, nr_shots);
        assert_eq!(s.apply_unary_gate_all(&X::new()), Ok(()));
        assert_eq!(s.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(s.apply_gate(&H::new(), &[2]), Ok(()));
        let mut result = ndarray::Array1::zeros(nr_shots);
        assert_eq!(s.peek_all_into(&[0, 1, 2], &mut result, &mut rng), Ok(()));
        // Ensure quantum state is preserved
        assert_eq!(s.counts, vec![nr_shots]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("-XII\n-IIX\n-IZI"));
        // Ensure measurement is correct
        assert_eq!(result.shape(), [nr_shots]);
        let mut n = vec![0; nr_bits];
        for &bits in result.iter()
        {
            for i in 0..nr_bits
            {
                n[i] += (bits as usize >> i) & 1;
            }
        }
        assert!(crate::stats::measurement_ok(n[0], nr_shots, 0.5, tol));
        assert_eq!(n[1], nr_shots);
        assert!(crate::stats::measurement_ok(n[0], nr_shots, 0.5, tol));
    }

    #[test]
    fn test_reset()
    {
        let nr_runs = 10;

        let mut rng = rand::thread_rng();

        let mut s = StabilizerState::new(1, nr_runs);
        assert_eq!(s.reset(0, &mut rng), Ok(()));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+Z"));

        let mut s = StabilizerState::new(1, nr_runs);
        assert_eq!(s.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(s.reset(0, &mut rng), Ok(()));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+Z"));

        let mut s = StabilizerState::new(2, nr_runs);
        assert_eq!(s.apply_unary_gate_all(&X::new()), Ok(()));
        assert_eq!(s.reset(0, &mut rng), Ok(()));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+ZI\n-IZ"));

        let mut s = StabilizerState::new(2, nr_runs);
        assert_eq!(s.apply_unary_gate_all(&X::new()), Ok(()));
        assert_eq!(s.reset(1, &mut rng), Ok(()));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("-ZI\n+IZ"));

        let mut s = StabilizerState::new(2, nr_runs);
        assert_eq!(s.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(s.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(s.reset(0, &mut rng), Ok(()));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("+ZI\n+IZ"));

        let mut s = StabilizerState::new(2, nr_runs);
        assert_eq!(s.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(s.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(s.reset(1, &mut rng), Ok(()));
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from("-XI\n+IZ"));
    }

    #[test]
    fn test_reset_all()
    {
        let nr_bits = 5;
        let nr_runs = 100;

        let mut s = StabilizerState::new(nr_bits, nr_runs);
        assert_eq!(s.apply_gate(&H::new(), &[2]), Ok(()));
        assert_eq!(s.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(s.apply_gate(&H::new(), &[4]), Ok(()));

        s.reset_all();
        assert_eq!(s.counts, vec![nr_runs]);
        let ss = s.tableaus.iter().map(|t| format!("{}", t)).collect::<Vec<String>>().join("\n\n");
        assert_eq!(ss, String::from(
r#"+ZIIII
+IZIII
+IIZII
+IIIZI
+IIIIZ"#));
    }

}
