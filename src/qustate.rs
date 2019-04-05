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

use rand::distributions::Distribution;

/// Quantum state.
///
/// Struct Qustate represents the quantum experiment. It consists of a series of
/// quantum states, combined with the number of times this state occurs in the
/// experiment. Each quantum state is a (normalized) superposition of basis states,
/// ∑<sub>i</sub>a<sub>i</sub>|i〉, where each basis function |i〉 is a Kronecker
/// product of quantum bits.
#[derive(Debug)]
pub struct QuState
{
    /// The number of qubits in this state
    nr_bits: usize,
    /// The number of separate runs for evolving this state
    nr_shots: usize,
    /// Run counts for each quantum state
    counts: Vec<usize>,
    /// The quantum states themselves
    states: crate::cmatrix::CMatrix,
}

impl QuState
{
    /// Create a new qustate of `nr_bits` qubits, all initialized to |0〉, which
    /// will be measured `nr_shots` times.
    pub fn new(nr_bits: usize, nr_shots: usize) -> Self
    {
        let mut states = crate::cmatrix::CMatrix::zeros((1 << nr_bits, 1));
        states[(0, 0)] = crate::cmatrix::COMPLEX_ONE;

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            counts: vec![nr_shots],
            states: states
        }
    }

    /// Create a new qustate from qubit coefficients.
    ///
    /// Create a new qustate as a direct product of qubits, where the
    /// coefficients of the |0〉 and |1〉 states in the qubits are given in
    /// `bit_coefs`. This array must be of size `2*n`, where `n` is the number
    /// of qubits in the system. The state will be evaluated in `nr_shots`
    /// separate runs.
    pub fn from_qubit_coefs(bit_coefs: &[num_complex::Complex64], nr_shots: usize) -> Self
    {
        assert!(bit_coefs.len() % 2 == 0, "Length of coefficient array is not even");

        let nr_bits = bit_coefs.len() / 2;

        let mut states = crate::cmatrix::CMatrix::ones((1, 1));
        for c in bit_coefs.chunks(2)
        {
            let norm = (c[0].norm_sqr() + c[1].norm_sqr()).sqrt();
            let bit = array![[c[0]/norm], [c[1]/norm]];
            states = crate::cmatrix::kron_mat(&states, &bit);
        }

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            counts: vec![nr_shots],
            states: states
        }
    }

    /// Return the number of qubits in this state
    pub fn nr_bits(&self) -> usize
    {
        self.nr_bits
    }

    /// Apply a n-ary quantum gate `gate` on the qubits from `bits` in this state.
    pub fn apply_gate<G>(&mut self, gate: &G, bits: &[usize])
    where G: crate::gates::Gate + ?Sized
    {
        let gate_bits = gate.nr_affected_bits();
        assert!(gate_bits == bits.len(),
            "The number of bits affected by the {} gate should be {}, but {} bits were provided.",
            gate.description(), gate_bits, bits.len()
        );

        crate::gates::apply_gate_mat_slice(self.states.view_mut(), gate, bits, self.nr_bits);
    }

    pub fn apply_unary_gate_all<G>(&mut self, gate: &G)
    where G: crate::gates::Gate
    {
        // XXX FIXME: this can be done smarter
        for bit in 0..self.nr_bits
        {
            self.apply_gate(gate, &[bit]);
        }
    }

    /// Apply a conditional n-ary quantum gate `gate`, controlled by classical
    /// bit `control`, on the qubits from `bits` in this state.
    pub fn apply_conditional_gate<G>(&mut self, control: &[bool], gate: &G,
        bits: &[usize])
    where G: crate::gates::Gate + ?Sized
    {
        assert!(control.len() == self.nr_shots
            "The number of control bits does not match the number of runs");

        let gate_bits = gate.nr_affected_bits();
        assert!(gate_bits == bits.len(),
            "The number of bits affected by the {} gate should be {}, but {} bits were provided.",
            gate.description(), gate_bits, bits.len()
        );

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

        let mut new_states = crate::cmatrix::CMatrix::zeros((1 << self.nr_bits, ranges.len()));
        for (new_icol, &(icol, _, apply)) in ranges.iter().enumerate()
        {
            new_states.column_mut(new_icol).assign(&self.states.column(icol));
            if apply
            {
                crate::gates::apply_gate_slice(new_states.column_mut(new_icol),
                    gate, bits, self.nr_bits);
            }
        }

        self.states = new_states;
        self.counts = ranges.iter().map(|t| t.1).collect();
    }

    fn collapse(mut coefs: crate::cmatrix::CVecSliceMut, block_size: usize, nr_blocks: usize,
        offset: usize, norm_sq: f64)
    {
        // Set coefficients for other bit to zero
        let mut off = offset;
        for _ in 0..nr_blocks
        {
            coefs.slice_mut(s![off..off+block_size]).fill(crate::cmatrix::COMPLEX_ZERO);
            off += 2 * block_size;
        }

        // Renormalize
        coefs *= num_complex::Complex::new(1.0 / norm_sq.sqrt(), 0.0);
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `qbit` in the state. Measurement is done
    /// in the `z`-basis. The random number generator `rng` is used for sampling.
    /// The result is returned as an array containing the measurement result for
    /// each run.
    pub fn measure<R: rand::Rng>(&mut self, qbit: usize, rng: &mut R) -> ndarray::Array1<u64>
    {
        let mut res = ndarray::Array1::zeros(self.nr_shots);
        self.measure_into(qbit, 0, &mut res, rng);
        res
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `qbit` in the state to classical bit
    /// `cbit` in `res`, which should be an array of sufficient length to store
    /// results for the total number of runs in the state. Measurement is done
    /// in the `z`-basis. The random number generator `rng` is used for sampling.
    pub fn measure_into<R: rand::Rng>(&mut self, qbit: usize, cbit: usize,
        res: &mut ndarray::Array1<u64>, rng: &mut R)
    {
        assert!(qbit < self.nr_bits, "Invalid bit index");
        assert!(res.len() >= self.nr_shots, "Not enough space to store the results");

        let block_size = 1 << (self.nr_bits - qbit - 1);
        let nr_blocks = 1 << qbit;

        // Compute chances of measuring 0
        let mut w0s = ndarray::Array1::zeros(self.states.cols());
        let mut off = 0;
        for _ in 0..nr_blocks
        {
            w0s += &self.states.slice(s![off..off+block_size, ..])
                .mapv(|c| c.norm_sqr())
                .sum_axis(ndarray::Axis(0));
            off += 2 * block_size;
        }

        // Compute how many times we measure 0
        let mut n0s = vec![];
        let mut new_nr_states = 0;
        for (&w0, &c) in w0s.iter().zip(self.counts.iter())
        {
            // Sometimes, the sum may add up to slightly more than 1, due to
            // numerical inaccuracies. This causes the Binomial distribution to
            // panic, so cap w0.
            let distribution = rand::distributions::Binomial::new(c as u64, w0.min(1.0));
            let n0 = rng.sample(distribution) as usize;
            n0s.push(n0);
            new_nr_states += if n0 == 0 || n0 == c { 1 } else { 2 };
        }

        let one_mask = 1 << cbit;
        let zero_mask = !one_mask;

        let mut new_states = crate::cmatrix::CMatrix::zeros((1 << self.nr_bits, new_nr_states));
        let mut new_counts = vec![];
        let mut new_idx = 0;
        let mut res_start = 0;
        for idx in 0..self.states.cols()
        {
            let w0 = w0s[[idx]];
            let n0 = n0s[idx];
            let count = self.counts[idx];

            // Store the result.
            res.slice_mut(s![res_start..res_start+n0]).map_inplace(
                |b| *b &= zero_mask
            );
            res.slice_mut(s![res_start+n0..res_start+count]).map_inplace(
                |b| *b |= one_mask
            );
            res_start += count;

            // Collapse the wave function
            new_states.column_mut(new_idx).assign(&self.states.column(idx));
            if n0 == count
            {
                Self::collapse(new_states.column_mut(new_idx), block_size, nr_blocks, block_size, w0);
                new_counts.push(count);
            }
            else if n0 == 0
            {
                Self::collapse(new_states.column_mut(new_idx), block_size, nr_blocks, 0, 1.0 - w0);
                new_counts.push(count);
            }
            else
            {
                Self::collapse(new_states.column_mut(new_idx), block_size, nr_blocks, block_size, w0);
                new_counts.push(n0);
                new_idx += 1;

                new_states.column_mut(new_idx).assign(&self.states.column(idx));
                Self::collapse(new_states.column_mut(new_idx), block_size, nr_blocks, 0, 1.0 - w0);
                new_counts.push(count - n0);
            }

            new_idx += 1;
        }

        self.counts = new_counts;
        self.states = new_states;
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `qbit` in the state to bit `cbit` in res,
    /// without affecting the quantum state, i.e. the wave function is not
    /// collapsed. The output array `res` should be of sufficient length to store
    /// results for the total number of runs in the state. Measurement is done
    /// in the `z`-basis. The random number generator `rng` is used for sampling.
    /// NOTE: this is not a physical process, and impossible to reproduce on
    /// a real quantum computer.
    pub fn peek_into<R: rand::Rng>(&mut self, qbit: usize, cbit: usize,
        res: &mut ndarray::Array1<u64>, rng: &mut R)
    {
        assert!(qbit < self.nr_bits, "Invalid bit index");
        assert!(res.len() >= self.nr_shots, "Not enough space to store the results");

        let block_size = 1 << (self.nr_bits - qbit - 1);
        let nr_blocks = 1 << qbit;

        // Compute chances of measuring 0
        let mut w0s = ndarray::Array1::zeros(self.states.cols());
        let mut off = 0;
        for _ in 0..nr_blocks
        {
            w0s += &self.states.slice(s![off..off+block_size, ..])
                .mapv(|c| c.norm_sqr())
                .sum_axis(ndarray::Axis(0));
            off += 2 * block_size;
        }

        let one_mask = 1 << cbit;
        let zero_mask = !one_mask;

        let mut res_start = 0;
        for (&w0, &c) in w0s.iter().zip(self.counts.iter())
        {
            // Compute how many times we measure 0
            let distribution = rand::distributions::Binomial::new(c as u64, w0.min(1.0));
            let n0 = rng.sample(distribution) as usize;

            // Store the result.
            res.slice_mut(s![res_start..res_start+n0]).map_inplace(
                |b| *b &= zero_mask
            );
            res.slice_mut(s![res_start+n0..res_start+c]).map_inplace(
                |b| *b |= one_mask
            );
            res_start += c;
        }
    }

    /// Measure all qubits
    ///
    /// Measure all qubits in this state, and return the results. The random
    /// number generator `rng` is used for sampling.
    pub fn measure_all<R: rand::Rng>(&mut self, rng: &mut R) -> ndarray::Array1<u64>
    {
        let mut res = ndarray::Array1::zeros(self.nr_shots);
        let cbits: Vec<usize> = (0..self.nr_bits).collect();
        self.measure_all_into(&cbits, &mut res, rng);
        res
    }

    fn measure_all_into_helper<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, collapse: bool, rng: &mut R)
    {
        assert!(res.len() >= self.nr_shots, "Not enough space to store the results");
        assert!(cbits.len() == self.nr_bits,
            "Number of measurement bits does not match number of quantum bits");

        let mut state_counts = vec![];
        for col_idx in 0..self.states.cols()
        {
            let mut count_map = ::std::collections::HashMap::new();
            let distr = rand::distributions::WeightedIndex::new(
                self.states.column(col_idx).iter().map(|c| c.norm_sqr())
            ).unwrap();
            for idx in distr.sample_iter(rng).take(self.counts[col_idx])
            {
                let entry = count_map.entry(idx).or_insert(0);
                *entry += 1;
            }

            state_counts.extend(count_map.into_iter());
        }

        let mask = !cbits.iter().fold(0u64, |m, b| m | (1u64 << b));
        let mut res_off = 0;
        for &(idx, count) in state_counts.iter()
        {
            // For each unique measurement, store n copies of it in the result,
            // where n is the number of times it was measured.
            let rev_idx = crate::support::reverse_bits(idx as u64, self.nr_bits);
            let perm_idx = crate::support::shuffle_bits(rev_idx, cbits);
            res.slice_mut(s![res_off..res_off+count]).map_inplace(
                |bits| *bits = (*bits & mask) | perm_idx
            );

            res_off += count;
        }

        if collapse
        {
            self.states = crate::cmatrix::CMatrix::zeros((1 << self.nr_bits, state_counts.len()));
            for (col_idx, &(idx, _)) in state_counts.iter().enumerate()
            {
                self.states[(idx, col_idx)] = crate::cmatrix::COMPLEX_ONE;
            }
            self.counts = state_counts.iter().map(|t| t.1).collect();

        }
    }
    /// Measure all qubits
    ///
    /// Measure all qubits in this state, and store the results in `res`,
    /// which should be of sufficient length to hold results for the number of
    /// runs in this state.  The first qubit measured is stored at the bit
    /// position indicated by the first element of `cbits`, and so on. The random
    /// number generator `rng` is used for sampling.
    pub fn measure_all_into<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, rng: &mut R)
    {
        self.measure_all_into_helper(cbits, res, true, rng);
    }

    /// Measure all qubits
    ///
    /// Measure all qubits in this state without affecting the quantum state,
    /// i.e. without collapsing the wave function. The measurement results
    /// are stored in `res`, which must be of sufficient length to hold results
    /// for the total number of runs in the state. The first qubit measured
    /// is stored at the bit position indicated by the first element of `cbits`,
    /// and so on. The random number generator `rng` is used for sampling.
    /// NOTE: this is not a physical process, and impossible to reproduce on
    /// a real quantum computer.
    pub fn peek_all_into<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, rng: &mut R)
    {
        self.measure_all_into_helper(cbits, res, false, rng);
    }

    /// Reset a qubit
    ///
    /// Reset the qubit with index `bit` to zero. This is done by measuring the
    /// bit, and rotating it back to zero if the result is 1. The random
    /// number generator `rng` is used for sampling in the measurement.
    pub fn reset<R: rand::Rng>(&mut self, bit: usize, rng: &mut R)
    {
        let measurement = self.measure(bit, rng);
        let control: Vec<bool> = measurement.iter().map(|&b| b != 0).collect();
        self.apply_conditional_gate(&control, &crate::gates::X::new(), &[bit]);
    }

    /// Reset all qubits
    ///
    /// Reset all qubits in this experiment, returning the state to |00...0〉
    /// for all runs.
    pub fn reset_all(&mut self)
    {
        self.states = crate::cmatrix::CMatrix::zeros((1 << self.nr_bits, 1));
        self.states[[0, 0]] = crate::cmatrix::COMPLEX_ONE;
        self.counts = vec![self.nr_shots];
    }
}

#[cfg(test)]
mod tests
{
    use super::QuState;
    use crate::gates::{CCX, CX, H, Kron, X, Y};

    #[test]
    fn test_new()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;

        let s = QuState::new(1, 1024);
        assert_eq!(s.nr_bits, 1);
        assert_eq!(s.nr_shots, 1024);
        assert_eq!(s.counts, vec![1024]);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z]]);

        let s = QuState::new(3, 1500);
        assert_eq!(s.nr_bits, 3);
        assert_eq!(s.nr_shots, 1500);
        assert_eq!(s.counts, vec![1500]);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z], [z], [z], [z], [z], [z], [z]]);
    }

    #[test]
    fn test_from_qubit_coefs()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let i = crate::cmatrix::COMPLEX_I;

        // |0〉⊗|1〉
        let s = QuState::from_qubit_coefs(&[o, z, z, o], 1);
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 1);
        assert_eq!(s.counts, vec![1]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [o], [z], [z]]);
        // |1〉⊗|0〉
        let s = QuState::from_qubit_coefs(&[z, o, o, z], 13);
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 13);
        assert_eq!(s.counts, vec![13]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [o], [z]]);
        // (H|0〉)⊗(Y|1〉), unnormalized
        let s = QuState::from_qubit_coefs(&[o, o, -i, z], 9);
        let x = ::std::f64::consts::FRAC_1_SQRT_2 * i;
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 9);
        assert_eq!(s.counts, vec![9]);
        assert_complex_matrix_eq!(&s.states, &array![[-x], [z], [-x], [z]]);
    }

    #[test]
    fn test_apply_conditional_gate()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[false, false, true, true, false], &X::new(), &[1]);
        assert_eq!(s.counts, vec![2, 2, 1]);
        assert_complex_matrix_eq!(&s.states,
            &array![[o, z, o], [z, o, z], [z, z, z], [z, z, z]]);

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[false, false, true, true, true], &X::new(), &[0]);
        assert_eq!(s.counts, vec![2, 3]);
        assert_complex_matrix_eq!(&s.states, &array![[o, z], [z, z], [z, o], [z, z]]);

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[true, false, true, true, false], &H::new(), &[1]);
        assert_eq!(s.counts, vec![1, 1, 2, 1]);
        assert_complex_matrix_eq!(&s.states,
            &array![[x, o, x, o], [x, z, x, z], [z, z, z, z], [z, z, z, z]]);

        let mut s = QuState::from_qubit_coefs(&[o, z, x, x], 5);
        s.apply_conditional_gate(&[true, false, true, true, false], &CX::new(), &[1, 0]);
        assert_eq!(s.counts, vec![1, 1, 2, 1]);
        assert_complex_matrix_eq!(&s.states,
            &array![[x, x, x, x], [z, x, z, x], [z, z, z, z], [x, z, x, z]]);

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[true, true, true, false, false], &H::new(), &[0]);
        assert_eq!(s.counts, vec![3, 2]);
        assert_complex_matrix_eq!(&s.states, &array![[x, o], [z, z], [x, z], [z, z]]);
        s.apply_conditional_gate(&[false, false, true, true, true], &H::new(), &[0]);
        assert_eq!(s.counts, vec![2, 1, 2]);
        assert_complex_matrix_eq!(&s.states, &array![[x, o, x], [z, z, z], [x, z, x], [z, z, z]]);
    }

    #[test]
    fn test_measure()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;

        let mut rng = rand::thread_rng();

        // |0〉
        let mut s = QuState::new(1, 3);
        let m = s.measure(0, &mut rng);
        assert_eq!(m, array![0, 0, 0]);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z]]);

        // |0〉⊗|0〉
        let mut s = QuState::from_qubit_coefs(&[o, z, o, z], 3);
        let m = s.measure(1, &mut rng);
        assert_eq!(m, ndarray::Array1::zeros(3));
        assert_complex_matrix_eq!(&s.states, &array![[o], [z], [z], [z]]);
        let m = s.measure(0, &mut rng);
        assert_eq!(m, ndarray::Array1::zeros(3));
        assert_complex_matrix_eq!(&s.states, &array![[o], [z], [z], [z]]);

        // (H|0〉)⊗(H|0〉), unnormalized
        let mut s = QuState::from_qubit_coefs(&[o, o, o, o], 1024);
        let m0 = s.measure(0, &mut rng);
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
                0 => assert_complex_vector_eq!(&s.states.column(sc_idx),
                    &array![x, x, z, z]),
                1 => assert_complex_vector_eq!(&s.states.column(sc_idx),
                    &array![z, z, x, x]),
                // LCOV_EXCL_START
                _ => panic!("Invalid value {} for bit", b)
                // LCOV_EXCL_STOP
            }
        }

        // After collapse, a new measurement should yield the same result
        let m0b = s.measure(0, &mut rng);
        assert_eq!(m0b, m0);

        // Measure second bit
        let m1 = s.measure(1, &mut rng);
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
            let coefs = &s.states.column(sc_idx);
            match (b0, b1)
            {
                (0, 0) => assert_complex_vector_eq!(coefs, &array![o, z, z, z]),
                (0, 1) => assert_complex_vector_eq!(coefs, &array![z, o, z, z]),
                (1, 0) => assert_complex_vector_eq!(coefs, &array![z, z, o, z]),
                (1, 1) => assert_complex_vector_eq!(coefs, &array![z, z, z, o]),
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

        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * o;

        let mut rng = rand::thread_rng();

        // |0〉
        let mut s = QuState::new(1, nr_shots);
        s.peek_into(0, 0, &mut measurements, &mut rng);
        assert!(measurements.iter().all(|&bits| bits == 0));
        assert_complex_matrix_eq!(&s.states, &array![[o], [z]]);

        // H|0〉
        let mut s = QuState::from_qubit_coefs(&[o, o], nr_shots);
        s.peek_into(0, 0, &mut measurements, &mut rng);
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        assert_complex_matrix_eq!(&s.states, &array![[x], [x]]);

        // H|0〉⊗ H|0〉
        let mut s = QuState::from_qubit_coefs(&[o, o, o, o], nr_shots);
        s.peek_into(0, 0, &mut measurements, &mut rng);
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        measurements.fill(0);
        s.peek_into(1, 0, &mut measurements, &mut rng);
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        assert_complex_matrix_eq!(&s.states, &array![[h], [h], [h], [h]]);

        // H|0〉⊗ |1〉
        let mut s = QuState::from_qubit_coefs(&[x, x, z, o], nr_shots);
        s.peek_into(0, 0, &mut measurements, &mut rng);
        assert!(crate::stats::measurement_ok(measurements.sum() as usize, nr_shots,
            0.5, 1.0e-5));
        measurements.fill(0);
        s.peek_into(1, 0, &mut measurements, &mut rng);
        assert_eq!(measurements.sum() as usize, nr_shots);
        assert_complex_matrix_eq!(&s.states, &array![[z], [x], [z], [x]]);
    }

    #[test]
    fn test_apply_unary_gate()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let i = crate::cmatrix::COMPLEX_I;

        let mut s = QuState::new(3, 1);
        s.apply_gate(&H::new(), &[0]);
        assert_complex_matrix_eq!(&s.states, &array![[x], [z], [z], [z], [x], [z], [z], [z]]);

        let mut s = QuState::new(3, 1);
        s.apply_gate(&H::new(), &[1]);
        assert_complex_matrix_eq!(&s.states, &array![[x], [z], [x], [z], [z], [z], [z], [z]]);

        let mut s = QuState::new(3, 1);
        s.apply_gate(&Y::new(), &[2]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [i], [z], [z], [z], [z], [z], [z]]);
    }

    #[test]
    fn test_apply_binary_gate()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;

        let mut s = QuState::new(3, 1);
        s.apply_gate(&CX::new(), &[0, 1]);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z], [z], [z], [z], [z], [z], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o, o, z, o, z], 1);
        s.apply_gate(&CX::new(), &[0, 1]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [z], [z], [z], [z], [o], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o, o, z, o, z], 1);
        s.apply_gate(&CX::new(), &[0, 2]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [z], [z], [z], [o], [z], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o, o, z, o, z], 1);
        let hh = Kron::new(H::new(), H::new());
        s.apply_gate(&hh, &[1, 2]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [z], [z], [h], [h], [h], [h]]);
    }

    #[test]
    fn test_apply_n_ary_gate()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;
        let hx = 0.5 * x;

        let mut s = QuState::new(3, 1);
        s.apply_gate(&CCX::new(), &[0, 1, 2]);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z], [z], [z], [z], [z], [z], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o, o, z], 1);
        s.apply_gate(&CCX::new(), &[0, 2, 1]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [z], [z], [z], [z], [o], [z]]);
        s.apply_gate(&CCX::new(), &[0, 1, 2]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [z], [z], [z], [z], [z], [o]]);

        let mut s = QuState::from_qubit_coefs(&[x, -x, x, -x, x, -x], 1);
        s.apply_gate(&CCX::new(), &[0, 2, 1]);
        assert_complex_matrix_eq!(&s.states,
            &array![[hx], [-hx], [-hx], [hx], [-hx], [-hx], [hx], [hx]]);
    }

    #[test]
    fn test_measure_all()
    {
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;

        let nr_bits = 3;
        let nr_shots = 5;

        let mut rng = rand::thread_rng();

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o, z, o], nr_shots);
        let result = s.measure_all(&mut rng);
        assert_eq!(result.shape(), [nr_shots]);
        assert!(result.iter().all(|&b| b == 0b111));

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o, o, z], nr_shots);
        let result = s.measure_all(&mut rng);
        assert_eq!(result.shape(), [nr_shots]);
        assert!(result.iter().all(|&b| b == 0b011));

        let mut s = QuState::new(nr_bits, nr_shots);
        s.apply_gate(&H::new(), &[2]);
        let result = s.measure_all(&mut rng);
        assert_eq!(result.shape(), [nr_shots]);
        assert!(result.iter().all(|&b| (b & 0b011) == 0));
    }

    #[test]
    fn test_peek_all()
    {
        let tol = 1.0e-5;
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;

        let nr_bits = 3;
        let nr_shots = 1024;

        let mut rng = rand::thread_rng();

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o, z, o], nr_shots);
        s.apply_gate(&H::new(), &[0]);
        s.apply_gate(&H::new(), &[2]);
        let mut result = ndarray::Array1::zeros(nr_shots);
        s.peek_all_into(&[0, 1, 2], &mut result, &mut rng);
        // Ensure quantum state is preserved
        assert_eq!(s.counts, vec![nr_shots]);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [h], [-h], [z], [z], [-h], [h]]);
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
        let z = crate::cmatrix::COMPLEX_ZERO;
        let o = crate::cmatrix::COMPLEX_ONE;
        let x = crate::cmatrix::COMPLEX_HSQRT2;

        let mut rng = rand::thread_rng();

        let mut s = QuState::from_qubit_coefs(&[o, z], nr_runs);
        s.reset(0, &mut rng);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o], nr_runs);
        s.reset(0, &mut rng);
        assert_complex_matrix_eq!(&s.states, &array![[o], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o], nr_runs);
        s.reset(0, &mut rng);
        assert_complex_matrix_eq!(&s.states, &array![[z], [o], [z], [z]]);

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o], nr_runs);
        s.reset(1, &mut rng);
        assert_complex_matrix_eq!(&s.states, &array![[z], [z], [o], [z]]);

        let mut s = QuState::from_qubit_coefs(&[x, -x, o, z], nr_runs);
        s.reset(0, &mut rng);
        match s.counts.len()
        {
            1 => { assert_complex_matrix_eq!(&s.states, &array![[o], [z], [z], [z]]); },
            2 => { assert_complex_matrix_eq!(&s.states, &array![[o, -o], [z, z], [z, z], [z, z]]); },
            _ => { assert!(false, "unexpected number of states ({})", s.counts.len()); }
        }

        let mut s = QuState::from_qubit_coefs(&[x, -x, o, z], nr_runs);
        s.reset(1, &mut rng);
        assert_complex_matrix_eq!(&s.states, &array![[x], [z], [-x], [z]]);
    }

    #[test]
    fn test_reset_all()
    {
        let nr_bits = 5;
        let nr_runs = 100;

        let mut s = QuState::new(nr_bits, nr_runs);
        s.apply_gate(&H::new(), &[2]);
        s.apply_gate(&X::new(), &[0]);
        s.apply_gate(&H::new(), &[4]);

        s.reset_all();
        assert_eq!(s.counts, vec![nr_runs]);
        let mut coefs = crate::cmatrix::CMatrix::zeros((1 << nr_bits, 1));
        coefs[[0, 0]] = crate::cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(&s.states, &coefs);
    }
}
