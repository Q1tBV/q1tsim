extern crate ndarray;
extern crate num_complex;
extern crate rand;

use cmatrix;
use gates;

use self::rand::distributions::Distribution;
use self::rand::Rng;

/// Single quantum state
///
/// The full state of the exoeriment is stored as a collection of quantum states,
/// combined with the number of times that this exact state occurs. Initially,
/// the experiment consists of one single state, which is the same for all
/// runs. Every time a measurement is made in the experiment, each state can split
/// into two or more different states, each with their own count.
#[derive(Clone, Debug, PartialEq)]
pub struct StateCount
{
    coefs: cmatrix::CVector,
    count: usize
}

impl StateCount
{
    fn new(coefs: cmatrix::CVector, count: usize) -> Self
    {
        StateCount { coefs: coefs, count: count }
    }
}

/// Quantum state.
///
/// Struct Qustate represents the quantum experiment. It consists of a series of
/// quantum states, combined with the number of times this state occurs in the
/// experiment. Each quantm state is a (normalized) superposition of basis states,
/// ∑<sub>i</sub>a<sub>i</sub>|i〉, where each basis function |i〉 is a Kronecker
/// product of quantum bits.
#[derive(Debug)]
pub struct QuState
{
    /// The number of qubits in this state
    nr_bits: usize,
    /// The number of separate runs for evolving this state
    nr_shots: usize,
    /// Coefficients of the basis states in this state
    state_counts: ::std::collections::LinkedList<StateCount>,
    /// Random number generator, for measurements
    rng: rand::rngs::ThreadRng
}

impl QuState
{
    /// Create a new qustate of `nr_bits` qubits, all initialized to |0〉, which
    /// will be measured `nr_shots` times.
    pub fn new(nr_bits: usize, nr_shots: usize) -> Self
    {
        let mut coefs = cmatrix::CVector::zeros(1 << nr_bits);
        coefs[0] = cmatrix::COMPLEX_ONE;

        let mut state_counts = ::std::collections::LinkedList::new();
        state_counts.push_back(StateCount::new(coefs, nr_shots));

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            state_counts: state_counts,
            rng: rand::thread_rng()
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

        let mut coefs = cmatrix::CVector::ones(1);
        for c in bit_coefs.chunks(2)
        {
            let norm = (c[0].norm_sqr() + c[1].norm_sqr()).sqrt();
            let bit = array![c[0]/norm, c[1]/norm];
            coefs = cmatrix::kron_vec(&coefs, &bit);
        }

        let mut state_counts = ::std::collections::LinkedList::new();
        state_counts.push_back(StateCount::new(coefs, nr_shots));

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            state_counts: state_counts,
            rng: rand::thread_rng()
        }
    }

    /// Return the number of qubits in this state
    pub fn nr_bits(&self) -> usize
    {
        self.nr_bits
    }

    /// Apply a n-ary quantum gate `gate` on the qubits from `bits` in this state.
    pub fn apply_gate<G>(&mut self, gate: &G, bits: &[usize])
    where G: gates::Gate + ?Sized
    {
        let gate_bits = gate.nr_affected_bits();
        assert!(gate_bits == bits.len(),
            "The number of bits affected by the gate does not match the provided number of bits.");

        if gate_bits == 1
        {
            let block_size = 1 << (self.nr_bits - bits[0]);
            let nr_blocks = 1 << bits[0];
            for m in self.state_counts.iter_mut()
            {
                for i in 0..nr_blocks
                {
                    gate.apply_slice(&mut m.coefs.slice_mut(s![i*block_size..(i+1)*block_size]));
                }
            }
        }
        else
        {
            let perm = gates::bit_permutation(self.nr_bits, bits);
            let inv_perm = perm.inverse();
            for m in self.state_counts.iter_mut()
            {
                inv_perm.apply_vec(&mut m.coefs);
                gate.apply(&mut m.coefs);
                perm.apply_vec(&mut m.coefs);
            }
        }
    }

    /// Apply a conditional n-ary quantum gate `gate`, controlled by classical
    /// bit `control`, on the qubits from `bits` in this state.
    pub fn apply_conditional_gate<G>(&mut self, control: &[bool], gate: &G,
        bits: &[usize])
    where G: gates::Gate + ?Sized
    {
        assert!(control.len() == self.nr_shots
            "The number of control bits does not match the number of runs");

        let gate_bits = gate.nr_affected_bits();
        assert!(gate_bits == bits.len(),
            "The number of bits affected by the gate does not match the provided number of bits.");

        // Create a list of ranges of runs that should have the gate applied
        let mut ranges = vec![];
        let mut begin = self.nr_shots;
        for (end, &cbit) in control.iter().enumerate()
        {
            if !cbit
            {
                if begin < self.nr_shots
                {
                    ranges.push((begin, end));
                    begin = self.nr_shots;
                }
            }
            else
            {
                if begin >= self.nr_shots
                {
                    begin = end;
                }
            }
        }
        if begin < self.nr_shots
        {
            ranges.push((begin, self.nr_shots));
        }

        if ranges.is_empty()
        {
            return;
        }

        // Now apply the gate for the runs in the collected ranges
        let mut new_state_counts = ::std::collections::LinkedList::new();
        let mut range_idx = 0;
        let mut off = 0;
        while !self.state_counts.is_empty()
        {
            let (begin, end) = ranges[range_idx];
            let mut m = self.state_counts.pop_front().unwrap();

            if begin >= off + m.count
            {
                // Next range starts after this collection of runs. Append
                // the runs to the result unchanged, and move to the next run.
                off += m.count;
                new_state_counts.push_back(m);
                continue;
            }

            if begin > off
            {
                // The current range starts in the middle of this collection of
                // runs. Append the first part of the runs to the result
                // unchanged, and continue with the remaining runs.
                let mut new_m = m.clone();
                new_m.count = begin - off;
                m.count -= new_m.count;
                off = begin;
                new_state_counts.push_back(new_m);
            }

            if end < off + m.count
            {
                // The current range ends in the middle of this collection of
                // runs. Push the last part of the runs back on the stack to
                // process later, and continue with the remaining runs.
                let mut new_m = m.clone();
                m.count = end - off;
                new_m.count -= m.count;
                self.state_counts.push_front(new_m);
            }

            // Now apply the gate to the runs in this range
            if gate_bits == 1
            {
                let block_size = 1 << (self.nr_bits - bits[0]);
                let nr_blocks = 1 << bits[0];
                for i in 0..nr_blocks
                {
                    gate.apply_slice(&mut m.coefs.slice_mut(s![i*block_size..(i+1)*block_size]));
                }
            }
            else
            {
                let perm = gates::bit_permutation(self.nr_bits, bits);
                let inv_perm = perm.inverse();
                inv_perm.apply_vec(&mut m.coefs);
                gate.apply(&mut m.coefs);
                perm.apply_vec(&mut m.coefs);
            }

            off += m.count;
            new_state_counts.push_back(m);

            if end > off
            {
                ranges[range_idx].0 = off;
            }
            else
            {
                range_idx += 1;
                if range_idx >= ranges.len()
                {
                    break;
                }
            }
        }

        new_state_counts.append(&mut self.state_counts);
        self.state_counts = new_state_counts;
    }

    fn collapse(coefs: &mut cmatrix::CVector, block_size: usize, nr_blocks: usize,
        offset: usize, norm_sq: f64)
    {
        // Set coefficients for other bit to zero
        let mut off = offset;
        for _ in 0..nr_blocks
        {
            coefs.slice_mut(s![off..off+block_size]).fill(cmatrix::COMPLEX_ZERO);
            off += 2 * block_size;
        }

        // Renormalize
        *coefs *= num_complex::Complex::new(1.0 / norm_sq.sqrt(), 0.0);
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `bit` in the state. Measurement is done
    /// in the `z`-basis.
    pub fn measure(&mut self, bit: usize) -> ::ndarray::Array1<u8>
    {
        assert!(bit < self.nr_bits, "Invalid bit index");

        let block_size = 1 << (self.nr_bits - bit - 1);
        let nr_blocks = 1 << bit;
        let mut res = ndarray::Array1::zeros(self.nr_shots);
        let mut res_start = 0;

        let mut new_state_counts = ::std::collections::LinkedList::new();
        while !self.state_counts.is_empty()
        {
            let mut m = self.state_counts.pop_front().unwrap();

            // Compute chance of measuring 0
            let mut w0 = 0.0f64;
            let mut off = 0;
            for _ in 0..nr_blocks
            {
                w0 += m.coefs.slice(s![off..off+block_size])
                    .iter().map(|c| c.norm_sqr())
                    .sum::<f64>();
                off += 2 * block_size;
            }

            // Sometimes, the sum may add up to slightly more than 1, due to
            // numerical inaccuracies. This causes the Binomial distribution to
            // panic, so cap w0.
            w0 = w0.min(1.0);

            // Compute how many times we measure 0
            let distribution = rand::distributions::Binomial::new(m.count as u64, w0);
            let n0 = self.rng.sample(distribution) as usize;

            // Store the result
            res.slice_mut(s![res_start..res_start+n0]).fill(0);
            res.slice_mut(s![res_start+n0..res_start+m.count]).fill(1);
            res_start += m.count;

            // Collapse the wave function
            if n0 == m.count
            {
                Self::collapse(&mut m.coefs, block_size, nr_blocks, block_size, w0);
                new_state_counts.push_back(m);
            }
            else if n0 == 0
            {
                Self::collapse(&mut m.coefs, block_size, nr_blocks, 0, 1.0 - w0);
                new_state_counts.push_back(m);
            }
            else
            {
                let mut m1 = m.clone();
                m1.count = m.count - n0;
                m.count = n0;
                Self::collapse(&mut m.coefs, block_size, nr_blocks, block_size, w0);
                Self::collapse(&mut m1.coefs, block_size, nr_blocks, 0, 1.0 - w0);
                new_state_counts.push_back(m);
                new_state_counts.push_back(m1);
            }
        }
        self.state_counts = new_state_counts;

        res
    }

    /// Measure all qubits
    ///
    /// Measure all qubits in this state, and return the results through
    /// `result`.
    pub fn measure_all(&mut self) -> ndarray::Array2<u8>
    {
        let nr_coefs = 1 << self.nr_bits;
        let mut counts = ::std::collections::HashMap::new();
        let mut result = ndarray::Array2::<u8>::zeros((self.nr_bits, self.nr_shots));
        let mut res_off = 0;
        let mut new_state_counts = ::std::collections::LinkedList::new();
        for m in &self.state_counts
        {
            // Collect which measurements are made for this collection of runs.
            // Store this in a hash map with the measurement result as key, and
            // the number of times it was measured as value.
            counts.clear();
            let distr = rand::distributions::WeightedIndex::new(
                m.coefs.iter().map(|c| c.norm_sqr())
            ).unwrap();
            for idx in distr.sample_iter(&mut self.rng).take(m.count)
            {
                let entry = counts.entry(idx).or_insert(0);
                *entry += 1;
            }

            // For each unique measurement, store n copies of it in the result,
            // where n is the number of times it was measured. Furthermore,
            // create a new set of n runs with the resulting state after
            // the measurement for the new state.
            for (&idx, &count) in counts.iter()
            {
                let mut coefs = cmatrix::CVector::zeros(nr_coefs);
                coefs[idx] = cmatrix::COMPLEX_ONE;
                new_state_counts.push_back(StateCount::new(coefs, count));

                for irow in 0..self.nr_bits
                {
                    if idx & (1 << irow) != 0
                    {
                        // MSB in idx is bit 0 in state
                        result.slice_mut(s![self.nr_bits-irow-1, res_off..res_off+count]).fill(1);
                    }
                }

                res_off += count;
            }
        }

        self.state_counts = new_state_counts;

        result
    }

    /// Reset a qubit
    ///
    /// Reset the qubit with index `bit` to zero. This is done by measuring the
    /// bit, and rotating it back to zero if the result is 1.
    pub fn reset(&mut self, bit: usize)
    {
        let measurement = self.measure(bit);
        let control: Vec<bool> = measurement.iter().map(|&b| b != 0).collect();
        self.apply_conditional_gate(&control, &gates::X::new(), &[bit]);
    }

    /// Reset all qubits
    ///
    /// Reset all qubits in this experiment, returning the state to |00...0〉
    /// for all runs.
    pub fn reset_all(&mut self)
    {
        let mut coefs = cmatrix::CVector::zeros(1 << self.nr_bits);
        coefs[0] = cmatrix::COMPLEX_ONE;

        self.state_counts.clear();
        self.state_counts.push_back(StateCount::new(coefs, self.nr_shots));
    }
}

#[cfg(test)]
mod tests
{
    extern crate ndarray;

    use cmatrix;
    use gates;
    use super::QuState;

    #[test]
    fn test_new()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let s = QuState::new(1, 1024);
        assert_eq!(s.nr_bits, 1);
        assert_eq!(s.nr_shots, 1024);
        assert_eq!(s.state_counts.front().unwrap().count, 1024);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![o, z]);
        let s = QuState::new(3, 1500);
        assert_eq!(s.nr_bits, 3);
        assert_eq!(s.nr_shots, 1500);
        assert_eq!(s.state_counts.front().unwrap().count, 1500);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![o, z, z, z, z, z, z, z]);
    }

    #[test]
    fn test_from_qubit_coefs()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;

        // |0〉⊗|1〉
        let s = QuState::from_qubit_coefs(&[o, z, z, o], 1);
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 1);
        assert_eq!(s.state_counts.front().unwrap().count, 1);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![z, o, z, z]);
        // |1〉⊗|0〉
        let s = QuState::from_qubit_coefs(&[z, o, o, z], 13);
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 13);
        assert_eq!(s.state_counts.front().unwrap().count, 13);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![z, z, o, z]);
        // (H|0〉)⊗(Y|1〉), unnormalized
        let s = QuState::from_qubit_coefs(&[o, o, -i, z], 9);
        let x = ::std::f64::consts::FRAC_1_SQRT_2 * i;
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 9);
        assert_eq!(s.state_counts.front().unwrap().count, 9);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![-x, z, -x, z]);
    }

    #[test]
    fn test_apply_conditional_gate()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[false, false, true, true, false], &gates::X::new(), &[1]);
        assert_eq!(s.state_counts.len(), 3);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 2);
        assert_complex_vector_eq!(&sc.coefs, &array![o, z, z, z]);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 2);
        assert_complex_vector_eq!(&sc.coefs, &array![z, o, z, z]);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 1);
        assert_complex_vector_eq!(&sc.coefs, &array![o, z, z, z]);

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[false, false, true, true, true], &gates::X::new(), &[0]);
        assert_eq!(s.state_counts.len(), 2);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 2);
        assert_complex_vector_eq!(&sc.coefs, &array![o, z, z, z]);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 3);
        assert_complex_vector_eq!(&sc.coefs, &array![z, z, o, z]);

        let mut s = QuState::new(2, 5);
        s.apply_conditional_gate(&[true, false, true, true, false], &gates::H::new(), &[1]);
        assert_eq!(s.state_counts.len(), 4);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 1);
        assert_complex_vector_eq!(&sc.coefs, &array![x, x, z, z]);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 1);
        assert_complex_vector_eq!(&sc.coefs, &array![o, z, z, z]);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 2);
        assert_complex_vector_eq!(&sc.coefs, &array![x, x, z, z]);
        let sc = s.state_counts.pop_front().unwrap();
        assert_eq!(sc.count, 1);
        assert_complex_vector_eq!(&sc.coefs, &array![o, z, z, z]);
    }

    #[test]
    fn test_measure()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        // |0〉
        let mut s = QuState::new(1, 3);
        let m = s.measure(0);
        assert_eq!(m, array![0, 0, 0]);
//         assert_close_l2!(s.coefs.as_ref(), array![o; z]);
//
//         // |0〉⊗|0〉
//         let mut s = QuState::from_qubit_coefs(&[o, z, o, z], 3);
//         let m = s.measure(1);
//         assert_vector_eq!(m, vector![0, 0, 0]);
//         assert_close_l2!(s.coefs.as_ref(), array![o, o, o; z, z, z; z, z, z; z, z, z]);
//         let m = s.measure(0);
//         assert_vector_eq!(m, vector![0, 0, 0]);
//         assert_close_l2!(s.coefs.as_ref(), array![o, o, o; z, z, z; z, z, z; z, z, z]);
//
//         // (H|0〉)⊗(H|0〉), unnormalized
//         let mut s = QuState::from_qubit_coefs(&[o, o, o, o], 1024);
//         let m0 = s.measure(0);
//         for j in 0..s.nr_shots
//         {
//             let b = m0[j];
//             match b
//             {
//                 0 => assert_close_l2!(s.coefs.as_ref().col(j), array![x; x; z; z]),
//                 1 => assert_close_l2!(s.coefs.as_ref().col(j), array![z; z; x; x]),
//                 _ => panic!("Invalid value {} for bit", b)
//             }
//         }
//         // After collapse, a new measurement should yield the same result
//         let m0b = s.measure(0);
//         assert_eq!(m0b, m0);
//         // Measure second bit
//         let m1 = s.measure(1);
//         for j in 0..s.nr_shots
//         {
//             let b0 = m0[j];
//             let b1 = m1[j];
//             assert!(b1 == 0 || b1 == 1);
//             match (b0, b1)
//             {
//                 (0, 0) => assert_close_l2!(s.coefs.as_ref().col(j), array![o; z; z; z]),
//                 (0, 1) => assert_close_l2!(s.coefs.as_ref().col(j), array![z; o; z; z]),
//                 (1, 0) => assert_close_l2!(s.coefs.as_ref().col(j), array![z; z; o; z]),
//                 (1, 1) => assert_close_l2!(s.coefs.as_ref().col(j), array![z; z; z; o], eps=1.0e-15),
//                 _      => panic!("Invalid value {:?} for bits", (b0, b1))
//             }
//         }
    }
//
//     #[test]
//     fn test_apply_unary_gate()
//     {
//         let z = cmatrix::COMPLEX_ZERO;
//         let x = cmatrix::COMPLEX_HSQRT2;
//         let i = cmatrix::COMPLEX_I;
//
//         let h = gates::H::new();
//         let y = gates::Y::new();
//
//         let mut s = QuState::new(3, 1);
//         s.apply_gate(&h, &[0]);
//         assert_close_l2!(s.coefs.as_ref(), array![x; z; z; z; x; z; z; z]);
//
//         let mut s = QuState::new(3, 1);
//         s.apply_gate(&h, &[1]);
//         assert_close_l2!(s.coefs.as_ref(), array![x; z; x; z; z; z; z; z]);
//
//         let mut s = QuState::new(3, 1);
//         s.apply_gate(&y, &[2]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; i; z; z; z; z; z; z]);
//     }
//
//     #[test]
//     fn test_apply_binary_gate()
//     {
//         let z = cmatrix::COMPLEX_ZERO;
//         let o = cmatrix::COMPLEX_ONE;
//         let h = 0.5 * o;
//
//         let mut s = QuState::new(3, 1);
//         let cx = gates::CX::new();
//         s.apply_gate(&cx, &[0, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![o; z; z; z; z; z; z; z]);
//
//         let mut s = QuState::from_qubit_coefs(&[z, o, o, z, o, z], 1);
//         s.apply_gate(&cx, &[0, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; z; o; z]);
//
//         let mut s = QuState::from_qubit_coefs(&[z, o, o, z, o, z], 1);
//         s.apply_gate(&cx, &[0, 2]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; o; z; z]);
//
//         let mut s = QuState::from_qubit_coefs(&[z, o, o, z, o, z], 1);
//         let hh = gates::Kron::new(gates::H::new(), gates::H::new());
//         s.apply_gate(&hh, &[1, 2]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; h; h; h; h]);
//     }
//
//     #[test]
//     fn test_apply_n_ary_gate()
//     {
//         let z = cmatrix::COMPLEX_ZERO;
//         let o = cmatrix::COMPLEX_ONE;
//         let x = cmatrix::COMPLEX_HSQRT2;
//         let hx = 0.5 * x;
//
//         let mut s = QuState::new(3, 1);
//         let ccx = gates::CCX::new();
//         s.apply_gate(&ccx, &[0, 1, 2]);
//         assert_close_l2!(s.coefs.as_ref(), array![o; z; z; z; z; z; z; z]);
//
//         let mut s = QuState::from_qubit_coefs(&[z, o, z, o, o, z], 1);
//         let ccx = gates::CCX::new();
//         s.apply_gate(&ccx, &[0, 2, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; z; o; z]);
//         s.apply_gate(&ccx, &[0, 1, 2]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; z; z; o]);
//
//         let mut s = QuState::from_qubit_coefs(&[x, -x, x, -x, x, -x], 1);
//         s.apply_gate(&ccx, &[0, 2, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![hx; -hx; -hx; hx; -hx; -hx; hx; hx]);
//     }

    #[test]
    fn test_measure_all()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let nr_bits = 3;
        let nr_shots = 5;

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o, z, o], nr_shots);
        let result = s.measure_all();
        assert_eq!(result.shape(), [nr_bits, nr_shots]);
        assert!(result.iter().all(|&b| b == 1));

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o, o, z], nr_shots);
        let result = s.measure_all();
        assert_eq!(result.shape(), [nr_bits, nr_shots]);
        assert!(result.slice(s![0..2, ..]).iter().all(|&b| b == 1));
        assert!(result.slice(s![2, ..]).iter().all(|&b| b == 0));

        let mut s = QuState::new(nr_bits, nr_shots);
        s.apply_gate(&gates::H::new(), &[2]);
        let result = s.measure_all();
        assert_eq!(result.shape(), [nr_bits, nr_shots]);
        assert!(result.slice(s![0..2, ..]).iter().all(|&b| b == 0));
    }

    #[test]
    fn test_reset()
    {
        let nr_runs = 10;

        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let mut s = QuState::from_qubit_coefs(&[o, z], nr_runs);
        s.reset(0);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![o, z]);

        let mut s = QuState::from_qubit_coefs(&[z, o], nr_runs);
        s.reset(0);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![o, z]);

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o], nr_runs);
        s.reset(0);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![z, o, z, z]);

        let mut s = QuState::from_qubit_coefs(&[z, o, z, o], nr_runs);
        s.reset(1);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![z, z, o, z]);

        let mut s = QuState::from_qubit_coefs(&[x, -x, o, z], nr_runs);
        s.reset(0);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![o, z, z, z]);

        let mut s = QuState::from_qubit_coefs(&[x, -x, o, z], nr_runs);
        s.reset(1);
        assert_complex_vector_eq!(&s.state_counts.front().unwrap().coefs, &array![x, z, -x, z]);
    }

    #[test]
    fn test_reset_all()
    {
        let nr_bits = 5;
        let nr_runs = 100;

        let mut s = QuState::new(nr_bits, nr_runs);
        s.apply_gate(&gates::H::new(), &[2]);
        s.apply_gate(&gates::X::new(), &[0]);
        s.apply_gate(&gates::H::new(), &[4]);

        s.reset_all();
        assert!(s.state_counts.len() == 1);
        let sc = s.state_counts.front().unwrap();
        assert_eq!(sc.count, nr_runs);
        let mut coefs = cmatrix::CVector::zeros(1 << nr_bits);
        coefs[0] = cmatrix::COMPLEX_ONE;
        assert_complex_vector_eq!(&sc.coefs, &coefs);
    }
}
