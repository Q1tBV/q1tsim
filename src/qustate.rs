extern crate ndarray;
extern crate num_complex;
extern crate rand;

use cmatrix;
use gates;

use self::rand::Rng;

/// Measured quantum state
///
/// The full state of the system is stored as a collection of measurements,
/// combined with the number of times that measurement was made.
#[derive(Clone, Debug, PartialEq)]
pub struct Measurement
{
    coefs: cmatrix::CVector,
    count: usize
}

impl Measurement
{
    fn new(coefs: cmatrix::CVector, count: usize) -> Self
    {
        Measurement { coefs: coefs, count: count }
    }
}

/// Quantum state.
///
/// Struct Qustate represents the quantum state of the system. It consists of a
/// (normalized) superposition of basis states, ∑<sub>i</sub>a<sub>i</sub>|i〉,
/// where each basis function |i〉 is a Kronecker product of quantum bits.
#[derive(Debug)]
pub struct QuState
{
    /// The number of qubits in this state
    nr_bits: usize,
    /// The number of separate runs for evolving this state
    nr_shots: usize,
    /// Coefficients of the basis states in this state
    measurements: ::std::collections::LinkedList<Measurement>,
    /// Random number generator, for measurements
    rng: rand::ThreadRng
}

impl QuState
{
    /// Create a new qustate of `nr_bits` qubits, all initialized to |0〉, which
    /// will be measured `nr_shots` times.
    pub fn new(nr_bits: usize, nr_shots: usize) -> Self
    {
        let mut coefs = cmatrix::CVector::zeros(1 << nr_bits);
        coefs[0] = cmatrix::COMPLEX_ONE;

        let mut measurements = ::std::collections::LinkedList::new();
        measurements.push_back(Measurement::new(coefs, nr_shots));

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            measurements: measurements,
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

        let mut measurements = ::std::collections::LinkedList::new();
        measurements.push_back(Measurement::new(coefs, nr_shots));

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            measurements: measurements,
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
            for m in self.measurements.iter_mut()
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
            for m in self.measurements.iter_mut()
            {
                inv_perm.apply_vec(&mut m.coefs);
                gate.apply(&mut m.coefs);
                perm.apply_vec(&mut m.coefs);
            }
        }
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
    /// Perform a measurement on qubit `bit` in the state.
    pub fn measure(&mut self, bit: usize) -> ::ndarray::Array1<u8>
    {
        assert!(bit < self.nr_bits, "Invalid bit index");

        let block_size = 1 << (self.nr_bits - bit - 1);
        let nr_blocks = 1 << bit;
        let mut res = ndarray::Array1::zeros(self.nr_shots);
        let mut res_start = 0;
        let mut msr_idx = 0;

        // The below approach with splitting the list and stitching it back
        // together again is awful, but unfortunately, it's the best we can do
        // until IterMut::insert_next() stabilizes.
        while msr_idx < self.measurements.len()
        {
            let mut tail = self.measurements.split_off(msr_idx);
            let mut m = tail.pop_front().unwrap();

            // Compute chance of measuring 0
            let mut w0 = 0.0f64;
            let mut off = 0;
            for _ in 0..nr_blocks
            {
                w0 += &m.coefs.slice(s![off..off+block_size]).iter().map(|c| c.norm_sqr()).sum();
                off += 2 * block_size;
            }

            // Sometimes, the sum may add up to slightly more than 1, die to
            // numerical inaccuracies. This causes the Binomial distribution to
            // panic, so cap w0.
            w0 = w0.min(1.0);

            // Compute how many times we measure 0
            let distribution = rand::distributions::Binomial::new(m.count as u64, w0);
            let n0 = self.rng.sample(distribution) as usize;
            //let mut n0 = 0;
            //for _ in 0..m.count
            //{
            //    n0 += self.rng.gen_bool(w0) as usize;
            //}

            // Store the result
            res.slice_mut(s![res_start..res_start+n0]).fill(0);
            res.slice_mut(s![res_start+n0..res_start+m.count]).fill(1);
            res_start += m.count;

            // Collapse the wave function
            if n0 == m.count
            {
                Self::collapse(&mut m.coefs, block_size, nr_blocks, block_size, w0);
                self.measurements.push_back(m);
                msr_idx += 1;
            }
            else if n0 == 0
            {
                Self::collapse(&mut m.coefs, block_size, nr_blocks, 0, 1.0 - w0);
                self.measurements.push_back(m);
                msr_idx += 1;
            }
            else
            {
                let mut m1 = m.clone();
                m1.count = m.count - n0;
                m.count = n0;

                Self::collapse(&mut m.coefs, block_size, nr_blocks, block_size, w0);
                self.measurements.push_back(m);
                msr_idx += 1;

                Self::collapse(&mut m1.coefs, block_size, nr_blocks, 0, 1.0 - w0);
                self.measurements.push_back(m1);
                msr_idx += 1;

            }

            // Re-append te unprocessed runs
            self.measurements.append(&mut tail);
        }

        res
    }
}

#[cfg(test)]
mod tests
{
    use cmatrix;
    use gates;
    use qustate::QuState;

    #[test]
    fn test_new()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let s = QuState::new(1, 1024);
        assert_eq!(s.nr_bits, 1);
        assert_eq!(s.nr_shots, 1024);
        assert_eq!(s.measurements.front().unwrap().count, 1024);
        assert_complex_vector_eq!(&s.measurements.front().unwrap().coefs, &array![o, z]);
        let s = QuState::new(3, 1500);
        assert_eq!(s.nr_bits, 3);
        assert_eq!(s.nr_shots, 1500);
        assert_eq!(s.measurements.front().unwrap().count, 1500);
        assert_complex_vector_eq!(&s.measurements.front().unwrap().coefs, &array![o, z, z, z, z, z, z, z]);
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
        assert_eq!(s.measurements.front().unwrap().count, 1);
        assert_complex_vector_eq!(&s.measurements.front().unwrap().coefs, &array![z, o, z, z]);
        // |1〉⊗|0〉
        let s = QuState::from_qubit_coefs(&[z, o, o, z], 13);
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 13);
        assert_eq!(s.measurements.front().unwrap().count, 13);
        assert_complex_vector_eq!(&s.measurements.front().unwrap().coefs, &array![z, z, o, z]);
        // (H|0〉)⊗(Y|1〉), unnormalized
        let s = QuState::from_qubit_coefs(&[o, o, -i, z], 9);
        let x = ::std::f64::consts::FRAC_1_SQRT_2 * i;
        assert_eq!(s.nr_bits, 2);
        assert_eq!(s.nr_shots, 9);
        assert_eq!(s.measurements.front().unwrap().count, 9);
        assert_complex_vector_eq!(&s.measurements.front().unwrap().coefs, &array![-x, z, -x, z]);
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
//         let mut s = QuState::from_qubit_coefs(&vec![z, o, o, z, o, z], 1);
//         s.apply_gate(&cx, &[0, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; z; o; z]);
//
//         let mut s = QuState::from_qubit_coefs(&vec![z, o, o, z, o, z], 1);
//         s.apply_gate(&cx, &[0, 2]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; o; z; z]);
//
//         let mut s = QuState::from_qubit_coefs(&vec![z, o, o, z, o, z], 1);
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
//         let mut s = QuState::from_qubit_coefs(&vec![z, o, z, o, o, z], 1);
//         let ccx = gates::CCX::new();
//         s.apply_gate(&ccx, &[0, 2, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; z; o; z]);
//         s.apply_gate(&ccx, &[0, 1, 2]);
//         assert_close_l2!(s.coefs.as_ref(), array![z; z; z; z; z; z; z; o]);
//
//         let mut s = QuState::from_qubit_coefs(&vec![x, -x, x, -x, x, -x], 1);
//         s.apply_gate(&ccx, &[0, 2, 1]);
//         assert_close_l2!(s.coefs.as_ref(), array![hx; -hx; -hx; hx; -hx; -hx; hx; hx]);
//     }
}
