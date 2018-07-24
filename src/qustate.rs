extern crate num_complex;
extern crate rand;
extern crate rulinalg;

use cmatrix;
use gates;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

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
    /// The number of measurement runs on this state
    nr_shots: usize,
    /// Coefficients of the basis states in this state
    coefs: cmatrix::CMatrix,
}

impl QuState
{
    /// Create a new qustate of `nr_bits` qubits, all initialized to |0〉, which
    /// will be measured `nr_shots` times.
    pub fn new(nr_bits: usize, nr_shots: usize) -> Self
    {
        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            coefs: cmatrix::CMatrix::eye(1 << nr_bits, 1)
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

        let mut coefs = cmatrix::CMatrix::eye(1, 1);
        for c in bit_coefs.chunks(2)
        {
            let norm = (c[0].norm_sqr() + c[1].norm_sqr()).sqrt();
            let bit = cmatrix::CMatrix::new(2, 1, vec![c[0]/norm, c[1]/norm]);
            coefs = coefs.kron(&bit);
        }

        QuState
        {
            nr_bits: nr_bits,
            nr_shots: nr_shots,
            coefs: coefs
        }
    }

    /// Apply a unary quantum gate `gate` on qubit `bit` in this state.
    pub fn apply_unary_gate<G>(&mut self, gate: &G, bit: usize)
    where G: gates::UnaryGate + ?Sized
    {
        let nr_measurements = self.coefs.as_ref().cols();
        let block_size = 1 << (self.nr_bits - bit);
        let nr_blocks = 1 << bit;
        for i in 0..nr_blocks
        {
            gate.apply_unary_slice(&mut self.coefs.as_mut()
                .sub_slice_mut([i*block_size, 0], block_size, nr_measurements));
        }
    }

    /// Get sort key of a row.
    ///
    /// When applying multi-bit gates, the rows in the state are shuffled
    /// such that:
    /// * The first half of the rows correspond to components with the first
    ///   affected bit being 0, the second half to those with this bit being 1.
    /// * Within each of these two blocks, the first half corresponds to
    ///   components with the second bit 0, the second half to those with the
    ///   second bit 1.
    /// * And so on, for each affected bit.
    ///
    /// This function generates the permutation that creates this order.
    fn get_sort_key(&self, idx: usize, bits: &[usize]) -> usize
    {
        let mut res = 0;
        for b in bits
        {
            let s = self.nr_bits - b - 1;
            res = (res << 1) | ((idx >> s) & 1);
        }
        res
    }

    /// Apply a binary quantum gate `gate` on qubits `bit0` and `bit1` in this
    /// state.
    pub fn apply_binary_gate<G>(&mut self, gate: &G, bit0: usize, bit1: usize)
    where G: gates::BinaryGate + ?Sized
    {
        let mut idxs = (0..(1 << self.nr_bits)).collect::<Vec<usize>>();
        idxs.sort_by_key(|&i| self.get_sort_key(i, &[bit0, bit1]));
        let perm = ::rulinalg::matrix::PermutationMatrix::from_array(idxs).unwrap();
        let inv_perm = perm.inverse();

        inv_perm.permute_rows_in_place(self.coefs.as_mut());
        gate.apply_binary(self.coefs.as_mut());
        perm.permute_rows_in_place(self.coefs.as_mut());
    }

    /// Apply a n-ary quantum gate `gate` on the qubits from `bits` in this state.
    pub fn apply_n_ary_gate<G>(&mut self, gate: &G, bits: &[usize])
    where G: gates::NaryGate + ?Sized
    {
        assert!(gate.nr_affected_bits() == bits.len(),
            "The number of bits affected by the gate does not match the provided number of bits.");

        let mut idxs = (0..(1 << self.nr_bits)).collect::<Vec<usize>>();
        idxs.sort_by_key(|&i| self.get_sort_key(i, bits));
        let perm = ::rulinalg::matrix::PermutationMatrix::from_array(idxs).unwrap();
        let inv_perm = perm.inverse();

        inv_perm.permute_rows_in_place(self.coefs.as_mut());
        gate.apply_n_ary(self.coefs.as_mut());
        perm.permute_rows_in_place(self.coefs.as_mut());
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `bit` in the state.
    pub fn measure(&mut self, bit: usize) -> ::rulinalg::vector::Vector<u8>
    {
        assert!(bit < self.nr_bits, "Invalid bit index");

        if self.coefs.as_ref().cols() == 1 && self.nr_shots > 1
        {
            // state has not been measured yet
            let new_coefs = self.coefs.as_ref() * rulinalg::matrix::Matrix::ones(1, self.nr_shots);
            self.coefs = cmatrix::CMatrix::from_matrix(new_coefs);
        }

        let block_size = 1 << (self.nr_bits - bit - 1);
        let nr_blocks = 1 << bit;
        let mut w0 = rulinalg::vector::Vector::zeros(self.nr_shots);
        let mut off = 0;
        for _ in 0..nr_blocks
        {
            let elems: Vec<f64> = self.coefs.as_ref()
                .sub_slice([off, 0], block_size, self.nr_shots)
                .iter().map(|c| c.norm_sqr()).collect();
            let block = rulinalg::matrix::Matrix::new(block_size, self.nr_shots, elems);
            w0 += block.sum_rows();
            off += 2 * block_size;
        }

        let mut res = ::rulinalg::vector::Vector::zeros(self.nr_shots);
        for (j, &w) in w0.iter().enumerate()
        {
            let mut off;
            let norm_sq;
            // Check if the measurement yields 0 or 1
            if rand::random::<f64>() >= w
            {
                res[j] = 1;
                off = 0;
                norm_sq = 1.0 - w;
            }
            else
            {
                off = block_size;
                norm_sq = w;
            }

            // Collapse the measured state
            for _ in 0..nr_blocks
            {
                self.coefs.as_mut().sub_slice_mut([off, j], block_size, 1)
                    .apply(&|_| cmatrix::COMPLEX_ZERO);
                off += 2 * block_size;
            }

            // Renormalize
            let inorm = 1.0 / norm_sq.sqrt();
            self.coefs.as_mut().col_mut(j).apply(&|c| c * inorm);
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
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_new()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        let s = QuState::new(1, 1024);
        assert_eq!(s.nr_bits, 1);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o; z]);
        let s = QuState::new(3, 1500);
        assert_eq!(s.nr_bits, 3);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o; z; z; z; z; z; z; z]);
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
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; o; z; z]);
        // |1〉⊗|0〉
        let s = QuState::from_qubit_coefs(&[z, o, o, z], 1);
        assert_eq!(s.nr_bits, 2);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; z; o; z]);
        // (H|0〉)⊗(Y|1〉), unnormalized
        let s = QuState::from_qubit_coefs(&[o, o, -i, z], 1);
        let x = ::std::f64::consts::FRAC_1_SQRT_2 * i;
        assert_eq!(s.nr_bits, 2);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![-x; z; -x; z]);
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
        assert_vector_eq!(m, vector![0, 0, 0]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o, o, o; z, z, z]);

        // |0〉⊗|0〉
        let mut s = QuState::from_qubit_coefs(&[o, z, o, z], 3);
        let m = s.measure(1);
        assert_vector_eq!(m, vector![0, 0, 0]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o, o, o; z, z, z; z, z, z; z, z, z]);
        let m = s.measure(0);
        assert_vector_eq!(m, vector![0, 0, 0]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o, o, o; z, z, z; z, z, z; z, z, z]);

        // (H|0〉)⊗(H|0〉), unnormalized
        let mut s = QuState::from_qubit_coefs(&[o, o, o, o], 1024);
        let m0 = s.measure(0);
        for j in 0..s.nr_shots
        {
            let b = m0[j];
            match b
            {
                0 => assert_complex_matrix_eq!(s.coefs.as_ref().col(j), matrix![x; x; z; z]),
                1 => assert_complex_matrix_eq!(s.coefs.as_ref().col(j), matrix![z; z; x; x]),
                _ => panic!("Invalid value {} for bit", b)
            }
        }
        // After collapse, a new measurement should yield the same result
        let m0b = s.measure(0);
        assert_eq!(m0b, m0);
        // Measure second bit
        let m1 = s.measure(1);
        for j in 0..s.nr_shots
        {
            let b0 = m0[j];
            let b1 = m1[j];
            assert!(b1 == 0 || b1 == 1);
            match (b0, b1)
            {
                (0, 0) => assert_complex_matrix_eq!(s.coefs.as_ref().col(j), matrix![o; z; z; z]),
                (0, 1) => assert_complex_matrix_eq!(s.coefs.as_ref().col(j), matrix![z; o; z; z]),
                (1, 0) => assert_complex_matrix_eq!(s.coefs.as_ref().col(j), matrix![z; z; o; z]),
                (1, 1) => assert_complex_matrix_eq!(s.coefs.as_ref().col(j), matrix![z; z; z; o], eps=1.0e-15),
                _      => panic!("Invalid value {:?} for bits", (b0, b1))
            }
        }

    }

    #[test]
    fn test_apply_unary_gate()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;

        let h = gates::Hadamard::new();
        let y = gates::Y::new();

        let mut s = QuState::new(3, 1);
        s.apply_unary_gate(&h, 0);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![x; z; z; z; x; z; z; z]);

        let mut s = QuState::new(3, 1);
        s.apply_unary_gate(&h, 1);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![x; z; x; z; z; z; z; z]);

        let mut s = QuState::new(3, 1);
        s.apply_unary_gate(&y, 2);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; i; z; z; z; z; z; z]);
    }

    #[test]
    fn test_apply_binary_gate()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let h = 0.5 * o;

        let mut s = QuState::new(3, 1);
        let cx = gates::CX::new();
        s.apply_binary_gate(&cx, 0, 1);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o; z; z; z; z; z; z; z]);

        let mut s = QuState::from_qubit_coefs(&vec![z, o, o, z, o, z], 1);
        s.apply_binary_gate(&cx, 0, 1);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; z; z; z; z; z; o; z]);

        let mut s = QuState::from_qubit_coefs(&vec![z, o, o, z, o, z], 1);
        s.apply_binary_gate(&cx, 0, 2);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; z; z; z; z; o; z; z]);

        let mut s = QuState::from_qubit_coefs(&vec![z, o, o, z, o, z], 1);
        let hh = gates::Kron::new(gates::Hadamard::new(), gates::Hadamard::new());
        s.apply_binary_gate(&hh, 1, 2);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; z; z; z; h; h; h; h]);
    }

    #[test]
    fn test_apply_n_ary_gate()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let hx = 0.5 * x;

        let mut s = QuState::new(3, 1);
        let ccx = gates::CCX::new();
        s.apply_n_ary_gate(&ccx, &[0, 1, 2]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![o; z; z; z; z; z; z; z]);

        let mut s = QuState::from_qubit_coefs(&vec![z, o, z, o, o, z], 1);
        let ccx = gates::CCX::new();
        s.apply_n_ary_gate(&ccx, &[0, 2, 1]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; z; z; z; z; z; o; z]);
        s.apply_n_ary_gate(&ccx, &[0, 1, 2]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![z; z; z; z; z; z; z; o]);

        let mut s = QuState::from_qubit_coefs(&vec![x, -x, x, -x, x, -x], 1);
        s.apply_n_ary_gate(&ccx, &[0, 2, 1]);
        assert_complex_matrix_eq!(s.coefs.as_ref(), matrix![hx; -hx; -hx; hx; -hx; -hx; hx; hx]);
    }
}
