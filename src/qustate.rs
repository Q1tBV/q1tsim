extern crate num_complex;
extern crate rand;
extern crate rulinalg;

use cmatrix;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

/// Quantum state.
///
/// Struct Qustate represents the quantum state of the system. It consists of a
/// (normalized) superposition of basis states, ∑<sub>i</sub>a<sub>i</sub>|i&rang;,
/// where each basis function |i&rang; is a Kronecker product of quantum bits.
pub struct QuState
{
    /// The number of qubits in this state
    nr_bits: usize,
    /// Coefficients of the basis states in this state
    coefs: cmatrix::CMatrix
}

impl QuState
{
    /// Create a new qustate of `n` qubits, all initialized to |0&rang;.
    pub fn new(n: usize) -> Self
    {
        QuState { nr_bits: n, coefs: cmatrix::CMatrix::eye(1 << n, 1) }
    }

    /// Create a new qustate from qubit coefficients.
    ///
    /// Create a new qustate as a direct product of qubits, where the
    /// coefficients of the |0〉 and |1〉 states in the qubits are given in
    // `bit_coefs`. This array must be of size `2*n`, where `n` is the number
    /// of qubits in the system.
    pub fn from_qubit_coefs(bit_coefs: &[num_complex::Complex64]) -> Self
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

        QuState { nr_bits: nr_bits, coefs: coefs }
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on one of the qubits in the state.
    pub fn measure(&mut self, i: usize) -> ::rulinalg::vector::Vector<u8>
    {
        assert!(i < self.nr_bits, "Invalid bit index");

        let nr_measurements = self.coefs.cols();

        let block_size = 1 << (self.nr_bits - i - 1);
        let nr_blocks = 1 << i;
        let mut w0 = ::rulinalg::vector::Vector::zeros(nr_measurements);
        let mut off = 0;
        for _ in 0..nr_blocks
        {
            let elems: Vec<f64> = self.coefs
                .sub_slice([off, 0], block_size, nr_measurements)
                .iter().map(|c| c.norm_sqr()).collect();
            let block = rulinalg::matrix::Matrix::new(block_size, nr_measurements,
                elems);
            w0 += block.sum_rows();
            off += 2 * block_size;
        }

        let mut res = ::rulinalg::vector::Vector::zeros(nr_measurements);
        for (i, &w) in w0.iter().enumerate()
        {
            let mut off;
            let norm_sq;
            if rand::random::<f64>() >= w
            {
                res[i] = 1;
                off = 0;
                norm_sq = 1.0 - w;
            }
            else
            {
                off = block_size;
                norm_sq = w;
            }
            for _ in 0..nr_blocks
            {
                self.coefs.col_mut(i).sub_slice_mut([off, 0], block_size, 1)
                    .apply(&|_| cmatrix::COMPLEX_ZERO);
                off += 2 * block_size;
            }

            let inorm = 1.0 / norm_sq.sqrt();
            self.coefs.col_mut(i).apply(&|c| c * inorm);
        }

        res
    }
}

#[cfg(test)]
mod tests
{
    use cmatrix;
    use qustate::QuState;

    #[test]
    fn test_new()
    {
        let s = QuState::new(1);
        assert_eq!(s.nr_bits, 1);
        assert_matrix_eq!(s.coefs.real(), matrix![1.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0], comp=float);
        let s = QuState::new(3);
        assert_eq!(s.nr_bits, 3);
        assert_matrix_eq!(s.coefs.real(), matrix![1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0], comp=float);
    }

    #[test]
    fn test_from_qubit_coefs()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;

        // |0〉⊗|1〉
        let s = QuState::from_qubit_coefs(&[o, z, z, o]);
        assert_matrix_eq!(s.coefs.real(), matrix![0.0; 1.0; 0.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        // |1〉⊗|0〉
        let s = QuState::from_qubit_coefs(&[z, o, o, z]);
        assert_matrix_eq!(s.coefs.real(), matrix![0.0; 0.0; 1.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        // (H|0〉)⊗(Y|1〉), unnormalized
        let s = QuState::from_qubit_coefs(&[o, o, -i, z]);
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        assert_matrix_eq!(s.coefs.real(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![ -x; 0.0;  -x; 0.0], comp=float);
    }

    #[test]
    fn test_measure()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;

        // |0>
        let mut s = QuState::new(1);
        assert_eq!(s.measure(0)[0], 0);
        assert_matrix_eq!(s.coefs.real(), matrix![1.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0], comp=float);

        // (H|0>)⊗|0>, unnormalized
        let mut s = QuState::from_qubit_coefs(&[o, z, o, z]);
        let m = s.measure(0);
        assert!(m[0] == 0 || m[0] == 1);
        if m[0] == 0
        {
            assert_matrix_eq!(s.coefs.real(), matrix![1.0; 0.0; 0.0; 0.0], comp=float);
        }
        else
        {
            assert_matrix_eq!(s.coefs.real(), matrix![0.0; 0.0; 1.0; 0.0], comp=float);
        }
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        let m1 = s.measure(0);
        assert_eq!(m1, m);
        let m = s.measure(1);
        assert_eq!(m[0], 0);

        // (H|0>)⊗(H|0>), unnormalized
        let mut s = QuState::from_qubit_coefs(&[o, o, o, o]);
        let m = s.measure(0);
        assert!(m[0] == 0 || m[0] == 1);
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        if m[0] == 0
        {
            assert_matrix_eq!(s.coefs.real(), matrix![x; x; 0.0; 0.0], comp=float);
        }
        else
        {
            assert_matrix_eq!(s.coefs.real(), matrix![0.0; 0.0; x; x], comp=float);
        }
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        let m1 = s.measure(0);
        assert_eq!(m1, m);
    }
}
