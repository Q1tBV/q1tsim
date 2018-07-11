extern crate rand;

use cmatrix;

use rulinalg::matrix::BaseMatrix;

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
    /// Create a new qustate of `n` qubits, all qinitialized to |0&rang;.
    pub fn new(n: usize) -> Self
    {
        QuState { nr_bits: n, coefs: cmatrix::CMatrix::eye(1 << n, 1) }
    }

    /// Create a new qustate from qubit coefficients.
    ///
    /// Create a new qustate as a direct product of qubits, where the real
    /// parts of the coefficients of the |0&rang; and |1&rang; states in the
    /// qubits are given in `r`, and the imaginary parts in `i`. Both arrays
    /// must be of size 2`n`, where `n` is the number of qubits in the system.
    pub fn from_qubit_coefs(r: &[f64], i: &[f64]) -> Self
    {
        assert!(r.len() == i.len(), "Length of real and imaginary part is different");
        assert!(r.len() % 2 == 0, "Length of coefficient array is not even");

        let nr_bits = r.len() / 2;

        let mut coefs = cmatrix::CMatrix::eye(1, 1);
        for (r, i) in r.chunks(2).zip(i.chunks(2))
        {
            let bit = if i[0] == 0.0 && i[1] == 0.0
            {
                let norm = r[0].hypot(r[1]);
                cmatrix::CMatrix::new_real(2, 1, vec![r[0]/norm, r[1]/norm])
            }
            else if r[0] == 0.0 && r[1] == 0.0
            {
                let norm = i[0].hypot(i[1]);
                cmatrix::CMatrix::new_imag(2, 1, vec![i[0]/norm, i[1]/norm])
            }
            else
            {
                let norm = (r[0]*r[0] + i[0]*i[0] + r[1]*r[1] + i[1]*i[1]).sqrt();
                cmatrix::CMatrix::new_complex(2, 1, vec![r[0]/norm, r[1]/norm], vec![i[0]/norm, i[1]/norm])
            };
            coefs = coefs.kron(&bit);
        }

        QuState { nr_bits: nr_bits, coefs: coefs }
    }

    /// Measure a qubit.
    ///
    /// Perform a measurement on one of the qubits in the state.
    pub fn measure(&mut self, i: usize) -> u8
    {
        assert!(i < self.nr_bits, "Invalid bit index");

        let weights = self.coefs.abs_sq();

        let block_size = 1 << (self.nr_bits - i - 1);
        let nr_blocks = 1 << i;
        let mut w0 = 0.0;
        let mut off = 0;
        for _ in 0..nr_blocks
        {
            w0 += weights.sub_slice([off, 0], block_size, 1).sum();
            off += 2 * block_size;
        }

        let ret1 = rand::random::<f64>() >= w0;

        let mut off = if ret1 { 0 } else { block_size };
        for _ in 0..nr_blocks
        {
            self.coefs.clear(off, 0, block_size, 1);
            off += 2 * block_size;
        }

        let norm_sq = if ret1 { 1.0 - w0 } else { w0 };
        self.coefs *= 1.0 / norm_sq.sqrt();

        ret1 as u8
    }

}

#[cfg(test)]
mod tests
{
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
        // |0>⊗|1>
        let r = [1.0, 0.0, 0.0, 1.0];
        let i = [0.0; 4];
        let s = QuState::from_qubit_coefs(&r, &i);
        assert_matrix_eq!(s.coefs.real(), matrix![0.0; 1.0; 0.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        // |1>⊗|0>
        let r = [0.0, 1.0, 1.0, 0.0];
        let i = [0.0; 4];
        let s = QuState::from_qubit_coefs(&r, &i);
        assert_matrix_eq!(s.coefs.real(), matrix![0.0; 0.0; 1.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        // (H|0>)⊗(Y|1>), unnormalized
        let r = [1.0, 1.0,  0.0, 0.0];
        let i = [0.0, 0.0, -1.0, 0.0];
        let s = QuState::from_qubit_coefs(&r, &i);
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        assert_matrix_eq!(s.coefs.real(), matrix![0.0; 0.0; 0.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![ -x; 0.0;  -x; 0.0], comp=float);
    }

    #[test]
    fn test_measure()
    {
        // |0>
        let mut s = QuState::new(1);
        assert_eq!(s.measure(0), 0);
        assert_matrix_eq!(s.coefs.real(), matrix![1.0; 0.0], comp=float);
        assert_matrix_eq!(s.coefs.imag(), matrix![0.0; 0.0], comp=float);

        // (H|0>)⊗|0>, unnormalized
        let r = [1.0, 0.0, 1.0, 0.0];
        let i = [0.0, 0.0, 0.0, 0.0];
        let mut s = QuState::from_qubit_coefs(&r, &i);
        let m = s.measure(0);
        assert!(m == 0 || m == 1);
        if m == 0
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
        assert_eq!(m, 0);

        // (H|0>)⊗(H|0>), unnormalized
        let r = [1.0, 1.0, 1.0, 1.0];
        let i = [0.0, 0.0, 0.0, 0.0];
        let mut s = QuState::from_qubit_coefs(&r, &i);
        let m = s.measure(0);
        assert!(m == 0 || m == 1);
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        if m == 0
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
