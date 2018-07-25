extern crate num_complex;
extern crate rulinalg;

use cmatrix;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

mod cx;
mod ccx;
mod custom;
mod hadamard;
mod identity;
mod kron;
mod x;
mod y;
mod z;

/// Generates the new row number for the row that initially
/// was at number `idx`, in a system of `nr_bits` qubits, in which a
/// gate is operating on the qubits in `affected_bits`.
fn get_sort_key(idx: usize, nr_bits: usize, affected_bits: &[usize]) -> usize
{
    let mut res = 0;
    for b in affected_bits
    {
        let s = nr_bits - b - 1;
        res = (res << 1) | ((idx >> s) & 1);
    }
    res
}

/// Reorder bits.
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
/// This function returns a permutation matrix `P`, such that the matrix
/// `P (G ⊗ I ⊗ ... ⊗ I) P`<sup>`T`</sup> describes the effect of operating with a
/// gate `G` on bits `affected_bits` in a `nr_bits`-sized system.
pub fn bit_permutation(nr_bits: usize, affected_bits: &[usize])
    -> rulinalg::matrix::PermutationMatrix<num_complex::Complex64>
{
    let mut idxs = (0..(1 << nr_bits)).collect::<Vec<usize>>();
    idxs.sort_by_key(|&i| get_sort_key(i, nr_bits, affected_bits));
    rulinalg::matrix::PermutationMatrix::from_array(idxs).unwrap()
}

pub trait Gate
{
    /// Return a short description of the gate. This may be the name of the
    /// gate (e.g. `"H"`, `"CX"`), or the way the gate was constructed (like
    /// `"I⊗Z"`)
    fn description(&self) -> &str;

    /// Return a matrix describing the unitary transformation that the gate
    /// provides
    fn matrix(&self) -> cmatrix::CMatrix;
}

pub trait UnaryGate: Gate
{
    /// Expanded matrix.
    ///
    /// Return the matrix describing the operation of this gate on qubit `bit`
    /// in a system of `nr_bits  qubits.
    fn expanded_matrix(&self, bit: usize, nr_bits: usize) -> cmatrix::CMatrix
    {
        assert!(bit < nr_bits, "Invalid bit index.");

        if nr_bits == 1
        {
            self.matrix()
        }
        else if bit == 0
        {
            self.matrix().kron(&cmatrix::CMatrix::eye_sq(1 << (nr_bits-1)))
        }
        else if bit == nr_bits-1
        {
            cmatrix::CMatrix::eye_sq(1 << bit).kron(&self.matrix())
        }
        else
        {
            cmatrix::CMatrix::eye_sq(1 << bit).kron(&self.matrix())
                .kron(&cmatrix::CMatrix::eye_sq(1 << (nr_bits-bit-1)))
        }
    }

    /// Apply a unary gate (working on a single qubit) to quantum state `state`.
    /// The number of rows in `state` must be even, with the first half
    /// corresponding to qustates with basis state |0〉 for the affected
    /// qubit, and the second half to states with basis state |1〉.
    fn apply_unary(&self, state: &mut rulinalg::matrix::Matrix<num_complex::Complex64>)
    {
        let (n, m) = (state.rows(), state.cols());
        self.apply_unary_slice(&mut state.sub_slice_mut([0, 0], n, m));
    }

    /// Apply a unary gate (working on a single qubit) to quantum state `state`.
    /// The number of rows in `state` must be even, with the first half
    /// corresponding to qustates with basis state |0〉 for the affected
    /// qubit, and the second half to states with basis state |1〉.
    fn apply_unary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let mat = self.matrix();

        let n = state.rows() / 2;
        let m = state.cols();
        let s0 = state.sub_slice([0, 0], n, m).into_matrix();
        let s1 = state.sub_slice([n, 0], n, m).into_matrix();

        state.sub_slice_mut([0, 0], n, m).set_to(&s0*mat[[0, 0]] + &s1*mat[[0, 1]]);
        state.sub_slice_mut([n, 0], n, m).set_to(&s0*mat[[1, 0]] + &s1*mat[[1, 1]]);
    }
}

pub trait BinaryGate: Gate
{
    /// Expanded matrix.
    ///
    /// Return the matrix describing the operation of this gate on qubits `bit0`
    /// and `bit1` in a system of `nr_bits  qubits.
    fn expanded_matrix(&self, bit0: usize, bit1: usize, nr_bits: usize) -> cmatrix::CMatrix
    {
        assert!(bit0 < nr_bits, "Invalid bit index.");
        assert!(bit1 < nr_bits, "Invalid bit index.");

        let p = bit_permutation(nr_bits, &[bit0, bit1]);
        let g = if nr_bits == 2
            {
                self.matrix()
            }
            else
            {
                self.matrix().kron(&cmatrix::CMatrix::eye_sq(1 << (nr_bits-2)))
            };
        // XXX WARNING!
        // Bug/inconsistency in rulinalg: right multiplying by a permutation
        // matrix actually results in a multiplication by its inverse. So instead
        // of the expected
        //
        // let invp = p.inverse();
        // cmatrix::CMatrix::from_matrix(p * g.as_ref() * invp)
        //
        // we do
        let invp = p.clone();
        cmatrix::CMatrix::from_matrix(p.as_matrix() * g.as_ref() * invp)
    }

    /// Apply a binary gate (working on two qubits) to quantum state `state`.
    /// The number of rows `r` in `state` must be a multiple of four, with the
    /// first block of `r/4` rows corresponding to qustates with basis states
    /// |00〉 for the affected qubits, the second block to |01〉, the
    /// third to |10〉 and the last to |11〉.
    fn apply_binary(&self, state: &mut rulinalg::matrix::Matrix<num_complex::Complex64>)
    {
        let (n, m) = (state.rows(), state.cols());
        self.apply_binary_slice(&mut state.sub_slice_mut([0, 0], n, m));
    }

    /// Apply a binary gate (working on two qubits) to quantum state `state`.
    /// The number of rows `r` in `state` must be a multiple of four, with the
    /// first block of `r/4` rows corresponding to qustates with basis states
    /// |00〉 for the affected qubits, the second block to |01〉, the
    /// third to |10〉 and the last to |11〉.
    fn apply_binary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        assert!(state.rows() % 4 == 0, "Number of rows is not a multiple of four.");

        let mat = self.matrix();

        let n = state.rows() / 4;
        let m = state.cols();
        let s0 = state.sub_slice([0, 0], n, m).into_matrix();
        let s1 = state.sub_slice([n, 0], n, m).into_matrix();
        let s2 = state.sub_slice([2*n, 0], n, m).into_matrix();
        let s3 = state.sub_slice([3*n, 0], n, m).into_matrix();

        state.sub_slice_mut([0, 0], n, m).set_to(
            &s0*mat[[0, 0]] + &s1*mat[[0, 1]] + &s2*mat[[0, 2]] + &s3*mat[[0, 3]]
        );
        state.sub_slice_mut([n, 0], n, m).set_to(
            &s0*mat[[1, 0]] + &s1*mat[[1, 1]] + &s2*mat[[1, 2]] + &s3*mat[[1, 3]]
        );
        state.sub_slice_mut([2*n, 0], n, m).set_to(
            &s0*mat[[2, 0]] + &s1*mat[[2, 1]] + &s2*mat[[2, 2]] + &s3*mat[[2, 3]]
        );
        state.sub_slice_mut([3*n, 0], n, m).set_to(
            &s0*mat[[3, 0]] + &s1*mat[[3, 1]] + &s2*mat[[3, 2]] + &s3*mat[[3, 3]]
        );
    }
}

pub trait NaryGate: Gate
{
    /// The number of qubits affected by this gate.
    fn nr_affected_bits(&self) -> usize;

    /// Expanded matrix.
    ///
    /// Return the matrix describing the operation of this gate on the qubits
    /// in `bits`, in a system of `nr_bits  qubits.
    fn expanded_matrix(&self, bits: &[usize], nr_bits: usize) -> cmatrix::CMatrix
    {
        let n = bits.len();
        let p = bit_permutation(nr_bits, bits);
        let g = if nr_bits == n
            {
                self.matrix()
            }
            else
            {
                self.matrix().kron(&cmatrix::CMatrix::eye_sq(1 << (nr_bits-n)))
            };
        // XXX WARNING!
        // Bug/inconsistency in rulinalg: right multiplying by a permutation
        // matrix actually results in a multiplication by its inverse. So instead
        // of the expected
        //
        // let invp = p.inverse();
        // cmatrix::CMatrix::from_matrix(p * g.as_ref() * invp)
        //
        // we do
        let invp = p.clone();
        cmatrix::CMatrix::from_matrix(p * g.as_ref() * invp)
    }

    /// Apply a `n`-ary gate (working on multiple qubits) to quantum state `state`.
    /// The number of rows `r` in `state` must be a multiple of 2<sup>`n`</sup>,
    /// with the first block of `r`/2<sup>`n`</sup> rows corresponding to qustates
    /// with basis states |00...0〉 for the affected qubits, the second block to
    /// |00...1〉, etc. up until |11...1〉.
    fn apply_n_ary(&self, state: &mut rulinalg::matrix::Matrix<num_complex::Complex64>)
    {
        let (nr, nc) = (state.rows(), state.cols());
        self.apply_n_ary_slice(&mut state.sub_slice_mut([0, 0], nr, nc));
    }

    /// Apply a `n`-ary gate (working on multiple qubits) to quantum state `state`.
    /// The number of rows `r` in `state` must be a multiple of 2<sup>`n`</sup>,
    /// with the first block of `r`/2<sup>`n`</sup> rows corresponding to qustates
    /// with basis states |00...0〉 for the affected qubits, the second block to
    /// |00...1〉, etc. up until |11...1〉.
    fn apply_n_ary_slice(&self, state: &mut rulinalg::matrix::MatrixSliceMut<num_complex::Complex64>)
    {
        let n = self.nr_affected_bits();
        assert!(state.rows() % (1 << n) == 0, "Number of rows is not a multiple of 2^{}.", n);

        let mat = self.matrix();
        let mut res = rulinalg::matrix::Matrix::zeros(state.rows(), state.cols());

        let nr = state.rows() / (1 << n);
        let nc = state.cols();
        for i in 0..(1 << n)
        {
            let mut slice = res.sub_slice_mut([i*nr, 0], nr, nc);
            for j in 0..(1 << n)
            {
                slice += state.sub_slice([j*nr, 0], nr, nc) * mat[[i,j]];
            }
        }

        state.sub_slice_mut([0, 0], nr << n, nc).set_to(res);
    }
}

pub use gates::cx::CX;
pub use gates::ccx::CCX;
pub use gates::custom::Custom;
pub use gates::hadamard::Hadamard;
pub use gates::identity::Identity;
pub use gates::kron::Kron;
pub use gates::x::X;
pub use gates::y::Y;
pub use gates::z::Z;
