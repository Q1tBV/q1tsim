extern crate num_complex;
extern crate rulinalg;

use cmatrix;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

mod cx;
mod hadamard;
mod identity;
mod kron;
mod x;
mod y;

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
    /// Apply a binary gate (working on two qubits) to quantum state `state`.
    /// The number of rows `n` in `state` must be a multiple of four, with the
    /// first block of `n/4` rows corresponding to qustates with basis states
    /// |00〉 for the affected qubits, the second block to |01〉, the
    /// third to |10〉 and the last to |11〉.
    fn apply_binary(&self, state: &mut rulinalg::matrix::Matrix<num_complex::Complex64>)
    {
        let (n, m) = (state.rows(), state.cols());
        self.apply_binary_slice(&mut state.sub_slice_mut([0, 0], n, m));
    }

    /// Apply a binary gate (working on two qubits) to quantum state `state`.
    /// The number of rows `n` in `state` must be a multiple of four, with the
    /// first block of `n/4` rows corresponding to qustates with basis states
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

pub use gates::cx::CX;
pub use gates::hadamard::Hadamard;
pub use gates::identity::Identity;
pub use gates::kron::Kron;
pub use gates::x::X;
pub use gates::y::Y;
