extern crate rulinalg;

mod hadamard;
mod identity;
mod kron;

use cmatrix;

pub trait Gate
{
    /// Return a short description of the gate. This may be the name of the
    /// gate (e.g. `"H"`, `"C_X"`), or the way the gate was constructed (like
    /// `"kron(I, Z)"`)
    fn description(&self) -> &str;

    /// Return a matrix describing the unitary transformation that the gate
    /// provides
    fn matrix(&self) -> cmatrix::CMatrix;
}

pub trait UnaryGate: Gate
{
    /// Apply a unary gate (working on a single qubit) to quantum state `state`.
    /// The number of rows in `state` must be even, with the first half
    /// corresponding to qustates with basis state |0&rang; for the affected
    /// qubit, and the second half to states with basis state |1&rang.
    fn apply_unary(&self, state: &mut cmatrix::CMatrix)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let mat = self.matrix();

        let n = state.rows() / 2;
        let m = state.cols();
        let s0 = state.sub_slice([0, 0], n, m).into_matrix();
        let s1 = state.sub_slice([n, 0], n, m).into_matrix();

        state.sub_slice_mut([0, 0], n, m).set_to(&s0*mat.at([0, 0]) + &s1*mat.at([0, 1]));
        state.sub_slice_mut([n, 0], n, m).set_to(&s0*mat.at([1, 0]) + &s1*mat.at([1, 1]));
    }
}

pub use gates::hadamard::Hadamard;
pub use gates::identity::Identity;
pub use gates::kron::Kron;
