extern crate rulinalg;

mod hadamard;
mod identity;

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

pub use gates::hadamard::Hadamard;
pub use gates::identity::Identity;
