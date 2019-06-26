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

mod controlled;
mod composite;
mod cx;
mod cy;
mod cz;
mod hadamard;
mod identity;
mod kron;
mod parameter;
mod rx;
mod ry;
mod rz;
mod s;
mod staticloop;
mod swap;
mod t;
mod v;
mod u1;
mod u2;
mod u3;
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
pub fn bit_permutation(nr_bits: usize, affected_bits: &[usize]) -> crate::permutation::Permutation
{
    let mut idxs: Vec<usize> = (0..(1 << nr_bits)).collect();
    idxs.sort_by_key(|&i| get_sort_key(i, nr_bits, affected_bits));
    crate::permutation::Permutation::new(idxs).unwrap()
}

/// Apply a gate
///
/// Apply gate `gate` operating on the bits in `bits` to a vector `vec`. The
/// number of elements in `vec` must be 2^`nr_bits`.
pub fn apply_gate_slice<G>(mut vec: crate::cmatrix::CVecSliceMut, gate: &G,
    bits: &[usize], nr_bits: usize)
where G: Gate + ?Sized
{
    let gate_bits = gate.nr_affected_bits();
    assert!(gate_bits == bits.len(),
        "The number of bits affected by the {} gate should be {}, but {} bits were provided.",
        gate.description(), gate_bits, bits.len()
    );
    assert!(vec.len() == 1 << nr_bits,
        "The number of bits in the state does not match the total number of bits.");

    if let &[bit] = bits
    {
        let block_size = 1 << (nr_bits - bit);
        let nr_blocks = 1 << bit;
        for i in 0..nr_blocks
        {
            gate.apply_slice(vec.slice_mut(s![i*block_size..(i+1)*block_size]));
        }
    }
    else
    {
        let perm = bit_permutation(nr_bits, bits);
        let mut work = crate::cmatrix::CVector::zeros(1 << nr_bits);
        perm.apply_vec_into(vec.view(), work.view_mut());
        gate.apply_slice(work.view_mut());
        perm.apply_inverse_vec_into(work.view(), vec);
    }
}

/// Apply a gate
///
/// Apply gate `gate` operating on the bits in `bits` to a matrix `matrix`. The
/// number of rows in `matrix` must be 2^`nr_bits`.
pub fn apply_gate_mat_slice<G>(mut matrix: crate::cmatrix::CMatSliceMut, gate: &G,
    bits: &[usize], nr_bits: usize)
where G: Gate + ?Sized
{
    let gate_bits = gate.nr_affected_bits();
    assert!(gate_bits == bits.len(),
        "The number of bits affected by the {} gate should be {}, but {} bits were provided.",
        gate.description(), gate_bits, bits.len()
    );
    assert!(matrix.rows() == 1 << nr_bits,
        "The number of bits in the state does not match the total number of bits.");

    if let &[bit] = bits
    {
        let block_size = 1 << (nr_bits - bit);
        let nr_blocks = 1 << bit;
        for i in 0..nr_blocks
        {
            gate.apply_mat_slice(matrix.slice_mut(s![i*block_size..(i+1)*block_size, ..]));
        }
    }
    else
    {
        let perm = bit_permutation(nr_bits, bits);
        let mut work = crate::cmatrix::CMatrix::zeros((matrix.rows(), matrix.cols()));
        ndarray::Zip::from(matrix.gencolumns()).and(work.gencolumns_mut())
            .apply(|s, d| perm.apply_vec_into(s, d));
        gate.apply_mat_slice(work.view_mut());
        ndarray::Zip::from(matrix.gencolumns_mut()).and(work.gencolumns())
            .apply(|d, s| perm.apply_inverse_vec_into(s, d));
    }
}

pub trait Gate
{
    /// Cost of this gate.
    ///
    /// An estimate of the cost of using this gate. This can be used e.g. in
    /// optimizing circuits, or when trying to decompose a transformation
    /// in an optimal gate sequence. The default implementation returns
    /// `std::f64::INFINITY`.
    fn cost(&self) -> f64 { ::std::f64::INFINITY }

    /// Return a short description of the gate. This may be the name of the
    /// gate (e.g. `"H"`, `"CX"`), or the way the gate was constructed (like
    /// `"I⊗Z"`)
    fn description(&self) -> &str;

    /// The number of qubits affected by this gate.
    fn nr_affected_bits(&self) -> usize;

    /// Return a matrix describing the unitary transformation that the gate
    /// provides
    fn matrix(&self) -> crate::cmatrix::CMatrix;

    /// Apply a gate.
    ///
    /// Apply a gate to quantum state `state`. The number of rows `r` in `state`
    /// must be a multiple of 2<sup>`n`</sup>, where `n` is the number of qubits
    /// this gate acts upon. The rows must be ordered, such that the first block
    /// of `r`/2<sup>`n`</sup> rows corresponds to qustates with basis states
    /// |00...0⟩ for the affected qubits, the second block to |00...1⟩, etc.,
    /// up until |11...1⟩.
    fn apply(&self, state: &mut crate::cmatrix::CVector)
    {
        self.apply_slice(state.view_mut());
    }

    /// Apply a gate.
    ///
    /// Apply a gate to quantum state `state`. The number of rows `r` in `state`
    /// must be a multiple of 2<sup>`n`</sup>, where `n` is the number of qubits
    /// this gate acts upon. The rows must be ordered, such that the first block
    /// of `r`/2<sup>`n`</sup> rows corresponds to qustates with basis states
    /// |00...0⟩ for the affected qubits, the second block to |00...1⟩, etc.,
    /// up until |11...1⟩.
    fn apply_mat(&self, state: &mut crate::cmatrix::CMatrix)
    {
        self.apply_mat_slice(state.view_mut());
    }

    /// Apply a gate.
    ///
    /// Apply a gate to quantum state `state`. The number of rows `r` in `state`
    /// must be a multiple of 2<sup>`n`</sup>, where `n` is the number of qubits
    /// this gate acts upon. The rows must be ordered, such that the first block
    /// of `r`/2<sup>`n`</sup> rows corresponds to qustates with basis states
    /// |00...0⟩ for the affected qubits, the second block to |00...1⟩, etc.,
    /// up until |11...1⟩.
    fn apply_slice(&self, mut state: crate::cmatrix::CVecSliceMut)
    {
        let nr_bits = self.nr_affected_bits();
        assert!(state.len() % (1 << nr_bits) == 0,
            "The number of rows in the state is {}, which is not valid for a {}-bit gate.",
            state.len(), nr_bits);

        let mat = self.matrix();
        let n = state.len() >> nr_bits;
        if nr_bits == 1
        {
            let s0 = state.slice(s![..n]).to_owned();
            let s1 = state.slice(s![n..]).to_owned();

            state.slice_mut(s![..n]).assign(&(&s0*mat[[0, 0]] + &s1*mat[[0, 1]]));
            state.slice_mut(s![n..]).assign(&(&s0*mat[[1, 0]] + &s1*mat[[1, 1]]));
        }
        else if nr_bits == 2
        {
            let s0 = state.slice(s![     ..n]).to_owned();
            let s1 = state.slice(s![  n..2*n]).to_owned();
            let s2 = state.slice(s![2*n..3*n]).to_owned();
            let s3 = state.slice(s![3*n..   ]).to_owned();

            state.slice_mut(s![..n]).assign(
                &(&s0*mat[[0, 0]] + &s1*mat[[0, 1]] + &s2*mat[[0, 2]] + &s3*mat[[0, 3]])
            );
            state.slice_mut(s![n..2*n]).assign(
                &(&s0*mat[[1, 0]] + &s1*mat[[1, 1]] + &s2*mat[[1, 2]] + &s3*mat[[1, 3]])
            );
            state.slice_mut(s![2*n..3*n]).assign(
                &(&s0*mat[[2, 0]] + &s1*mat[[2, 1]] + &s2*mat[[2, 2]] + &s3*mat[[2, 3]])
            );
            state.slice_mut(s![3*n..]).assign(
                &(&s0*mat[[3, 0]] + &s1*mat[[3, 1]] + &s2*mat[[3, 2]] + &s3*mat[[3, 3]])
            );
        }
        else
        {
            let mut res = crate::cmatrix::CVector::zeros(state.len());

            for i in 0..(1 << nr_bits)
            {
                let mut slice = res.slice_mut(s![i*n..(i+1)*n]);
                for j in 0..(1 << nr_bits)
                {
                    let x = &state.slice(s![j*n..(j+1)*n]).to_owned() * mat[[i,j]];
                    slice += &x;
                }
            }

            state.assign(&res);
        }
    }

    /// Apply a gate.
    ///
    /// Apply a gate to quantum state `state`. The number of rows `r` in `state`
    /// must be a multiple of 2<sup>`n`</sup>, where `n` is the number of qubits
    /// this gate acts upon. The rows must be ordered, such that the first block
    /// of `r`/2<sup>`n`</sup> rows corresponds to qustates with basis states
    /// |00...0⟩ for the affected qubits, the second block to |00...1⟩, etc.,
    /// up until |11...1⟩.
    fn apply_mat_slice(&self, mut state: crate::cmatrix::CMatSliceMut)
    {
        let nr_bits = self.nr_affected_bits();
        assert!(state.len() % (1 << nr_bits) == 0,
            "The number of rows in the state is {}, which is not valid for a {}-bit gate.",
            state.len(), nr_bits);

        let mat = self.matrix();
        let n = state.rows() >> nr_bits;
        if nr_bits == 1
        {
            let s0 = state.slice(s![..n, ..]).to_owned();
            let s1 = state.slice(s![n.., ..]).to_owned();

            state.slice_mut(s![..n, ..]).assign(&(&s0*mat[[0, 0]] + &s1*mat[[0, 1]]));
            state.slice_mut(s![n.., ..]).assign(&(&s0*mat[[1, 0]] + &s1*mat[[1, 1]]));
        }
        else if nr_bits == 2
        {
            let s0 = state.slice(s![     ..n, ..]).to_owned();
            let s1 = state.slice(s![  n..2*n, ..]).to_owned();
            let s2 = state.slice(s![2*n..3*n, ..]).to_owned();
            let s3 = state.slice(s![3*n..   , ..]).to_owned();

            state.slice_mut(s![..n, ..]).assign(
                &(&s0*mat[[0, 0]] + &s1*mat[[0, 1]] + &s2*mat[[0, 2]] + &s3*mat[[0, 3]])
            );
            state.slice_mut(s![n..2*n, ..]).assign(
                &(&s0*mat[[1, 0]] + &s1*mat[[1, 1]] + &s2*mat[[1, 2]] + &s3*mat[[1, 3]])
            );
            state.slice_mut(s![2*n..3*n, ..]).assign(
                &(&s0*mat[[2, 0]] + &s1*mat[[2, 1]] + &s2*mat[[2, 2]] + &s3*mat[[2, 3]])
            );
            state.slice_mut(s![3*n.., ..]).assign(
                &(&s0*mat[[3, 0]] + &s1*mat[[3, 1]] + &s2*mat[[3, 2]] + &s3*mat[[3, 3]])
            );
        }
        else
        {
            let mut res = crate::cmatrix::CMatrix::zeros((state.rows(), state.cols()));

            for i in 0..(1 << nr_bits)
            {
                let mut slice = res.slice_mut(s![i*n..(i+1)*n, ..]);
                for j in 0..(1 << nr_bits)
                {
                    let x = &state.slice(s![j*n..(j+1)*n, ..]).to_owned() * mat[[i,j]];
                    slice += &x;
                }
            }

            state.assign(&res);
        }
    }

    /// Check the number of bits
    ///
    /// Check if the number of bit indices `n` is equal to the number
    /// of bits this gate operates on. If not, return an InvalidNrBits error.
    fn check_nr_bits(&self, n: usize) -> crate::error::Result<()>
    {
        if self.nr_affected_bits() != n
        {
            Err(crate::error::Error::InvalidNrBits(n, self.nr_affected_bits(),
                String::from(self.description())))
        }
        else
        {
            Ok(())
        }
    }

    /// Whether this gate is a stabilizer gate
    ///
    /// Return `true` if this gate is a stabilizer gate, i.e. if conjugating
    /// a Pauli operator (or tensor product thereof for multi-bit gates) with
    /// this gate, again returns a Pauli operator. Circuits consisting of only
    /// these types of gates can be simulated more efficiently. The default
    /// implementation returns `false`.
    fn is_stabilizer(&self) -> bool
    {
        false
    }

    /// Conjugate Pauli operators
    ///
    /// Conjugate the (tensor product of) Pauli operator `ops` with this gate,
    /// and return the sign of the resulting operator. I.e., if this gate is
    /// called `G`, replace the operator `O` described by ops with
    /// `GOG`<sup>`†`</sup>. This of course only works if the result can be
    /// expressed in Pauli operators, so only if this gate can be written
    /// in terms of `H`, `S` and `CX` gates. The default implementation of
    /// this function returns a `NotAStabilizer` error.
    fn conjugate(&self, _ops: &mut [crate::stabilizer::PauliOp]) -> crate::error::Result<bool>
    {
        Err(crate::error::Error::NotAStabilizer(String::from(self.description())))
    }
}

#[cfg(test)]
fn gate_test<G>(gate: G, state: &mut crate::cmatrix::CMatrix, result: &crate::cmatrix::CMatrix)
where G: Gate
{
    for i in 0..state.cols()
    {
        gate.apply_slice(state.column_mut(i));
    }
    assert_complex_matrix_eq!(&*state, result);
}

pub use self::parameter::Parameter;

pub use self::controlled::{C, CH, CRX, CRY, CRZ, CS, CSdg, CT, CTdg, CU1, CU2, CU3, CV, CVdg, CCRX, CCRY, CCRZ, CCX, CCZ};
pub use self::composite::Composite;
pub use self::cx::CX;
pub use self::cy::CY;
pub use self::cz::CZ;
pub use self::hadamard::H;
pub use self::identity::I;
pub use self::kron::Kron;
pub use self::rx::RX;
pub use self::ry::RY;
pub use self::rz::RZ;
pub use self::s::{S, Sdg};
pub use self::staticloop::Loop;
pub use self::t::{T, Tdg};
pub use self::swap::Swap;
pub use self::u1::U1;
pub use self::u2::U2;
pub use self::u3::U3;
pub use self::v::{V, Vdg};
pub use self::x::X;
pub use self::y::Y;
pub use self::z::Z;
