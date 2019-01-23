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


extern crate ndarray;
extern crate num_complex;

use cmatrix;
use permutation;

#[macro_use] mod controlled;
mod composite;
mod hadamard;
mod identity;
mod kron;
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
pub fn bit_permutation(nr_bits: usize, affected_bits: &[usize]) -> permutation::Permutation
{
    let mut idxs: Vec<usize> = (0..(1 << nr_bits)).collect();
    idxs.sort_by_key(|&i| get_sort_key(i, nr_bits, affected_bits));
    permutation::Permutation::new(idxs).inverse()
}

/// Apply a gate
///
/// Apply gate `gate` operating on the bits in `bits` to a matrix `matrix`. The
/// number of rows in `matrix` must be 2^`nr_bits`.
fn apply_gate<G>(matrix: &mut cmatrix::CMatrix, gate: &G, bits: &[usize], nr_bits: usize)
where G: Gate + ?Sized
{
    let gate_bits = gate.nr_affected_bits();
    assert!(gate_bits == bits.len(),
        "The number of bits affected by the gate does not match the provided number of bits.");

    if gate_bits == 1
    {
        let block_size = 1 << (nr_bits - bits[0]);
        let nr_blocks = 1 << bits[0];
        for i in 0..nr_blocks
        {
            gate.apply_mat_slice(&mut matrix.slice_mut(s![i*block_size..(i+1)*block_size, ..]));
        }
    }
    else
    {
        let perm = bit_permutation(nr_bits, bits);
        let inv_perm = perm.inverse();

        // FIXME: see if we can find a cleaner way to do this
        *matrix = matrix.select(ndarray::Axis(0), inv_perm.indices());
        gate.apply_mat(matrix);
        *matrix = matrix.select(ndarray::Axis(0), perm.indices());
    }
}

pub trait Gate
{
    /// An estimate of the cost of using this gate
    fn cost(&self) -> f64;

    /// Return a short description of the gate. This may be the name of the
    /// gate (e.g. `"H"`, `"CX"`), or the way the gate was constructed (like
    /// `"I⊗Z"`)
    fn description(&self) -> &str;

    /// The number of qubits affected by this gate.
    fn nr_affected_bits(&self) -> usize;

    /// Return a matrix describing the unitary transformation that the gate
    /// provides
    fn matrix(&self) -> cmatrix::CMatrix;

    /// Expanded matrix.
    ///
    /// Return the matrix describing the operation of this gate on the qubits
    /// in `bits`, in a system of `nr_bits  qubits.
    fn expanded_matrix(&self, bits: &[usize], nr_bits: usize) -> cmatrix::CMatrix
    {
        let gate_bits = self.nr_affected_bits();

        assert_eq!(bits.len(), gate_bits,
            "The number of bit indices provided does not match the number of bits affected by this gate.");
        for &bit in bits
        {
            assert!(bit < nr_bits, "Invalid bit index {} for {}-bit system.", bit, nr_bits);
        }

        if gate_bits == 1
        {
            cmatrix::kron_mat(&cmatrix::CMatrix::eye(1 << bits[0]),
                &cmatrix::kron_mat(&self.matrix(),
                    &cmatrix::CMatrix::eye(1 << (nr_bits-bits[0]-1))))
        }
        else
        {
            let gate_mat = cmatrix::kron_mat(&self.matrix(),
                &cmatrix::CMatrix::eye(1 << (nr_bits-gate_bits)));
            bit_permutation(nr_bits, bits).transform(&gate_mat)
        }
    }

    /// Apply a gate.
    ///
    /// Apply a gate to quantum state `state`. The number of rows `r` in `state`
    /// must be a multiple of 2<sup>`n`</sup>, where `n` is the number of qubits
    /// this gate acts upon. The rows must be ordered, such that the first block
    /// of `r`/2<sup>`n`</sup> rows corresponds to qustates with basis states
    /// |00...0〉 for the affected qubits, the second block to |00...1〉, etc.,
    /// up until |11...1〉.
    fn apply(&self, state: &mut cmatrix::CVector)
    {
        self.apply_slice(&mut state.view_mut());
    }

    fn apply_mat(&self, state: &mut cmatrix::CMatrix)
    {
        self.apply_mat_slice(&mut state.view_mut());
    }

    /// Apply a gate.
    ///
    /// Apply a gate to quantum state `state`. The number of rows `r` in `state`
    /// must be a multiple of 2<sup>`n`</sup>, where `n` is the number of qubits
    /// this gate acts upon. The rows must be ordered, such that the first block
    /// of `r`/2<sup>`n`</sup> rows corresponds to qustates with basis states
    /// |00...0〉 for the affected qubits, the second block to |00...1〉, etc.,
    /// up until |11...1〉.
    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
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
            let mut res = cmatrix::CVector::zeros(state.len());

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

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
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
            let mut res = cmatrix::CMatrix::zeros((state.rows(), state.cols()));

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
}

#[cfg(test)]
fn gate_test<G>(gate: G, state: &mut cmatrix::CMatrix, result: &cmatrix::CMatrix)
where G: Gate
{
    for i in 0..state.cols()
    {
        gate.apply_slice(&mut state.column_mut(i));
    }
    assert_complex_matrix_eq!(&*state, result);
}

pub use self::controlled::{C, CH, CRX, CRY, CRZ, CS, CSdg, CT, CTdg, CU1, CU2, CU3, CV, CVdg, CX, CY, CZ, CCX, CCZ};
pub use self::composite::Composite;
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
