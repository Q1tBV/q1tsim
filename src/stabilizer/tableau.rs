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

use crate::stabilizer::PauliOp;

/// Structure describing the expected outcome of a measurement
#[derive(Debug, PartialEq)]
pub enum MeasurementInfo
{
    /// Measurement yields zero if the associated parameter is `false`, one otherwise
    Deterministic(bool),
    /// Measurement yields zero or one with equal probability. The associated
    /// parameter is the index of the row that will need collapsing
    Random(usize)
}

/// Structure describing a single stabilizer state
#[derive(Clone)]
pub struct StabilizerTableau
{
    /// The number of qubits in the state
    nr_bits: usize,
    /// Generators of the stabilizer group for this state
    xz: Vec<u64>,
    /// Signs of the generators
    signs: Vec<u64>
}

impl StabilizerTableau
{
    pub fn new(nr_bits: usize) -> Self
    {
        let xz_size = (2*nr_bits*nr_bits + 0x3f) >> 6;
        let signs_size = (nr_bits + 0x3f) >> 6;
        let mut res = StabilizerTableau {
            nr_bits: nr_bits,
            xz: vec![0; xz_size],
            signs: vec![0; signs_size]
        };

        for i in 0..nr_bits
        {
            res.set(i, i, PauliOp::Z);
        }

        res
    }

    #[inline(always)]
    fn bit_indices(&self, i: usize, j: usize) -> (usize, usize)
    {
        let idx = 2*(i*self.nr_bits + j);
        (idx >> 6, idx & 0x3f)
    }

    fn get_bits(&self, i: usize, j: usize) -> u64
    {
        let (byte_idx, bit_idx) = self.bit_indices(i, j);
        (self.xz[byte_idx] >> bit_idx) & 0x03
    }

    fn get(&self, i: usize, j: usize) -> PauliOp
    {
        PauliOp::from_bits(self.get_bits(i, j))
    }

    fn get_z(&self, i: usize, j: usize) -> bool
    {
        let (byte_idx, bit_idx) = self.bit_indices(i, j);
        self.xz[byte_idx] & (1 << bit_idx) != 0
    }

    fn get_x(&self, i: usize, j: usize) -> bool
    {
        let (byte_idx, bit_idx) = self.bit_indices(i, j);
        self.xz[byte_idx] & (2 << bit_idx) != 0
    }

    fn set_bits(&mut self, i: usize, j: usize, op: u64)
    {
        let (byte_idx, bit_idx) = self.bit_indices(i, j);
        self.xz[byte_idx] = (self.xz[byte_idx] & !(0x03 << bit_idx)) | ((op & 0x03) << bit_idx);
    }

    fn set(&mut self, i: usize, j: usize, op: PauliOp)
    {
        self.set_bits(i, j, op.to_bits());
    }

    fn get_sign(&self, i: usize) -> bool
    {
        let (byte_idx, bit_idx) = (i >> 6, i & 0x3f);
        (self.signs[byte_idx] & (1 << bit_idx)) != 0
    }

    fn set_sign(&mut self, i: usize, sign: bool)
    {
        let (byte_idx, bit_idx) = (i >> 6, i & 0x3f);
        self.signs[byte_idx] = (self.signs[byte_idx] & !(1 << bit_idx))
            | ((sign as u64) << bit_idx);
    }

    fn xor_sign(&mut self, i: usize, sign: bool)
    {
        if sign
        {
            let (byte_idx, bit_idx) = (i >> 6, i & 0x3f);
            self.signs[byte_idx] ^= 1 << bit_idx;
        }
    }

    fn swap_rows(&mut self, i0: usize, i1: usize)
    {
        for j in 0..self.nr_bits
        {
            let b = self.get_bits(i0, j);
            self.set_bits(i0, j, self.get_bits(i1, j));
            self.set_bits(i1, j, b);
        }

        let sign = self.get_sign(i0);
        self.set_sign(i0, self.get_sign(i1));
        self.set_sign(i1, sign);
    }

    fn multiply_row(&mut self, i0: usize, i1: usize)
    {
        const PHASE_FACTORS: [u8; 16] = [
            0, 0, 0, 0,
            0, 0, 1, 3,
            0, 1, 0, 3,
            0, 1, 3, 0
        ];

        let mut i_pow = if self.get_sign(i1) { 2 } else { 0 };
        for j in 0..self.nr_bits
        {
            let xz0 = self.get_bits(i0, j);
            let xz1 = self.get_bits(i1, j);
            self.set_bits(i0, j, xz0 ^ xz1);

            let idx = (xz0 << 2 | xz1) as usize;
            i_pow = (i_pow + PHASE_FACTORS[idx]) & 0x03;
        }

        assert!(i_pow == 0 || i_pow == 2);
        self.xor_sign(i0, i_pow == 2);
    }

    fn normalize(&mut self)
    {
        let n = self.nr_bits;

        let mut i = 0;
        for j in 0..n
        {
            if let Some(k) = (i..n).filter(|&k| self.get_x(k, j)).next()
            {
                self.swap_rows(i, k);
                for m in 0..n
                {
                    if m != i && self.get_x(m, j)
                    {
                        self.multiply_row(m, i);
                    }
                }
                i += 1;
            }
        }

        for j in 0..n
        {
            if let Some(k) = (i..n).filter(|&k| self.get_z(k, j)).next()
            {
                self.swap_rows(i, k);
                for m in 0..n
                {
                    if m != i && self.get_z(m, j)
                    {
                        self.multiply_row(m, i);
                    }
                }
                i += 1;
            }
        }
    }

    /// Apply a n-ary quantum gate `gate` on the qubits from `bits` in this state.
    pub fn apply_gate<G>(&mut self, gate: &G, bits: &[usize]) -> crate::error::Result<()>
    where G: crate::gates::Gate + ?Sized
    {
        let mut ops = vec![];
        for i in 0..self.nr_bits
        {
            ops.clear();
            ops.extend(bits.iter().map(|&j| self.get(i, j)));
            let flip_sign = gate.conjugate(&mut ops)?;
            for (&j, &op) in bits.iter().zip(ops.iter())
            {
                self.set(i, j, op);
            }
            self.xor_sign(i, flip_sign);
        }

        self.normalize();

        Ok(())
    }

    /// Collapse a wave function after measurement.
    ///
    /// Collapse qbit `bit` to value `value` after a measurement hase been made.
    /// Calling this function is only done for measurements with random outcomes,
    /// and in this case there should be at least one row in the tableau with
    /// an X or Y operator at the `bit` position. The number `i` is the index
    /// of the last such row in the tableau.
    /// Use this function only in combination with `measure()`.
    pub fn collapse(&mut self, i: usize, bit: usize, value: bool)
    {
        for k in 0..i
        {
            if self.get_x(k, bit)
            {
                self.multiply_row(k, i);
            }
        }

        for j in 0..self.nr_bits
        {
            self.set(i, j, if j == bit { PauliOp::Z } else { PauliOp::I });
        }
        self.set_sign(i, value);

        self.normalize();
    }

    /// Measure a qubit
    ///
    /// Measure the qubit with index `bit` in this tableau. The outcome of the
    /// measurement is either zero or one, in which case a `Deterministic`
    /// result is returned, or an equal superposition of these states, in
    /// which case a `Random` result is returned.
    pub fn measure(&self, bit: usize) -> MeasurementInfo
    {
        if let Some(i) = (0..self.nr_bits).rev().filter(|&i| self.get_x(i, bit)).next()
        {
            MeasurementInfo::Random(i)
        }
        else
        {
            let i = (0..self.nr_bits).rev()
                .filter(|&i| self.get(i, bit) == PauliOp::Z)
                .next().unwrap();
            MeasurementInfo::Deterministic(self.get_sign(i))
        }
    }

    /// Reset a qubit
    ///
    /// Reset the qubit with index `bit` to zero.
    pub fn reset(&mut self, bit: usize)
    {
        if let Some(i) = (0..self.nr_bits).rev().filter(|&i| self.get_x(i, bit)).next()
        {
            // Bit is in superposition. Collapse it to zero.
            self.collapse(i, bit, false);
        }
        else
        {
            // Bit is either zero or one. Force it to zero.
            let i = (0..self.nr_bits).rev()
                .filter(|&i| self.get(i, bit) == PauliOp::Z)
                .next().unwrap();
            self.set_sign(i, false);
        }
    }
}


impl ::std::fmt::Display for StabilizerTableau
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        let n = self.nr_bits;
        for i in 0..n
        {
            write!(f, "{}", if self.get_sign(i) { '-' } else { '+' })?;
            for j in 0..n
            {
                write!(f, "{}", self.get(i, j))?;
            }
            if i < n-1
            {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests
{
    use super::{MeasurementInfo, StabilizerTableau};
    use crate::gates::{CX, CY, CZ, H, S, Sdg, V, Vdg, X, Y, Z};
    use crate::stabilizer::PauliOp;

    use ::std::fmt::Write;

    #[test]
    fn test_new()
    {
        let m = StabilizerTableau::new(3);
        assert_eq!(m.nr_bits, 3);
        assert_eq!(m.xz, vec![0x0000000000010101]);
        assert_eq!(m.signs, vec![0]);

        let m = StabilizerTableau::new(6);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0x0100040010004001, 0x0000000000000040]);
        assert_eq!(m.signs, vec![0]);
    }

    #[test]
    fn test_set()
    {
        let mut m = StabilizerTableau::new(6);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0x0100040010004001, 0x0000000000000040]);

        m.set(0, 0, PauliOp::X);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0x0100040010004002, 0x0000000000000040]);

        m.set(1, 2, PauliOp::Y);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0x0100040010034002, 0x0000000000000040]);

        m.set(0, 0, PauliOp::Z);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0x0100040010034001, 0x0000000000000040]);

        m.set(5, 1, PauliOp::Y);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0xc100040010034001, 0x0000000000000040]);

        m.set(5, 3, PauliOp::X);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0xc100040010034001, 0x0000000000000048]);

        m.set(5, 5, PauliOp::X);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0xc100040010034001, 0x0000000000000088]);
    }

    #[test]
    fn test_get()
    {
        let nr_bits = 15;
        let mut m = StabilizerTableau::new(nr_bits);
        for i in 0..nr_bits
        {
            for j in 0..nr_bits
            {
                assert_eq!(m.get(i, j), if i == j { PauliOp::Z } else { PauliOp:: I });
            }
        }

        m.set(1, 2, PauliOp::X);
        m.set(3, 14, PauliOp::Y);
        m.set(8, 7, PauliOp::Z);
        for i in 0..nr_bits
        {
            for j in 0..nr_bits
            {
                let op = m.get(i, j);
                match (i, j)
                {
                    (1, 2) => assert_eq!(op, PauliOp::X),
                    (3, 14) => assert_eq!(op, PauliOp::Y),
                    (8, 7) => assert_eq!(op, PauliOp::Z),
                    (a, b) if a == b => assert_eq!(op, PauliOp::Z),
                    _ => assert_eq!(op, PauliOp::I)
                }
            }
        }
    }

    #[test]
    fn test_display()
    {
        let mut m = StabilizerTableau::new(13);
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+ZIIIIIIIIIIII
+IZIIIIIIIIIII
+IIZIIIIIIIIII
+IIIZIIIIIIIII
+IIIIZIIIIIIII
+IIIIIZIIIIIII
+IIIIIIZIIIIII
+IIIIIIIZIIIII
+IIIIIIIIZIIII
+IIIIIIIIIZIII
+IIIIIIIIIIZII
+IIIIIIIIIIIZI
+IIIIIIIIIIIIZ"));

        m.set(1, 3, PauliOp::X);
        m.set(12, 11, PauliOp::Y);
        s.clear();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+ZIIIIIIIIIIII
+IZIXIIIIIIIII
+IIZIIIIIIIIII
+IIIZIIIIIIIII
+IIIIZIIIIIIII
+IIIIIZIIIIIII
+IIIIIIZIIIIII
+IIIIIIIZIIIII
+IIIIIIIIZIIII
+IIIIIIIIIZIII
+IIIIIIIIIIZII
+IIIIIIIIIIIZI
+IIIIIIIIIIIYZ"));
    }

    #[test]
    fn test_apply()
    {
        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[2]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IIX
+ZII
+IZI"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&X::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+ZII
-IZI
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&Y::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+ZII
-IZI
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&Z::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"-IXI
+ZII
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&S::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IYI
+ZII
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&Sdg::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"-IYI
+ZII
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&V::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"-IYI
+ZII
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&Vdg::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IYI
+ZII
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(m.apply_gate(&CX::new(), &[0, 2]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+XIX
+ZIZ
+IZI"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&CY::new(), &[0, 1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IXI
+ZII
+IIZ"));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&H::new(), &[2]), Ok(()));
        assert_eq!(m.apply_gate(&CZ::new(), &[1, 2]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IXZ
+IZX
+ZII"));
    }

    #[test]
    fn test_measure_deterministic()
    {
        let m = StabilizerTableau::new(3);
        assert_eq!(m.measure(0), MeasurementInfo::Deterministic(false));
        assert_eq!(m.measure(1), MeasurementInfo::Deterministic(false));
        assert_eq!(m.measure(2), MeasurementInfo::Deterministic(false));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&X::new(), &[1]), Ok(()));
        assert_eq!(m.measure(0), MeasurementInfo::Deterministic(false));
        assert_eq!(m.measure(1), MeasurementInfo::Deterministic(true));
        assert_eq!(m.measure(2), MeasurementInfo::Deterministic(false));

        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(m.apply_gate(&H::new(), &[2]), Ok(()));
        assert_eq!(m.measure(0), MeasurementInfo::Deterministic(true));
    }

    #[test]
    fn test_measure_random()
    {
        let mut m = StabilizerTableau::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.measure(1), MeasurementInfo::Random(0));

        let mut m = StabilizerTableau::new(2);
        assert_eq!(m.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(m.apply_gate(&CX::new(), &[0, 1]), Ok(()));
        assert_eq!(m.measure(0), MeasurementInfo::Random(0));
    }
}
