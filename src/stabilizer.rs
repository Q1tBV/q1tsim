/// Enumeration for the Pauli spin operations
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum PauliOp
{
    /// Identity
    I,
    /// Pauli Z
    Z,
    /// Pauli X
    X,
    /// Pauli Y
    Y
}

impl PauliOp
{
    fn from_bits(bits: u64) -> Self
    {
        match bits & 0x03
        {
            0 => PauliOp::I,
            1 => PauliOp::Z,
            2 => PauliOp::X,
            3 => PauliOp::Y,
            _ => unreachable!()
        }
    }

    fn to_bits(self) -> u64
    {
        match self
        {
            PauliOp::I => 0,
            PauliOp::Z => 1,
            PauliOp::X => 2,
            PauliOp::Y => 3
        }
    }
}

impl ::std::fmt::Display for PauliOp
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        match *self
        {
            PauliOp::I => write!(f, "I"),
            PauliOp::Z => write!(f, "Z"),
            PauliOp::X => write!(f, "X"),
            PauliOp::Y => write!(f, "Y"),
        }
    }
}

/// Structure describing a single stabilizer state
struct StabilizerMatrix
{
    /// The number of qubits in the state
    nr_bits: usize,
    /// Generators of the stabilizer group for this state
    xz: Vec<u64>,
    /// Signs of the generators
    signs: Vec<u64>
}

impl StabilizerMatrix
{
    fn new(nr_bits: usize) -> Self
    {
        let xz_size = (2*nr_bits*nr_bits + 0x3f) >> 6;
        let signs_size = (nr_bits + 0x3f) >> 6;
        let mut res = StabilizerMatrix {
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

    fn measure_random<R: rand::Rng>(&mut self, i: usize, bit: usize, rng: &mut R) -> u8
    {
        let res: bool = rng.gen();

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
        self.set_sign(i, res);

        self.normalize();

        res as u8
    }

    fn measure_deterministic(&mut self, bit: usize) -> u8
    {
        let i = (0..self.nr_bits).rev()
            .filter(|&i| self.get(i, bit) == PauliOp::Z)
            .next().unwrap();
        return self.get_sign(i) as u8;
    }

    pub fn measure<R: rand::Rng>(&mut self, bit: usize, rng: &mut R) -> u8
    {
        if let Some(i) = (0..self.nr_bits).rev().filter(|&i| self.get_x(i, bit)).next()
        {
            self.measure_random(i, bit, rng)
        }
        else
        {
            self.measure_deterministic(bit)
        }
    }
}


impl ::std::fmt::Display for StabilizerMatrix
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
    use super::{PauliOp, StabilizerMatrix};
    use crate::gates::{CX, CY, CZ, H, S, Sdg, V, Vdg, X, Y, Z};

    use ::std::fmt::Write;

    #[test]
    fn test_new()
    {
        let m = StabilizerMatrix::new(3);
        assert_eq!(m.nr_bits, 3);
        assert_eq!(m.xz, vec![0x0000000000010101]);
        assert_eq!(m.signs, vec![0]);

        let m = StabilizerMatrix::new(6);
        assert_eq!(m.nr_bits, 6);
        assert_eq!(m.xz, vec![0x0100040010004001, 0x0000000000000040]);
        assert_eq!(m.signs, vec![0]);
    }

    #[test]
    fn test_set()
    {
        let mut m = StabilizerMatrix::new(6);
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
        let mut m = StabilizerMatrix::new(nr_bits);
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
        let mut m = StabilizerMatrix::new(13);
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
        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[2]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IIX
+ZII
+IZI"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&X::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+ZII
-IZI
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&Y::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+ZII
-IZI
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&Z::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"-IXI
+ZII
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&S::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IYI
+ZII
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&Sdg::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"-IYI
+ZII
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&V::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"-IYI
+ZII
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&Vdg::new(), &[1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IYI
+ZII
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[0]), Ok(()));
        assert_eq!(m.apply_gate(&CX::new(), &[0, 2]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+XIX
+ZIZ
+IZI"));

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&CY::new(), &[0, 1]), Ok(()));
        let mut s = String::new();
        assert_eq!(write!(s, "{}", m), Ok(()));
        assert_eq!(s, String::from(
r"+IXI
+ZII
+IIZ"));

        let mut m = StabilizerMatrix::new(3);
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
        let mut rng = rand::thread_rng();
        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.measure(0, &mut rng), 0);
        assert_eq!(m.measure(1, &mut rng), 0);
        assert_eq!(m.measure(2, &mut rng), 0);

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&X::new(), &[1]), Ok(()));
        assert_eq!(m.measure(0, &mut rng), 0);
        assert_eq!(m.measure(1, &mut rng), 1);
        assert_eq!(m.measure(2, &mut rng), 0);

        let mut m = StabilizerMatrix::new(3);
        assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
        assert_eq!(m.apply_gate(&X::new(), &[0]), Ok(()));
        assert_eq!(m.apply_gate(&H::new(), &[2]), Ok(()));
        assert_eq!(m.measure(0, &mut rng), 1);
    }

    #[test]
    fn test_measure_random()
    {
        let nr_runs = 1024;
        let tol = 1.0e-5;
        let mut rng = rand::thread_rng();

        let mut n1 = 0;
        let mut s = String::new();
        for _ in 0..nr_runs
        {
            let mut m = StabilizerMatrix::new(3);
            assert_eq!(m.apply_gate(&H::new(), &[1]), Ok(()));
            let res = m.measure(1, &mut rng);
            n1 += res as usize;

            s.clear();
            assert_eq!(write!(s, "{}", m), Ok(()));
            if res == 0
            {
                assert_eq!(s, "+ZII\n+IZI\n+IIZ");
            }
            else
            {
                assert_eq!(s, "+ZII\n-IZI\n+IIZ");
            }
        }
        assert!(crate::stats::measurement_ok(n1, nr_runs, 0.5, tol));


        let mut n0 = 0;
        let mut s = String::new();
        for _ in 0..nr_runs
        {
            let mut m = StabilizerMatrix::new(2);
            assert_eq!(m.apply_gate(&H::new(), &[0]), Ok(()));
            assert_eq!(m.apply_gate(&CX::new(), &[0, 1]), Ok(()));
            let res = m.measure(0, &mut rng);
            n0 += res as usize;

            s.clear();
            assert_eq!(write!(s, "{}", m), Ok(()));
            if res == 0
            {
                assert_eq!(s, "+ZI\n+IZ");
            }
            else
            {
                assert_eq!(s, "-ZI\n-IZ");
            }
        }
        assert!(crate::stats::measurement_ok(n1, nr_runs, 0.5, tol));
    }
}
