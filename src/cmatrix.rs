extern crate rulinalg;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

type RMatrix = rulinalg::matrix::Matrix<f64>;

fn kron(a0: &RMatrix, a1: &RMatrix) -> RMatrix
{
    let (n0, m0, n1, m1) = (a0.rows(), a0.cols(), a1.rows(), a1.cols());
    let mut res = RMatrix::zeros(n0*n1, m0*m1);
    for (i0, r0) in a0.row_iter().enumerate()
    {
        for (j0, &x0) in r0.iter().enumerate()
        {
            res.sub_slice_mut([n1*i0, m1*j0], n1, m1).set_to(a1 * x0);
        }
    }
    res
}

fn abs_sq(a: &RMatrix) -> RMatrix
{
    let mut res = RMatrix::zeros(a.rows(), a.cols());
    for (x, &y) in res.iter_mut().zip(a.iter())
    {
        *x = y * y;
    }
    res
}

/// Type for a matrix of complex numbers.
///
/// Type CMatrix represents a matrix of complex numbers, where both the real
/// and the imaginary part are encoded as double precision floating point
/// numbers (`f64`).
pub enum CMatrix
{
    /// A real-valued matrix.
    Real(RMatrix),
    /// A fully imaginary matrix.
    Imag(RMatrix),
    /// A complex-valued matrix
    Complex(RMatrix, RMatrix)
}

impl CMatrix
{
    /// The identity matrix.
    ///
    /// Return an identity matrix of size `n` × `m`. If `n` &ne; `m`, the
    /// remaining rows or columns are filled with zeros.
    pub fn eye(n: usize, m: usize) -> Self
    {
        if n == m
        {
            CMatrix::Real(RMatrix::identity(n))
        }
        else
        {
            let mut res: RMatrix = RMatrix::zeros(n, m);
            for x in res.diag_iter_mut(rulinalg::matrix::DiagOffset::Main)
            {
                *x = 1.0;
            }
            CMatrix::Real(res)
        }
    }

    /// Return a new real-valued matrix, with the coefficients from `data`.
    pub fn new_real<U>(rows: usize, cols: usize, data: U) -> CMatrix
    where U: Into<Vec<f64>>
    {
        CMatrix::Real(RMatrix::new(rows, cols, data))
    }

    /// Return a new imaginary matrix, with the coefficients from `data`.
    pub fn new_imag<U>(rows: usize, cols: usize, data: U) -> CMatrix
    where U: Into<Vec<f64>>
    {
        CMatrix::Imag(RMatrix::new(rows, cols, data))
    }

    /// Return a new complex matrix, with the real parts of the coefficients
    /// from `rdata`, and the imaginary part from `idata`.
    pub fn new_complex<U>(rows: usize, cols: usize, rdata: U, idata: U) -> CMatrix
    where U: Into<Vec<f64>>
    {
        CMatrix::Complex(RMatrix::new(rows, cols, rdata), RMatrix::new(rows, cols, idata))
    }

    /// Return the number of rows in this matrix.
    pub fn rows(&self) -> usize
    {
        match *self
        {
            CMatrix::Real(ref mat)       => mat.rows(),
            CMatrix::Imag(ref mat)       => mat.rows(),
            CMatrix::Complex(ref mat, _) => mat.rows()
        }
    }

    /// Return the number of columns in this matrix.
    pub fn cols(&self) -> usize
    {
        match *self
        {
            CMatrix::Real(ref mat)       => mat.cols(),
            CMatrix::Imag(ref mat)       => mat.cols(),
            CMatrix::Complex(ref mat, _) => mat.cols()
        }
    }

    /// Return the real part of this matrix.
    pub fn real(&self) -> ::std::borrow::Cow<RMatrix>
    {
        match *self
        {
            CMatrix::Real(ref mat)       => ::std::borrow::Cow::Borrowed(mat),
            CMatrix::Imag(_)             => ::std::borrow::Cow::Owned(RMatrix::zeros(self.rows(), self.cols())),
            CMatrix::Complex(ref mat, _) => ::std::borrow::Cow::Borrowed(mat)
        }
    }

    /// Return the imaginary part of this matrix.
    pub fn imag(&self) -> ::std::borrow::Cow<RMatrix>
    {
        match *self
        {
            CMatrix::Real(_)             => ::std::borrow::Cow::Owned(RMatrix::zeros(self.rows(), self.cols())),
            CMatrix::Imag(ref mat)       => ::std::borrow::Cow::Borrowed(mat),
            CMatrix::Complex(_, ref mat) => ::std::borrow::Cow::Borrowed(mat)
        }
    }

    /// Squared absolute value of this matrix.
    ///
    /// Return the matrix with filled with the squared absolute values of the
    /// elements in this matrix.
    pub fn abs_sq(&self) -> RMatrix
    {
        match *self
        {
            CMatrix::Real(ref r)           => abs_sq(r),
            CMatrix::Imag(ref i)           => abs_sq(i),
            CMatrix::Complex(ref r, ref i) => abs_sq(r) + abs_sq(i)
        }
    }

    /// Compute the Kronecker product of this matrix and `m`.
    pub fn kron(&self, m: &Self) -> Self
    {
        match *self
        {
            CMatrix::Real(ref r0)            => {
                match *m
                {
                    CMatrix::Real(ref r1)            => { CMatrix::Real(kron(r0, r1)) },
                    CMatrix::Imag(ref i1)            => { CMatrix::Imag(kron(r0, i1)) },
                    CMatrix::Complex(ref r1, ref i1) => {
                        CMatrix::Complex(kron(r0, r1), kron(r0, i1))
                    }
                }
            },
            CMatrix::Imag(ref i0)            => {
                match m
                {
                    CMatrix::Real(ref r1)            => { CMatrix::Imag(kron(i0, r1)) },
                    CMatrix::Imag(ref i1)            => { CMatrix::Real(-kron(i0, i1)) },
                    CMatrix::Complex(ref r1, ref i1) => {
                        CMatrix::Complex(-kron(i0, i1), kron(i0, r1))
                    }
                }
            },
            CMatrix::Complex(ref r0, ref i0) => {
                match m
                {
                    CMatrix::Real(ref r1)            => {
                        CMatrix::Complex(kron(r0, r1), kron(i0, r1))
                    },
                    CMatrix::Imag(ref i1)            => {
                        CMatrix::Complex(-kron(i0, i1), kron(r0, i1))
                    },
                    CMatrix::Complex(ref r1, ref i1) => {
                        CMatrix::Complex(kron(r0, r1) - kron(i0, i1), kron(r0, i1) + kron(i0, r1))
                    }
                }
            },
        }
    }

    /// Set all elements in the range [`i`..`i`+`ni`, `j`..`j`+`nj`] to zero.
    pub fn clear(&mut self, i: usize, j: usize, ni: usize, nj: usize)
    {
        match *self
        {
            CMatrix::Real(ref mut rl)                => {
                rl.sub_slice_mut([i, j], ni, nj).set_to(RMatrix::zeros(ni, nj));
            },
            CMatrix::Imag(ref mut im)                => {
                im.sub_slice_mut([i, j], ni, nj).set_to(RMatrix::zeros(ni, nj));
            }
            CMatrix::Complex(ref mut rl, ref mut im) => {
                rl.sub_slice_mut([i, j], ni, nj).set_to(RMatrix::zeros(ni, nj));
                im.sub_slice_mut([i, j], ni, nj).set_to(RMatrix::zeros(ni, nj));
            }
        }
    }
}

impl ::std::ops::MulAssign<f64> for CMatrix
{
    fn mul_assign(&mut self, x: f64)
    {
        match *self
        {
            CMatrix::Real(ref mut r)               => *r *= x,
            CMatrix::Imag(ref mut i)               => *i *= x,
            CMatrix::Complex(ref mut r, ref mut i) => { *r *= x; *i *= x }
        }
    }
}
