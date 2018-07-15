extern crate rulinalg;

use cmatrix;
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

/// Type for a single complex number
#[derive(Clone, Copy)]
pub enum CValue
{
    /// A real number.
    Real(f64),
    /// An imaginary number.
    Imag(f64),
    /// A complex number.
    Complex(f64, f64)
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

    /// Return whether this matrix is real-valued
    pub fn is_real(&self) -> bool
    {
        match *self
        {
            CMatrix::Real(_) => true,
            _                => false
        }
    }

    /// Return whether this matrix is fully imaginary
    pub fn is_imaginary(&self) -> bool
    {
        match *self
        {
            CMatrix::Imag(_) => true,
            _                => false
        }
    }

    /// Return the element at position `idx`.
    pub fn at(&self, idx: [usize; 2]) -> CValue
    {
        match *self
        {
            CMatrix::Real(ref r)           => CValue::Real(r[idx]),
            CMatrix::Imag(ref i)           => CValue::Imag(i[idx]),
            CMatrix::Complex(ref r, ref i) => CValue::Complex(r[idx], i[idx])
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

    /// Take an immutable slice of size `rows` × `cols` from this matrix,
    /// starting at position `offset`.
    pub fn sub_slice(&self, offset: [usize; 2], rows: usize, cols: usize)
        -> cmatrix::CMatrixSlice
    {
        cmatrix::CMatrixSlice::new(self, offset, rows, cols)
    }

    /// Take an mutable slice of size `rows` × `cols` from this matrix,
    /// starting at position `offset`.
    pub fn sub_slice_mut(&mut self, offset: [usize; 2], rows: usize, cols: usize)
        -> cmatrix::CMatrixSliceMut
    {
        cmatrix::CMatrixSliceMut::new(self, offset, rows, cols)
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

    /// Expand a real valued or imaginary matrix into a fully complex one,
    /// with zeros for the previously omitted part
    pub fn expand_to_complex(&mut self)
    {
        let (n, m) = (self.rows(), self.cols());
        let r;
        let i;
        match *self
        {
            CMatrix::Real(ref mut cr)                => {
                r = ::std::mem::replace(cr, RMatrix::zeros(0, 0));
                i = RMatrix::zeros(n, m);
            },
            CMatrix::Imag(ref mut ci)                => {
                r = RMatrix::zeros(n, m);
                i = ::std::mem::replace(ci, RMatrix::zeros(0, 0));
            },
            CMatrix::Complex(ref mut cr, ref mut ci) => {
                r = ::std::mem::replace(cr, RMatrix::zeros(0, 0));
                i = ::std::mem::replace(ci, RMatrix::zeros(0, 0));
            }
        }
        *self = CMatrix::Complex(r, i);
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

impl ::std::ops::Add<CMatrix> for CMatrix
{
    type Output = CMatrix;

    fn add(self, m: CMatrix) -> Self::Output
    {
        &self + &m
    }
}

impl<'a, 'b> ::std::ops::Add<&'b CMatrix> for &'a CMatrix
{
    type Output = CMatrix;

    fn add(self, m: &CMatrix) -> Self::Output
    {
        match *self
        {
            CMatrix::Real(ref r)           => {
                match *m
                {
                    CMatrix::Real(ref mr)            => CMatrix::Real(r + mr),
                    CMatrix::Imag(ref mi)            => CMatrix::Complex(r.clone(), mi.clone()),
                    CMatrix::Complex(ref mr, ref mi) => CMatrix::Complex(r + mr, mi.clone())
                }
            },
            CMatrix::Imag(ref i)           => {
                match *m
                {
                    CMatrix::Real(ref mr)            => CMatrix::Complex(mr.clone(), i.clone()),
                    CMatrix::Imag(ref mi)            => CMatrix::Imag(i + mi),
                    CMatrix::Complex(ref mr, ref mi) => CMatrix::Complex(mr.clone(), i + mi)
                }
            },
            CMatrix::Complex(ref r, ref i) => {
                match *m
                {
                    CMatrix::Real(ref mr)            => CMatrix::Complex(r + mr, i.clone()),
                    CMatrix::Imag(ref mi)            => CMatrix::Complex(r.clone(), i + mi),
                    CMatrix::Complex(ref mr, ref mi) => CMatrix::Complex(r + mr, i + mi)
                }
            }
        }
    }
}

impl ::std::ops::Sub<CMatrix> for CMatrix
{
    type Output = CMatrix;

    fn sub(self, m: CMatrix) -> Self::Output
    {
        &self - &m
    }
}

impl<'a, 'b> ::std::ops::Sub<&'b CMatrix> for &'a CMatrix
{
    type Output = CMatrix;

    fn sub(self, m: &CMatrix) -> Self::Output
    {
        match *self
        {
            CMatrix::Real(ref r)           => {
                match *m
                {
                    CMatrix::Real(ref mr)            => CMatrix::Real(r - mr),
                    CMatrix::Imag(ref mi)            => CMatrix::Complex(r.clone(), -mi),
                    CMatrix::Complex(ref mr, ref mi) => CMatrix::Complex(r - mr, -mi)
                }
            },
            CMatrix::Imag(ref i)           => {
                match *m
                {
                    CMatrix::Real(ref mr)            => CMatrix::Complex(-mr, i.clone()),
                    CMatrix::Imag(ref mi)            => CMatrix::Imag(i - mi),
                    CMatrix::Complex(ref mr, ref mi) => CMatrix::Complex(-mr, i - mi)
                }
            },
            CMatrix::Complex(ref r, ref i) => {
                match *m
                {
                    CMatrix::Real(ref mr)            => CMatrix::Complex(r - mr, i.clone()),
                    CMatrix::Imag(ref mi)            => CMatrix::Complex(r.clone(), i - mi),
                    CMatrix::Complex(ref mr, ref mi) => CMatrix::Complex(r - mr, i - mi)
                }
            }
        }
    }
}

impl<'a> ::std::ops::Mul<CValue> for &'a CMatrix
{
    type Output = CMatrix;

    fn mul(self, x: CValue) -> Self::Output
    {
        match *self
        {
            CMatrix::Real(ref mr)            => {
                match x
                {
                    CValue::Real(ref xr)            => CMatrix::Real(mr*xr),
                    CValue::Imag(ref xi)            => CMatrix::Imag(mr*xi),
                    CValue::Complex(ref xr, ref xi) => CMatrix::Complex(mr*xr, mr*xi)
                }
            },
            CMatrix::Imag(ref mi)            => {
                match x
                {
                    CValue::Real(ref xr)            => CMatrix::Imag(mi*xr),
                    CValue::Imag(ref xi)            => CMatrix::Real(mi*(-xi)),
                    CValue::Complex(ref xr, ref xi) => CMatrix::Complex(mi*(-xi), mi*xr)
                }
            },
            CMatrix::Complex(ref mr, ref mi) => {
                match x
                {
                    CValue::Real(xr)        => CMatrix::Complex(mr*xr, mi*xr),
                    CValue::Imag(xi)        => CMatrix::Complex(mi*(-xi), mr*xi),
                    CValue::Complex(xr, xi) => CMatrix::Complex(mr*xr - mi*xi, mr*xi + mi*xr)
                }
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
