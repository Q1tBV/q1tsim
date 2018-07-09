extern crate rulinalg;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};
use rulinalg::matrix::Matrix;

fn kron(a0: &Matrix<f64>, a1: &Matrix<f64>) -> Matrix<f64>
{
    let (n0, m0, n1, m1) = (a0.rows(), a0.cols(), a1.rows(), a1.cols());
    let mut res = Matrix::zeros(n0*n1, m0*m1);
    for (i0, r0) in a0.row_iter().enumerate()
    {
        for (j0, &x0) in r0.iter().enumerate()
        {
            res.sub_slice_mut([n1*i0, m1*j0], n1, m1).set_to(a1 * x0);
        }
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
    Real(Matrix<f64>),
    /// A fully imaginary matrix.
    Imag(Matrix<f64>),
    /// A complex-valued matrix
    Complex(Matrix<f64>, Matrix<f64>)
}

impl CMatrix
{
    /// Return an identity matrix of size `n` × `n`.
    pub fn eye(n: usize) -> Self
    {
        CMatrix::Real(Matrix::identity(n))
    }

    /// Return a new real-valued matrix, with the coefficients from `mat`.
    pub fn new_real(mat: Matrix<f64>) -> Self
    {
        CMatrix::Real(mat)
    }

    /// Return the number of rows in this matrix.
    pub fn rows(&self) -> usize
    {
        match self
        {
            &CMatrix::Real(ref mat)       => mat.rows(),
            &CMatrix::Imag(ref mat)       => mat.rows(),
            &CMatrix::Complex(ref mat, _) => mat.rows()
        }
    }

    /// Return the number of columns in this matrix.
    pub fn cols(&self) -> usize
    {
        match self
        {
            &CMatrix::Real(ref mat)       => mat.cols(),
            &CMatrix::Imag(ref mat)       => mat.cols(),
            &CMatrix::Complex(ref mat, _) => mat.cols()
        }
    }

    /// Return the real part of this matrix.
    pub fn real(&self) -> ::std::borrow::Cow<Matrix<f64>>
    {
        match self
        {
            &CMatrix::Real(ref mat)       => ::std::borrow::Cow::Borrowed(mat),
            &CMatrix::Imag(_)             => ::std::borrow::Cow::Owned(Matrix::zeros(self.rows(), self.cols())),
            &CMatrix::Complex(ref mat, _) => ::std::borrow::Cow::Borrowed(mat)
        }
    }

    /// Return the imaginary part of this matrix. If the matrix is purely real,
    /// `None` is returned.
    pub fn imag(&self) -> ::std::borrow::Cow<Matrix<f64>>
    {
        match self
        {
            &CMatrix::Real(_)             => ::std::borrow::Cow::Owned(Matrix::zeros(self.rows(), self.cols())),
            &CMatrix::Imag(ref mat)       => ::std::borrow::Cow::Borrowed(mat),
            &CMatrix::Complex(_, ref mat) => ::std::borrow::Cow::Borrowed(mat)
        }
    }

    /// Compute the Kronecker product of this matrix and `m`.
    pub fn kron(&self, m: &Self) -> Self
    {
        match self
        {
            &CMatrix::Real(ref r0)            => {
                match m
                {
                    &CMatrix::Real(ref r1)            => { CMatrix::Real(kron(r0, r1)) },
                    &CMatrix::Imag(ref i1)            => { CMatrix::Imag(kron(r0, i1)) },
                    &CMatrix::Complex(ref r1, ref i1) => {
                        CMatrix::Complex(kron(r0, r1), kron(r0, i1))
                    }
                }
            },
            &CMatrix::Imag(ref i0)            => {
                match m
                {
                    &CMatrix::Real(ref r1)            => { CMatrix::Imag(kron(i0, r1)) },
                    &CMatrix::Imag(ref i1)            => { CMatrix::Real(-kron(i0, i1)) },
                    &CMatrix::Complex(ref r1, ref i1) => {
                        CMatrix::Complex(-kron(i0, i1), kron(i0, r1))
                    }
                }
            },
            &CMatrix::Complex(ref r0, ref i0) => {
                match m
                {
                    &CMatrix::Real(ref r1)            => {
                        CMatrix::Complex(kron(r0, r1), kron(i0, r1))
                    },
                    &CMatrix::Imag(ref i1)            => {
                        CMatrix::Complex(-kron(i0, i1), kron(r0, i1))
                    },
                    &CMatrix::Complex(ref r1, ref i1) => {
                        CMatrix::Complex(kron(r0, r1) - kron(i0, i1), kron(r0, i1) + kron(i0, r1))
                    }
                }
            },
        }
    }
}
