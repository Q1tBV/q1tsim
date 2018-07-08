extern crate rulinalg;

use rulinalg::matrix::BaseMatrix;

/// Type for a matrix of complex numbers.
///
/// Type CMatrix represents a matrix of complex numbers, where both the real
/// and the imaginary part are encoded as double precision floating point
/// numbers (`f64`).
pub struct CMatrix
{
    /// The real part of the matrix, `None` if the matrix is completely imaginary.
    real: Option<rulinalg::matrix::Matrix<f64>>,
    /// The imaginary part of the matrix, `None` if the matrix is real.
    imag: Option<rulinalg::matrix::Matrix<f64>>
}

impl CMatrix
{
    /// Return an identity matrix of size `n` × `n`.
    pub fn eye(n: usize) -> Self
    {
        CMatrix { real: Some(rulinalg::matrix::Matrix::identity(n)), imag: None }
    }

    /// Return a new real-valued matrix, with the coefficients from `mat`.
    pub fn new_real(mat: rulinalg::matrix::Matrix<f64>) -> Self
    {
        CMatrix { real: Some(mat), imag: None }
    }

    /// Return the number of rows in this matrix.
    pub fn rows(&self) -> usize
    {
        match self.real
        {
            Some(ref mat) => mat.rows(),
            None          => self.imag().unwrap().rows()
        }
    }

    /// Return the number of columns in this matrix.
    pub fn cols(&self) -> usize
    {
        match self.real
        {
            Some(ref mat) => mat.cols(),
            None          => self.imag().unwrap().cols()
        }
    }

    /// Return the real part of this matrix. If the matrix is purely imaginary,
    /// `None` is returned.
    pub fn real(&self) -> Option<&rulinalg::matrix::Matrix<f64>>
    {
        self.real.as_ref()
    }

    /// Return the imaginary part of this matrix. If the matrix is purely real,
    /// `None` is returned.
    pub fn imag(&self) -> Option<&rulinalg::matrix::Matrix<f64>>
    {
        self.imag.as_ref()
    }

    /// Return the real part of this matrix.
    pub fn real_expand(&self) -> ::std::borrow::Cow<rulinalg::matrix::Matrix<f64>>
    {
        match &self.real
        {
            &Some(ref mat) => ::std::borrow::Cow::Borrowed(mat),
            &None          => ::std::borrow::Cow::Owned(rulinalg::matrix::Matrix::zeros(self.rows(), self.cols()))
        }
    }

    /// Return the imaginary part of this matrix.
    pub fn imag_expand(&self) -> ::std::borrow::Cow<rulinalg::matrix::Matrix<f64>>
    {
        match &self.imag
        {
            &Some(ref mat) => ::std::borrow::Cow::Borrowed(mat),
            &None          => ::std::borrow::Cow::Owned(rulinalg::matrix::Matrix::zeros(self.rows(), self.cols()))
        }
    }
}
