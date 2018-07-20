extern crate num_complex;
pub extern crate rulinalg;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

pub const COMPLEX_ZERO:   num_complex::Complex64 = num_complex::Complex { re: 0.0, im: 0.0 };
pub const COMPLEX_ONE:    num_complex::Complex64 = num_complex::Complex { re: 1.0, im: 0.0 };
pub const COMPLEX_HSQRT2: num_complex::Complex64 = num_complex::Complex { re: ::std::f64::consts::FRAC_1_SQRT_2, im: 0.0 };
pub const COMPLEX_I:      num_complex::Complex64 = num_complex::Complex { re: 0.0, im: 1.0 };

pub struct CMatrix
{
    matrix: rulinalg::matrix::Matrix<num_complex::Complex64>
}

impl CMatrix
{
    /// Create a new complex matrix of size `rows` × `cols`, with elements from
    /// `data`. The elements should be stored in row-major order.
    pub fn new<U>(rows: usize, cols: usize, data: U) -> CMatrix
    where U: Into<Vec<num_complex::Complex64>>
    {
        CMatrix { matrix: rulinalg::matrix::Matrix::new(rows, cols, data) }
    }

    /// Create an identity matrix of size `rows` × `cols`. If `rows` ≠ `cols`,
    /// the remaining rows or columns are filled with zeros.
    pub fn eye(rows: usize, cols: usize) -> Self
    {
        let mut matrix;
        if rows == cols
        {
            matrix = rulinalg::matrix::Matrix::identity(rows);
        }
        else
        {
            matrix = rulinalg::matrix::Matrix::zeros(rows, cols);
            for x in matrix.diag_iter_mut(rulinalg::matrix::DiagOffset::Main)
            {
                *x = COMPLEX_ONE;
            }
        }
        CMatrix { matrix: matrix }
    }

    /// Compute the Kronecker product `self` ⊗ `m`.
    pub fn kron(&self, m: &Self) -> Self
    {
        let (n0, m0, n1, m1) = (self.rows(), self.cols(), m.rows(), m.cols());
        let mut data = rulinalg::matrix::Matrix::zeros(n0*n1, m0*m1);
        for i in 0..n0
        {
            for j in 0..m0
            {
                data.sub_slice_mut([i*n1, j*m1], n1, m1).set_to(&m.matrix * self[[i,j]]);
            }
        }
        CMatrix { matrix: data }
    }

    /// Return the real part of this matrix
    pub fn real(&self) -> rulinalg::matrix::Matrix<f64>
    {
        rulinalg::matrix::Matrix::new(self.rows(), self.cols(), self.iter()
            .map(|c| c.re).collect::<Vec<f64>>())
    }

    /// Return the imaginary part of this matrix
    pub fn imag(&self) -> rulinalg::matrix::Matrix<f64>
    {
        rulinalg::matrix::Matrix::new(self.rows(), self.cols(), self.iter()
            .map(|c| c.im).collect::<Vec<f64>>())
    }
}

impl ::std::ops::Deref for CMatrix
{
    type Target = rulinalg::matrix::Matrix<num_complex::Complex<f64>>;

    fn deref(&self) -> &Self::Target
    {
        &self.matrix
    }
}

impl ::std::ops::DerefMut for CMatrix
{
    fn deref_mut(&mut self) -> &mut Self::Target
    {
        &mut self.matrix
    }
}

impl ::std::fmt::Display for CMatrix
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        write!(f, "{}", self.matrix)
    }
}

impl ::std::ops::MulAssign<f64> for CMatrix
{
    fn mul_assign(&mut self, x: f64)
    {
        // XXX TODO: Check which is better
        self.matrix *= num_complex::Complex::new(x, 0.0);
//         self.matrix.apply(&|c| c * x);
    }
}

#[cfg(test)]
macro_rules! assert_complex_matrix_eq
{
    ($x:expr, $y:expr) => {
        let rx = $crate::rulinalg::matrix::Matrix::new($x.rows(), $x.cols(), $x.iter()
            .map(|c| c.re).collect::<Vec<f64>>());
        let ry = $crate::rulinalg::matrix::Matrix::new($y.rows(), $y.cols(), $y.iter()
            .map(|c| c.re).collect::<Vec<f64>>());
        assert_matrix_eq!(rx, ry, comp=float);

        let ix = $crate::rulinalg::matrix::Matrix::new($x.rows(), $x.cols(), $x.iter()
            .map(|c| c.im).collect::<Vec<f64>>());
        let iy = $crate::rulinalg::matrix::Matrix::new($y.rows(), $y.cols(), $y.iter()
            .map(|c| c.im).collect::<Vec<f64>>());
        assert_matrix_eq!(ix, iy, comp=float);
    }
}
