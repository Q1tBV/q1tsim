extern crate num_complex;
extern crate ndarray;

pub const COMPLEX_ZERO:   num_complex::Complex64 = num_complex::Complex { re: 0.0, im: 0.0 };
pub const COMPLEX_ONE:    num_complex::Complex64 = num_complex::Complex { re: 1.0, im: 0.0 };
pub const COMPLEX_HSQRT2: num_complex::Complex64 = num_complex::Complex { re: ::std::f64::consts::FRAC_1_SQRT_2, im: 0.0 };
pub const COMPLEX_I:      num_complex::Complex64 = num_complex::Complex { re: 0.0, im: 1.0 };

pub type CVector = ndarray::Array1<num_complex::Complex64>;
pub type CMatrix = ndarray::Array2<num_complex::Complex64>;
pub type CVecSliceMut<'a> = ndarray::ArrayViewMut1<'a, num_complex::Complex64>;
pub type CSliceMut<'a> = ndarray::ArrayViewMut2<'a, num_complex::Complex64>;

/// Compute the Kronecker product `v0` ⊗ `v1`.
pub fn kron_vec(v0: &CVector, v1: &CVector) -> CVector
{
    let (n0, n1) = (v0.len(), v1.len());
    let mut res = CVector::zeros(n0*n1);
    for i in 0..n0
    {
        res.slice_mut(s![i*n1..(i+1)*n1]).assign(&(v1 * v0[i]));
    }
    res
}

/// Compute the Kronecker product `m0` ⊗ `m1`.
pub fn kron_mat(a0: &CMatrix, a1: &CMatrix) -> CMatrix
{
    let (n0, m0, n1, m1) = (a0.rows(), a0.cols(), a1.rows(), a1.cols());
    let mut res = CMatrix::zeros((n0*n1, m0*m1));
    for i in 0..n0
    {
        for j in 0..m0
        {
            res.slice_mut(s![i*n1..(i+1)*n1, j*m1..(j+1)*m1]).assign(&(a1 * a0[[i, j]]));
        }
    }
    res
}

#[cfg(test)]
macro_rules! assert_complex_matrix_eq
{
    ($x:expr, $y:expr) => {
        {
            let rx = $crate::rulinalg::matrix::Matrix::new($x.rows(), $x.cols(),
                $x.iter().map(|c| c.re).collect::<Vec<f64>>());
            let ry = $crate::rulinalg::matrix::Matrix::new($y.rows(), $y.cols(),
                $y.iter().map(|c| c.re).collect::<Vec<f64>>());
            assert_matrix_eq!(rx, ry, comp=float);

            let ix = $crate::rulinalg::matrix::Matrix::new($x.rows(), $x.cols(),
                $x.iter().map(|c| c.im).collect::<Vec<f64>>());
            let iy = $crate::rulinalg::matrix::Matrix::new($y.rows(), $y.cols(),
                $y.iter().map(|c| c.im).collect::<Vec<f64>>());
            assert_matrix_eq!(ix, iy, comp=float);
        }
    };
    ($x:expr, $y:expr, eps=$t:expr) => {
        {
            let rx = $crate::rulinalg::matrix::Matrix::new($x.rows(), $x.cols(),
                $x.iter().map(|c| c.re).collect::<Vec<f64>>());
            let ry = $crate::rulinalg::matrix::Matrix::new($y.rows(), $y.cols(),
                $y.iter().map(|c| c.re).collect::<Vec<f64>>());
            assert_matrix_eq!(rx, ry, comp=float, eps=$t);

            let ix = $crate::rulinalg::matrix::Matrix::new($x.rows(), $x.cols(),
                $x.iter().map(|c| c.im).collect::<Vec<f64>>());
            let iy = $crate::rulinalg::matrix::Matrix::new($y.rows(), $y.cols(),
                $y.iter().map(|c| c.im).collect::<Vec<f64>>());
            assert_matrix_eq!(ix, iy, comp=float, eps=$t);
        }
    }
}
