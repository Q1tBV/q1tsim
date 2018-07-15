use cmatrix::CMatrix;
use rulinalg;

use rulinalg::matrix::{BaseMatrix, BaseMatrixMut};

pub struct CMatrixSlice<'a>
{
    matrix: &'a CMatrix,
    offset: [usize; 2],
    rows: usize,
    cols: usize
}

impl<'a> CMatrixSlice<'a>
{
    pub fn new(matrix: &'a CMatrix, offset: [usize; 2], rows: usize, cols: usize) -> Self
    {
        CMatrixSlice { matrix: matrix, offset: offset, rows: rows, cols: cols }
    }

    pub fn into_matrix(&self) -> CMatrix
    {
        match *self.matrix
        {
            CMatrix::Real(ref r)           => {
                CMatrix::Real(r.sub_slice(self.offset, self.rows, self.cols).into_matrix())
            },
            CMatrix::Imag(ref i)           => {
                CMatrix::Imag(i.sub_slice(self.offset, self.rows, self.cols).into_matrix())
            },
            CMatrix::Complex(ref r, ref i) => {
                CMatrix::Complex(
                    r.sub_slice(self.offset, self.rows, self.cols).into_matrix(),
                    i.sub_slice(self.offset, self.rows, self.cols).into_matrix()
                )
            },
        }
    }
}

pub struct CMatrixSliceMut<'a>
{
    matrix: &'a mut CMatrix,
    offset: [usize; 2],
    rows: usize,
    cols: usize
}

impl<'a> CMatrixSliceMut<'a>
{
    pub fn new(matrix: &'a mut CMatrix, offset: [usize; 2], rows: usize, cols: usize) -> Self
    {
        CMatrixSliceMut { matrix: matrix, offset: offset, rows: rows, cols: cols }
    }

    pub fn set_to(&mut self, m: CMatrix)
    {
        if (self.matrix.is_real() && !m.is_real()) || (self.matrix.is_imaginary() && !m.is_imaginary())
        {
            self.matrix.expand_to_complex();
        }

        match *self.matrix
        {
            CMatrix::Real(ref mut r)               => {
                match m
                {
                    CMatrix::Real(mr) => {
                        r.sub_slice_mut(self.offset, self.rows, self.cols).set_to(mr);
                    },
                    _                 => {
                        panic!("Trying to assign non-real matrix to real slice.");
                    }
                }
            },
            CMatrix::Imag(ref mut i)               => {
                match m
                {
                    CMatrix::Imag(mi) => {
                        i.sub_slice_mut(self.offset, self.rows, self.cols).set_to(mi);
                    },
                    _                 => {
                        panic!("Trying to assign non-imaginary matrix to imaginary slice.");
                    }
                }
            },
            CMatrix::Complex(ref mut r, ref mut i) => {
                match m
                {
                    CMatrix::Real(mr) => {
                        r.sub_slice_mut(self.offset, self.rows, self.cols).set_to(mr);
                        i.sub_slice_mut(self.offset, self.rows, self.cols).set_to(rulinalg::matrix::Matrix::zeros(self.rows, self.cols));
                    },
                    CMatrix::Imag(mi) => {
                        r.sub_slice_mut(self.offset, self.rows, self.cols).set_to(rulinalg::matrix::Matrix::zeros(self.rows, self.cols));
                        i.sub_slice_mut(self.offset, self.rows, self.cols).set_to(mi);
                    },
                    CMatrix::Complex(mr, mi) => {
                        r.sub_slice_mut(self.offset, self.rows, self.cols).set_to(mr);
                        i.sub_slice_mut(self.offset, self.rows, self.cols).set_to(mi);
                    }
                }
            },
        }
    }
}

impl<'a> ::std::ops::MulAssign<f64> for CMatrixSliceMut<'a>
{
    fn mul_assign(&mut self, x: f64)
    {
        match self.matrix
        {
            CMatrix::Real(ref mut r)               => {
                let mut slice = r.sub_slice_mut(self.offset, self.rows, self.cols);
                slice *= x;
            },
            CMatrix::Imag(ref mut i)               => {
                let mut slice = i.sub_slice_mut(self.offset, self.rows, self.cols);
                slice *= x;
            },
            CMatrix::Complex(ref mut r, ref mut i) => {
                let mut rslice = r.sub_slice_mut(self.offset, self.rows, self.cols);
                let mut islice = i.sub_slice_mut(self.offset, self.rows, self.cols);
                rslice *= x;
                islice *= x;
            },
        }
    }
}
