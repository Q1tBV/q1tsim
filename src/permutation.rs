extern crate ndarray;

use cmatrix;

pub struct Permutation
{
    idxs: Vec<usize>
}

impl Permutation
{
    pub fn new(idxs: Vec<usize>) -> Self
    {
        Permutation { idxs: idxs }
    }

    pub fn size(&self) -> usize
    {
        self.idxs.len()
    }

    pub fn inverse(&self) -> Self
    {
        let mut invidxs = vec![0; self.size()];
        for (i, &pi) in self.idxs.iter().enumerate()
        {
            invidxs[pi] = i;
        }
        Permutation::new(invidxs)
    }

    pub fn apply<F>(&self, mut swap: F)
    where F: FnMut(usize, usize)
    {
        let mut perm = self.idxs.clone();
        for i in 0..self.size()
        {
            let mut pi = perm[i];
            while pi != i
            {
                let next = perm[pi];
                swap(i, pi);
                perm[pi] = pi;
                pi = next;
            }
        }
    }

    pub fn apply_vec(&self, v: &mut cmatrix::CVector)
    {
        self.apply(|i, j| v.swap(i, j));
    }

    pub fn transform(&self, a: &cmatrix::CMatrix) -> cmatrix::CMatrix
    {
        let invidxs = self.inverse().idxs;
        a.select(ndarray::Axis(0), &invidxs).select(ndarray::Axis(1), &invidxs)
    }
}
