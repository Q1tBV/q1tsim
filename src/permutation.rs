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


extern crate ndarray;
extern crate num_traits;

/// Struct for permutations
///
/// Struct Permutation is used to represents permutations. It can be used to
/// shuffle the elements in a vector, or rows and columns in a matrix.
pub struct Permutation
{
    /// The permuted indexes.
    idxs: Vec<usize>
}

impl Permutation
{
    /// Create a new permutation.
    ///
    /// Create a new permutation from `idxs`. If `n` is the length of `idxs`,
    /// the array should contain all elements in the range `[0, n-1]` (inclusive).
    /// When `idxs[i] = j`, application of the permutation on a vector means that
    /// the the `j`th element of the original vector is taken to position `i` in
    /// the permuted vector.
    pub fn new(idxs: Vec<usize>) -> crate::error::Result<Self>
    {
        if idxs.is_empty()
        {
            return Err(crate::error::Error::EmptyPermutation);
        }

        let n = idxs.len();
        let max = *idxs.iter().max().unwrap();
        if max >= n
        {
            return Err(crate::error::Error::InvalidPermutationElement(max, n));
        }

        let mut seen = vec![false; n];
        for &pi in idxs.iter()
        {
            if seen[pi]
            {
                return Err(crate::error::Error::DoublePermutationElement(pi));
            }
            seen[pi] = true;
        }

        Ok(Permutation { idxs: idxs })
    }

    /// Return the permutation index vector
    pub fn indices(&self) -> &[usize]
    {
        &self.idxs
    }

    /// Return the size of the permutation
    pub fn size(&self) -> usize
    {
        self.idxs.len()
    }

    /// The inverse permutation
    ///
    /// Compute the inverse of this permutation. If this permutation maps
    /// element `i` to `j`, then then in the inverse `j` is mapped to `i`.
    pub fn inverse(&self) -> Self
    {
        let mut invidxs = vec![0; self.size()];
        for (i, &pi) in self.idxs.iter().enumerate()
        {
            invidxs[pi] = i;
        }
        Permutation::new(invidxs).unwrap()
    }

    /// Matrix representation.
    ///
    /// Return a matrix representation of the permutation. For a matrix
    /// representation `P` of permutation `perm`, the result of `P.dot(v)` for a
    /// vector `v` is the same as `perm.apply_vec(&mut v)`.
    pub fn matrix<T>(&self) -> ndarray::Array2<T>
    where T: Clone + num_traits::Zero + num_traits::One
    {
        let n = self.size();
        let mut res = ndarray::Array::zeros((n, n));
        for (i, &pi) in self.idxs.iter().enumerate()
        {
            res[[i, pi]] = T::one();
        }
        res
    }

    fn apply<F>(&self, mut swap: F)
    where F: FnMut(usize, usize)
    {
        let mut in_place = vec![false; self.size()];
        for i in 0..self.size()
        {
            if !in_place[i]
            {
                let mut j = i;
                let mut pj = self.idxs[j];
                while pj != i
                {
                    swap(j, pj);
                    in_place[j] = true;
                    j = pj;
                    pj = self.idxs[pj];
                }
                in_place[j] = true;
            }
        }
    }

    /// Permute a vector
    ///
    /// Apply this permutation to vector `v`.
    pub fn apply_vec_in_place<T>(&self, v: &mut ndarray::Array1<T>)
    {
        self.apply(|i, j| v.swap(i, j));
    }

    /// Permute a vector
    ///
    /// Apply this permutation to vector `v`, and store the result in `perm_v`.
    pub fn apply_vec_into<T>(&self, v: ndarray::ArrayView1<T>,
        mut perm_v: ndarray::ArrayViewMut1<T>)
    where T: Copy
    {
        for (new_idx, &old_idx) in self.idxs.iter().enumerate()
        {
            perm_v[new_idx] = v[old_idx];
        }
    }

    /// Permute a vector
    ///
    /// Apply the inverse of this permutation to vector `v`, and store the
    /// result in `inv_perm_v`.
    pub fn apply_inverse_vec_into<T>(&self, v: ndarray::ArrayView1<T>,
        mut inv_perm_v: ndarray::ArrayViewMut1<T>)
    where T: Copy
    {
        for (new_idx, &old_idx) in self.idxs.iter().enumerate()
        {
            inv_perm_v[old_idx] = v[new_idx];
        }
    }

    /// Transform a matrix
    ///
    /// Transform matrix `a`, permuting both its rows and columns by this
    /// permutation.
    pub fn transform<T>(&self, a: &ndarray::Array2<T>) -> ndarray::Array2<T>
    where T: Copy
    {
        a.select(ndarray::Axis(0), &self.idxs).select(ndarray::Axis(1), &self.idxs)
    }
}

#[cfg(test)]
mod tests
{
    use super::Permutation;

    #[test]
    fn test_new()
    {
        let perm = Permutation::new(vec![0, 1, 2, 3]).unwrap();
        assert_eq!(perm.indices(), &[0, 1, 2, 3]);
        let perm = Permutation::new(vec![3, 1, 0, 2]).unwrap();
        assert_eq!(perm.indices(), &[3, 1, 0, 2]);

        assert!(matches!(Permutation::new(vec![]), Err(crate::error::Error::EmptyPermutation)));
        assert!(matches!(Permutation::new(vec![0, 2, 3, 4]), Err(crate::error::Error::InvalidPermutationElement(4, 4))));
        assert!(matches!(Permutation::new(vec![0, 2, 3, 2]), Err(crate::error::Error::DoublePermutationElement(2))));
    }

    #[test]
    fn test_inverse()
    {
        let iperm = Permutation::new(vec![0, 1, 2, 3]).unwrap().inverse();
        assert_eq!(iperm.indices(), &[0, 1, 2, 3]);
        let iperm = Permutation::new(vec![3, 0, 1, 2]).unwrap().inverse();
        assert_eq!(iperm.indices(), &[1, 2, 3, 0]);
        let iperm = Permutation::new(vec![6, 2, 7, 1, 4, 3, 0, 5]).unwrap().inverse();
        assert_eq!(iperm.indices(), &[6, 3, 1, 5, 4, 7, 0, 2]);
    }

    #[test]
    fn test_matrix()
    {
        let perm = Permutation::new(vec![1, 3, 0, 2]).unwrap();
        assert_eq!(perm.matrix::<usize>(), array![
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [1, 0, 0, 0],
            [0, 0, 1, 0]
        ]);
        let perm = Permutation::new(vec![1, 3, 5, 2, 7, 0, 4, 6]).unwrap();
        assert_eq!(perm.matrix::<f64>(), array![
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        ]);

        let perm = Permutation::new(vec![1, 3, 0, 2]).unwrap();
        let v = array![1, 3, 5, 7];
        assert_eq!(perm.matrix().dot(&v), array![3, 7, 1, 5]);

        let perm = Permutation::new(vec![1, 3, 5, 2, 7, 0, 4, 6]).unwrap();
        let v = array![1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9];
        assert_eq!(perm.matrix().dot(&v), array![2.3, 4.5, 6.7, 3.4, 8.9, 1.2, 5.6, 7.8]);
    }

    #[test]
    fn test_apply_vec_in_place()
    {
        let perm = Permutation::new(vec![1, 3, 0, 2]).unwrap();
        let mut v = array![1, 3, 5, 7];
        perm.apply_vec_in_place(&mut v);
        assert_eq!(v, array![3, 7, 1, 5]);

        let perm = Permutation::new(vec![1, 3, 5, 2, 7, 0, 4, 6]).unwrap();
        let mut v = array![1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9];
        perm.apply_vec_in_place(&mut v);
        assert_eq!(v, array![2.3, 4.5, 6.7, 3.4, 8.9, 1.2, 5.6, 7.8]);
    }

    #[test]
    fn test_apply_vec()
    {
        let perm = Permutation::new(vec![1, 3, 0, 2]).unwrap();
        let v = array![1, 3, 5, 7];
        let mut work = array![0, 0, 0, 0];
        perm.apply_vec_into(v.view(), work.view_mut());
        assert_eq!(work, array![3, 7, 1, 5]);

        let perm = Permutation::new(vec![1, 3, 5, 2, 7, 0, 4, 6]).unwrap();
        let v = array![1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9];
        let mut work = array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        perm.apply_vec_into(v.view(), work.view_mut());
        assert_eq!(work, array![2.3, 4.5, 6.7, 3.4, 8.9, 1.2, 5.6, 7.8]);
    }

    #[test]
    fn test_transform()
    {
        let a = array![[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let perm = Permutation::new(vec![2, 0, 1]).unwrap();
        let b = perm.transform(&a);
        assert_eq!(b, array![[9, 7, 8], [3, 1, 2], [6, 4, 5]]);
    }
}
