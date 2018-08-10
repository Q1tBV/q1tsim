extern crate num_complex;

use cmatrix;
use gates;

/// The `Swap` gate
///
/// The `Swap` gate swap two qubits.
pub struct Swap
{
}

impl Swap
{
    /// Create a `Swap` gate
    pub fn new() -> Self
    {
        Swap { }
    }

    pub fn transform(state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 4 == 0, "Number of rows is not a mutiple of 4.");

        let n = state.len() / 4;
        for i in n..2*n
        {
            state.swap(i, i+n);
        }
    }

    pub fn transform_mat(state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.len() % 4 == 0, "Number of rows is not a multiple of 4.");

        let n = state.rows() / 4;
        let m = state.cols();
        for i in n..2*n
        {
            for j in 0..m
            {
                state.swap((i, j), (i+n, j));
            }
        }
    }
}

impl gates::Gate for Swap
{
    fn cost(&self) -> f64
    {
        3.0 * gates::CX::cost()
    }

    fn description(&self) -> &str
    {
        "Swap"
    }

    fn nr_affected_bits(&self) -> usize
    {
        2
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        array![
            [o, z, z, z],
            [z, z, o, z],
            [z, o, z, z],
            [z, z, z, o]
        ]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        Self::transform(state);
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        Self::transform_mat(state);
    }
}

#[cfg(test)]
mod tests
{
    use gates::{Gate, Swap};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = Swap::new();
        assert_eq!(gate.description(), "Swap");
    }

    #[test]
    fn test_matrix()
    {
        let gate = Swap::new();
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(gate.matrix(), array![
            [o, z, z, z],
            [z, z, o, z],
            [z, o, z, z],
            [z, z, z, o]
        ]);
    }

    #[test]
    fn test_apply_mat()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * o;
        let mut state = array![
            [o, z, x,  x, x,  h, z],
            [z, o, x, -x, z, -h, z],
            [z, z, z,  z, x,  h, o],
            [z, z, z,  z, z, -h, z]
        ];
        Swap::new().apply_mat(&mut state);
        let result = array![
            [o, z, x,  x, x,  h, z],
            [z, z, z,  z, x,  h, o],
            [z, o, x, -x, z, -h, z],
            [z, z, z,  z, z, -h, z]
        ];
        assert_complex_matrix_eq!(&state, &result);
    }
}
