extern crate num_complex;

use cmatrix;
use gates;

/// The Pauli Z gate.
///
/// The Z gate rotates the state over π radians around the `z` axis of
/// the Bloch sphere, i.e. it flips the sign of the |1〉 components of the qubit.
pub struct Z
{
}

impl Z
{
    /// Create a new Pauli Z gate.
    pub fn new() -> Self
    {
        Z { }
    }
}

impl gates::Gate for Z
{
    fn cost(&self) -> f64
    {
        gates::U1::cost()
    }

    fn description(&self) -> &str
    {
        "Z"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        array![[o, z], [z, -o]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        state.slice_mut(s![n..]).mapv_inplace(|c| -c);
    }

    fn apply_mat_slice(&self, state: &mut cmatrix::CMatSliceMut)
    {
        assert!(state.rows() % 2 == 0, "Number of rows is not even.");

        let n = state.rows() / 2;
        state.slice_mut(s![n.., ..]).mapv_inplace(|c| -c);
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, Z};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let z = Z::new();
        assert_eq!(z.description(), "Z");
    }

    #[test]
    fn test_matrix()
    {
        let gate = Z::new();
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, -o]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [ z, -o, -x, x]];
        gate_test(Z::new(), &mut state, &result);
    }
}
