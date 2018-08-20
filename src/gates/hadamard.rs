extern crate num_complex;

use cmatrix;
use gates;

/// The Hadamard gate.
///
/// The Hadamard gate maps the zero state |0&rang; to the symmetric combination
/// of |0&rang; and |1&rang;, and the |1&rang; state to the anti-symmetric
/// combination.
pub struct H
{
}

impl H
{
    /// Create a new Hadamard gate.
    pub fn new() -> Self
    {
        H { }
    }

    pub fn transform(state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let s0 = state.slice(s![..n]).to_owned() * cmatrix::COMPLEX_HSQRT2;
        let s1 = state.slice(s![n..]).to_owned() * cmatrix::COMPLEX_HSQRT2;

        state.slice_mut(s![..n]).assign(&(&s0 + &s1));
        state.slice_mut(s![n..]).assign(&(s0 - s1));
    }

    pub fn cost() -> f64
    {
        gates::U2::cost()
    }

}

impl gates::Gate for H
{
    fn cost(&self) -> f64
    {
        Self::cost()
    }

    fn description(&self) -> &str
    {
        "H"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let x = cmatrix::COMPLEX_HSQRT2;
        array![[x, x], [x, -x]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        Self::transform(state);
    }

    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("h {}", bit_names[bits[0]])
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, H};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let h = H::new();
        assert_eq!(h.description(), "H");
    }

    #[test]
    fn test_matrix()
    {
        let h = H::new();
        let s = cmatrix::COMPLEX_HSQRT2;
        assert_complex_matrix_eq!(h.matrix(), array![[s, s], [s, -s]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[x, x, o, z], [x, -x, z, o]];
        gate_test(H::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = H::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "h qb");
    }
}
