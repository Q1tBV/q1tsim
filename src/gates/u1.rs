extern crate num_complex;

use cmatrix;
use gates;

/// Phase gate.
///
/// The `U`<sub>`1`</sub>`(λ)` gate shifts the phase of the |1〉 component of a
/// qubit over an angle `λ`.
pub struct U1
{
    lambda: f64,
    desc: String
}

impl U1
{
    /// Create a new `U`<sub>`1`</sub> gate.
    pub fn new(lambda: f64) -> Self
    {
        U1 { lambda: lambda, desc: format!("U1({:.4})", lambda) }
    }

    pub fn cost() -> f64
    {
        7.0
    }
}

impl gates::Gate for U1
{
    fn cost(&self) -> f64
    {
        Self::cost()
    }

    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let p = num_complex::Complex::from_polar(&1.0, &self.lambda);
        array![[o, z], [z, p]]
    }

    fn apply_slice(&self, state: &mut cmatrix::CVecSliceMut)
    {
        assert!(state.len() % 2 == 0, "Number of rows is not even.");

        let n = state.len() / 2;
        let mut slice = state.slice_mut(s![n..]);
        slice *= num_complex::Complex::from_polar(&1.0, &self.lambda);
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, U1};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        assert_eq!(gate.description(), "U1(0.7854)");
    }

    #[test]
    fn test_matrix()
    {
        let gate = U1::new(::std::f64::consts::FRAC_PI_2);
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        assert_complex_matrix_eq!(gate.matrix(), array![[o, z], [z, i]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![[o, z, x, x], [z, o, x, -x]];
        let result = array![[o, z, x, x], [z, x*(o+i), 0.5*(o+i), -0.5*(o+i)]];
        let gate = U1::new(::std::f64::consts::FRAC_PI_4);
        gate_test(gate, &mut state, &result);
    }
}
