extern crate num_complex;

use cmatrix;
use gates;

/// Rotation around `y` axis.
///
/// The `R`<sub>`Y`</sub>`(λ)` gate rotates the qubit around the `y` axis of the
/// Bloch sphere over an angle `theta`.
pub struct RY
{
    half_theta: f64,
    desc: String
}

impl RY
{
    /// Create a new `R`<sub>`Y`</sub> gate.
    pub fn new(theta: f64) -> Self
    {
        RY { half_theta: 0.5 * theta, desc: format!("RY({:.4})", theta) }
    }
}

impl gates::Gate for RY
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
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
        let c = num_complex::Complex::new(self.half_theta.cos(), 0.0);
        let s = num_complex::Complex::new(self.half_theta.sin(), 0.0);
        array![[c, -s], [s, c]]
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, RY};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = RY::new(0.21675627161);
        assert_eq!(gate.description(), "RY(0.2168)");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;

        let gate = RY::new(::std::f64::consts::FRAC_PI_2);
        assert_complex_matrix_eq!(gate.matrix(), array![[x, -x], [x, x]]);

        let gate = RY::new(::std::f64::consts::PI);
        assert_complex_matrix_eq!(gate.matrix(), array![[z, -o], [o, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [x, -x, z, o],
            [x,  x, o, z]
        ];
        let gate = RY::new(::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }
}
