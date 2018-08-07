extern crate num_complex;

use cmatrix;
use gates;

/// Rotation around `x` axis.
///
/// The `R`<sub>`X`</sub>`(θ)` gate rotates the qubit around the `x` axis of the
/// Bloch sphere over an angle `θ`. The associated matrix is
/// ```text
/// ┌                     ┐
/// │  cos(θ/2) -isin(θ/2)│
/// │                     │
/// │-isin(θ/2)   cos(θ/2)│
/// └                     ┘
/// ```
pub struct RX
{
    half_theta: f64,
    desc: String
}

impl RX
{
    /// Create a new `R`<sub>`X`</sub> gate.
    pub fn new(theta: f64) -> Self
    {
        RX { half_theta: 0.5 * theta, desc: format!("RX({:.4})", theta) }
    }
}

impl gates::Gate for RX
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
        let c  = num_complex::Complex::new(self.half_theta.cos(), 0.0);
        let si = num_complex::Complex::new(0.0, self.half_theta.sin());
        array![[c, -si], [-si, c]]
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, RX};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = RX::new(::std::f64::consts::PI);
        assert_eq!(gate.description(), "RX(3.1416)");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;

        let gate = RX::new(::std::f64::consts::FRAC_PI_2);
        assert_complex_matrix_eq!(gate.matrix(), array![[x, -x*i], [-x*i, x]]);

        let gate = RX::new(::std::f64::consts::PI);
        assert_complex_matrix_eq!(gate.matrix(), array![[z, -i], [-i, z]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [   x, -x*i, 0.5*(o-i),  0.5*(o+i)],
            [-x*i,    x, 0.5*(o-i), -0.5*(o+i)]
        ];
        let gate = RX::new(::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }
}
