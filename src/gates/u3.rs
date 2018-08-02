extern crate num_complex;

use cmatrix;
use gates;

/// U<sub>3</sub> gate.
///
/// The `U`<sub>`3`</sub>`(θ, ϕ, λ)` gate is the univeral single-qubit
/// transformation, all unary gates can be written in terms of `U`<sub>`3`</sub>.
/// It transforms a qubit by the matrix
/// ```text
/// ┌                                     ┐
/// │       cos(θ/2)      -exp(iλ)sin(θ/2)│
/// │                                     │
/// │exp(iϕ)sin(θ/2)   exp(i(λ+ϕ))cos(θ/2)│
/// └                                     ┘
/// ```
pub struct U3
{
    half_theta: f64,
    phi: f64,
    lambda: f64,
    desc: String
}

impl U3
{
    /// Create a new `U`<sub>`3`</sub> gate.
    pub fn new(theta: f64, phi: f64, lambda: f64) -> Self
    {
        U3
        {
            half_theta: 0.5 * theta,
            phi: phi,
            lambda: lambda,
            desc: format!("U3({:.4}, {:.4}, {:.4})", theta, phi, lambda)
        }
    }

    pub fn cost() -> f64
    {
        201.0
    }
}

impl gates::Gate for U3
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
        let (c, s) = (self.half_theta.cos(), self.half_theta.sin());
        array![[ num_complex::Complex::new(c, 0.0),
                -num_complex::Complex::from_polar(&s, &self.lambda)],
               [ num_complex::Complex::from_polar(&s, &self.phi),
                 num_complex::Complex::from_polar(&c, &(self.phi+self.lambda))]]
    }
}

#[cfg(test)]
mod tests
{
    extern crate num_complex;

    use gates::{gate_test, Gate, U3};
    use cmatrix;
    use self::num_complex::Complex;

    #[test]
    fn test_description()
    {
        let gate = U3::new(0.17, ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_eq!(gate.description(), "U3(0.1700, 0.7854, 0.6931)");
    }

    #[test]
    fn test_matrix()
    {
        let gate = U3::new(0.32, ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [Complex::new(0.9872272833756269,                0.0), Complex::new( -0.1225537622232209, -0.1017981646382380)],
            [Complex::new(0.1126549842634128, 0.1126549842634128), Complex::new(0.09094356700076842,  0.9830294892130130)]
        ]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![[x, x], [z, -x], [x, z], [z, z]];
        let result = array![
            [  -o-i,       -o],
            [     z,        o],
            [ 2.0*x,  (o+i)*x],
            [     z, -(o+i)*x]
        ] * (0.5 * o);
        let gate = U3::new(3.0*::std::f64::consts::FRAC_PI_2,
            ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }
}
